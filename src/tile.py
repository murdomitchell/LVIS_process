'''lvis batch processing and tiling.
Contains the Tile class for processing multiple LVIS files in a folder and writing a merged DEM geotiff.'''

import glob
import logging
import pathlib
import numpy as np
import xarray as xr
from pyproj import Transformer

from . import lvis, utils

class Tile():
    '''Tiling LVIS data from a folder'''

    def __init__(self, folderloc, tilesize = 0.1, outEPSG = 3031, max_bounds = None):
        '''Initialise the Tile object by reading the bounds of all LVIS files in the folder and setting up tiling parameters.
        Process the tiles with process(res) and write a merged DEM geotiff with writeDEM(filename).
        args:
            folderloc: location of the folder containing LVIS files
            tilesize: size of each tile in degrees
            outEPSG: EPSG code for the output coordinate system
            max_bounds: optional bounds to limit processing [minX, minY, maxX, maxY]
        returns:
            None, but sets up the Tile object for processing'''

        self.merged_da = None
        self.tilesize = tilesize
        self.outEPSG = outEPSG
        self.max_bounds = np.array(max_bounds)
        if pathlib.Path(folderloc).is_dir():
            self.fileList = glob.glob(folderloc + '/*.h5')
        else:
            raise ValueError(f'{folderloc} is not a valid folder path')

        self._get_all_bounds()

    ### Private Methods ###

    def _get_all_bounds(self):
        '''Get the bounds of all LVIS files in the folder and save to self.all_bounds. 
        Also does the tile info setup'''

        all_bounds =  np.array([360.0, 180.0, 0.0, -180.0])
        logging.debug(all_bounds)
        for i in self.fileList:

            logging.debug(i)
            # get bounds of every file
            clvis = lvis.LVIS(i, onlyBounds = True)
            logging.debug(clvis.bounds)

            all_bounds[:2] = np.minimum(all_bounds[:2], clvis.bounds[:2])
            all_bounds[2:] = np.maximum(all_bounds[2:], clvis.bounds[2:])

        if self.max_bounds.ndim > 0:
            all_bounds[:2] = np.maximum(all_bounds[:2], self.max_bounds[:2])
            all_bounds[2:] = np.minimum(all_bounds[2:], self.max_bounds[2:])

        logging.info('The bounds are: %s', all_bounds)
        # max bounds
        self.all_bounds = all_bounds

        tileamount_x = int((all_bounds[2] - all_bounds[0]) / self.tilesize) + 1
        tileamount_y = int((all_bounds[3] - all_bounds[1]) / self.tilesize) + 1

        # important tiling info 
        logging.warning(f'There are {tileamount_x * tileamount_y} tiles to process,' + 
            fr' {tileamount_x} x by {tileamount_y} y, {self.tilesize} degrees in size.')

    def _reproject(self, bounds, tolatlon = False):
        '''Reproject bounds to or from lat/lon. 
        If tolatlon is True, reprojects from self.outEPSG to lat/lon. 
        If False, reprojects from lat/lon to self.outEPSG. Returns reprojected bounds.'''

        latlon = 'EPSG:4326'
        if tolatlon:
            ba = bounds[2]; bb = bounds[0]
            inEPSG = f'EPSG:{self.outEPSG}'; outEPSG = latlon
            
        else:
            ba = bounds[0]; bb = bounds[2]
            inEPSG = latlon; outEPSG = f'EPSG:{self.outEPSG}'

        transformer = Transformer.from_crs(inEPSG, outEPSG, always_xy = True)
        x0, y0 = transformer.transform(ba, bounds[1])
        x1, y1 = transformer.transform(bb, bounds[3])

        if tolatlon:
            return [x0 % 360, y0, x1 % 360, y1] # back to 0,360 like the data
        return [min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)] # make sure in correct order

    def _iterate(self, res):
        '''Iterate over all tiles, process each tile, and accumulate binned raster data.
        This avoids storing all raw points in memory.'''

        self.res = res
        
        # Pre-compute raster dimensions from all_bounds (reprojected coords)
        bounds_reproj = self._reproject([self.all_bounds[0], self.all_bounds[1], 
                                         self.all_bounds[2], self.all_bounds[3]])
        self.x_min = bounds_reproj[0]
        self.y_max = bounds_reproj[3]
        self.n_x = int(round((bounds_reproj[2] - bounds_reproj[0]) / res)) + 1
        self.n_y = int(round((bounds_reproj[3] - bounds_reproj[1]) / res)) + 1
        n_pixels = self.n_y * self.n_x
        
        logging.info(f'Raster grid: {self.n_x} x {self.n_y} pixels at {res} m resolution')
        
        # Initialize accumulation arrays (use int64 to match np.bincount return type)
        self.sum_flat = np.zeros(n_pixels, dtype = np.float32)
        self.count_flat = np.zeros(n_pixels, dtype = np.int32)
        
        #loop setup
        x_iter = np.arange(self.all_bounds[0], self.all_bounds[2], self.tilesize)
        y_iter = np.arange(self.all_bounds[1], self.all_bounds[3], self.tilesize)

        for x0 in utils.get_tqdm(x_iter, desc = '  x tiles', position = 0):
            x1 = x0 + self.tilesize
            for y0 in utils.get_tqdm(y_iter, desc = '  y tiles', leave = False, position = 1):
                y1 = y0 + self.tilesize

                logging.debug(f'Processing tile between {x0}, {y0} to {x1}, {y1}')
                self._process_tile([x0, y0, x1, y1])
                logging.info(f'Tile between {x0}, {y0} to {x1}, {y1} processed and binned')
        
        if self.count_flat.sum() > 0:
            logging.info(f'Accumulated {self.count_flat.sum():.0f} binned points across raster')
        else:
            self.sum_flat = None
            self.count_flat = None

    def _process_tile(self, bounds):
        '''Single tile processing function. 
        Reads all LVIS files that intersect the tile, processes them, 
        and bins them directly into the raster accumulation arrays.'''

        zG_list = []
        x_list = []
        y_list = []

        for i in self.fileList:

            clvis = lvis.LVIS(i, minX = bounds[0], minY = bounds[1], maxX = bounds[2], maxY = bounds[3], tiling = True)
            
            if clvis.nWaves == 0:
                logging.debug(i + 'does not have data in tile ' + str(bounds))
                continue

            logging.debug(i + 'has data in tile ' + str(bounds))

            clvis.process()
            data, x, y = clvis.getDemData(outEPSG = self.outEPSG)

            # collect data
            zG_list.append(data)
            x_list.append(x)
            y_list.append(y)

        if not x_list:
            logging.debug(f'No data for tile {bounds}')
            return

        # concatenate all data
        zG_all = np.concatenate(zG_list)
        x_all = np.concatenate(x_list)
        y_all = np.concatenate(y_list)
        
        # compute integer pixel indices for each point and ignore out-of-bounds or nodata values
        x_inds = np.floor((x_all - self.x_min) / self.res).astype(int)
        y_inds = np.floor((self.y_max - y_all) / self.res).astype(int)
        
        valid = (
            (x_inds >= 0) & (x_inds < self.n_x)
            & (y_inds >= 0) & (y_inds < self.n_y)
            & np.isfinite(zG_all)
            & (zG_all >= 0.0)
        )
        
        # Flatten 2D pixel indices to 1D, accumulate with bincount
        flat_inds = y_inds[valid] * self.n_x + x_inds[valid]
        n_pixels = self.n_y * self.n_x
        
        sum_tile = np.bincount(flat_inds, weights = zG_all[valid], minlength = n_pixels)
        count_tile = np.bincount(flat_inds, minlength = n_pixels)
        
        # accumulate into class arrays
        self.sum_flat += sum_tile
        self.count_flat += count_tile
        
        logging.debug(f'Binned {valid.sum()} valid points for tile {bounds}')

    ### Public Interface ###

    def process(self, res):
        '''Process all tiles and accumulate binned raster data.
        args:
            res: resolution of the output DEM in meters
        returns:
            None, but accumulates binned data ready for writeDEM()'''

        self._iterate(res)

    def writeDEM(self, filename, interp_dist = 0, method = 'linear'):
        '''Write accumulated binned LVIS data to a geotiff.
        args:
            filename: name of the output geotiff file
            interp_dist: distance for interpolation to fill holes in the DEM in meters. Default 0.
            method: interpolation method ('linear', 'nearest', 'cubic', 'cKDTree', or 'fillnodata')
        returns:
            None, but saves a geotiff of the merged DEM to filename.'''

        if self.sum_flat is None or self.count_flat is None:
            logging.error('No data to write. Call process(res) first. (writeDEM)')
            return
        
        nodata = np.nan
        
        filename = str(filename)
        filename += '.tif' if not filename.endswith('.tif') else ''
        
        # average per pixel
        count_arr = self.count_flat.reshape(self.n_y, self.n_x)
        sum_arr = self.sum_flat.reshape(self.n_y, self.n_x)
        image_arr = np.where(count_arr > 0, sum_arr / count_arr, nodata)

        # interpolation
        image_arr = lvis.interpDEM(image_arr, interp_dist, self.res, method = method)
        
        # uild x/y coordinates for the raster grid
        x_coords = self.x_min + np.arange(self.n_x) * self.res
        y_coords = self.y_max - np.arange(self.n_y) * self.res
        
        # Create xarray and write to geotiff
        image_da = xr.DataArray(
            image_arr,
            dims = ('y', 'x'),
            coords = {'y': ('y', y_coords), 'x': ('x', x_coords)},
            attrs = {'nodata': nodata, 'method': method, 'interp_dist': interp_dist, 'resolution': self.res},
        )
        image_da = image_da.astype('float32')
        logging.info(f'Image array created with shape {image_arr.shape}')
        logging.debug(f'Image has resolution of {image_da.rio.resolution()} and bounds of {image_da.rio.bounds()}')
        
        image_da.rio.write_crs(self.outEPSG, inplace = True)
        image_da.rio.write_nodata(nodata, inplace = True)
        logging.debug(f'CRS set to EPSG:{self.outEPSG} and nodata value set to {nodata}')
        image_da.rio.to_raster(filename, driver = 'GTiff')
        
        logging.info('Geotiff written to %s', filename)
