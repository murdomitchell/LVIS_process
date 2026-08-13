'''lvis data handler.
Contains the LVIS class, Post class, and interpDEM function.'''

import logging
import pathlib
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from pyproj import Transformer
from scipy.ndimage.filters import gaussian_filter1d 
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
from rasterio.fill import fillnodata

from . import utils

class LVIS():
    '''
    LVIS data handler
    '''

    def __init__(self, fileloc, minX = -np.inf, maxX = np.inf, tiling = False,
        minY = -np.inf, maxY = np.inf, onlyBounds = False):
        '''
        Read LVIS data from fileloc and save to self. Only reads data within bounds. 
        Futher processing can then be done with process() to estimate ground elevation, 
        and writeDEM() to write a geotiff of the ground elevation estimates.
        getDemData() can be used to get the x, y, zG arrays of the estimated ground points for further processing or custom plotting.
        plotWave() can be used to plot a single waveform for inspection.
        args:
            fileloc: location of the input LVIS file
            minX, maxX, minY, maxY: bounds to read data within (in lat/lon) If not specified, reads all data
            tiling: whether this LVIS object is being used for tiling. If True, will disable some logging and tqdm bars for cleaner output. Default False.
            onlyBounds: whether to only save the bounds of the data then exit. Default False.
        '''
        
        self.tiling = tiling
        
        if not pathlib.Path(fileloc).is_file():
            raise ValueError(f'{fileloc} is not a valid file path')

        # call the file reader and load in to the self
        self._readLVIS(fileloc, minX, minY, maxX, maxY, onlyBounds)
        
        if not onlyBounds:
            if self.nWaves > 0:
                logging.info(f'LVIS data read from {fileloc}')
            else:
                logging.debug(f'LVIS data read from {fileloc}')

    ### Private Methods ###

    def _readLVIS(self, fileloc, minX, minY, maxX, maxY, onlyBounds):
        '''
        Read LVIS data from file
        '''

        with xr.open_dataset(fileloc, engine = 'netcdf4') as ds:
            wave_dim, bin_dim = ds['RXWAVE'].dims
            self.nBins = ds['RXWAVE'].sizes[bin_dim]
            last_idx = self.nBins - 1

            ds = ds.assign(
                tempLon = (ds['LON0'] + ds[f'LON{last_idx}']) / 2.0,
                tempLat = (ds['LAT0'] + ds[f'LAT{last_idx}']) / 2.0,
            )

            if onlyBounds:
                self._saveBounds(ds['tempLon'], ds['tempLat'])
                return

            mask = (
                (ds['tempLon'] >= minX)
                & (ds['tempLon'] < maxX)
                & (ds['tempLat'] >= minY)
                & (ds['tempLat'] < maxY)
            )
            valid = np.nonzero(mask.values)[0]

            if valid.size == 0:
                logging.debug('No data contained in that region')
                self.nWaves = 0
                return

            subset = ds.isel({wave_dim: valid})
            self.nWaves = subset.sizes[wave_dim]

            self._saveBounds(subset['tempLon'], subset['tempLat'])
            self.lfid = subset['LFID'].values
            self.lShot = subset['SHOTNUMBER'].values
            self.waves = subset['RXWAVE'].rename({wave_dim: 'wave', bin_dim: 'bin'})
            self.nBins = self.waves.sizes['bin']
            self.lZ0 = subset['Z0'].values
            self.lZN = subset[f'Z{last_idx}'].values

    def _saveBounds(self, tlon, tlat):
        '''save lat lon and bounds '''

        self.lon = np.asarray(tlon)
        self.lat = np.asarray(tlat)
        self.bounds = [np.min(self.lon), np.min(self.lat),
                       np.max(self.lon), np.max(self.lat)]
        
    def _estimateGround(self):
        '''
        Processes waveforms to estimate ground
        Only works for bare Earth. DO NOT USE IN TREES
        '''

        # find resolution
        res = (self.lZ0[0] - self.lZN[0]) / self.nBins
        bin_index = np.arange(self.nBins)
        self.zG = np.full((self.nWaves), -999.0)

        if self.tiling:
            itert = utils.get_tqdm(range(self.nWaves), desc = '  Estimating zG', leave = False, position = 2)
        else:
            itert = utils.get_tqdm(range(self.nWaves), desc = '  Estimating zG')
        
        for i in itert:

            #denoise waveform data
            row = self.waves.isel(wave = i).values.astype(float)
            row -= self.meanNoise[i]
            row[row < self.threshold[i]] = 0.0

            # minimum acceptable width: remove interior positive bins that are not consecutive
            valid_bins = np.nonzero(row > 0.0)[0]
            if valid_bins.shape[0] > 2:
                interior = valid_bins[1:-1]
                prev_bins = valid_bins[:-2]
                next_bins = valid_bins[2:]
                bad = (interior != prev_bins + 1) | (interior != next_bins - 1)
                row[interior[bad]] = 0.0

            #smooth
            row = gaussian_filter1d(row, self.sWidth / res)

            #cofg
            if row.sum() > 0.0:  # avoid empty waveforms (clouds etc)
                # Decodes LVIS's RAM efficient elevation on the fly
                z0 = self.lZ0[i]; zn = self.lZN[i]
                z = z0 + (zn - z0) * bin_index / (self.nBins - 1)
                self.zG[i] = np.average(z, weights = row) # centre of gravity
            # zG: array of ground elevation estimates

        logging.info(f'Denoised {self.nWaves} waves')

    def _findStats(self):
        '''Finds standard deviation and mean of noise'''

        # determine number of bins to calculate stats over
        res = float((self.lZ0[0] - self.lZN[0])) / self.nBins    # range resolution
        noiseBins = int(self.statsLen / res)   # number of bins within 'self.statsLen'

        # no loop needed now!
        noise_slice = self.waves.isel(bin = slice(0, noiseBins))
        self.meanNoise = noise_slice.mean(dim = 'bin').values
        self.stdevNoise = noise_slice.std(dim = 'bin').values

        # Set a noise threshold
        self.threshold = self.meanNoise + self.sigThresh * self.stdevNoise

        if not self.tiling:
            logging.info('Noise stats calculated')

    def _reprojectLVIS(self, dobounds = True):
        '''A method to reproject the footprint coordinates'''

        transformer = Transformer.from_crs('EPSG:4326', f'EPSG:{self.outEPSG}', always_xy = True)
        self.x, self.y = transformer.transform(self.lon, self.lat)

        if dobounds:
            self.bounds = [self.x.min(), self.y.min(), self.x.max(), self.y.max()]

        if not self.tiling:
            logging.info('Coordinates reprojected to EPSG:' + str(self.outEPSG))

    def _createTiff(self, res, filename, outEPSG, interp_dist, method = 'linear'):
        '''Make a geotiff from an array of points using xarray.
        args:
            res: resolution of the output geotiff in meters
            filename: name of the output geotiff file
            outEPSG: EPSG code to reproject coordinates to
            interp_dist: maximum distance to interpolate missing values (m)
            method: interpolation method ('linear', 'nearest', 'cubic', 'cKDTree', or 'fillnodata')
        returns:
            None, but writes a geotiff to filename'''

        nodata = np.nan

        filename += '.tif' if not filename.endswith('.tif') else ''

        data = np.asarray(self.zG)
        x = np.asarray(self.x)
        y = np.asarray(self.y)

        # determine bounds and image size
        n_x = int(round((x.max() - x.min()) / res)) + 1
        n_y = int(round((y.max() - y.min()) / res)) + 1
        logging.debug(f'Raster grid: {n_x} x {n_y} pixels at {res} m resolution')

        # compute integer pixel indices for each point and ignore out-of-bounds or nodata values
        x_inds = np.floor((x - x.min()) / res).astype(int)
        y_inds = np.floor((y.max() - y) / res).astype(int)

        valid = (
            (x_inds >= 0) & (x_inds < n_x)
            & (y_inds >= 0) & (y_inds < n_y)
            & np.isfinite(data)
            & (data >= 0.0)
        )
        logging.debug(f'{valid.sum()} valid points out of {len(data)} total points')

        if valid.sum() == 0:
            raise ValueError('No valid points')

        # Flatten 2D pixel indices to 1D, accumulate with bincount
        flat_inds = y_inds[valid] * n_x + x_inds[valid]
        n_pixels = n_y * n_x

        sum_flat  = np.bincount(flat_inds, weights = data[valid], minlength = n_pixels)
        count_flat = np.bincount(flat_inds, minlength = n_pixels)

        # average per pixel
        count_arr = count_flat.reshape(n_y, n_x)
        sum_arr   = sum_flat.reshape(n_y, n_x)
        image_arr = np.where(count_arr > 0, sum_arr / count_arr, nodata)

        image_arr = interpDEM(image_arr, interp_dist, res, method = method)

        # build x/y coordinates for the raster grid
        x_coords = x.min() + np.arange(n_x) * res
        y_coords = y.max() - np.arange(n_y) * res

        # Explicitly assign the dimension names in tuples to get coordinate binding
        image_da = xr.DataArray(
            image_arr,
            dims = ('y', 'x'),
            coords = {'y': ('y', y_coords), 'x': ('x', x_coords)},
            attrs = {'nodata': nodata, 'method': method, 'interp_dist': interp_dist, 'resolution': res},
        )
        image_da = image_da.astype('float32')
        logging.info(f'Image array created with shape {image_arr.shape}')
        logging.debug(f'Image has resolution of {image_da.rio.resolution()} and bounds of {image_da.rio.bounds()}')

        image_da.rio.write_crs(outEPSG, inplace = True)
        image_da.rio.write_nodata(nodata, inplace = True)
        logging.debug(f'CRS set to EPSG:{outEPSG} and nodata value set to {nodata}')
        image_da.rio.to_raster(filename, driver = 'GTiff')

        logging.info('Geotiff written to %s', filename)

    ### Public Interface ###
    
    def process(self, sigThresh = 5, statsLen = 10, sWidth = 0.5):
        '''Does the processing steps to reproject the data and then estimates ground elevation
        args:
            sigThresh: number of stdev above mean to set noise threshold
            statsLen: length of the noise statistics window
            sWidth: width of the Gaussian filter
        returns: 
            None, but sets self.zG to the estimated ground elevation'''
        
        self.sigThresh = sigThresh
        self.statsLen = statsLen
        self.sWidth = sWidth

        self._findStats()
        self._estimateGround()

    def getDemData(self, outEPSG = 3031):
        '''Returns x, y, zG arrays for the estimated ground points.
        args:
            outEPSG: EPSG code to reproject coordinates to 
        returns:
            data: array of ground elevation estimates
            x: array of x coordinates
            y: array of y coordinates'''

        if not hasattr(self, 'zG'):
            logging.error('Ground elevation has not been estimated yet. Call process() first. (getDemData)')
            return None, None, None

        self.outEPSG = outEPSG
        self._reprojectLVIS()

        data = np.asarray(self.zG)
        x = np.asarray(self.x)
        y = np.asarray(self.y)

        return data, x, y

    def writeDEM(self, res, filename = 'temp', outEPSG = 3031, interp_dist = 0, method = 'linear'):
        '''Write LVIS ground elevation data to a geotiff.
        args:
            res: resolution of the output geotiff in meters
            filename: name of the output geotiff file
            outEPSG: EPSG code to reproject coordinates to
            interp_dist: distance for interpolation to fill gaps in the output geotiff
        returns:
            None, but writes a geotiff to filename'''

        if not hasattr(self, 'zG'):
            logging.error('Ground elevation has not been estimated yet. Call process() first. (writeDEM)')
            return
        
        self.outEPSG = outEPSG
        self._reprojectLVIS()

        self._createTiff(res, filename = filename, outEPSG = outEPSG, interp_dist = interp_dist, method = method)

    def plotWave(self, i, outName = 'waveform'):
        '''A method to plot a single waveform
        args:
            i: index of the waveform to plot
            outName: name of the output image file (without extension)
        returns:
            None, but saves a plot of the waveform to outName.i.png'''

        utils.pltsetup()
        outName = outName + '-' + str(i) + '.png'

        wave = self.waves.isel(wave = i).values
        z0 = self.lZ0[i]
        zn = self.lZN[i]
        elev = np.linspace(z0, zn, self.nBins)

        plt.plot(wave, elev)
        plt.xlabel('Waveform return')
        plt.ylabel('Elevation (m)')
        plt.savefig(outName)
        plt.close()
        logging.info(f'Plotted wave {i} to {outName}')

class Post:

    def __init__(self, fileloc):
        '''A class to handle post-processing of DEMs, including interpolation and plotting.
        args:
            fileloc: location of the input DEM file or an xarray DataArray'''

        da = utils.openDEM(fileloc)

        self.da = da

    ### Public Interface ###

    def writeTiff(self, filename, outEPSG = 3031):
        '''Write the DEM to a geotiff.
        args:
            res: resolution of the output geotiff in meters
            filename: name of the output geotiff file
            outEPSG: EPSG code to reproject coordinates to
        returns:
            None, but writes a geotiff to filename'''
        
        filename = str(filename)
        filename += '.tif' if not filename.endswith('.tif') else ''

        self.da.rio.write_crs(outEPSG, inplace = True)
        logging.debug(f'CRS set to EPSG:{outEPSG} and nodata value set to {self.da.rio.nodata}')
        self.da.rio.to_raster(filename, driver = 'GTiff')

        logging.info('Geotiff written to %s', filename)

    def interpolate(self, pixel_dist, method = 'linear'):
        '''Interpolate to fill holes in the DEM.
        args:
            pixel_dist: distance for interpolation to fill holes in the DEM in pixels
            method: interpolation method ('linear', 'nearest', 'cubic', 'cKDTree', or 'fillnodata')
        returns:
            None, but sets self.da to the interpolated DEM'''

        before = np.isnan(self.da.values).sum()
        self.da.values = interpDEM(self.da.values, pixel_dist = pixel_dist, method = method)
        after = np.isnan(self.da.values).sum()
        logging.info(f'Gap filling applied: {before - after} pixels filled')

    def plot(self, outName, mapda = True, title = 'DEM', cmap = 'plasma'):
        '''A method to plot the post-processed DEM
        args:
            outName: name of the output image file (without extension)
            mapda: whether to plot the data with coastlines
            title: title of the plot
            cmap: colormap to use for the plot
        returns:
            None, but saves a plot of the DEM to outName.png'''

        utils.pltsetup()

        if mapda:
            try:
                epsg_code = self.da.rio.crs.to_epsg()
                if epsg_code == 3031 or epsg_code is None:
                    proj = ccrs.Stereographic(central_longitude = 0, central_latitude = -90, true_scale_latitude = -71)    
                else:
                    proj = ccrs.epsg(epsg_code)
                logging.debug(f'Plotting with EPSG code: {epsg_code}')
            except:
                proj = ccrs.Stereographic(central_longitude = 0, central_latitude = -90, true_scale_latitude = -71)    
                logging.debug(f'Plotting without EPSG code')

            ax = plt.axes(projection = proj)
            self.da.plot.imshow(ax = ax, add_colorbar = True, cmap = cmap, transform = proj)
            ax.coastlines()
            ax.gridlines(draw_labels = True, dms = True, x_inline = False, y_inline = False)
        else:
            self.da.plot.imshow(add_colorbar = True, cmap = cmap)
            plt.axis('off')

        plt.title(title)
        plt.savefig(outName)
        logging.info(f'Saved as {outName}.png')
        plt.close('all')

def interpDEM(image_arr, interp_dist = 0, res = 1, pixel_dist = None, method = 'linear'):
    '''Interpolate to fill holes in the DEM.
    args:
        image_arr: 2D array of elevation values with NaNs for holes
        interp_dist: distance for interpolation to fill holes in the DEM in meters
        res: resolution of the DEM in meters per pixel (used to convert interp_dist to pixels)
        pixel_dist: distance for interpolation in pixels (overrides interp_dist/res if provided)
        method: interpolation method ('linear', 'nearest', 'cubic', 'cKDTree', or 'fillnodata')'''

    if interp_dist <= 0 and pixel_dist is None:
        logging.debug(f'No gap filling applied to the image with interp_dist = {interp_dist}')
        return image_arr
        
    # Calculate search distance in pixels (dist / resolution)
    if pixel_dist:
        dist_threshold_px = pixel_dist
    else:
        dist_threshold_px = interp_dist / res

    # Create a boolean mask where True/1 is valid data to interpolate FROM
    valid = ~np.isnan(image_arr)
    rows, cols = np.indices(image_arr.shape)

    if method == 'cKDTree':
    
        power = 2
        k = 64
        dist_from_valid = distance_transform_edt(~valid)
        to_fill = (~valid) & (dist_from_valid <= dist_threshold_px)

        valid_coords = np.column_stack([rows[valid], cols[valid]])
        valid_values = image_arr[valid]
        tree = cKDTree(valid_coords)

        fill_coords = np.column_stack([rows[to_fill], cols[to_fill]])

        # Query k nearest neighbours for every fill pixel at once: fully vectorised
        dists, idxs = tree.query(fill_coords, k = k, workers = -1) # all workers

        # Zero out neighbours beyond the radius
        dists = np.where(dists <= dist_threshold_px, dists, np.nan)
        dists = np.maximum(dists, 1e-10)

        weights = np.where(np.isnan(dists), 0.0, 1.0 / dists ** power)
        weight_sums = weights.sum(axis = 1, keepdims = True)

        neighbour_vals = valid_values[idxs]
        neighbour_vals = np.where(np.isnan(dists), 0.0, neighbour_vals)

        filled_values = (weights * neighbour_vals).sum(axis = 1) / weight_sums.squeeze()

        result = image_arr.copy()
        result[to_fill] = filled_values

    elif method == 'fillnodata':
        
        result = fillnodata(
            image_arr,
            mask = valid,
            max_search_distance = dist_threshold_px,
            smoothing_iterations = 2 # blend
        )

    elif method in ['linear', 'nearest', 'cubic']:

        filled = griddata(
            points = (rows[valid], cols[valid]),
            values = image_arr[valid],
            xi = (rows, cols),
            method = method,
            fill_value = np.nan
        )

        # Mask out any filled pixel that's too far from a known point
        dist_from_valid = distance_transform_edt(~valid)
        filled[dist_from_valid > dist_threshold_px] = np.nan
        result = filled

    else:
        logging.error('interpolation method should be one of "cKDTree", ' +
            f'"fillnodata", "linear", "nearest", or "cubic". Got {method}')
        raise ValueError(f'Invalid interpolation method: {method}')

    logging.info(f'Gap filling to {dist_threshold_px:.0f} pixels away ({interp_dist} m)')

    return result