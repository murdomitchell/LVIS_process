'''Module for comparing two DEMs and calculating elevation differences and volume changes.'''

import logging
from . import utils, lvis

class Ice:

    def __init__(self, demloc1, demloc2):
        '''Initialise the Ice object by loading two DEMs and calculating their 
        difference and volume change in the overlapping region (second minus first).
        args:
            demloc1 (str or xr.DataArray): Path to the first DEM file.
            demloc2 (str or xr.DataArray): Path to the second DEM file.
        returns:
            None, but sets up the Ice object with diff and vol_da attributes for difference and volume'''

        self.da1 = utils.openDEM(demloc1)
        self.da2 = utils.openDEM(demloc2)

        self._diff()
        self._volume()

    ### Private Methods ###

    def _diff(self):
        '''Calculate the elevation difference between two DEMs in their overlapping region.'''

        # Get bounds of both DEMs
        bounds1 = self.da1.rio.bounds()
        bounds2 = self.da2.rio.bounds()
        
        # Find intersection bounds (left, bottom, right, top)
        left = max(bounds1[0], bounds2[0])
        bottom = max(bounds1[1], bounds2[1])
        right = min(bounds1[2], bounds2[2])
        top = min(bounds1[3], bounds2[3])
        
        # Check if there is overlap
        if left >= right or bottom >= top:
            logging.warning('No overlap between DEMs')
            return None
        
        logging.debug(f'Overlapping region bounds: [{left}, {bottom}, {right}, {top}]')
        
        # Clip both DEMs to overlapping region
        da1_clip = self.da1.rio.clip_box(left, bottom, right, top)
        da2_clip = self.da2.rio.clip_box(left, bottom, right, top)
        
        # Ensure both have same grid by reindexing to nearest neighbor
        da2_clip = da2_clip.reindex_like(da1_clip, method = 'nearest')
        
        # Calculate elevation difference (da2 - da1)
        self.diff = da2_clip - da1_clip
        self.da1_clip = da1_clip
        self.da2_clip = da2_clip

    def _volume(self):
        '''Calculate volume change between two DEMs in overlapping region only.
        Returns volume_per_pixel_da: DataArray with volume change in m^3 for each pixel'''
        
        # Get pixel area from resolution
        res = self.da1_clip.rio.resolution()
        pixel_area = abs(res[0] * res[1])
        
        logging.debug(f'Pixel area: {pixel_area:.2f} m^2')
        logging.debug(f'DEM1 shape: {self.da1_clip.shape}, DEM2 shape: {self.da2_clip.shape}')
        
        # Calculate volume per pixel (difference * pixel area)
        vol_per_pixel = self.diff * pixel_area
        try:
            vol_per_pixel = vol_per_pixel.drop_attrs()
        except:
            pass
        vol_per_pixel.attrs['units'] = 'm^3'
        vol_per_pixel.attrs['description'] = 'Volume change per pixel'
        
        # Calculate total volume (sum of differences * pixel area)
        volume_change = float((vol_per_pixel).sum(skipna = True).values)
        volume_desc = f'Volume change: {volume_change:.2e} m^3 ({volume_change / 1e9:.3f} km^3)'
        
        logging.info(volume_desc)
        logging.info(f'Volume per pixel DataArray created with shape {vol_per_pixel.shape}')

        vol_per_pixel.attrs['volume_desc'] = volume_desc
        vol_per_pixel.attrs['total_change'] = volume_change / 1e9
        self.vol_da = vol_per_pixel

    ### Public Interface ###

    def diff_array(self):
        '''Returns the elevation difference DataArray.'''

        return self.diff

    def volume_array(self):
        '''Returns the volume change DataArray.'''

        return self.vol_da
    
    def save_diff(self, filename = 'diff'):
        '''Save the elevation difference DataArray as a geotiff.
        args:
            filename: name of the output geotiff file'''
        
        lvis.Post(self.diff).writeTiff(filename)

    def save_volume(self, filename = 'vol'):
        '''Save the volume change DataArray as a geotiff.
        args:
            filename: name of the output geotiff file'''
        
        lvis.Post(self.vol_da).writeTiff(filename)
    
    def plot_diff(self, mapda = True, outName = 'diff', title = 'Elevation Difference'):
        '''Plot the elevation difference between two DEMs in their overlapping region.
        args:
            mapda: whether to plot the data with coastlines
            outName: name of the output image file (without extension)
            title: title of the plot'''

        lvis.Post(self.diff).plot(outName = outName, mapda = mapda, title = title, cmap = 'RdBu_r')
    
    def plot_volume(self, mapda = True, outName = 'vol', title = None):
        '''Plot the volume change between two DEMs in their overlapping region.
        args:
            mapda: whether to plot the data with coastlines
            outName: name of the output image file (without extension)
            title: title of the plot'''

        thetitle = title or fr'{self.vol_da.attrs["total_change"]:.3f} km^3'
        lvis.Post(self.vol_da).plot(outName = outName, mapda = mapda, title = thetitle, cmap = 'RdBu_r')