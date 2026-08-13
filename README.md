# LVIS Processing Pipeline 2026

A modular Python project designed to process **NASA Land, Vegetation, and Ice Sensor (LVIS)** L1B data (HDF5/NetCDF). This pipeline estimates ground elevation using waveform analysis, generates Digital Elevation Models (DEMs) as GeoTIFFs, and provides tools for spatial tiling and temporal comparison.

## Background

NASA's LVIS instrument returns raw airborne lidar waveforms, not ready-to-use elevation data — turning that raw sensor output into a usable elevation grid normally means a lot of manual, file-by-file processing. This project automates that: it reads raw LVIS waveforms, estimates ground elevation, and stitches large batches of files into a single gridded DEM, with tools to compare elevation change between two time periods (e.g. ice volume change over a glacier). Built as an independent project during my Earth Physics Master's at the University of Edinburgh (Mar–Apr 2026).

<img width="1443" height="1177" alt="big15" src="https://github.com/user-attachments/assets/f8e36b05-7365-4ba1-bce6-9cad9abd3d00" />

## Features

* **Waveform Processing**: Automated denoising and ground elevation estimation using a Center of Gravity (COFG) algorithm.
* **Spatial Tiling**: Seamlessly process and merge large folders of LVIS files into a single, high-resolution DEM.
* **Gap Filling**: Multiple interpolation methods (`linear`, `nearest`, `cubic`, `cKDTree`, `fillnodata`) to handle data gaps.
* **DEM Comparison**: Dedicated `Ice` class to calculate elevation differences and volume changes between two datasets.
* **Visualization**: High-quality plotting tools for individual waveforms and spatial DEM maps with optional coastline overlays.
* **Diagnostics**: Integrated hardware resource monitoring (RAM/CPU/Disk) and smart progress bars.

## Quick Start

1.  **Locate your data**: Ensure your LVIS `.h5` files are in a reachable directory.
2.  **Run the pipeline**: Execute `run.py`. The script will automatically create the `dems/` and `plots/` directories if they do not exist.

```bash
# Process all files in the default directory
python3 run.py
```

## Requirements

The project requires **Python 3.8+** and the following libraries:

```bash
pip install numpy xarray scipy rioxarray rasterio pyproj matplotlib cartopy psutil tqdm netCDF4
```

If running with the command line a `dems/` and `plots/` folder are created automatically.

## File Structure

```
lvis/
  run.py           — File to run that does the processing
  task.py          — Hardcoded tasks for ease of running
  dems/            — Stores the DEM geotiffs
  plots/           — Stores output plots
  src/
    __init__.py        — Makes it a module
    lvis.py            — Read the raw data in and then produces a DEM out of it
    tile.py            — Process multiple data files and then tile a large area
    compare.py         — Comparing two DEMs and calculating differences
    utils.py           — Utility functions
```

<details>
<summary><strong>Command Line Options</strong> (click to expand)</summary>

### General Settings
* **`-h, --help`**: Show help message and exit.
* **`-db, --debug`**: Enable `DEBUG` logging level for verbose output.
* **`-tw, --testwave`**: Plots a sample waveform from the first available file and exits.

### Input/Output Configuration
* **`-d, --directory <PATH>`**: Folder containing LVIS files. *Default:* `./data/lvis/`.
* **`-f, --filename <NAME>`**: Process a specific file within the directory.
* **`-dn, --demname <NAME>`**: Set the output filename (`all` saves to `dems/all.tif`).
* **`-pn, --plotname <NAME>`**: Set the output plot filename (saves to `plots/<NAME>.png`).

### Processing Parameters
* **`-r, --resolution <float>`**: DEM resolution in meters. *Default:* `100`.
* **`-t, --tilesize <float>`**: Spatial tile size in degrees for batch processing. *Default:* `0.1`.
* **`-gf, --gapfill <int>`**: Max distance (meters) to interpolate across gaps. *Default:* `0`.
* **`-pg, --pineglacier`**: Constrains processing to Pine Island Glacier bounds.
* **`-np, --noplotdem`**: Skip the generation of a final DEM map plot.
* **`-sp, --simpleplot`**: Disables coastline overlays and coordinate projection on maps.

</details>

<details>
<summary><strong>Usage Examples</strong> (click to expand)</summary>

### Command Line Usage

1. Reading in an LVIS L1B file and plotting a graph of an arbitrary waveform
    - For the first file in the default folder:
    ```bash
    python3 run.py --testwave
    ```
    - For a specific file in a specific folder:
    ```bash
    python3 run.py --testwave -d <dir_path> -f <file_name>
    ```
2. Process a LVIS file in to a DEM
    - With the default folder and other defaults:
    ```bash
    python3 run.py -f <file_name>
    ```
    - With a specific folder, resolution of 30, and simple plots, and setting the name of the dem to "res30":
    ```bash
    python3 run.py -d <dir_path> -f <file_name> -r 30 -sp -dn res30
    ```
3. Process a whole folder of LVIS files into a singular DEM
    - With the default folder and other defaults:
    ```bash
    python3 run.py
    ```
    - With a specific folder, a subset of Pine Island Glacier, a tilesize of 0.5 degrees, and no plotting:
    ```bash
    python3 run.py -d <dir_path> -pg -t 0.5 -np
    ```
4. Use a gap filling algorithm for the output dem
    - From raw data, same as 2. and 3. except with `-gf <distance>`. So here 500 meters:
    ```bash
    python3 run.py -gf 500
    ```
    - From a previously generated DEM, with a specific interpolation method, or for reading processed DEMs and analysing the difference: these are post-processing steps, done with Python directly (see below) rather than `run.py`, which is for initial processing only.

Or of course `task.py` can be used to run all tasks like
```bash
python3 task.py
```
Or a specific task like
```bash
python3 task.py -t <task_no>
```

### Python Usage

With importing the package the Classes themselves can be used for greater control:
```python
from src import lvis, tile, compare
```
Having a plots and dems folder is not hardcoded in `/src` so it is not necessary here, but these examples stick to that convention anyway.

1. Reading in an LVIS L1B file and plotting a graph of an arbitrary waveform
    - For a specific file in a specific folder (can change the waveform by changing 0):
    ```python
    current_lvis = lvis.LVIS(full_file_location)
    current_lvis.plotWave(0, outName = 'plots/testwave')
    ```
2. Process a LVIS file in to a DEM
    - With filling holes, a resolution of 30, and gap filling to 200 meters:
    ```python
    one_file = lvis.LVIS(full_file_location)
    one_file.process()
    one_file.writeDEM(30, filename = 'dems/one', interp_dist = 200)
    ```
    - Plotting the DEM:
    ```python
    lvis.Post('dems/one').plot(outName = 'plots/one')
    ```
3. Process a whole folder of LVIS files into a singular DEM
    - Using a tilesize of .5 degrees, and a resolution of 100:
    ```python
    all_files = tile.Tile(full_folder_location, tilesize = .5)
    all_files.process(res = 100)
    all_files.writeDEM(filename = 'dems/all.tif')
    ```
    - Plotting the DEM with simpleplot:
    ```python
    lvis.Post('dems/all').plot(outName = 'plots/all', mapda = False)
    ```
4. Use a gap filling algorithm for the output dem
    - From raw data, do as just shown in 2.
    - From a previously generated DEM (note that the distance is raster pixels instead of meters):
    ```python
    dem = lvis.Post('dems/all')
    dem.interpolate(pixel_dist = 10)
    dem.plot(outName = 'plots/all', mapda = False)
    ```
    - It can then be saved to a new tiff:
        ```python
        dem.writeTiff('dems/all_interpolated')
        ```
    - With a specific interpolation method, one of ('linear', 'nearest', 'cubic', 'cKDTree', 'fillnodata'):
    ```python
    dem = lvis.Post('dems/all')
    dem.interpolate(pixel_dist = 10, method = 'linear')
    ```
5. Read processed DEMs and perform analysis on the difference
    - For getting the volume change (dem15 minus dem09):
    ```python
    ice = compare.Ice('dems/dem09', 'dems/dem15')
    print(ice.volume_array().attrs['total_change'])
    ```
    - For saving the two arrays as tiffs:
    ```python
    ice = compare.Ice('dems/dem09', 'dems/dem15')
    ice.save_diff(filename = 'dems/diff')
    ice.save_volume(filename = 'dems/volchng')
    ```
    - For plotting the two arrays appropriately:
    ```python
    ice = compare.Ice('dems/dem09', 'dems/dem15')
    ice.plot_diff(outName = 'plots/diff')
    ice.plot_volume(outName = 'plots/volchng')
    ```

</details>
