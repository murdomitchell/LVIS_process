'''Utility functions for LVIS. Contains get_tqdm, pltsetup, and resource_usage functions.'''

import os
import sys
import psutil
import time
import tqdm
import logging
import pathlib
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray

# Import resource module for Unix/Linux peak memory
try:
    import resource
except ImportError:
    resource = None

def get_tqdm(*iterable, **kwargs):
    '''
    Smart tqdm wrapper that disables the progress bar when logging level is below WARNING.
    Works with both positional and keyword arguments.
    '''

    level = logging.getLogger().getEffectiveLevel()
    disable = level < logging.WARNING
    
    # Default disable
    if 'disable' not in kwargs:
        kwargs['disable'] = disable
    
    # Default width
    if 'ncols' not in kwargs:
        kwargs['ncols'] = 100
    
    if iterable:
        return tqdm.tqdm(iterable[0], **kwargs)
    else:
        return tqdm.tqdm(**kwargs)
    
def pltsetup():
    '''my matplotlib settings'''

    plt.rcParams['axes.prop_cycle'] = plt.cycler(color = # specify by c = 'C0'
        ['#332288', '#117733', '#AA4499', '#999933', 
        '#44AA99', '#882255', '#88CCEE', '#CC6677', '#DDCC77'])
    # cmaps = 'plasma' # other choices: 'afmhot' 'ocean' 'nipy_spectral' 'plasma' 'RdBu_r'
    plt.rcParams['image.cmap'] = 'plasma'

    # normal latex doc is 6.69291" wide
    factor =  1.6 # .9 for double latex
    plt.rcParams['figure.figsize'] = (5 * factor, 3 * factor)
    plt.rcParams['figure.dpi'] = 500
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.pad_inches'] = 0
    plt.rcParams['savefig.bbox'] = 'tight'

def openDEM(fileloc):
    '''Open a DEM file as an xarray DataArray or passes through an existing DataArray.
    args:
        fileloc: location of the input DEM file or an xarray DataArray'''

    if isinstance(fileloc, str) or isinstance(fileloc, pathlib.Path):

        fileloc = str(fileloc)
        fileloc += '.tif' if not fileloc.endswith('.tif') else ''

        if not pathlib.Path(fileloc).is_file():
            raise ValueError(f'{fileloc} is not a valid file path')
        if sys.platform == 'win32':
            # with modern rioxarray
            da = rioxarray.open_rasterio(fileloc, masked = True)
        else:
            # have to do this because the linux cluster has an old rioxarray
            da = xr.open_rasterio(fileloc)

        # Remove the band dimension
        da = da.squeeze()

    else:
        da = fileloc

    return da

def resource_usage(start_time):
    '''Prints a summary of the current process's resource usage, 
    including peak and current RAM usage, CPU time, CPU percentage, and disk stats.'''

    # Get the current process
    process = psutil.Process(os.getpid())

    # Memory info handling based on platform
    mem_info = process.memory_info()
    current_ram_mb = mem_info.rss / (1024 * 1024)
    
    peak_ram_mb = None
    if sys.platform == 'win32':
        # On Windows with modern psutil, peak_wset is available
        if hasattr(mem_info, 'peak_wset'):
            peak_ram_mb = mem_info.peak_wset / (1024 * 1024)
    else:
        # On Linux/Unix, use resource module for peak memory
        if resource is not None:
            usage = resource.getrusage(resource.RUSAGE_SELF)
            # ru_maxrss is in kilobytes on Linux
            peak_ram_mb = usage.ru_maxrss / 1024

    # ── Gather ALL processes (main + children) ──
    try:
        children = process.children(recursive = True)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        children = []

    all_procs = [process] + children

    # ── CPU times (sum across main + all children) ──
    total_user = 0.0
    total_system = 0.0

    for p in all_procs:
        try:
            cpu = p.cpu_times()
            total_user += cpu.user
            total_system += cpu.system
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

    total_cpu = total_user + total_system

    # ── CPU % (average over the whole run) ──
    # We calculate it manually from total CPU time + wall time
    wall_time = time.time() - start_time
    cpu_percent_avg = (total_cpu / wall_time * 100) if wall_time > 0 else 0.0

    # Disk stats (current working directory)
    disk_info = psutil.disk_usage('.')
    disk_total_gb = disk_info.total / (1024 * 1024 * 1024)
    disk_used_gb = disk_info.used / (1024 * 1024 * 1024)
    disk_free_gb = disk_info.free / (1024 * 1024 * 1024)
    disk_percent = disk_info.percent

    print('\n' + '=' * 60)
    print('RESOURCE USAGE SUMMARY')
    print('=' * 60)
    if peak_ram_mb is not None:
        print(f'Peak RAM usage                  : {peak_ram_mb:8.2f} MB')
    print(f'Current RAM usage (RSS)         : {current_ram_mb:8.2f} MB')
    print(f'Total CPU time                  : {total_cpu:8.2f} seconds')
    print(f'   └─ User time                 : {total_user:8.2f} seconds')
    print(f'   └─ System time               : {total_system:8.2f} seconds')
    if wall_time is not None:
        print(f'Wall-clock time                 : {wall_time:8.2f} seconds')
        print(f'Average CPU utilisation         : {cpu_percent_avg:7.1f}%')
    print(f'Disk total                      : {disk_total_gb:8.2f} GB')
    print(f'Disk used                       : {disk_used_gb:8.2f} GB')
    print(f'Disk free                       : {disk_free_gb:8.2f} GB')
    print(f'Disk usage percentage           : {disk_percent:7.1f}%')
    print('=' * 60)
