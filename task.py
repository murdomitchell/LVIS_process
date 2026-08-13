'''Main script to run the assignment tasks. 
del and gc.collect() are to lower the memory artifically below 2GB'''

import sys
import time
import pathlib
import argparse
import gc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from src import lvis, tile, utils, compare

# setup
loc = '/geos/netdata/oosa/assignment/lvis/2015'
boxbig = [258, -82, 275, -73]
boxsmall = [258, -77, 265, -75]

def parse_arguments():
    '''Argument parser'''

    parser = argparse.ArgumentParser(prog = 'task.py',
                        description = 'Do the assignment tasks',)
    parser.add_argument('-t', '--task', type = int, default = 0, metavar = 'int',
                        help = 'Set the task number to execute (default: 0, which runs all tasks)')
    parser.add_argument('-p', '--plots', action = 'store_true',
                        help = 'Just do the plots for the tasks (default: False)')
    args = parser.parse_args()

    if args.task not in [0, 1, 2, 3, 4, 5]:
        parser.error('task must be a number between 0 and 5')

    return args

def task1():
    print('Running Task 1')

    all_files = tile.Tile(loc)
    current_lvis = lvis.LVIS(all_files.fileList[0])
    current_lvis.plotWave(0, outName = 'plots/testwave')

    del all_files, current_lvis
    gc.collect()

def task2():
    print('Running Task 2')

    all_files = tile.Tile(loc)
    one_file = lvis.LVIS(all_files.fileList[0])
    one_file.process()
    one_file.writeDEM(30, filename = 'dems/one', interp_dist = 200)
    lvis.Post('dems/one').plot(outName = 'plots/one', title = all_files.fileList[0][-25:-3])

    del all_files, one_file
    gc.collect()

def task3():
    print('Running Task 3')

    all_files = tile.Tile(loc, 
        tilesize = .2, max_bounds = boxbig)
    all_files.process(100)
    all_files.writeDEM(filename = 'dems/all.tif')
    lvis.Post('dems/all').plot(outName = 'plots/all', title = str(boxbig))

    del all_files
    gc.collect()

def task4():
    print('Running Task 4')

    for i in ['dems/big09', 'dems/big15']:
        cloc = '/geos/netdata/oosa/assignment/lvis/20' + i[-2:]
        all_files = tile.Tile(cloc, tilesize = .2, max_bounds = boxsmall)
        all_files.process(100)
        all_files.writeDEM(filename = i, interp_dist = 1000)
        lvis.Post(i).plot(outName = i.replace('dems', 'plots'), title = '')

        del all_files
        gc.collect()

def task5():
    print('Running Task 5')

    if not pathlib.Path('dems/big09.tif').is_file() or not pathlib.Path('dems/big15.tif').is_file():
        print('Doing all the processing for the big files')
        task4()

    ice = compare.Ice('dems/big09', 'dems/big15')
    print(ice.volume_array().attrs['total_change'])
    ice.plot_diff(outName = 'plots/diff')
    ice.plot_volume(outName = 'plots/volchng')

    del ice
    gc.collect()

def plots():

    # comparing interpolation methods
    if not pathlib.Path('dems/raw09.tif').is_file():
        all_files = tile.Tile(loc, tilesize = .2, max_bounds = boxsmall)
        all_files.process(100)
        all_files.writeDEM(filename = 'dems/raw09.tif')
        del all_files
        gc.collect()

    linear = lvis.Post('dems/raw09')
    linear.interpolate(pixel_dist = 50)
    linear.plot(outName = 'plots/linear', title = 'linear')
    methods = ['nearest', 'cubic', 'cKDTree', 'fillnodata']
    for i in methods:
        dem = lvis.Post('dems/raw09')
        dem.interpolate(pixel_dist = 50, method = i)
        dem.plot(outName = f'plots/{i}', title = i)
        compare.Ice(dem.da, linear.da).plot_diff(outName = f'plots/{i}_d', title = i)
        
        del dem
    del linear
    gc.collect()

    one = lvis.Post('dems/one')
    one.interpolate(50, 'fillnodata')
    one.plot(outName = 'plots/onefill', title = 'fillnodata')
    del one
    gc.collect()

    # task 1
    all_files = tile.Tile(loc)
    current_lvis = lvis.LVIS(all_files.fileList[0])
    current_lvis.plotWave(0, outName = 'plots/testwave')
    del current_lvis
    gc.collect()

    # task 2
    lvis.Post('dems/one').plot(outName = 'plots/one', title = all_files.fileList[0][-25:-3])
    gc.collect()

    # task 3
    lvis.Post('dems/all').plot(outName = 'plots/all', title = str(boxbig))
    gc.collect()

    # task 4
    for i in ['dems/big09', 'dems/big15']:
        lvis.Post(i).plot(outName = i.replace('dems', 'plots'), title = '')
        gc.collect()

    # task 5
    ice = compare.Ice('dems/big09', 'dems/big15')
    print(ice.volume_array().attrs['total_change'])
    ice.plot_diff(outName = 'plots/diff', title = '')
    ice.plot_volume(outName = 'plots/volchng')
    gc.collect()

def main():

    args = parse_arguments()

    if args.plots:
        plots()
        return
    
    match args.task:
        case 0:
            task1()
            task2()
            task3()
            task4()
            task5()
        case 1:
            task1()
        case 2:
            task2()
        case 3:
            task3()
        case 4:
            task4()
        case 5:
            task5()

if __name__ == '__main__':

    start_time = time.time()
    main()
    utils.resource_usage(start_time)
