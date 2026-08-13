'''Main script to run the LVIS processing pipeline.'''

import sys
import time
import pathlib
import logging
import argparse

# Make src importable
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from src import lvis, tile, utils, compare

def setup_logging(debug = False, info = False):
    '''Configure logging with appropriate level and format.'''

    level = logging.DEBUG if debug else logging.INFO if info else logging.WARNING
    logging.basicConfig(
        level = level,
        format = '%(asctime)s - %(levelname)s - %(message)s',
        datefmt = '%H:%M:%S',
    )

def parse_arguments():
    '''Argument parser'''

    parser = argparse.ArgumentParser(prog = 'run.py',
                        description = 'Run the LVIS processing pipeline',)
    parser.add_argument('-db', '--debug', action = 'store_true', 
                        help = 'Enable debug mode')
    parser.add_argument('-tw', '--testwave', action = 'store_true',
                        help = 'just plots a basis waveform then exits')
    parser.add_argument('-d', '--directory', default = '/geos/netdata/oosa/assignment/lvis/2015',
                        help = 'Directory containing LVIS files to process', metavar = 'NAME')
    parser.add_argument('-f', '--filename', metavar = 'NAME', 
                        help = 'Select specific file within directory to process (default all files)')
    parser.add_argument('-r', '--resolution', type = float, default = 100, metavar = 'float',
                        help = 'Set the output DEM resolution in meters (default: 100)')
    parser.add_argument('-gf', '--gapfill', type = int, default = 0, metavar = 'DIST',
                        help = 'Distance for gap filling in DEM in meters (default: 0)')
    parser.add_argument('-dn', '--demname', metavar = 'NAME',
                        help = 'Set the output DEM filename (no prefix or suffix)')
    parser.add_argument('-pn', '--plotname', metavar = 'NAME',
                        help = 'Set the output plot filename (no prefix or suffix)')
    parser.add_argument('-t', '--tilesize', type = float, default = 0.1, metavar = 'float',
                        help = 'Set the tile size in degrees for processing (default: 0.1)')
    parser.add_argument('-np', '--noplotdem', action = 'store_true',
                        help = 'Does not Plot the output DEM after processing (default: False)')
    parser.add_argument('-pg', '--pineglacier', action = 'store_true',
                        help = 'Constrains the data to the rough bounds of pine island glacier (default: False)')
    parser.add_argument('-sp', '--simpleplot', action = 'store_true',
                        help = 'Removes the coastlines and equal projection of the dem map (default: False)')
    
    args = parser.parse_args()

    if args.resolution <= 0:
        parser.error('--resolution must be a positive number')
    if args.tilesize <= 0:
        parser.error('--tilesize must be a positive number')

    return args

def main():

    args = parse_arguments()
    
    setup_logging(args.debug)

    # Make sure dems and plots is there
    pathlib.Path('dems').mkdir(exist_ok = True)
    pathlib.Path('plots').mkdir(exist_ok = True)
    demname = f'dems/{args.demname or "test"}'
    plotname = f'plots/{args.plotname or "test"}'

    folderloc = pathlib.Path(args.directory).resolve()
    if not folderloc.is_dir():
        raise ValueError(f'{args.directory} is not a valid directory path')
    
    max_bounds = [258, -80, 270, -75] if args.pineglacier else None
    
    logging.debug(f'Input directory: {str(folderloc)}')
    logging.debug(f'Output DEM name: {demname}, Output plot name: {plotname}')
    logging.debug(f'Resolution: {args.resolution} m | Gapfill: {args.gapfill} | Tilesize: {args.tilesize}°')
    
    if args.filename:
        
        fileloc =  folderloc / args.filename
        if not fileloc.is_file():
            raise ValueError(f'{fileloc} is not a valid file path')
        
        one_file = lvis.LVIS(str(fileloc))
        if args.testwave:
            one_file.plotWave(0, outName = 'plots/testwave')
            return
        
        one_file.process()
        one_file.writeDEM(args.resolution, demname, interp_dist = args.gapfill)
        if not args.noplotdem:
            lvis.Post(demname).plot(outName = plotname, mapda = not args.simpleplot)

    else:

        all_files = tile.Tile(str(folderloc), tilesize = args.tilesize, max_bounds = max_bounds) # for testwave to get filelist
        if args.testwave:
            clvis = lvis.LVIS(all_files.fileList[0])
            clvis.plotWave(0, outName = 'plots/testwave')
            return
        
        all_files.process(args.resolution)
        all_files.writeDEM(demname, interp_dist = args.gapfill)
        if not args.noplotdem:
            lvis.Post(demname).plot(outName = plotname, mapda = not args.simpleplot)

    logging.info('Processing complete')

if __name__ == '__main__':

    start_time = time.time()
    main()
    utils.resource_usage(start_time)
