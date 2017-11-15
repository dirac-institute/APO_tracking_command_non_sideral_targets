
import numpy as np
import scipy
import glob
import os
import re
import sys
import warnings
import pandas as pd
import astropy.io.fits as pyfits
import pyslalib.slalib as sla
import argparse
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/CFHT_observing/scripts/')
from observing_functions import *
import glob

'''
sample execution: 

ipython -i -- center_data.py -dd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/ -cd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/centered_frames/ -sf object_stars_positions



file info:

gethead FILTER *fits |grep asteroid | grep "SDSS r" | awk '{print $'1'}' | xargs -n 1 gethead DATE-OBS
gethead FILTER *fits |grep asteroid | grep "SDSS r" | awk '{print $'1'}' | xargs -n 1 gethead EXPTIME

'''

parser = argparse.ArgumentParser()
parser.add_argument("-dd", "--data_directory", help="directory for the storage of data.", nargs='*')
parser.add_argument("-cd", "--center_directory", help="directory for the storage of centered data", nargs='*')
#parser.add_argument("-af","--argument_file", type=open, action=LoadFromFile)
parser.add_argument("-sf", "--stack_file", help="text file containing file name, object x,y, ref star x,y", nargs='*')
args = parser.parse_args()

data_directory = args.data_directory[0]
center_directory = args.center_directory[0]
stack_file = args.stack_file[0]

files = np.loadtxt(stack_file,usecols=(0,),dtype='string')
image_data_frame = pd.DataFrame(files,columns=['fname'])

obj_x_y_coords_star_x_y_coords = np.loadtxt(stack_file,usecols=(1,2,3,4),dtype='string')

image_data_frame['ast_x'] = pd.Series("", index=image_data_frame.index)
image_data_frame['ast_x'] = obj_x_y_coords_star_x_y_coords[:,0]

image_data_frame['ast_y'] = pd.Series("", index=image_data_frame.index)
image_data_frame['ast_y'] = obj_x_y_coords_star_x_y_coords[:,1]

image_data_frame['star_x'] = pd.Series("", index=image_data_frame.index)
image_data_frame['star_x'] = obj_x_y_coords_star_x_y_coords[:,2]

image_data_frame['star_y'] = pd.Series("", index=image_data_frame.index)
image_data_frame['star_y'] = obj_x_y_coords_star_x_y_coords[:,3]

image_data_frame['exptime_s'] = pd.Series("", index=image_data_frame.index)
image_data_frame['filter'] = pd.Series("", index=image_data_frame.index)
image_data_frame['date_obs'] = pd.Series("", index=image_data_frame.index)


#uncomment to center images
for ff,fname in enumerate(files):
    #if ff >59 and ff < 70:
    fits_file_name = data_directory+fname
    image_data_frame['exptime_s'][ff] = pyfits.open(fits_file_name)[0].header['EXPTIME']
    image_data_frame['filter'][ff] = pyfits.open(fits_file_name)[0].header['FILTER']
    image_data_frame['date_obs'][ff] = cal_date_fits_format_to_mjd(pyfits.open(fits_file_name)[0].header['DATE-OBS'])
    #create centered_star_frames

    x_pos, y_pos = image_data_frame.ix[ff,['star_x', 'star_y']]
    centered_frame, dat_head = create_frame_for_stacking(fits_file_name, int(x_pos), int(y_pos))
    centered_name = fits_file_name.replace('.fits','_centered_stars.fits')
    pyfits.writeto(centered_name,centered_frame.astype(np.float32),overwrite=True,header=dat_head)
    os.system('mv '+ centered_name + ' ' + center_directory)
    #create centered asteroid frames
    x_pos, y_pos = image_data_frame.ix[ff,['ast_x', 'ast_y']]
    centered_frame, dat_head = create_frame_for_stacking(fits_file_name, int(x_pos), int(y_pos))
    centered_name = fits_file_name.replace('.fits','_centered_asteroid.fits')
    pyfits.writeto(centered_name,centered_frame.astype(np.float32),overwrite=True,header=dat_head)
    os.system('mv '+ centered_name + ' ' + center_directory)
