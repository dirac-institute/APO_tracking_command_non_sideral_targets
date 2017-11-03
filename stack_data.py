#data for stacking APO data from the night of 2017-10-29

import numpy as np
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

'''
sample execution: 

ipython -i -- stack_data.py -dd reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/ -od stack_output/AK17U010/2017_10_29/ -sf object_stars_positions

'''

parser = argparse.ArgumentParser()
parser.add_argument("-dd", "--data_directory", help="directory for the storage of data.", nargs='*')
parser.add_argument("-od", "--output_directory", help="directory for the storage of stacked data", nargs='*')
parser.add_argument("-sf", "--stack_file", help="text file containing file name, object x,y, ref star x,y", nargs='*')
#parser.add_argument("-af","--argument_file", type=open, action=LoadFromFile)
args = parser.parse_args()

data_directory = args.data_directory[0]
output_directory = args.output_directory[0]
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



for ff,fname in enumerate(files):
    image_data_frame['exptime_s'][ff] = pyfits.open(data_directory+fname)[0].header['EXPTIME']
    image_data_frame['filter'][ff] = pyfits.open(data_directory+fname)[0].header['FILTER']
    image_data_frame['date_obs'][ff] = cal_date_fits_format_to_mjd(pyfits.open(data_directory+fname)[0].header['DATE-OBS'])


#image_data_frame[image_data_frame['filter']=='SDSS i']

fits_file, x_pos, y_pos = 'reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/reduced_2017U1.0072.fits', 990, 1214

centered_frame, dat_head = create_frame_for_stacking(fits_file, x_pos, y_pos)

centered_name = fits_file.replace('.fits','_centered.fits')

pyfits.writeto(centered_name,centered_frame.astype(np.float32),overwrite=True,header=dat_head)

'''
datfile = pyfits.getdata(f, header=True)
dat_raw = datfile[0]
dat_head = datfile[1]

amp = pyfits.open(f)[0].header['READAMPS']

biases = np.array([trim_image(df['fname'][n])[0] for n in bias_idx])

datfile = trim_image(df['fname'][n])
dat_raw = datfile[0]
dat_head = datfile[1]
time = df['exp'][n]

pyfits.writeto(name,dat.astype(np.float32),overwrite=True,header=dat_head)
'''