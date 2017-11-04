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
import glob

'''
sample execution: 

ipython -i -- stack_data.py -dd reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/ -od reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/stacked_frames/ -sf object_stars_positions -cd reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/centered_frames/

'''

parser = argparse.ArgumentParser()
parser.add_argument("-dd", "--data_directory", help="directory for the storage of data.", nargs='*')
parser.add_argument("-od", "--output_directory", help="directory for the storage of stacked data", nargs='*')
parser.add_argument("-cd", "--center_directory", help="directory for the storage of centered data", nargs='*')
parser.add_argument("-sf", "--stack_file", help="text file containing file name, object x,y, ref star x,y", nargs='*')
#parser.add_argument("-af","--argument_file", type=open, action=LoadFromFile)
args = parser.parse_args()

data_directory = args.data_directory[0]
center_directory = args.center_directory[0]
center_directory = args.center_directory[0]
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



#uncomment to center images
for ff,fname in enumerate(files):
    #if ff >59 and ff < 70:
    fits_file_name = data_directory+fname
    image_data_frame['exptime_s'][ff] = pyfits.open(fits_file_name)[0].header['EXPTIME']
    image_data_frame['filter'][ff] = pyfits.open(fits_file_name)[0].header['FILTER']
    image_data_frame['date_obs'][ff] = cal_date_fits_format_to_mjd(pyfits.open(fits_file_name)[0].header['DATE-OBS'])
    #create centered_star_frames
    '''
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
    '''


#i frames

fname =  np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS i']['fname'].tolist())
time_s = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS i']['exptime_s'].tolist())
dates_mjd = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS i']['date_obs'].tolist())


fits_file_name = center_directory+fname[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))

filter_name = pyfits.open(centered_name_asteroid)[0].header['FILTER'][pyfits.open(centered_name_asteroid)[0].header['FILTER'].find('SDSS ')+5:]

for i in range(0, len(fname)):
    fits_file_name = center_directory+fname[i]
    centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
    datfile = pyfits.getdata(centered_name_asteroid, header=True)
    dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
    dat_head = datfile[1]
    stack_array[:,:,i] = dat_raw

frame_interval = '_frames_'+fname[0][fname[0].find('00'):].replace('.fits','') + '_to_' + fname[-1][fname[-1].find('00'):].replace('.fits','')
#stack_array -= np.median(stack_array)
#stack_array /= fname.shape[0]
fits_file_name = output_directory + fname[i]
stacked_name_asteroid = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_asteroid.fits')
dat_head['EXPTIME'] = time_s.sum()
dat_head['DATE-OBS'] = np.mean(dates_mjd)
pyfits.writeto(stacked_name_stars,np.median(stack_array.astype(np.float32),axis=2),overwrite=True,header=dat_head)


fits_file_name = center_directory+fname[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))

#stack i frames for stars
for i in range(0, len(fname)):
    fits_file_name = center_directory+fname[i]
    centered_name_star = fits_file_name.replace('.fits','_centered_star.fits')
    datfile = pyfits.getdata(centered_name_asteroid, header=True)
    dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
    dat_head = datfile[1]
    stack_array[:,:,i] = dat_raw

fits_file_name = output_directory + fname[i]
#stack_array -= np.median(stack_array)
#stack_array /= fname.shape[0]

stacked_name_stars = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_stars.fits')
dat_head['EXPTIME'] = time_s.sum()
dat_head['DATE-OBS'] = np.mean(dates_mjd)
pyfits.writeto(stacked_name_stars,np.median(stack_array.astype(np.float32),axis=2),overwrite=True,header=dat_head)


#g frames

fname =  np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS g']['fname'].tolist())
time_s = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS g']['exptime_s'].tolist())
dates_mjd = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS g']['date_obs'].tolist())

fits_file_name = center_directory+fname[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))

filter_name = pyfits.open(centered_name_asteroid)[0].header['FILTER'][pyfits.open(centered_name_asteroid)[0].header['FILTER'].find('SDSS ')+5:]

for i in range(0, len(fname)):
    fits_file_name = center_directory+fname[i]
    centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
    datfile = pyfits.getdata(centered_name_asteroid, header=True)
    dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
    dat_head = datfile[1]
    stack_array[:,:,i] = dat_raw

#stack_array -= np.median(stack_array)
#stack_array /= fname.shape[0]
frame_interval = '_frames_'+fname[0][fname[0].find('00'):].replace('.fits','') + '_to_' + fname[-1][fname[-1].find('00'):].replace('.fits','')

fits_file_name = output_directory + fname[i]
stacked_name_asteroid = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_asteroid.fits')
dat_head['EXPTIME'] = time_s.sum()
dat_head['DATE-OBS'] = np.mean(dates_mjd)
pyfits.writeto(stacked_name_stars,np.median(stack_array.astype(np.float32),axis=2),overwrite=True,header=dat_head)


fits_file_name = center_directory+fname[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))

#stack i frames for stars
for i in range(0, len(fname)):
    fits_file_name = center_directory+fname[i]
    centered_name_star = fits_file_name.replace('.fits','_centered_star.fits')
    datfile = pyfits.getdata(centered_name_asteroid, header=True)
    dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
    dat_head = datfile[1]
    stack_array[:,:,i] = dat_raw

#stack_array -= np.median(stack_array)
#stack_array /= fname.shape[0]
fits_file_name = output_directory + fname[i]
stacked_name_stars = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_stars.fits')
dat_head['EXPTIME'] = time_s.sum()
dat_head['DATE-OBS'] = np.mean(dates_mjd)
pyfits.writeto(stacked_name_stars,np.median(stack_array.astype(np.float32),axis=2),overwrite=True,header=dat_head)



#r frames

fname_temp =  np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS r']['fname'].tolist())
time_s_temp = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS r']['exptime_s'].tolist())
dates_mjd_temp = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS r']['date_obs'].tolist())

fits_file_name = center_directory+fname_temp[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape)

filter_name = pyfits.open(centered_name_asteroid)[0].header['FILTER'][pyfits.open(centered_name_asteroid)[0].header['FILTER'].find('SDSS ')+5:]
#start, stop = 0, 15


start_stop = np.array([[51,46],[46,41],[41,39],[39,37],[37,35],[35,34],[34,32],[32,30],[30,28],[28,26],[26,23],[23,19],[19,15],[15,0]])
for qq in range(0,len(start_stop)):
    stop, start = start_stop[qq]
    fname = fname_temp[start:stop]
    time_s = time_s_temp[start:stop]
    dates_mjd = dates_mjd_temp[start:stop]
    stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))
    for i in range(0, len(fname)):
        fits_file_name = center_directory+fname[i]
        centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
        datfile = pyfits.getdata(centered_name_asteroid, header=True)
        dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
        dat_head = datfile[1]
        stack_array[:,:,i] = dat_raw

    #stack_array -= np.median(stack_array)
    #stack_array /= fname_temp[start:stop].shape[0]
    frame_interval = '_frames_'+fname[0][fname[0].find('00'):].replace('.fits','') + '_to_' + fname[-1][fname[-1].find('00'):].replace('.fits','')

    fits_file_name = output_directory + fname[i]
    stacked_name_asteroid = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_asteroid.fits')
    dat_head['EXPTIME'] = time_s.sum()
    dat_head['DATE-OBS'] = np.mean(dates_mjd)
    pyfits.writeto(stacked_name_stars,np.median(stack_array.astype(np.float32),axis=2),overwrite=True,header=dat_head)


    fits_file_name = center_directory+fname[0]#stack i frames for asteroids
    centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
    datfile = pyfits.getdata(centered_name_asteroid, header=True)
    dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
    dat_head = datfile[1]
    stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))

    #stack i frames for stars
    for i in range(0, len(fname)):
        fits_file_name = center_directory+fname[i]
        centered_name_star = fits_file_name.replace('.fits','_centered_star.fits')
        datfile = pyfits.getdata(centered_name_asteroid, header=True)
        dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
        dat_head = datfile[1]
        stack_array[:,:,i] = dat_raw

    #stack_array -= np.median(stack_array)
    #stack_array /= fname_temp[start:stop].shape[0]
    fits_file_name = output_directory + fname[i]
    stacked_name_stars = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_stars.fits')
    dat_head['EXPTIME'] = time_s.sum()
    dat_head['DATE-OBS'] = np.mean(dates_mjd)
    pyfits.writeto(stacked_name_stars,np.median(stack_array.astype(np.float32),axis=2),overwrite=True,header=dat_head)



'''
#fix date obs

files1 = glob.glob(os.path.join(output_directory, "*r_stacked_asteroid*.fits"))
files2 = glob.glob(os.path.join(output_directory, "*r_stacked_stars*.fits"))


for i in range (0, len(files1)):
    print files1[i], files2
    datfile1 = pyfits.getdata(files1[i], header=True)
    dat_raw1 = datfile1[0]#[::-1,:] #must flip data then flip back
    dat_head1 = datfile1[1]
    datfile2 = pyfits.getdata(files2[i], header=True)
    dat_raw2 = datfile2[0]#[::-1,:] #must flip data then flip back
    dat_head2 = datfile2[1]
    dat_head1['DATE-OBS'] = dat_head2['DATE-OBS']
    pyfits.writeto(files1[i],dat_raw1.astype(np.float32),overwrite=True,header=dat_head1)

'''
