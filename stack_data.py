#data for stacking APO data from the night of 2017-10-29

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
import scipy.ndimage

'''
sample execution: 

#i
ipython -i -- stack_data.py -dd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/ -od reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/stacked_frames/ -sf object_stars_positions -cd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/centered_frames/ -sf object_stars_positions -m i -nt 4 -op median_average -cut 0.16 -stst 8 0

#g
ipython -i -- stack_data.py -dd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/ -od reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/stacked_frames/ -sf object_stars_positions -cd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/centered_frames/ -sf object_stars_positions -m g -nt 4 -op median_average -cut 0.16 -stst 26 7

#r

ipython -i -- stack_data.py -dd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/ -od reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/stacked_frames/ -sf object_stars_positions -cd /Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/data/centered_frames/ -sf object_stars_positions -m r -nt 4 -op median_average -cut 0.16 -stst 45 31


file info:

gethead FILTER *fits |grep asteroid | grep "SDSS r" | awk '{print $'1'}' | xargs -n 1 gethead DATE-OBS
gethead FILTER *fits |grep asteroid | grep "SDSS r" | awk '{print $'1'}' | xargs -n 1 gethead EXPTIME

'''

parser = argparse.ArgumentParser()
parser.add_argument("-dd", "--data_directory", help="directory for the storage of data.", nargs='*')
parser.add_argument("-od", "--output_directory", help="directory for the storage of stacked data", nargs='*')
parser.add_argument("-cd", "--center_directory", help="directory for the storage of centered data", nargs='*')
parser.add_argument("-sf", "--stack_file", help="text file containing file name, object x,y, ref star x,y", nargs='*')
parser.add_argument("-m", "--mode", help="g, r or i.", nargs='*')
parser.add_argument("-nt", "--number_mean_median_stack_trials", help="g, r or i.", nargs='*')
parser.add_argument("-op", "--operation", help="mean or median",nargs='*')
parser.add_argument("-cut", "--cut", help="robust average cut",nargs='*')
parser.add_argument("-stst", "--start_stop", help="start and stopping point of indices, stop goes first, eg 8 0",nargs='*')

#parser.add_argument("-af","--argument_file", type=open, action=LoadFromFile)
args = parser.parse_args()

data_directory = args.data_directory[0]
center_directory = args.center_directory[0]
center_directory = args.center_directory[0]
output_directory = args.output_directory[0]
stack_file = args.stack_file[0]
mode = args.mode[0]
number_mean_median_stack_trials = int(args.number_mean_median_stack_trials[0])
operation = args.operation[0]
cut = float(args.cut[0])
stop, start = int(args.start_stop[0]), int(args.start_stop[1])

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




#i mode
#cut = 0.16 is good
#display reduced_2017U1.0034_frames_0005_to_0034_filter_i_stacked_asteroid_mean.fits 1 zr- zs- z1=2490 z2=2580
#test
#display reduced_2017U1.0034_frames_0005_to_0034_filter_i_stacked_asteroid_median_average_4.fits 1 zr- zs- z1=2490 z2=2580

#g frames
#Yan's is best version  display
#display stack-g-15frames-45min-cl2.fits 2 zr- zs- z1=0.01 z2=0.
#my g frames
#display reduced_2017U1.0037_frames_0019_to_0037_filter_g_stacked_asteroid_median.fits 1 zr- zs- z1=2400 z2=2550

#test
#display reduced_2017U1.0037_frames_0019_to_0037_filter_g_stacked_asteroid_median_average_4.fits 1 zr- zs- z1=2400 z2=2550

#r frames
#cut = 0.4 works
#display reduced_2017U1.0067_frames_0054_to_0067_filter_r_stacked_asteroid_mean.fits 1 zr- zs- z1=1700 z2=1800

#test
#display reduced_2017U1.0067_frames_0054_to_0067_filter_r_stacked_asteroid_mean.fits 1 zr- zs- z1=1700 z2=1800


fname_temp =  np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS ' + mode]['fname'].tolist())
time_s_temp = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS ' + mode]['exptime_s'].tolist())
dates_mjd_temp = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS ' + mode]['date_obs'].tolist())

fits_file_name = center_directory+fname_temp[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape)

filter_name = pyfits.open(centered_name_asteroid)[0].header['FILTER'][pyfits.open(centered_name_asteroid)[0].header['FILTER'].find('SDSS ')+5:]
#start, stop = 0, 15

fname_temp =  np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS ' + mode]['fname'].tolist())
time_s_temp = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS ' + mode]['exptime_s'].tolist())
dates_mjd_temp = np.asarray(image_data_frame.ix[image_data_frame['filter']=='SDSS ' + mode]['date_obs'].tolist())

fits_file_name = center_directory+fname_temp[0]#stack i frames for asteroids
centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
datfile = pyfits.getdata(centered_name_asteroid, header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head = datfile[1]
stack_array = np.zeros(dat_raw.shape)

filter_name = pyfits.open(centered_name_asteroid)[0].header['FILTER'][pyfits.open(centered_name_asteroid)[0].header['FILTER'].find('SDSS ' + mode)+5:]

for qq in range(0,1):
#for qq in range(0,1):
    fname = fname_temp[start:stop]
    time_s = time_s_temp[start:stop]
    dates_mjd = dates_mjd_temp[start:stop]
#experimental rotation
    number_trials = number_mean_median_stack_trials
    stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1] *number_trials*2).reshape(dat_raw.shape[0], dat_raw.shape[1], number_trials,2)


for mm in range(0,number_trials):
    print ("trial stage",mm)
    per_trial_stack_array = np.zeros(dat_raw.shape[0] * dat_raw.shape[1]* len(fname)).reshape(dat_raw.shape[0], dat_raw.shape[1], len(fname))
    for i in range(0, len(fname)):
        print ("file number",i)
        fits_file_name = center_directory+fname[i]
        centered_name_asteroid = fits_file_name.replace('.fits','_centered_asteroid.fits')
        datfile = pyfits.getdata(centered_name_asteroid, header=True)
        dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
        dat_head = datfile[1]
        dat_raw[np.where(dat_raw == np.inf)] = 0.0
        per_trial_stack_array[:,:,i] = scipy.ndimage.interpolation.rotate(dat_raw,np.random.randint(0,359), reshape = False)
        #stack_array[:,:,i+(len(fname)*mm)] = np.rot90(dat_raw,np.random.randint(1,5))
        #stack_array[:,:,i+(len(fname)*mm)] = dat_raw
    if operation == "mean":
        stack_array[:,:,mm,0] = scipy.stats.trim_mean(per_trial_stack_array[::-1,:].astype(np.float32),cut,axis=2)
    if operation == "median":
        stack_array[:,:,mm,0] = np.median(per_trial_stack_array[::-1,:].astype(np.float32),axis=2)
    if operation == "both":
        stack_array[:,:,mm,0] = scipy.stats.trim_mean(per_trial_stack_array[::-1,:].astype(np.float32),cut,axis=2)
        stack_array[:,:,mm,1] = np.median(per_trial_stack_array[::-1,:].astype(np.float32),axis=2)
    if operation == "median_average":
        stack_array[:,:,mm,0] = scipy.stats.trim_mean(per_trial_stack_array[::-1,:].astype(np.float32),cut,axis=2)


frame_interval = '_frames_'+fname[0][fname[0].find('00'):].replace('.fits','') + '_to_' + fname[-1][fname[-1].find('00'):].replace('.fits','')
#stack_array -= scipy.stats.trim_mean(stack_array)
#stack_array /= fname.shape[0]
fits_file_name = output_directory + fname[i]
stacked_name_asteroid = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_asteroid.fits')
dat_head['EXPTIME'] = time_s.sum()
dat_head['DATE-OBS'] = np.mean(dates_mjd)
if operation == "mean":
    write_stack_array = np.mean(stack_array[:,:,:,0],axis=2)
    pyfits.writeto(stacked_name_asteroid.replace('.fits','_'+operation+'_' +args.number_mean_median_stack_trials[0]+ '.fits'),write_stack_array,overwrite=True,header=dat_head)
if operation == "median":
    write_stack_array = np.median(stack_array[:,:,:,1],axis=2)
    pyfits.writeto(stacked_name_asteroid.replace('.fits','_'+operation+'_' +args.number_mean_median_stack_trials[0]+ '.fits'),write_stack_array,overwrite=True,header=dat_head)
if operation == "both":
    write_stack_array = np.mean(stack_array[:,:,:,0],axis=2)
    pyfits.writeto(stacked_name_asteroid.replace('.fits','_mean_' +args.number_mean_median_stack_trials[0]+ '.fits'),write_stack_array,overwrite=True,header=dat_head)
    write_stack_array = np.median(stack_array[:,:,:,1],axis=2)
    pyfits.writeto(stacked_name_asteroid.replace('.fits','_median_' +args.number_mean_median_stack_trials[0]+ '.fits'),write_stack_array,overwrite=True,header=dat_head)
if operation == "median_average":
    write_stack_array = np.median(stack_array[:,:,:,0],axis=2)
    pyfits.writeto(stacked_name_asteroid.replace('.fits','_median_average_' +args.number_mean_median_stack_trials[0]+ '.fits'),write_stack_array,overwrite=True,header=dat_head)


'''
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
#stack_array -= scipy.stats.trim_mean(stack_array)
#stack_array /= fname.shape[0]

stacked_name_stars = fits_file_name.replace('.fits',frame_interval+'_filter_' + filter_name + '_stacked_stars.fits')
dat_head['EXPTIME'] = time_s.sum()
dat_head['DATE-OBS'] = np.mean(dates_mjd)
pyfits.writeto(stacked_name_stars,scipy.stats.trim_mean(stack_array[::-1,:].astype(np.float32),cut,axis=2),overwrite=True,header=dat_head)
'''


'''
g = 23.51
g_unc = 0.22
r = 23.1
r_unc = 0.09
i = 22.87
i_unc = 0.23
a_slope(g, g_unc, r, r_unc, i, i_unc)


Fix yan's flats:

#r
datfile = pyfits.getdata('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/skyflat_r.0099.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head1 = datfile[1]

datfile = pyfits.getdata('/Users/bolin/Dropbox/Interstellar/yans-flats/illum-2017U1-r-frm51to72.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head2 = dat_head1

# ll, ul, lr, ur
quads = ['DSEC11', 'DSEC21', 'DSEC12', 'DSEC22']

dat = [[],[],[],[]]
for i,quad in enumerate(quads):
    idx_string = pyfits.open('/Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q4DD04/UT171029/skyflat_r.0099.fits')[0].header[quad]
    idx = re.split('[: ,]',idx_string.rstrip(']').lstrip('['))
    dat[i] = dat_raw[int(idx[2])-1:int(idx[3]),int(idx[0])-1:int(idx[1])]

sci_lo = np.concatenate((dat[2], dat[3]), axis = 1)
sci_up = np.concatenate((dat[0], dat[1]), axis = 1)
sci = np.concatenate((sci_up, sci_lo), axis = 0)

pyfits.writeto('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/cals/2master_flat_r.fits',sci.astype(np.float32),overwrite=True,header=dat_head1)

datfile = pyfits.getdata('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/skyflat_r.0099.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head1 = datfile[1]

datfile = pyfits.getdata('/Users/bolin/Dropbox/Interstellar/yans-flats/flat-r.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head2 = dat_head1

# ll, ul, lr, ur
quads = ['DSEC11', 'DSEC21', 'DSEC12', 'DSEC22']

dat = [[],[],[],[]]
for i,quad in enumerate(quads):
    idx_string = pyfits.open('/Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q4DD04/UT171029/skyflat_r.0099.fits')[0].header[quad]
    idx = re.split('[: ,]',idx_string.rstrip(']').lstrip('['))
    dat[i] = dat_raw[int(idx[2])-1:int(idx[3]),int(idx[0])-1:int(idx[1])]

sci_lo = np.concatenate((dat[2], dat[3]), axis = 1)
sci_up = np.concatenate((dat[0], dat[1]), axis = 1)
sci = np.concatenate((sci_up, sci_lo), axis = 0)

pyfits.writeto('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/cals/master_flat_r.fits',sci.astype(np.float32),overwrite=True,header=dat_head1)

#i
datfile = pyfits.getdata('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/skyflat_i.0080.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head1 = datfile[1]

datfile = pyfits.getdata('/Users/bolin/Dropbox/Interstellar/yans-flats/illum-2017U1-i-allframes.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head2 = dat_head1

# ll, ul, lr, ur
quads = ['DSEC11', 'DSEC21', 'DSEC12', 'DSEC22']

dat = [[],[],[],[]]
for i,quad in enumerate(quads):
    idx_string = pyfits.open('/Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q4DD04/UT171029/skyflat_r.0099.fits')[0].header[quad]
    idx = re.split('[: ,]',idx_string.rstrip(']').lstrip('['))
    dat[i] = dat_raw[int(idx[2])-1:int(idx[3]),int(idx[0])-1:int(idx[1])]

sci_lo = np.concatenate((dat[2], dat[3]), axis = 1)
sci_up = np.concatenate((dat[0], dat[1]), axis = 1)
sci = np.concatenate((sci_up, sci_lo), axis = 0)

pyfits.writeto('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/cals/2master_flat_i.fits',sci.astype(np.float32),overwrite=True,header=dat_head1)

datfile = pyfits.getdata('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/skyflat_i.0080.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head1 = datfile[1]

datfile = pyfits.getdata('/Users/bolin/Dropbox/Interstellar/yans-flats/flat-i.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head2 = dat_head1

# ll, ul, lr, ur
quads = ['DSEC11', 'DSEC21', 'DSEC12', 'DSEC22']

dat = [[],[],[],[]]
for i,quad in enumerate(quads):
    idx_string = pyfits.open('/Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q4DD04/UT171029/skyflat_r.0099.fits')[0].header[quad]
    idx = re.split('[: ,]',idx_string.rstrip(']').lstrip('['))
    dat[i] = dat_raw[int(idx[2])-1:int(idx[3]),int(idx[0])-1:int(idx[1])]

sci_lo = np.concatenate((dat[2], dat[3]), axis = 1)
sci_up = np.concatenate((dat[0], dat[1]), axis = 1)
sci = np.concatenate((sci_up, sci_lo), axis = 0)

pyfits.writeto('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/cals/master_flat_i.fits',sci.astype(np.float32),overwrite=True,header=dat_head1)

#g
datfile = pyfits.getdata('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/skyflat_g.0089.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head1 = datfile[1]

datfile = pyfits.getdata('/Users/bolin/Dropbox/Interstellar/yans-flats/illum-2017U1-g-allframes.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head2 = dat_head1

# ll, ul, lr, ur
quads = ['DSEC11', 'DSEC21', 'DSEC12', 'DSEC22']

dat = [[],[],[],[]]
for i,quad in enumerate(quads):
    idx_string = pyfits.open('/Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q4DD04/UT171029/skyflat_r.0099.fits')[0].header[quad]
    idx = re.split('[: ,]',idx_string.rstrip(']').lstrip('['))
    dat[i] = dat_raw[int(idx[2])-1:int(idx[3]),int(idx[0])-1:int(idx[1])]

sci_lo = np.concatenate((dat[2], dat[3]), axis = 1)
sci_up = np.concatenate((dat[0], dat[1]), axis = 1)
sci = np.concatenate((sci_up, sci_lo), axis = 0)

pyfits.writeto('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/cals/2master_flat_g.fits',sci.astype(np.float32),overwrite=True,header=dat_head1)

datfile = pyfits.getdata('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/skyflat_g.0089.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head1 = datfile[1]

datfile = pyfits.getdata('/Users/bolin/Dropbox/Interstellar/yans-flats/flat-g.fits', header=True)
dat_raw = datfile[0]#[::-1,:] #must flip data then flip back
dat_head2 = dat_head1

# ll, ul, lr, ur
quads = ['DSEC11', 'DSEC21', 'DSEC12', 'DSEC22']

dat = [[],[],[],[]]
for i,quad in enumerate(quads):
    idx_string = pyfits.open('/Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q4DD04/UT171029/skyflat_r.0099.fits')[0].header[quad]
    idx = re.split('[: ,]',idx_string.rstrip(']').lstrip('['))
    dat[i] = dat_raw[int(idx[2])-1:int(idx[3]),int(idx[0])-1:int(idx[1])]

sci_lo = np.concatenate((dat[2], dat[3]), axis = 1)
sci_up = np.concatenate((dat[0], dat[1]), axis = 1)
sci = np.concatenate((sci_up, sci_lo), axis = 0)

pyfits.writeto('/Users/bolin/NEO/Follow_up/APO_observing/reduced_data/AK17U010/2017_10_29/rawdata/reduced/cals/master_flat_g.fits',sci.astype(np.float32),overwrite=True,header=dat_head1)

'''
