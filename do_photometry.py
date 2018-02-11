#Bryce Bolin

import sys
#special numpy
sys.path.insert(1, '/usr/local/Cellar/numpy/1.14.0/lib/python2.7/site-packages')
import numpy as np
sys.path.append('/Users/bolin/anaconda3/envs/py27/lib/python2.7/site-packages')
import argparse
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/APO_observing/')
from apo_observing_functions import *
import zscale as z

#astropy photometry stuff
from astropy.io import fits
import astropy.io.fits as pyfits
from astropy.time import Time
from astropy.utils.console import ProgressBar
from astropy.modeling import models, fitting
#from photutils.morphology import centroid_com
from photutils import centroid_com, centroid_1dg, centroid_2dg
from photutils import CircularAperture, CircularAnnulus, aperture_photometry





'''
reads in a list of x,y positions to access a list of reduced ccd image fits files, cuts out a square of size specified by the user around the x, y coords, perform centroiding and then aperture photometry using median background subtraction, then spits out the instrumental magnitudes with the uncertainties

example execution:

ipython -i -- do_photometry.py -dd /Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q1UW07/UT180120/ARCTIC_2018_01_20_UTC/reduced/data/ -fl /Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q1UW07/UT180120/ARCTIC_2018_01_20_UTC/reduced/data/2018_AV2_filenames -cl /Users/bolin/NEO/Follow_up/APO_observing/rawdata/Q1UW07/UT180120/ARCTIC_2018_01_20_UTC/reduced/data/2018_AV2_X_Y_positions -apc 9 15 20 -shl 10 -cm com -ofn 2018_AV2_APO_2018_01_20

'''

parser = argparse.ArgumentParser()
parser.add_argument("-dd", "--data_directory", help="directory for the storage of data.", nargs='*')
parser.add_argument("-fl", "--file_list", help="list of files to reduce", nargs='*')
parser.add_argument("-cl", "--center_list", help="list of target positions in x,y coordiantes, two columns", nargs='*')
parser.add_argument("-apc", "--aperture_components", help="aperture components in units of pixes: aperture_radius sky_inner_radius sky_outer_radius, e.g., 6 15 20", nargs='*')
parser.add_argument("-shl", "--square_half_length", help="half the length of the square centered on the center position of the psf. This is used for the centroid step. Must be large enough to accomodate the sky rings in the perture photometry step.", nargs='*')
parser.add_argument("-cm", "--centroid_mode", help="set to 1d, 2d gaussian fitting or center of mass for centroiding, e.g., 1d, 2d, com", nargs='*')
parser.add_argument("-ofn", "--output_filename", help="name of the pyc filename for output of the photometry measurements", nargs='*')
args = parser.parse_args()

data_directory = args.data_directory[0]
file_list = args.file_list[0]
center_file = args.center_list[0]
aperture_radius_pixels, inner_sky_ring_pixels, outer_sky_ring_pixels = string_seperated_to_array_spaces(args.aperture_components,'float')
square_half_length_pixels =  int(args.square_half_length[0])
output_filename = args.output_filename[0]
centroid_mode = args.centroid_mode[0]

files = np.loadtxt(file_list,dtype='string')
center_list_X_Y = np.loadtxt(center_file,dtype='int')
center_X, Center_Y = center_list_X_Y[:,0], center_list_X_Y[:,1]

#to properly orient data, rotate by 270 degrees and flip vertical
#np.rot90(datfile[0][::-1,:],k=3), then flip x and y coords from list.
true_x_coords, true_y_coords = Center_Y, center_X


#output array all strings to contain filename mjd exptime filter mag mag_unc airmass
number_of_entries_per_row = 7
max_char_size = 20
output_array_filename_mjd_exptime_filter_mag_mag_unc_airmass = np.chararray((len(files),number_of_entries_per_row),itemsize=max_char_size)
output_array[:] = "aaaaaaaaaaaaaaaaaaaa"

print "#filename_0 mjd_1 exptime_s_2 filter_3 mag_4 mag_unc_5 airmass_6
#for i in range(0,len(files)):
for i in range(14,15):
    fits_file_name = data_directory+files[i]
    exp_time_s = pyfits.open(fits_file_name)[0].header['EXPTIME']
    filter = pyfits.open(fits_file_name)[0].header['FILTER']
    airmass = pyfits.open(fits_file_name)[0].header['AIRMASS']
    date_mjd = cal_date_fits_format_to_mjd(pyfits.open(fits_file_name)[0].header['DATE-OBS'])
    datfile = pyfits.getdata(fits_file_name, header=True)
    #to properly orient data, rotate by 270 degrees and flip vertical
    #np.rot90(datfile[0][::-1,:],k=3), then flip x and y coords from list.
    dat_raw = np.rot90(datfile[0][::-1,:],k=3)
    dat_head = datfile[1]

    #centroid stamp
    image_stamp = dat_raw[true_y_coords[i]-square_half_length_pixels:true_y_coords[i]+square_half_length_pixels,true_x_coords[i]-square_half_length_pixels:true_x_coords[i]+square_half_length_pixels]

    if centroid_mode == "1d":
        x_stamp_centroid, y_stamp_centroid = centroid_1dg(image_stamp)
    if centroid_mode == "2d":
        x_stamp_centroid, y_stamp_centroid = centroid_2dg(image_stamp)
    if centroid_mode == "com":
        x_stamp_centroid, y_stamp_centroid = centroid_com(image_stamp)

    #x_centroid, y_centroid true location
    #x_centroid, y_centroid, swap the x and y
    x_centroid = int(y_stamp_centroid + true_x_coords[i] - square_half_length_pixels)
    y_centroid = int(x_stamp_centroid + true_y_coords[i] - square_half_length_pixels)

    target_apertures = CircularAperture((x_centroid, y_centroid), aperture_radius_pixels)
    background_annuli = CircularAnnulus((x_centroid, y_centroid), r_in=inner_sky_ring_pixels, r_out= outer_sky_ring_pixels)
    flux_in_annuli = aperture_photometry(dat_raw,
                                         background_annuli)['aperture_sum'].data
    background = flux_in_annuli/background_annuli.area()
    flux = aperture_photometry(dat_raw,
                               target_apertures)['aperture_sum'].data
    background_subtracted_flux = (flux - background *
                                  target_apertures.area())
    error_flux = np.sqrt(flux)

    snr  = background_subtracted_flux / error_flux

    mag = magnitude_calc_no_background(background_subtracted_flux, exp_time_s)

    mag_unc = 1./snr


    print files[i], str(date_mjd), str(exp_time), filter, str(mag), str(mag_unc), str(airmass)
    output_array_filename_mjd_exptime_filter_mag_mag_unc_airmass[i] = files[i], str(date_mjd), str(exp_time), filter, str(mag), str(mag_unc), str(airmass)


output_array_filename_mjd_exptime_filter_mag_mag_unc_airmass.dump(output_filename + ".pyc")

'''
#view of entire image
vmin,vmax = z.zscale(dat_raw)
plt.ion()
plt.figure()
plt.imshow(dat_raw,vmin=vmin,vmax=vmax,cmap='gray')
plt.scatter(x=[x_centroid], y=[y_centroid], c='r', s=40)

#stamp view of centroid of source
i=0
vmin,vmax = z.zscale(dat_raw)
plt.ion()
plt.figure()
plt.imshow(dat_raw[y_centroid-square_half_length_pixels:y_centroid+square_half_length_pixels,x_centroid-square_half_length_pixels:x_centroid+square_half_length_pixels],vmin=vmin,vmax=vmax,cmap='gray')
plt.scatter(x=[x_stamp_centroid], y=[y_stamp_centroid], c='r', s=40)

'''

#print out filename mjd exptime filter mag mag_unc

