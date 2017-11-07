#Bryce Bolin

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import string
import random
import matplotlib.ticker
from threading import Thread
import matplotlib.ticker
import scipy.ndimage
from time import clock, time
import sys
import io
import os
import argparse
#import functions from /Users/bolin/Families/scripts/family_functions.py
sys.path.insert(0, '/Users/bolin/Families/scripts/')
from family_functions import *
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/CFHT_observing/scripts/')
from observing_functions import *
import warnings
import glob
from matplotlib.ticker import FuncFormatter
import argparse
import scipy.signal as signal
import pyslalib.slalib as sla

plt.ion()

lightcurve_mjd_mag_mag_unc = np.loadtxt('A_2017_U1_APO_2017_10_29_lightcurve')



time_mjd = lightcurve_mjd_mag_mag_unc[:,0]
apo_mag = lightcurve_mjd_mag_mag_unc[:,1]
apo_mag_unc = lightcurve_mjd_mag_mag_unc[:,2]

#DCT lightcurve
lightcurve_hours_mag = np.loadtxt('A_2017_U1_DCT_2017_10_30_lightcurve')

dct_date_MJD = cal_date_to_mjd(2017,10,30)
DCT_2017_10_30_hours = lightcurve_hours_mag[:,0]
dct_fulldate_mjd = dct_date_MJD + (DCT_2017_10_30_hours/24.)
dct_fulldate_seconds = (dct_fulldate_mjd - time_mjd[0]) * 3600*24
dct_mag = lightcurve_hours_mag[:,1]
dct_mag_unc = lightcurve_hours_mag[:,2]

print ("Date (MJD) Magnitude r' Magnitude uncertainty")
for i in range(0,len(dct_fulldate_mjd)):
    print (dct_fulldate_mjd[i], dct_mag[i], dct_mag_unc[i])

mag = np.append(apo_mag+0.2,dct_mag)#normalized APO mags to 10/30
mag_unc = np.append(apo_mag_unc,dct_mag_unc)
full_times_mjd = np.append(time_mjd, dct_fulldate_mjd)

print ("Date (MJD) Magnitude r' Magnitude uncertainty")
for i in range(0,len(mag)):
    print (full_times_mjd[i], mag[i], mag_unc[i])

from astropy.table import Table
t = Table([full_times_mjd, mag, mag], names=('time', 'mag', 'unc'))
#t.write('A_2017_U1_APO_DCT_2017_10_29_30_combined_lightcurve.fits', format='fits')

time_seconds = (full_times_mjd - full_times_mjd[0]) * 24 * 3600.
#longest_period_s = (dct_fulldate_mjd[0]-time_mjd[0])*24*3600.
#shortest_period_s = 1200.
#frequencies = np.arange(1.0/(longest_period_s/(2.0)),1.0/(shortest_period_s/(2.0)), 0.000001)
#periodigram = signal.lombscargle(time_seconds, mag, frequencies)

from astropy.stats import LombScargle
frequency, power = LombScargle(time_seconds, mag).autopower(samples_per_peak=100)


line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
plt.errorbar(time_mjd, apo_mag, yerr=apo_mag_unc, ecolor='black',capsize=3,capthick=1.25,markeredgecolor='black',markeredgewidth=1.2, linestyle='none')
plt.xlabel(r'$\mathrm{\mathrm{Time \; of \; observation \; (MJD)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.show()
plt.savefig('APO_only_lightcurve_2017_10_29_to_30.png')

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
plt.errorbar(dct_fulldate_mjd, dct_mag, yerr=dct_mag_unc, ecolor='black',capsize=3,capthick=1.25,markeredgecolor='black',markeredgewidth=1.2, linestyle='none')
plt.xlabel(r'$\mathrm{\mathrm{Time \; of \; observation \; (MJD)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.show()
plt.savefig('DCT_only_lightcurve_2017_10_29_to_30.png')

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
plt.errorbar(full_times_mjd, mag, yerr=mag_unc, ecolor='black',capsize=3,capthick=1.25,markeredgecolor='black',markeredgewidth=1.2, linestyle='none')
plt.xlabel(r'$\mathrm{\mathrm{Time \; of \; observation \; (MJD)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.show()
plt.savefig('APO_DCT_combined_lightcurve_2017_10_29_to_30.png')

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
#plt.plot(((1.0/frequencies) * 2.0 * np.pi)/3600., periodigram)
#plt.plot(frequencies, periodigram)
plt.plot(((1.0/frequency))/3600.,power)
plt.ylabel(r'$\mathrm{Power}$')
plt.xlabel(r'$\mathrm{Period \; (h)}$')
plt.xlim(3.14202585*.33333,23)
plt.show()
plt.savefig('APO_DCT_combined_power_spectrum_2017_10_29_to_30.png')

#MPC lightcurve

mpc_year_month_day_mag = np.loadtxt('MPC_phot_date_mag_filt', usecols=(0,1, 2, 3))
year = mpc_year_month_day_mag[:,0]
month = mpc_year_month_day_mag[:,1]
day = np.modf(mpc_year_month_day_mag[:,2])[1]
mjd_date = np.array(map(cal_date_to_mjd,year,month,day))
fract_day = np.modf(mpc_year_month_day_mag[:,2])[0]
mjd_date_with_frac =  mjd_date + fract_day

MPCmags = mpc_year_month_day_mag[:,3]
MPCmags_unc = 0.1

mpc_filt_obs_code = np.loadtxt('MPC_phot_date_mag_filt', usecols=(4,5),dtype='string')

#normalize magnitude to 2017-10-30 UTC and r'
oct_day_geo_brightness = np.loadtxt('HORIZONS_U1_brightness_GEO_10_14_10_30')

mag_correction = np.zeros(len(day))
for i in range(0,len(mag_correction)):
    brightness_correction = oct_day_geo_brightness[-1,1] - oct_day_geo_brightness[np.where(oct_day_geo_brightness[:,0]==day[i])][0,1]
    mag_correction[i] = brightness_correction

#fix mags for secular change in brightness:

MPC_fixed_mags = MPCmags + mag_correction

#V-G for sun is 0.212, V-R for sun is 0.363, V-w is 0.16, V-r is 0.214

#fix colors

#Everything to
new_r_mags = np.zeros(len(mpc_filt_obs_code[:,0]))
for i in range(0,len(mpc_filt_obs_code[:,0])):
    new_V = (magnitude_loss_filter_mean_S_C_taxonomy_ugrizywGRI_lsst(mpc_filt_obs_code[i,0]))
    new_r_mags[i] = new_V-0.214

MPCmagsfixed_brightness_color = MPC_fixed_mags + new_r_mags


for i in range(0,len(MPCmagsfixed_brightness_color)):
    print (mjd_date_with_frac[i], MPCmagsfixed_brightness_color[i], MPCmags_unc, mpc_filt_obs_code[:,1][i])


#add apo + dct stuff to mpc stuff

APO_2017_10_29_mjd_mags_mag_unc = np.loadtxt('A_2017_U1_APO_2017_10_29_lightcurve.txt')
#fix APO mags
APO_2017_10_29_mjd_mags_mag_unc[:,1] + 0.2
DCT_2017_10_30_mjd_mags_mag_unc = np.loadtxt('A_2017_U1_DCT_2017_10_30_lightcurve.txt')


#APO
for i in range(0, len(APO_2017_10_29_mjd_mags_mag_unc)):
    print (APO_2017_10_29_mjd_mags_mag_unc[:,0][i],APO_2017_10_29_mjd_mags_mag_unc[:,1][i],np.round(APO_2017_10_29_mjd_mags_mag_unc[:,2][i],3), "705")

#DCT
for i in range(0, len(DCT_2017_10_30_mjd_mags_mag_unc)):
    print (DCT_2017_10_30_mjd_mags_mag_unc[:,0][i],DCT_2017_10_30_mjd_mags_mag_unc[:,1][i],DCT_2017_10_30_mjd_mags_mag_unc[:,2][i],"G37")

#combine everything
APO_string = np.loadtxt('APO_A2017U1_2017_10_29_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
APO_date = np.loadtxt('APO_A2017U1_2017_10_29_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0))
DCT_string = np.loadtxt('DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
DCT_date = np.loadtxt('DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0))
MPC_string = np.loadtxt('MPC_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
MPC_date = np.loadtxt('MPC_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0))

#DCT + APO

combined_string = np.append(APO_string, DCT_string,axis=0)
combined_date = np.append(APO_date,DCT_date,axis=0)

combined_string_sort = combined_string[np.argsort(combined_date)]
for i in range(0,len(combined_string_sort)):
    if combined_string_sort[i,3] == '705': #normalize brightness to 2017-10-30
        print(combined_string_sort[i,0], str(float(combined_string_sort[i,1])+0.2), combined_string_sort[i,2], combined_string_sort[i,3])
    if combined_string_sort[i,3] != '705':
        print(combined_string_sort[i,0], combined_string_sort[i,1], combined_string_sort[i,2], combined_string_sort[i,3])

#DCT + APO + MPC
combined_string = np.append(np.append(APO_string, DCT_string,axis=0),MPC_string, axis=0)
combined_date = np.append(np.append(APO_date,DCT_date,axis=0), MPC_date,axis=0)

combined_string_sort = combined_string[np.argsort(combined_date)]
for i in range(0,len(combined_string_sort)):
    if combined_string_sort[i,3] == '705': #normalize brightness to 2017-10-30
        print(combined_string_sort[i,0], str(float(combined_string_sort[i,1])+0.2), combined_string_sort[i,2], combined_string_sort[i,3])
    if combined_string_sort[i,3] != '705':
        print(combined_string_sort[i,0], combined_string_sort[i,1], combined_string_sort[i,2], combined_string_sort[i,3])


#periodogram analysis

from astropy.stats import LombScargle

#Test periods combine with periodogram

num_peak = 2.0

#combine DCT + APO
DCTAPO_date_MJD_mag_mag_unc = np.loadtxt('APO_DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2))
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

#plt.ion()
test_period_1 = 6.826/24.
test_period_2 = 9.97/24

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 9.5*1.75
margin = 0.5

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(2,1,1)

minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 2.0
best_frequency = frequency[np.argmax(power)]/num_peak
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)

best_frequency = 1/test_period_1
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)

best_frequency = 1/test_period_2
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)

t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)*0.995
ax1.plot(t, y, color='black')
ax1.invert_yaxis()
ax1.set(xlabel=r'$\mathrm{Phase}$', ylabel=r'$r \; \mathrm{Magnitude}$')
#plt.title(r'$\mathrm{Phased \; data \;  at \; period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$')
plt.gca().invert_yaxis()
ax1.set_xlim(0.0,1.0)

ax1 = fig.add_subplot(2,1,2)

minimum_frequency = 1.0
maximum_frequency=80.0
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

best_frequency = frequency[np.argmax(power)]/num_peak

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
#plt.plot(((1.0/frequencies) * 2.0 * np.pi)/3600., periodigram)
#plt.plot(frequencies, periodigram)
#ax1.semilogx((1.0/(frequency/num_peak)) *24.,(power*50)+8965,color='grey')
ax1.semilogx((1.0/(frequency/num_peak)) *24.,(power),color='grey')
best_frequency1 = 1/test_period_1
ax1.axvline((1.0/(best_frequency1)) *24., color='green', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency1)*24,2))+'\;  h}$',linewidth=2.2)
ax1.axvline((1.0/(best_frequency)) *24., color='blue', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$',linewidth=2.2)
best_frequency2 = 1/test_period_2
ax1.axvline((1.0/(best_frequency2)) *24., color='orange', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency2)*24,2))+'\;  h}$',linewidth=2.2)

#ax1.set(ylabel=r'$\mathrm{Power \; level}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set(ylabel=r'$\mathrm{Power}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set_xlim(0.99,40.0)
ax1.legend(loc='upper right',prop={'size':19})

plt.savefig('APO_DCT_combined_phased_data_2017_10_29_to_30_three_periods.eps')
plt.savefig('APO_DCT_combined_phased_data_2017_10_29_to_30_three_periods.png')


#Time vs r' mag

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 9.5*1.75
margin = 0.5
#plt.ion()

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(2,1,1)
minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)
num_peak = 2.0
best_frequency = frequency[np.argmax(power)]/num_peak
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                   frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.plot(t, y, color='black')
ax1.invert_yaxis()
ax1.set(xlabel=r'$\mathrm{Phase}$', ylabel=r'$r\; \mathrm{Magnitude}$')
plt.title(r'$\mathrm{Phased \; data \;  at \; period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$')
plt.gca().invert_yaxis()
plt.xlim(0.0,1.0)
ax1 = fig.add_subplot(2,1,2)
num_peaks = 2.0
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

t = np.linspace(-2, (DCTAPO_date_MJD[-1]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2)))*1.5,10000.)*24.0
Amplitude = 2
offset = -1.0 * np.pi *1.25
y = (Amplitude * 0.5* np.sin((2*np.pi*t*(best_frequency/24.)*num_peaks)+offset)) +np.median(DCTAPO_mag)*0.997

plt.plot(t, y,alpha=0.55, color="blue",linewidth=5.0)
plt.errorbar(((DCTAPO_date_MJD-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2))))*24., DCTAPO_mag, yerr=DCTAPO_mag_unc, ecolor='black',capsize=3,capthick=1.25,markeredgecolor='black',markeredgewidth=1.2, linestyle='none')
plt.xlabel(r'$\mathrm{\mathrm{Time \; from \; MJD \;'+ str(np.round(DCTAPO_date_MJD[0],2))+' \; (hr)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.xlim(-2.5,30.0)
plt.show()
plt.savefig('APO_DCT_phase_combined_lightcurve_2017_10_29_to_30.eps')
plt.savefig('APO_DCT_phase_combined_lightcurve_2017_10_29_to_30.png')


#time vs mag different periods

num_peak = 2.0

#combine DCT + APO
DCTAPO_date_MJD_mag_mag_unc = np.loadtxt('APO_DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2))
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

#plt.ion()
test_period_1 = 6.826/24.
test_period_2 = 9.97/24

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 9.5*1.75
margin = 0.5

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(2,1,1)

minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 2.0
best_frequency = frequency[np.argmax(power)]/num_peak
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5,color='blue')

best_frequency = 1/test_period_1
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5,color='green')

best_frequency = 1/test_period_2
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5,color='orange')

t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)*0.995
ax1.plot(t, y, color='black')
ax1.invert_yaxis()
ax1.set(xlabel=r'$\mathrm{Phase}$', ylabel=r'$r \; \mathrm{Magnitude}$')
#plt.title(r'$\mathrm{Phased \; data \;  at \; period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$')
plt.gca().invert_yaxis()
ax1.set_xlim(0.0,1.0)

ax1 = fig.add_subplot(2,1,2)
num_peaks = 2.0
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 2.0
best_frequency = frequency[np.argmax(power)]/num_peak

t = np.linspace(-2, (DCTAPO_date_MJD[-1]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2)))*1.5,10000.)*24.0
Amplitude = 2
offset = -np.pi*1.25
y = (Amplitude * 0.5* np.sin((2*np.pi*t*(best_frequency/24.)*num_peaks)+offset)) +np.median(DCTAPO_mag)*0.997

plt.plot(t, y,alpha=0.55, color="blue",linewidth=5.0)
best_frequency1 = 1/test_period_1
t = np.linspace(-2, (DCTAPO_date_MJD[-1]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2)))*1.5,10000.)*24.0
Amplitude = 2
offset = -np.pi*1.25
y = (Amplitude * 0.5* np.sin((2*np.pi*t*(best_frequency1/24.)*num_peaks)+offset)) +np.median(DCTAPO_mag)*0.997
plt.plot(t, y,alpha=0.55, color="green",linewidth=5.0)
best_frequency2 = 1/test_period_2
t = np.linspace(-2, (DCTAPO_date_MJD[-1]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2)))*1.5,10000.)*24.0
Amplitude = 2
offset = -np.pi*1.25
y = (Amplitude * 0.5* np.sin((2*np.pi*t*(best_frequency2/24.)*num_peaks)+offset)) +np.median(DCTAPO_mag)*0.997
plt.plot(t, y,alpha=0.55, color="orange",linewidth=5.0)
plt.errorbar(((DCTAPO_date_MJD-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2))))*24., DCTAPO_mag, yerr=DCTAPO_mag_unc, ecolor='black',capsize=3,capthick=1.25,markeredgecolor='black',markeredgewidth=1.2, linestyle='none')
plt.xlabel(r'$\mathrm{\mathrm{Time \; from \; MJD \;'+ str(np.round(DCTAPO_date_MJD[0],2))+' \; (hr)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.xlim(-2.5,30.0)
plt.show()
plt.savefig('APO_DCT_phase_combined_lightcurve_2017_10_29_to_30_three_curves.eps')
plt.savefig('APO_DCT_phase_combined_lightcurve_2017_10_29_to_30_three_curves.png')


#matplotlib.pyplot.close("all")


#combine DCT + APO + MPC LombScargle


#double peak
num_peak = 2.0


DCTAPOMPC_date_MJD_mag_mag_unc = np.loadtxt('DCT_APO_MPC_combined_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2))
DCTAPOMPC_date_MJD = DCTAPOMPC_date_MJD_mag_mag_unc[:,0]
DCTAPOMPC_mag = DCTAPOMPC_date_MJD_mag_mag_unc[:,1]
DCTAPOMPC_mag_unc = DCTAPOMPC_date_MJD_mag_mag_unc[:,2]

from astropy.table import Table
t = Table([DCTAPOMPC_date_MJD, DCTAPOMPC_mag, DCTAPOMPC_mag_unc], names=('time', 'mag', 'unc'))
#t.write('DCT_APO_MPC_combined_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.fits', format='fits')

from astropy.stats import LombScargle

minimum_frequency = 0.5
maximum_frequency=20.
frequency, power = LombScargle(DCTAPOMPC_date_MJD, DCTAPOMPC_mag, DCTAPOMPC_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
#plt.plot(((1.0/frequencies) * 2.0 * np.pi)/3600., periodigram)
#plt.plot(frequencies, periodigram)
plt.semilogx((1.0/(frequency/num_peak)) *24.,power)
plt.ylabel(r'$\mathrm{Power}$')
plt.xlabel(r'$\mathrm{Period \; (h)}$')
plt.xlim((1/(maximum_frequency/num_peak))*24,(1/(minimum_frequency/num_peak))*24)
plt.show()
plt.savefig('APO_DCT_MPC_combined_power_spectrum_2017_10_29_to_30.png')



#double peak
minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPOMPC_date_MJD, DCTAPOMPC_mag, DCTAPOMPC_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 1.0
best_frequency = frequency[np.argmax(power)]/num_peak
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPOMPC_date_MJD, DCTAPOMPC_mag, DCTAPOMPC_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPOMPC_date_MJD * best_frequency) % 1

t = np.linspace(0, 1.0,1000.)
Amplitude = 2
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak)) +np.median(DCTAPOMPC_mag)

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
#plt.ion()
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(111)
ax1.errorbar(phase,  DCTAPOMPC_mag, DCTAPOMPC_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)
ax1.plot(phase_fit/num_peak, y_fit, color='black')
#ax1.plot(t, y, color='black')
ax1.invert_yaxis()
ax1.set(xlabel=r'$\mathrm{Phase}$', ylabel=r'$\mathrm{Magnitude}$')
plt.title(r'$\mathrm{Phased \; data \;  at \; period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$')
plt.gca().invert_yaxis()
plt.show()
plt.savefig('APO_DCT_MPC_combined_phased_data_2017_10_29_to_30.png')
