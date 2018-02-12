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
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/APO_observing/')
from apo_observing_functions import *
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
lightcurve_hours_mag = np.loadtxt('2017u1_lc.txt')

dct_date_MJD = cal_date_to_mjd(2017,10,30)
DCT_2017_10_30_hours = lightcurve_hours_mag[:,0]
dct_fulldate_mjd = dct_date_MJD + (DCT_2017_10_30_hours/24.)
dct_fulldate_seconds = (dct_fulldate_mjd - time_mjd[0]) * 3600*24
dct_mag = lightcurve_hours_mag[:,1]
dct_mag_unc_upper = lightcurve_hours_mag[:,2]
dct_mag_unc_lower = lightcurve_hours_mag[:,3]

dct_mag_unc = (dct_mag_unc_upper + dct_mag_unc_lower)/2

print ("#Date (MJD) Magnitude r' Magnitude uncertainty")
for i in range(0,len(dct_fulldate_mjd)):
    print (dct_fulldate_mjd[i], dct_mag[i], dct_mag_unc[i])

mag = np.append(apo_mag+0.2,dct_mag)#normalized APO mags to 10/30
mag_unc = np.append(apo_mag_unc,dct_mag_unc)
full_times_mjd = np.append(time_mjd, dct_fulldate_mjd)

print ("#Date (MJD) Magnitude r' Magnitude uncertainty")
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
APO_date = np.loadtxt('APO_A2017U1_2017_10_29_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,))
DCT_string = np.loadtxt('DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
DCT_date = np.loadtxt('DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,))
MPC_string = np.loadtxt('MPC_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
MPC_date = np.loadtxt('MPC_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,))

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

num_peak = 1.0

#combine DCT + APO
DCTAPO_date_MJD_mag_mag_unc = np.loadtxt('APO_DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2))
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

#plt.ion()
test_period_1 = 4.06/24.
test_period_2 = 6.36/24

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 9.5*1.75
margin = 0.5

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
fig.subplots_adjust(hspace=.35)
ax1 = fig.add_subplot(2,1,1)

minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 1.0
best_frequency = frequency[np.argmax(power)]/num_peak
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)

#2x
num_peak = 1.0
best_frequency = (frequency[np.argmax(power)]/num_peak) *0.5
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5, color='red')

#0.5x
num_peak = 1.0
best_frequency = (frequency[np.argmax(power)]/num_peak) * 2
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5, color='grey')

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

ax1.axvline((0.5/(best_frequency)) *24., color='red', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency* 0.5)*24,2))+'\;  h}$',linewidth=2.2)

ax1.axvline((1.0/(best_frequency1)) *24., color='green', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency1)*24,2))+'\;  h}$',linewidth=2.2)

ax1.axvline((1.0/(best_frequency)) *24., color='blue', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$',linewidth=2.2)

best_frequency2 = 1/test_period_2
ax1.axvline((1.0/(best_frequency2)) *24., color='orange', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency2)*24,2))+'\;  h}$',linewidth=2.2)

ax1.axvline((2.0/(best_frequency)) *24., color='grey', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency* 2)*24,2))+'\;  h}$',linewidth=2.2)


#ax1.set(ylabel=r'$\mathrm{Power \; level}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set(ylabel=r'$\mathrm{Power}$', xlabel=r'$\mathrm{Lightcurve \; period \; (h)}$')
ax1.set_xlim(0.99,15.0)
ax1.legend(loc='upper left',prop={'size':19})

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
fig.subplots_adjust(hspace=.35)
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
ax1.errorbar(phase[:13],  DCTAPO_mag[:13]+6, DCTAPO_mag_unc[:13], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='blue')
ax1.errorbar(phase[13:],  DCTAPO_mag[13:]+6, DCTAPO_mag_unc[13:], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='orange')
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.plot(t, y+6, color="green",linewidth=3.0)
ax1.invert_yaxis()
ax1.set(ylabel=r'$r\; \mathrm{Magnitude}$')
plt.title(r'$\mathrm{Phased \; data \;  at \; period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$')
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)
num_peak = 2.0
best_frequency = (frequency[np.argmax(power)]/num_peak) * 2 # half period
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                   frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)
ax1.errorbar(phase[:13],  DCTAPO_mag[:13]+3, DCTAPO_mag_unc[:13], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='blue')
ax1.errorbar(phase[13:],  DCTAPO_mag[13:]+3, DCTAPO_mag_unc[13:], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='orange')
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.plot(t, y+3, color="red",linewidth=3.0)

frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)
num_peak = 2.0
best_frequency = frequency[np.argmax(power)]/num_peak * 0.5
phase_fit = np.linspace(0, num_peak*2)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                   frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)
ax1.plot(t, y, color="black",linewidth=3.0)
ax1.errorbar(phase[:13],  DCTAPO_mag[:13], DCTAPO_mag_unc[:13], fmt='o', mew=0, capsize=0, elinewidth=1.5,label=r'$\mathrm{2017-10-29 \; UTC \; APO \; 3.5 \; m}$',color='blue')
ax1.errorbar(phase[13:],  DCTAPO_mag[13:], DCTAPO_mag_unc[13:], fmt='o', mew=0, capsize=0, elinewidth=1.5,label=r'$\mathrm{2017-10-30 \; UTC \; DCT}$',color='orange')
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.invert_yaxis()
ax1.set(xlabel=r'$\mathrm{Phase}$', ylabel=r'$r\; \mathrm{Magnitude}$')
plt.xlim(0.0,1.0)
plt.ylim(20,35.3)
plt.legend(loc='upper right',prop={'size':19})
ax1 = fig.add_subplot(2,1,2)
num_peaks = 2.0
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

t = np.linspace(-2, (DCTAPO_date_MJD[-1]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2)))*1.5,10000.)*24.0
Amplitude = 1.7
offset = -1.0 * np.pi *1.3
y = (Amplitude * 0.5* np.sin((2*np.pi*t*(best_frequency/24.)*num_peaks)+offset)) +np.median(DCTAPO_mag)*0.99999999

plt.plot(t, y,alpha=0.55, color="grey",linewidth=3.0)
plt.errorbar(((DCTAPO_date_MJD[:12]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2))))*24., DCTAPO_mag[:12], yerr=DCTAPO_mag_unc[:12],fmt='o', mew=0, capsize=0, elinewidth=1.5)
plt.errorbar(((DCTAPO_date_MJD[12:]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2))))*24., DCTAPO_mag[12:], yerr=DCTAPO_mag_unc[12:],fmt='o', mew=0, capsize=0, elinewidth=1.5)
plt.xlabel(r'$\mathrm{\mathrm{Time \; from \; MJD \;'+ str(np.round(DCTAPO_date_MJD[0],2))+' \; (hr)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.xlim(-2.5,30.0)
plt.ylim(22,25.3)
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
fig.subplots_adjust(hspace=.35)
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


#extrapolation spectral slope thing

masiero_nm_ref = np.loadtxt('masiero_nm_ref')
nm = masiero_nm_ref[:,0]
ref = masiero_nm_ref[:,1]
nm_continuous = np.linspace(nm[0],nm[-1],1000.)

nm_function_lambda = lambda x, a, b: a*x + b
slope, intercept = 3.0, 1/1650.
slope_pm, intercept_pm = 1.5, 1/1650.
fit_ref = nm_function_lambda(nm_continuous,slope, intercept)

norm_factor = fit_ref[find_nearest(nm_continuous,500)[0]]
ref_norm = fit_ref/norm_factor

corellation_coefficient = scipy.stats.pearsonr(nm_continuous, ref_norm)[0]

#errors on each extrapolated y point
nm_continuous_vs_ref_norm_err = np.array(map(error_propagation_weighted_sum, np.ones(len(ref_norm)), np.ones(len(ref_norm)) * intercept_pm, ref_norm, np.ones(len(ref_norm)) * slope_pm, np.ones(len(ref_norm))*corellation_coefficient))

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(nm,ref,'.')
ax1.plot(nm_continuous,ref_norm)
ax1.errorbar(nm_continuous[::50],ref_norm[::50],yerr=nm_continuous_vs_ref_norm_err[::50])
ax1.set_ylim(-1.2,6.2)
ax1.set(xlabel=r'$\mathrm{wavelength \; (nm)}$', ylabel=r'$\mathrm{Normalized \; reflectance}$')

#daniela plot
#amp and amp unc
amp_frac = np.loadtxt('amplitude_pdf')

bins = amp_frac[:,0]
frac = amp_frac[:,1]


number_entries_per_bin = np.round(amp_frac[:,1]*10000000.).astype('int')

bin_width = 0.8
for i in range(0, len(number_entries_per_bin)):
    if i == 0:
        synthetic_entries = (np.ones(number_entries_per_bin[0])*bins[0])+ np.random.uniform(-1.0 * bin_width/2,bin_width/2,number_entries_per_bin[0])
        #print (number_entries_per_bin[i])
    if i > 0:
        synthetic_entries = np.append(synthetic_entries, np.ones(number_entries_per_bin[i])*bins[i]+ np.random.uniform(-1.0 * bin_width/2,bin_width/2,number_entries_per_bin[i]))
        #print (number_entries_per_bin[i])

plt.figure()
plt.hist(synthetic_entries,50)

shape, mean, unc =  skewnorm.fit(synthetic_entries[::100])


#gatspy periodogram
from gatspy import periodic

model = periodic.LombScargleFast()
model.fit(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc)
minimum_period = 0.01
maximum_period=2
periods = np.linspace(minimum_period, maximum_period, 10000)
scores = model.score(periods)
ax1 = fig.add_subplot(1,1,1)
line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
ax1.semilogx(periods,scores,color='blue')

minimum_frequency = 1.0
maximum_frequency=80.0
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

#plt.ion()
#plt.plot(((1.0/frequencies) * 2.0 * np.pi)/3600., periodigram)
#plt.plot(frequencies, periodigram)
#ax1.semilogx((1.0/(frequency/num_peak)) *24.,(power*50)+8965,color='grey')
ax1.semilogx((1.0/(frequency/num_peak)) *24.,(power),color='red')
ax1.set(ylabel=r'$\mathrm{Power}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set_xlim(0.01,2)
plt.savefig('APO_DCT_combined_phased_data_gatspy.eps')
plt.savefig('APO_DCT_combined_phased_data_gatspy.png')

#color plot showing a/b vs rho vs critical period hours

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 8.5*1.75
margin = 0.5
origin = 'lower'
n = 1000
#x = np.random.standard_normal(n)
#y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
y = np.linspace(1,7,n) #rho
x = np.linspace(3.5,10.3,n)  #a/b
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()
X, Y = np.meshgrid(x,y)
Z=critical_period_h = critical_period_axial_ratio_s_prolate(X,Y)/3600.
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
increment = 0.25
levels = np.arange(np.floor(np.min(critical_period_h)), np.ceil(np.max(critical_period_h))+increment,increment)
increment = 4
contour_levels = np.arange(np.floor(np.min(critical_period_h)), np.ceil(np.max(critical_period_h))+increment,increment)
contour_levels = np.delete(contour_levels,np.where(contour_levels==8.0)[0][0])
levels_white = np.array([8.0])
CS = plt.contourf(X, Y, Z, 10, cmap=parula_map, origin=origin, levels = levels)
CS4 = plt.contour(X, Y, Z, 10, origin=origin,linewidths=2, levels=contour_levels, alphas=0.8)
CS5 = plt.contour(X, Y, Z, 10, origin=origin,linewidths=4, levels=levels_white, alphas=0.8,colors='white')
#manual_locations_odd = np.array([(5.47,6.66), (5.47,4.01), (5.47,2.69), (5.47,1.85), (5.47,1.39), (5.47,1.1), (5.47,0.88), (5.47,0.73)])
#manual_locations_even = np.array([(6.0,3.95),(6.0,2.76),(6.0,2.02),(6.0,1.54),(6.0,1.22),(6.0,0.97),(6.0,0.8)])
#manual_locations = np.concatenate((manual_locations_odd, manual_locations_even),axis=0).tolist()
manual_locations = [(7.23,3.81),(8.26,2.88),(9,2.16), (9.5,1.7),(9.9,1.3),(10.0,1.1)]
manual_location_white = [(5.56,5.28)]
plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18,manual=manual_locations, colors='k')
plt.clabel(CS5, CS5.levels, inline=True, fmt='%1.1f', fontsize=18,manual=manual_location_white, colors='w')
#plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18, colors='k')
#plt.clabel(CS5, CS5.levels, inline=True, fmt='%1.1f', fontsize=18, colors='w')
ax1.set(xlabel=r'$b/a$', ylabel=r'$\rho \; \mathrm{(g cm^{-3})}$')
cb = plt.colorbar(CS)
tick_locs = np.array([4,6,8,10,12,14,16,18,20,22, 24,26, 28, 30, 32, 34, 36, 38, 40, 42, 44])[::2]
tick_labels = np.array([r'$4.0$',r'$6.0$',r'$8.0$', r'$10.0$', r'$12.0$', r'$14.0$', r'$16.0$', r'$18.0$',r'$20.0$',r'$22.0$',r'$24.0$',r'$26.0$',r'$28.0$',r'$30.0$',r'$32.0$',r'$34.0$',r'$36.0$',r'$38.0$',r'$40.0$',r'$42.0$',r'$44.0$'])[::2]
cb.locator     = matplotlib.ticker.FixedLocator(tick_locs)
cb.formatter   = matplotlib.ticker.FixedFormatter(tick_labels)
cb.update_ticks()
cb.set_label(r'$\mathrm{Critical \; period \; (h)}$')
plt.savefig('axial_ratio_vs_densit_vs_critical_period.pdf')
plt.savefig('axial_ratio_vs_densit_vs_critical_period.png')

#cohesion

radius_km = .130 #km assuming pv = 0.1
period_s = 8.1 * 3600.
line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 8.5*1.75
margin = 0.5
origin = 'lower'
n = 1000
#x = np.random.standard_normal(n)
#y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
y = np.linspace(1,7,n) #rho
x = np.linspace(3.5,10.3,n)  #a/b
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()
X, Y = np.meshgrid(x,y)
Z=cohseive_strenght = cohesive_strength_asteroid_axial_ratio_density_pascals(period_s, radius_km,Y,X)
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
increment = 1
#levels = np.arange(np.floor(np.min(Z)), np.ceil(np.max(Z))+increment,increment)
increment = 0.25
levels = np.arange(np.floor(np.min(Z)), np.ceil(np.max(Z))+increment,increment)
increment = 5
contour_levels = np.arange(1,61,5)-1
CS = plt.contourf(X, Y, Z, 10, cmap=parula_map, origin=origin, levels = levels)
CS4 = plt.contour(X, Y, Z, 10, origin=origin,linewidths=2, alphas=0.8, levels = contour_levels)
manual_locations = [(4.95,2.12),(5.43,2.65), (6.3,3.45),(6.78,4.1),(7.28,4.67), (7.7,4.9), (8.25,5.36), (8.8,5.87), (9.06,6.04), (9.43,6.36), (9.8,6.67)]
#fmt = {}
#strs = np.arange(3.0,33,3).astype('string').tolist()
#for l,s in zip( CS.levels, strs ):
#    fmt[l] = s
#plt.clabel(CS4, CS4.levels, inline=True, fmt=fmt, fontsize=18, colors='k',manual=manual_locations)
plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18, colors='k',manual=manual_locations)
#plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18, colors='k')
ax1.set(xlabel=r'$b/a$', ylabel=r'$\rho \; \mathrm{(g cm^{-3})}$')
cb = plt.colorbar(CS)
#tick_locs = np.array([0,3,6,9,12,15,18,21,24,27,30])
#tick_labels = np.array([r'$0.0$',r'$3.0$',r'$6.0$',r'$9.0$', r'$12.0$', r'$15.0$', r'$18.0$', r'$21.0$', r'$24.0$',r'$27.0$',r'$30.0$'])
tick_locs = np.arange(1,66,5)
tick_labels = np.array([r'$1.0$',r'$6.0$',r'$11.0$',r'$16.0$', r'$21.0$', r'$26.0$', r'$31.0$', r'$36.0$', r'$41.0$',r'$46.0$',r'$51.0$',r'$56.0$',r'$61.0$'])
cb.locator     = matplotlib.ticker.FixedLocator(tick_locs)
cb.formatter   = matplotlib.ticker.FixedFormatter(tick_labels)
cb.update_ticks()
cb.set_label(r'$\mathrm{Cohesive \; strength \; (Pa)}$')
plt.savefig('axial_ratio_vs_densit_vs_cohesive_strength.pdf')
plt.savefig('axial_ratio_vs_densit_vs_cohesive_strength.png')

#order combined data

#light time correction 2017 10 30
light_time_correction_2017_10_30 = 0.003260956790123456
date_mjd_mags_mags_unc = np.loadtxt('combined_unordered', usecols=(0,1,2))
date, mags,mags_unc = date_mjd_mags_mags_unc[:,0] , date_mjd_mags_mags_unc[:,1], date_mjd_mags_mags_unc[:,2]
obs_code = np.loadtxt('combined_unordered',usecols=(3,),dtype='string')
combined_string = np.loadtxt('combined_unordered',dtype='string')
combined_string_sort = combined_string[np.argsort(date)]
for i in range(0,len(combined_string_sort)):
    print(str(float(combined_string_sort[i,0]) - light_time_correction_2017_10_30), combined_string_sort[i,1], combined_string_sort[i,2], combined_string_sort[i,3])


#de-phasing Drahus data

#(phase + some_large_integer)/period
#drahus' data were taken on 58053.0 and 58054.0
#get phased data and add some_large_integer and divide sum by period until you get the sum to align with the observation dates taken from CDAC.


#2018 AV2 data

num_peak = 1.0

#combine DCT + APO
DCTAPO_date_MJD_mag_mag_unc = np.loadtxt('2018_AV2_APO_lightcurve_2018_01_20_uncorrected_mag_mag_unc_filter_MJD.txt',usecols=(3,0,1))
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]
DCTAPO_filter = np.loadtxt('2018_AV2_APO_lightcurve_2018_01_20_uncorrected_mag_mag_unc_filter_MJD.txt',usecols=(2,),dtype='string')

#plt.ion()
test_period_1 = 4.06/24.
test_period_2 #Bryce Bolin

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
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/APO_observing/')
from apo_observing_functions import *
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
lightcurve_hours_mag = np.loadtxt('2017u1_lc.txt')

dct_date_MJD = cal_date_to_mjd(2017,10,30)
DCT_2017_10_30_hours = lightcurve_hours_mag[:,0]
dct_fulldate_mjd = dct_date_MJD + (DCT_2017_10_30_hours/24.)
dct_fulldate_seconds = (dct_fulldate_mjd - time_mjd[0]) * 3600*24
dct_mag = lightcurve_hours_mag[:,1]
dct_mag_unc_upper = lightcurve_hours_mag[:,2]
dct_mag_unc_lower = lightcurve_hours_mag[:,3]

dct_mag_unc = (dct_mag_unc_upper + dct_mag_unc_lower)/2

print ("#Date (MJD) Magnitude r' Magnitude uncertainty")
for i in range(0,len(dct_fulldate_mjd)):
    print (dct_fulldate_mjd[i], dct_mag[i], dct_mag_unc[i])

mag = np.append(apo_mag+0.2,dct_mag)#normalized APO mags to 10/30
mag_unc = np.append(apo_mag_unc,dct_mag_unc)
full_times_mjd = np.append(time_mjd, dct_fulldate_mjd)

print ("#Date (MJD) Magnitude r' Magnitude uncertainty")
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
APO_date = np.loadtxt('APO_A2017U1_2017_10_29_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,))
DCT_string = np.loadtxt('DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
DCT_date = np.loadtxt('DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,))
MPC_string = np.loadtxt('MPC_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2,3),dtype='string')
MPC_date = np.loadtxt('MPC_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,))

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

num_peak = 1.0

#combine DCT + APO
DCTAPO_date_MJD_mag_mag_unc = np.loadtxt('APO_DCT_A2017U1_2017_10_30_date_mjd_mag_r_mag_unc_obs_code.txt',usecols=(0,1,2))
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

#plt.ion()
test_period_1 = 4.06/24.
test_period_2 = 6.36/24

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 9.5*1.75
margin = 0.5

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
fig.subplots_adjust(hspace=.35)
ax1 = fig.add_subplot(2,1,1)

minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 1.0
best_frequency = frequency[np.argmax(power)]/num_peak
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5)

#2x
num_peak = 1.0
best_frequency = (frequency[np.argmax(power)]/num_peak) *0.5
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5, color='red')

#0.5x
num_peak = 1.0
best_frequency = (frequency[np.argmax(power)]/num_peak) * 2
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                    frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
ax1.errorbar(phase,  DCTAPO_mag, DCTAPO_mag_unc, fmt='o', mew=0, capsize=0, elinewidth=1.5, color='grey')

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

ax1.axvline((0.5/(best_frequency)) *24., color='red', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency* 0.5)*24,2))+'\;  h}$',linewidth=2.2)

ax1.axvline((1.0/(best_frequency1)) *24., color='green', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency1)*24,2))+'\;  h}$',linewidth=2.2)

ax1.axvline((1.0/(best_frequency)) *24., color='blue', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$',linewidth=2.2)

best_frequency2 = 1/test_period_2
ax1.axvline((1.0/(best_frequency2)) *24., color='orange', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency2)*24,2))+'\;  h}$',linewidth=2.2)

ax1.axvline((2.0/(best_frequency)) *24., color='grey', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency* 2)*24,2))+'\;  h}$',linewidth=2.2)


#ax1.set(ylabel=r'$\mathrm{Power \; level}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set(ylabel=r'$\mathrm{Power}$', xlabel=r'$\mathrm{Lightcurve \; period \; (h)}$')
ax1.set_xlim(0.99,15.0)
ax1.legend(loc='upper left',prop={'size':19})

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
fig.subplots_adjust(hspace=.35)
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
ax1.errorbar(phase[:13],  DCTAPO_mag[:13]+6, DCTAPO_mag_unc[:13], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='blue')
ax1.errorbar(phase[13:],  DCTAPO_mag[13:]+6, DCTAPO_mag_unc[13:], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='orange')
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.plot(t, y+6, color="green",linewidth=3.0)
ax1.invert_yaxis()
ax1.set(ylabel=r'$r\; \mathrm{Magnitude}$')
plt.title(r'$\mathrm{Phased \; data \;  at \; period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$')
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)
num_peak = 2.0
best_frequency = (frequency[np.argmax(power)]/num_peak) * 2 # half period
phase_fit = np.linspace(0, num_peak)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                   frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)
ax1.errorbar(phase[:13],  DCTAPO_mag[:13]+3, DCTAPO_mag_unc[:13], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='blue')
ax1.errorbar(phase[13:],  DCTAPO_mag[13:]+3, DCTAPO_mag_unc[13:], fmt='o', mew=0, capsize=0, elinewidth=1.5,color='orange')
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.plot(t, y+3, color="red",linewidth=3.0)

frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)
num_peak = 2.0
best_frequency = frequency[np.argmax(power)]/num_peak * 0.5
phase_fit = np.linspace(0, num_peak*2)
y_fit = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).model(t=phase_fit / (best_frequency),
                                   frequency=best_frequency)
phase = (DCTAPO_date_MJD * best_frequency) % 1
t = np.linspace(0, 1.0,1000.)
Amplitude = 2
set_phase = np.pi*1.15
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak + set_phase)) +np.median(DCTAPO_mag)
ax1.plot(t, y, color="black",linewidth=3.0)
ax1.errorbar(phase[:13],  DCTAPO_mag[:13], DCTAPO_mag_unc[:13], fmt='o', mew=0, capsize=0, elinewidth=1.5,label=r'$\mathrm{2017-10-29 \; UTC \; APO \; 3.5 \; m}$',color='blue')
ax1.errorbar(phase[13:],  DCTAPO_mag[13:], DCTAPO_mag_unc[13:], fmt='o', mew=0, capsize=0, elinewidth=1.5,label=r'$\mathrm{2017-10-30 \; UTC \; DCT}$',color='orange')
#ax1.plot(phase_fit[::-1]/num_peak, y_fit, color='black')
ax1.invert_yaxis()
ax1.set(xlabel=r'$\mathrm{Phase}$', ylabel=r'$r\; \mathrm{Magnitude}$')
plt.xlim(0.0,1.0)
plt.ylim(20,35.3)
plt.legend(loc='upper right',prop={'size':19})
ax1 = fig.add_subplot(2,1,2)
num_peaks = 2.0
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]

t = np.linspace(-2, (DCTAPO_date_MJD[-1]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2)))*1.5,10000.)*24.0
Amplitude = 1.7
offset = -1.0 * np.pi *1.3
y = (Amplitude * 0.5* np.sin((2*np.pi*t*(best_frequency/24.)*num_peaks)+offset)) +np.median(DCTAPO_mag)*0.99999999

plt.plot(t, y,alpha=0.55, color="grey",linewidth=3.0)
plt.errorbar(((DCTAPO_date_MJD[:12]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2))))*24., DCTAPO_mag[:12], yerr=DCTAPO_mag_unc[:12],fmt='o', mew=0, capsize=0, elinewidth=1.5)
plt.errorbar(((DCTAPO_date_MJD[12:]-DCTAPO_date_MJD[0] + (DCTAPO_date_MJD[0]-np.round(DCTAPO_date_MJD[0],2))))*24., DCTAPO_mag[12:], yerr=DCTAPO_mag_unc[12:],fmt='o', mew=0, capsize=0, elinewidth=1.5)
plt.xlabel(r'$\mathrm{\mathrm{Time \; from \; MJD \;'+ str(np.round(DCTAPO_date_MJD[0],2))+' \; (hr)}}$')
plt.ylabel(r'$r \; \mathrm{Magnitude}$')
plt.xlim(-2.5,30.0)
plt.ylim(22,25.3)
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
fig.subplots_adjust(hspace=.35)
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


#extrapolation spectral slope thing

masiero_nm_ref = np.loadtxt('masiero_nm_ref')
nm = masiero_nm_ref[:,0]
ref = masiero_nm_ref[:,1]
nm_continuous = np.linspace(nm[0],nm[-1],1000.)

nm_function_lambda = lambda x, a, b: a*x + b
slope, intercept = 3.0, 1/1650.
slope_pm, intercept_pm = 1.5, 1/1650.
fit_ref = nm_function_lambda(nm_continuous,slope, intercept)

norm_factor = fit_ref[find_nearest(nm_continuous,500)[0]]
ref_norm = fit_ref/norm_factor

corellation_coefficient = scipy.stats.pearsonr(nm_continuous, ref_norm)[0]

#errors on each extrapolated y point
nm_continuous_vs_ref_norm_err = np.array(map(error_propagation_weighted_sum, np.ones(len(ref_norm)), np.ones(len(ref_norm)) * intercept_pm, ref_norm, np.ones(len(ref_norm)) * slope_pm, np.ones(len(ref_norm))*corellation_coefficient))

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(nm,ref,'.')
ax1.plot(nm_continuous,ref_norm)
ax1.errorbar(nm_continuous[::50],ref_norm[::50],yerr=nm_continuous_vs_ref_norm_err[::50])
ax1.set_ylim(-1.2,6.2)
ax1.set(xlabel=r'$\mathrm{wavelength \; (nm)}$', ylabel=r'$\mathrm{Normalized \; reflectance}$')

#daniela plot
#amp and amp unc
amp_frac = np.loadtxt('amplitude_pdf')

bins = amp_frac[:,0]
frac = amp_frac[:,1]


number_entries_per_bin = np.round(amp_frac[:,1]*10000000.).astype('int')

bin_width = 0.8
for i in range(0, len(number_entries_per_bin)):
    if i == 0:
        synthetic_entries = (np.ones(number_entries_per_bin[0])*bins[0])+ np.random.uniform(-1.0 * bin_width/2,bin_width/2,number_entries_per_bin[0])
        #print (number_entries_per_bin[i])
    if i > 0:
        synthetic_entries = np.append(synthetic_entries, np.ones(number_entries_per_bin[i])*bins[i]+ np.random.uniform(-1.0 * bin_width/2,bin_width/2,number_entries_per_bin[i]))
        #print (number_entries_per_bin[i])

plt.figure()
plt.hist(synthetic_entries,50)

shape, mean, unc =  skewnorm.fit(synthetic_entries[::100])


#gatspy periodogram
from gatspy import periodic

model = periodic.LombScargleFast()
model.fit(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc)
minimum_period = 0.01
maximum_period=2
periods = np.linspace(minimum_period, maximum_period, 10000)
scores = model.score(periods)
ax1 = fig.add_subplot(1,1,1)
line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
ax1.semilogx(periods,scores,color='blue')

minimum_frequency = 1.0
maximum_frequency=80.0
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

#plt.ion()
#plt.plot(((1.0/frequencies) * 2.0 * np.pi)/3600., periodigram)
#plt.plot(frequencies, periodigram)
#ax1.semilogx((1.0/(frequency/num_peak)) *24.,(power*50)+8965,color='grey')
ax1.semilogx((1.0/(frequency/num_peak)) *24.,(power),color='red')
ax1.set(ylabel=r'$\mathrm{Power}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set_xlim(0.01,2)
plt.savefig('APO_DCT_combined_phased_data_gatspy.eps')
plt.savefig('APO_DCT_combined_phased_data_gatspy.png')

#color plot showing a/b vs rho vs critical period hours

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 8.5*1.75
margin = 0.5
origin = 'lower'
n = 1000
#x = np.random.standard_normal(n)
#y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
y = np.linspace(1,7,n) #rho
x = np.linspace(3.5,10.3,n)  #a/b
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()
X, Y = np.meshgrid(x,y)
Z=critical_period_h = critical_period_axial_ratio_s_prolate(X,Y)/3600.
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
increment = 0.25
levels = np.arange(np.floor(np.min(critical_period_h)), np.ceil(np.max(critical_period_h))+increment,increment)
increment = 4
contour_levels = np.arange(np.floor(np.min(critical_period_h)), np.ceil(np.max(critical_period_h))+increment,increment)
contour_levels = np.delete(contour_levels,np.where(contour_levels==8.0)[0][0])
levels_white = np.array([8.0])
CS = plt.contourf(X, Y, Z, 10, cmap=parula_map, origin=origin, levels = levels)
CS4 = plt.contour(X, Y, Z, 10, origin=origin,linewidths=2, levels=contour_levels, alphas=0.8)
CS5 = plt.contour(X, Y, Z, 10, origin=origin,linewidths=4, levels=levels_white, alphas=0.8,colors='white')
#manual_locations_odd = np.array([(5.47,6.66), (5.47,4.01), (5.47,2.69), (5.47,1.85), (5.47,1.39), (5.47,1.1), (5.47,0.88), (5.47,0.73)])
#manual_locations_even = np.array([(6.0,3.95),(6.0,2.76),(6.0,2.02),(6.0,1.54),(6.0,1.22),(6.0,0.97),(6.0,0.8)])
#manual_locations = np.concatenate((manual_locations_odd, manual_locations_even),axis=0).tolist()
manual_locations = [(7.23,3.81),(8.26,2.88),(9,2.16), (9.5,1.7),(9.9,1.3),(10.0,1.1)]
manual_location_white = [(5.56,5.28)]
plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18,manual=manual_locations, colors='k')
plt.clabel(CS5, CS5.levels, inline=True, fmt='%1.1f', fontsize=18,manual=manual_location_white, colors='w')
#plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18, colors='k')
#plt.clabel(CS5, CS5.levels, inline=True, fmt='%1.1f', fontsize=18, colors='w')
ax1.set(xlabel=r'$b/a$', ylabel=r'$\rho \; \mathrm{(g cm^{-3})}$')
cb = plt.colorbar(CS)
tick_locs = np.array([4,6,8,10,12,14,16,18,20,22, 24,26, 28, 30, 32, 34, 36, 38, 40, 42, 44])[::2]
tick_labels = np.array([r'$4.0$',r'$6.0$',r'$8.0$', r'$10.0$', r'$12.0$', r'$14.0$', r'$16.0$', r'$18.0$',r'$20.0$',r'$22.0$',r'$24.0$',r'$26.0$',r'$28.0$',r'$30.0$',r'$32.0$',r'$34.0$',r'$36.0$',r'$38.0$',r'$40.0$',r'$42.0$',r'$44.0$'])[::2]
cb.locator     = matplotlib.ticker.FixedLocator(tick_locs)
cb.formatter   = matplotlib.ticker.FixedFormatter(tick_labels)
cb.update_ticks()
cb.set_label(r'$\mathrm{Critical \; period \; (h)}$')
plt.savefig('axial_ratio_vs_densit_vs_critical_period.pdf')
plt.savefig('axial_ratio_vs_densit_vs_critical_period.png')

#cohesion

radius_km = .130 #km assuming pv = 0.1
period_s = 8.1 * 3600.
line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 8.5*1.75
margin = 0.5
origin = 'lower'
n = 1000
#x = np.random.standard_normal(n)
#y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
y = np.linspace(1,7,n) #rho
x = np.linspace(3.5,10.3,n)  #a/b
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()
X, Y = np.meshgrid(x,y)
Z=cohseive_strenght = cohesive_strength_asteroid_axial_ratio_density_pascals(period_s, radius_km,Y,X)
fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
ax1 = fig.add_subplot(1,1,1)
increment = 1
#levels = np.arange(np.floor(np.min(Z)), np.ceil(np.max(Z))+increment,increment)
increment = 0.25
levels = np.arange(np.floor(np.min(Z)), np.ceil(np.max(Z))+increment,increment)
increment = 5
contour_levels = np.arange(1,61,5)-1
CS = plt.contourf(X, Y, Z, 10, cmap=parula_map, origin=origin, levels = levels)
CS4 = plt.contour(X, Y, Z, 10, origin=origin,linewidths=2, alphas=0.8, levels = contour_levels)
manual_locations = [(4.95,2.12),(5.43,2.65), (6.3,3.45),(6.78,4.1),(7.28,4.67), (7.7,4.9), (8.25,5.36), (8.8,5.87), (9.06,6.04), (9.43,6.36), (9.8,6.67)]
#fmt = {}
#strs = np.arange(3.0,33,3).astype('string').tolist()
#for l,s in zip( CS.levels, strs ):
#    fmt[l] = s
#plt.clabel(CS4, CS4.levels, inline=True, fmt=fmt, fontsize=18, colors='k',manual=manual_locations)
plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18, colors='k',manual=manual_locations)
#plt.clabel(CS4, CS4.levels, inline=True, fmt='%1.1f', fontsize=18, colors='k')
ax1.set(xlabel=r'$b/a$', ylabel=r'$\rho \; \mathrm{(g cm^{-3})}$')
cb = plt.colorbar(CS)
#tick_locs = np.array([0,3,6,9,12,15,18,21,24,27,30])
#tick_labels = np.array([r'$0.0$',r'$3.0$',r'$6.0$',r'$9.0$', r'$12.0$', r'$15.0$', r'$18.0$', r'$21.0$', r'$24.0$',r'$27.0$',r'$30.0$'])
tick_locs = np.arange(1,66,5)
tick_labels = np.array([r'$1.0$',r'$6.0$',r'$11.0$',r'$16.0$', r'$21.0$', r'$26.0$', r'$31.0$', r'$36.0$', r'$41.0$',r'$46.0$',r'$51.0$',r'$56.0$',r'$61.0$'])
cb.locator     = matplotlib.ticker.FixedLocator(tick_locs)
cb.formatter   = matplotlib.ticker.FixedFormatter(tick_labels)
cb.update_ticks()
cb.set_label(r'$\mathrm{Cohesive \; strength \; (Pa)}$')
plt.savefig('axial_ratio_vs_densit_vs_cohesive_strength.pdf')
plt.savefig('axial_ratio_vs_densit_vs_cohesive_strength.png')

#order combined data

#light time correction 2017 10 30
light_time_correction_2017_10_30 = 0.003260956790123456
date_mjd_mags_mags_unc = np.loadtxt('combined_unordered', usecols=(0,1,2))
date, mags,mags_unc = date_mjd_mags_mags_unc[:,0] , date_mjd_mags_mags_unc[:,1], date_mjd_mags_mags_unc[:,2]
obs_code = np.loadtxt('combined_unordered',usecols=(3,),dtype='string')
combined_string = np.loadtxt('combined_unordered',dtype='string')
combined_string_sort = combined_string[np.argsort(date)]
for i in range(0,len(combined_string_sort)):
    print(str(float(combined_string_sort[i,0]) - light_time_correction_2017_10_30), combined_string_sort[i,1], combined_string_sort[i,2], combined_string_sort[i,3])


#de-phasing Drahus data

#(phase + some_large_integer)/period
#drahus' data were taken on 58053.0 and 58054.0
#get phased data and add some_large_integer and divide sum by period until you get the sum to align with the observation dates taken from CDAC.


#2018 AV2 data
#2018 02 11
#put a decent dent in analysis by writing photometry script and doing basic photometry
#not able to measure a period with single 10s image data
#recommend stacking 2-4 recomputing photometry and redoing period analysis
#colors it is recommended to stack all data by filter before frame object_g.0201.fits and do colors

num_peak = 1.0

#combine DCT + APO
DCTAPO_date_MJD_mag_mag_unc = np.loadtxt('2018_AV2_APO_lightcurve_2018_01_20_uncorrected_mag_mag_unc_filter_MJD.txt',usecols=(3,0,1))
DCTAPO_date_MJD = DCTAPO_date_MJD_mag_mag_unc[:,0]
DCTAPO_mag = DCTAPO_date_MJD_mag_mag_unc[:,1]
DCTAPO_mag_unc = DCTAPO_date_MJD_mag_mag_unc[:,2]
DCTAPO_filter = np.loadtxt('2018_AV2_APO_lightcurve_2018_01_20_uncorrected_mag_mag_unc_filter_MJD.txt',usecols=(2,),dtype='string')

cuts = np.where((DCTAPO_filter=='r') & (DCTAPO_mag_unc<0.3)& (DCTAPO_mag>19.5))
#clean data r only
DCTAPO_date_MJD = DCTAPO_date_MJD[cuts]
DCTAPO_mag = DCTAPO_mag[cuts]
DCTAPO_mag_unc = DCTAPO_mag_unc[cuts]


#plt.ion()
test_period_1 = 4.06/24.

line_width = 2.5
mult = 1.2
paperheight = 6.5*1.75
paperwidth = 9.5*1.75
margin = 0.5

fig = plt.figure(figsize=(paperwidth - 2*margin, paperheight - 2*margin))
fig.subplots_adjust(hspace=.35)
ax1 = fig.add_subplot(2,1,1)

minimum_frequency = 1.0
maximum_frequency=40.
frequency, power = LombScargle(DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc).autopower(samples_per_peak=1000, minimum_frequency = minimum_frequency, maximum_frequency=maximum_frequency)

num_peak = 1.0
best_frequency = frequency[np.argmax(power)]/num_peak
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

ax1.axvline((1.0/(best_frequency)) *24., color='blue', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round((1/best_frequency)*24,2))+'\;  h}$',linewidth=2.2)

#ax1.set(ylabel=r'$\mathrm{Power \; level}$', xlabel=r'$\mathrm{Period \; (h)}$')
ax1.set(ylabel=r'$\mathrm{Power}$', xlabel=r'$\mathrm{Lightcurve \; period \; (h)}$')
ax1.set_xlim(0.99,15.0)
ax1.legend(loc='upper left',prop={'size':19})

plt.savefig('APO_2018_AV2_2018_01_20_phase.pdf')
plt.savefig('APO_2018_AV2_2018_01_20_phase.png')

