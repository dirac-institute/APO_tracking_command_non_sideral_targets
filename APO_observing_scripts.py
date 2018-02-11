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

def get_rates(rate, pa, dec_deg, dec_min, dec_sec): #rate is in "/min, pa is in degs
    #print (np.sign(dec_deg) * (dec_min/60.))
    RA = (rate * (1000./60.) * np.sin(np.radians(pa)))/np.cos(np.radians(dec_deg + ((np.sign(dec_deg) * dec_min)/60.) + ((np.sign(dec_deg) * dec_sec)/3600.)))
    DEC = (rate * (1000./60.)) * np.cos(np.radians(pa))
    return RA, DEC #mili arcsec per sec

def get_rates_no_cos_dec(rate, pa_deg): #rate is in "/min, pa is in degs USE FOR UH 88" when cos dec is turned on
    #print (np.sign(dec_deg) * (dec_min/60.))                                                                                                                                                              
    RA = (rate * (1000./60.) * np.sin(np.radians(pa_deg)))
    DEC = (rate * (1000./60.)) * np.cos(np.radians(pa_deg))
    return RA, DEC #mili arcsec per sec  

def get_rates_cos_dec(rate, pa, dec_deg, dec_min, dec_sec): #rate is in "/min, pa is in degs                                                                                                                    
    #print (np.sign(dec_deg) * (dec_min/60.))                                                                                                                                                           
    RA = (rate * (1000./60.) * np.sin(np.radians(pa)))*np.cos(np.radians(dec_deg + ((np.sign(dec_deg) * dec_min)/60.) + ((np.sign(dec_deg) * dec_sec)/3600.)))
    DEC = (rate * (1000./60.)) * np.cos(np.radians(pa))
    return RA, DEC #mili arcsec per sec           


def M_anom_and_mean_motion_to_time_of_peri(M,n,epoch_mjd_ut): #both M and n are in degrees and degrees/day.
    if M > 180:
        time_peri = epoch_mjd_ut - ((M - 360.) / n)
    if M <=180:
        time_peri = epoch_mjd_ut - (M / n)
    return time_peri
