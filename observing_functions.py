from matplotlib.ticker import FuncFormatter
import pylab as plt
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import string
import os
import numpy as np
import random
from decimal import *
import warnings
import re
import pylab
import math
from matplotlib.widgets import Slider
from bisect import bisect_left
from bisect import bisect_right
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy import stats
from sklearn.neighbors import KernelDensity
from statsmodels.nonparametric.kde import KDEUnivariate
from statsmodels.nonparametric.kernel_density import KDEMultivariate
#from family_functions import *
import pyslalib.slalib as sla
import subprocess
import astropy.io.fits as pyfits


sed_paretheses =  'awk -F "[()]" \'{ for (i=2; i<NF; i+=2) print $i }\''

hours_to_deg = 15.0
years_to_months  = 12.
years_to_days = 365.
days_to_hours = 24.
hours_to_minutes = 60.
minutes_to_seconds = 60.
deg_to_arcsec = 3600.
djcal_precision = 5
astrores_seconds_rounding_format_RA = 2
astrores_seconds_rounding_format_Dec = 1
number_array_entries_for_hms = 3
number_coordinates = 2 #1 for RA and 1 for Dec
pipe = '|'
colon = ':'
epoch_cfht = '2000.0'
epoch_2000 = 2000.0
j_2000_mjd = '51544.5'
one = '1'
plus = '+'
minus = '-'
zero = '0'
blank = ''
space = ' '
dot = '.'
underscore = '_'
newline = '\n'
object_name_max_length_cfht = 20
RA_s_length_decimal = 3#string length for astrores formating
Jy_to_mJy = 1000.0
numbered_asteroid_string_max_length = 5
UT_to_HST = (10.0 / days_to_hours)
G = 6.67408e-11 #m^3 kg e^-1 s^-2
mass_sun_kg = 1.989e30
au_to_meters = 1.496e11
seconds_to_days = 1.0 / (24.*3600.)

visir_constant = 36.0

'''
if canceling y motion, there is only motion in x:
you start offset from 15" from center top of chip 22, so you have the 186" of half the current gap, 4 chips and 4 gaps for a total of
186 + (4 * 372) + (4 * 14) = 1730".
'''
x_positive_direction_center_offset_arcsec = 1730.
x_positive_direction_center_offset_arcsec = 1730.

'''
cancel motion offset directions for block time calculation:
if canceling x motion, there is only motion in y:
you start offset from 15" from center top of chip 22, so you have the 11" gap in the middle of the mosaic, 800" of chip,
then 81" of gap then 800" of chip again for a total of 1707".
'''
y_positive_direction_center_offset_arcsec = 1692.
y_negative_direction_center_offset_arcsec = 1707.

def angular_separation(RA_1,Dec_1,RA_2,Dec_2):
    return sla.sla_sep(RA_1,Dec_1,RA_2,Dec_2)

def append_zeros_less_than_10_beginning(number):
    number_abs = np.abs(number)
    if number_abs<10 and number_abs>=1:
        number_str = "0" + str(number_abs)

    if number_abs<1:
        number_str = "0" + str(number_abs)

    if number_abs >= 10:
        number_str = str(number_abs)
    return number_str

def append_zeros_end(number_str,required_spaces):
    number_spaces = len(number_str[number_str.find('.')+1:])
    if number_spaces < required_spaces:
        for i in range(0, int(required_spaces-number_spaces)):
            number_str += zero
    return number_str

def cal_date_to_mjd(year,month,date):
    return sla.sla_caldj(year,month,date)[0]

def cal_date_fits_format_to_mjd(date_fits_header_string):
    #converts fits format of header DATE-OBS, e.g.,  2017-10-29T08:37:21.386925 to MJD
    calendar_date = date_fits_header_string[:date_fits_header_string.find('T')]
    year, month, day = calendar_date.split('-')
    date_mjd = sla.sla_caldj(int(year), int(month), int(day))[0]

    day_split = date_fits_header_string[date_fits_header_string.find('T')+1:].split(':')
    hour,minute,second = [float(i) for i in day_split]
    fraction_day = (hour/24.0) + (minute/(60*24.)) + (second / (24*3600.))

    return date_mjd + fraction_day

def convert_deg_to_hms_RA(deg):
    decimal_m, h = np.modf(deg/hours_to_deg )
    decimal_s, m = np.modf(np.abs(decimal_m) * hours_to_minutes)
    s = np.round(decimal_s*minutes_to_seconds,astrores_seconds_rounding_format_RA)
    return h,m,s

def convert_deg_to_hms_Dec(deg):
    decimal_m, h = np.modf(deg)
    decimal_s, m = np.modf(np.abs(decimal_m) * hours_to_minutes)
    s = np.round(decimal_s*minutes_to_seconds,astrores_seconds_rounding_format_Dec)
    return h,m,s

def convert_h_m_s_to_deg(h,m,s):
    return (h * np.sign(h)) + ((m* np.sign(h))/60.) + ((s* np.sign(h))/3600.)

def convert_h_m_s_to_deg_no_sign(h,m,s):
    return (h) + ((m)/60.) + ((s)/3600.)

def create_frame_for_stacking(fits_file, x_pos, y_pos):
    #must be square
    buffer = 10
    datfile = pyfits.getdata(fits_file, header=True)
    dat_raw = datfile[0][::-1,:] #must flip data then flip back
    dat_head = datfile[1]
    frame_size = pyfits.open(fits_file)[0].header['NAXIS1']
    standard_frame_size = int(np.ceil(np.sqrt(frame_size**2 + frame_size**2) + buffer))
    null_frame = np.zeros(standard_frame_size * standard_frame_size).reshape(standard_frame_size,standard_frame_size)
    center_y, center_x = np.array(null_frame.shape)/2
    add_array_edge_coord_y, add_array_edge_coord_x = center_y - y_pos, center_x - x_pos
    #add arrays on top of each others
    null_frame[add_array_edge_coord_y:add_array_edge_coord_y+frame_size,add_array_edge_coord_x:add_array_edge_coord_x+frame_size] += dat_raw
    return null_frame[::-1,:]



def css_efficiency(m,epsilon_0, m_lim, m_drop):
    return epsilon_0 / (1 + np.exp((m -m_lim)/m_drop))

def dictionary_month(numerical_month):
    month = dict([[1,'Jan'], [2,'Feb'], [3,'Mar'], [4,'Apr'], [5,'May'], [6,'Jun'], [7,'Jul'], [8,'Aug'], [9,'Sep'], [10,'Oct'], [11,'Nov'], [12,'Dec']])
    return month[numerical_month]

def get_rates_no_cos_dec(rate, pa_deg): #rate is in "/min, pa is in degs USE FOR UH 88" when cos dec is turned on
    RA = rate * np.sin(np.radians(pa_deg))
    DEC = rate  * np.cos(np.radians(pa_deg))
    return RA, DEC #mili arcsec per sec

def magnitude_calc(zero_point_magnitude, flux_circle_counts, aperture_radius_pixels, sky_flux_counts, exp_time_s):
    #flux is flux_circle - sky_value* area
    flux_counts = flux_circle_counts - (aperture_radius_pixels**2 * np.pi * sky_flux_counts)
    M = zero_point_magnitude - (2.5*np.log10(flux_counts)) - (2.5 * np.log10(exp_time_s))
    return M

def M_anom_and_mean_motion_to_time_of_peri(M,n,epoch_mjd_ut): #both M and n are in degrees and degrees/day.
    if M > 180:
        time_peri = epoch_mjd_ut - ((M - 360.) / n)
    if M <=180:
        time_peri = epoch_mjd_ut - (M / n)
    return time_peri

def mjd_to_mpc_date(time_ut_mjd):
    end_time_year, end_time_month, end_time_day, end_time_day_fraction, marker = sla.sla_djcl(time_ut_mjd)
    end_time_hours_frac, end_time_hour = np.modf(end_time_day_fraction*days_to_hours)
    end_time_minutes = np.round(end_time_hours_frac * hours_to_minutes, 0)
    end_time_ut = np.array([end_time_year, end_time_month, end_time_day, end_time_hour, end_time_minutes])
    return end_time_ut

def mjd_to_mpc_date_inc_seconds(time_ut_mjd):
    end_time_year, end_time_month, end_time_day, end_time_day_fraction, marker = sla.sla_djcl(time_ut_mjd)
    end_time_hours_frac, end_time_hour = np.modf(end_time_day_fraction*days_to_hours)
    end_time_min_frac, end_time_minutes = np.modf(end_time_hours_frac * hours_to_minutes)
    end_time_sec = np.round(end_time_min_frac * minutes_to_seconds,0)
    if end_time_sec == minutes_to_seconds:
        end_time_ut = np.array([end_time_year, end_time_month, end_time_day, end_time_hour, end_time_minutes+1,0.0])
    if end_time_sec !=minutes_to_seconds:
        end_time_ut = np.array([end_time_year, end_time_month, end_time_day, end_time_hour, end_time_minutes,end_time_sec])
    return end_time_ut

def semi_major_axis_massive_body_to_orbital_period_seconds(a_au,mass_kg):
    a_m = a_au *au_to_meters
    period_seconds = np.pi * 2.0 * np.sqrt(a_m**3 / (G*mass_kg))
    return period_seconds


def phaseFunction ( pa,G, ast_hel_dist, g_dist , system):#phase angle, G12 slope, r, Delta
# returns the phase function e.g. V = H - phaseFunction using G12 system, Muinonen et al. 2010
    '''


    '''
    #For HG1G2 basis functions
    knots12 = np.array([
         0.13089969389957470, 0.13089969389957470,
         0.52359877559829882, 0.52359877559829882,
         1.0471975511965976, 1.0471975511965976,
         1.5707963267948966, 1.5707963267948966,
         2.0943951023931953, 2.0943951023931953,
         2.6179938779914944, 2.6179938779914944,
         0.75000000000000000, 0.92500000000000004,
         0.33486016000000002, 0.62884169000000001,
         0.13410559999999999, 0.31755495000000000,
         5.11047560000000012E-002, 0.12716367000000001,
         2.14656870000000007E-002, 2.23739030000000005E-002,
         3.63969890000000011E-003, 1.65056890000000008E-004,
         -1.9098593171027440, -0.57295779513082323,
         -0.55463432402954305, -0.76705367127692980,
         -0.24404598605661942, -0.45665789025666331,
         -9.49804380669390519E-002, -0.28071808974438800,
         -2.14114236377018034E-002, -0.11173256930106369,
         -9.13286120000000035E-002, -8.65731380000000014E-008])\
        .reshape(3,6,2)

    knots3 = np.array([
         0.00000000000000000, 5.23598775598298812e-003,
         1.74532925199432955E-002, 3.49065850398865909E-002,
         6.98131700797731819E-002, 0.13962634015954636,
         0.20943951023931956, 0.34906585039886590,
         0.52359877559829882, 1.0000000000000000,
         0.83381185000000002, 0.57735424000000002,
         0.42144772000000003, 0.23174230000000001,
         0.10348178000000000, 6.17334729999999970E-002,
         1.61070060000000001E-002, 0.00000000000000000E+000,
         -0.10630096999999999, -41.180439484750558,
         -10.366914741841031, -7.5784614901018310,
         -3.6960950182737276, -0.78605651603897575,
         -0.46527011787268774, -0.20459544750797629,
         0.00000000000000000E+000])\
        .reshape(3,9)


    if system == 'HG12':
        G12 = G

        basis_function = np.zeros(3)
        if G12 < 0.2:
            G1 = (0.7527*G12) + 0.06164
            G2 = (-0.9612*G12) + 0.6270
        if G12 >= 0.2:
            G1 = (0.9529*G12) + 0.02162
            G2 = (-0.6125*G12) + 0.5572

        if pa < 0.1308996938995747: #7.5 degrees, basis functions 1 and 2
            for i in range(0,2):
                basis_function[i] = ((pa - knots12[i,0,1]) * knots12[2,0,i]) + knots12[1,0,i]

        if pa >=  0.1308996938995747:#7.5 degrees
            for i in range(0,2):#loop over basis functions 1 and 2
                for j in range(0, knots12.shape[1]-1): #loop over phase-angle bins
                    if knots12[0,j,i] <= pa and pa <= knots12[0,j+1,i]:
                        phase_angle_lo = knots12[0,j,i]
                        phase_angle_hi = knots12[0,j+1,i]
                        basis_function_lo = knots12[1,j,i]
                        basis_function_hi = knots12[1,j+1,i]
                        phase_angle_derivative_lo = knots12[2,j,i]
                        phase_angle_derivative_hi = knots12[2,j+1,i]
                        Delta_phase_angle = phase_angle_hi - phase_angle_lo
                        Delta_basis_function = basis_function_hi - basis_function_lo
                        a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
                        b = (-1.0*phase_angle_derivative_hi*Delta_phase_angle) + Delta_basis_function
                        t = (pa - phase_angle_lo) / Delta_phase_angle
                        basis_function[i] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))

        if pa > 2.6179938779914944:#150.0 degrees
            #2
            phase_angle_lo = knots12[0,5,0]
            phase_angle_hi = pa
            basis_function_lo = knots12[1,4,0]
            basis_function_hi = knots12[1,5,0]
            phase_angle_derivative_lo = knots12[2,4,0]
            phase_angle_derivative_hi = knots12[2,5,0]
            Delta_phase_angle = phase_angle_hi - phase_angle_lo
            Delta_basis_function = basis_function_hi - basis_function_lo
            a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
            b = (-1.0*phase_angle_derivative_hi*Delta_phase_angle) + Delta_basis_function
            t = (pa - phase_angle_lo) / Delta_phase_angle
            basis_function[0] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))
            #2
            phase_angle_lo = knots12[0,5,1]
            phase_angle_hi = pa
            basis_function_lo = knots12[1,4,1]
            basis_function_hi = knots12[1,5,1]
            phase_angle_derivative_lo = knots12[2,4,1]
            phase_angle_derivative_hi = knots12[2,5,1]
            Delta_phase_angle = phase_angle_hi - phase_angle_lo
            Delta_basis_function = basis_function_hi - basis_function_lo
            a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
            b = (-1.0*phase_angle_derivative_hi*Delta_phase_angle) + Delta_basis_function
            t = (pa - phase_angle_lo) / Delta_phase_angle
            basis_function[1] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))

        #basis function 3
        if pa > 0.52359877559829882:#30 degrees
            basis_function[2] = 0.0
        if pa <= 0.52359877559829882:
            for i in range(0, knots3.shape[1] -1):
                if knots3[0,i] <= pa and pa <= knots3[0,i+1]:
                    phase_angle_lo = knots3[0,i]
                    phase_angle_hi = knots3[0,i+1]
                    basis_function_lo = knots3[1,i]
                    basis_function_hi = knots3[1,i+1]
                    phase_angle_derivative_lo = knots3[2,i]
                    phase_angle_derivative_hi = knots3[2,i+1]
                    Delta_phase_angle = phase_angle_hi - phase_angle_lo
                    Delta_basis_function = basis_function_hi - basis_function_lo
                    a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
                    b = ((-1.0 *phase_angle_derivative_hi)*Delta_phase_angle) + Delta_basis_function
                    t = (pa - phase_angle_lo) / Delta_phase_angle
                    basis_function[2] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))

        appr_mag = -2.5 * np.log10( (G1 * basis_function[0]) + (G2 * basis_function[1]) + ((1.0 - G1 - G2) * basis_function[2])) + 5.0*np.log10(ast_hel_dist*g_dist)

    if system == 'HG':
        phi_1 =  (np.exp(-90.56*np.tan(0.5*pa)**2) * (1.0 - (0.986*np.sin(pa)) / ( 0.119 + (1.341 * np.sin(pa)) - (0.754*(np.sin(pa)**2)) ) ) ) + ( (1 - np.exp(-90.56*np.tan(0.5*pa)**2)) * ( np.exp(-3.332*np.tan(0.5*pa)**0.631) ) )
        phi_2 =  (np.exp(-90.56*np.tan(0.5*pa)**2) * (1.0 - (0.238*np.sin(pa)) / ( 0.119 + (1.341 * np.sin(pa)) - (0.754*(np.sin(pa)**2)) ) ) ) + ( (1 - np.exp(-90.56*np.tan(0.5*pa)**2)) * ( np.exp(-1.862*np.tan(0.5*pa)**1.218) ) )

        appr_mag = -2.5*np.log10( ((1.0 - G)*phi_1) + (G*phi_2) ) + 5.0*np.log10(ast_hel_dist*g_dist)

    return appr_mag

def phase_function_bowell_muinonen (system, phase_angle_radians, G_or_G_12, sun_to_asteroid_distance_au, observer_to_asteroid_distance_au):
    """

    Arguments:
        :param system (str): only two possibilities 'HG' or 'HG12'. The HG system uses the classic photometric model from Bowell et al. 1988 'Application of Photometric Models to Asteroids'. Asteroids IIThe HG12 system uses the photometric model from Muinonen et al. 'A three-parameter magnitude phase function for asteroids'.
        :param phase_angle_radians (float): the sun to asteroid to observer angle in radians.
        :param G_or_G_12 (float): The phase slope parameter. It is called 'G' in the Bowell HG model and it is called 'G12' in the Muinonen model. Details about how G and G12 vary with asteroid spectral type are found in Oszkiewicz, D.A. et al., 2012. 'Asteroid taxonomic signatures from photometric phase curves'. and Pravec, P. et al., 2012. 'Absolute magnitudes of asteroids and a revision of asteroid albedo estimates from WISE thermal observations'.
        :param sun_to_asteroid_distance_au (float): The sun to asteroid distance in au.
        :param observer_to_asteroid_distance_au (float): The observer to asteroid distance in au.


    Returns:
        float: V - H where V is the apparent magnitude of an asteroid, H is the absolute magnitude of an asteroid

    Example:
        From JPL HORIZONS:
            Ephemeris Type [change] : 	OBSERVER
            Target Body [change] : 	367789 (2011 AG5)
            Observer Location [change] : 	Geocentric [500]
            Time Span [change] : 	Start=2017-07-20, Stop=2017-07-21, Step=1 h
            Date and Time     App.Mag  helio dist. (au)  geo dist. (au)     phase angle (radians)
            2017-Jul-20 00:00 26.00    1.771234568438   2.29925674246612   0.43361309002397525

        With the HG system:
            system = 'HG'
            phase_angle_radians = 0.43361309002397525
            G_or_G_12 = 0.15 #average of S + C type
            sun_to_asteroid_distance_au = 1.771234568438
            observer_to_asteroid_distance_au = 2.29925674246612
            H_2011_AG5 =  21.8
            V_minus_H_mag = phase_function_bowell_muinonen (system, phase_angle_radians, G_or_G_12, sun_to_asteroid_distance_au, observer_to_asteroid_distance_au)
            print V_minus_H_mag
            4.19786444941
            print V_minus_H_mag + H_2011_AG5
            25.9978644494

        With the HG12 system:
            system = 'HG12'
            phase_angle_radians = 0.43361309002397525
            G_or_G_12 = 0.53 #average of S + C type
            sun_to_asteroid_distance_au = 1.771234568438
            observer_to_asteroid_distance_au = 2.29925674246612
            H_2011_AG5 =  21.8
            V_minus_H_mag = phase_function_bowell_muinonen (system, phase_angle_radians, G_or_G_12, sun_to_asteroid_distance_au, observer_to_asteroid_distance_au)
            print V_minus_H_mag
            4.13188162363
            print V_minus_H_mag + H_2011_AG5
            25.9318816236

    """

    #Cubic splines of the first and second basis functions in the HG12 system
    #The first pair of rows contain 6 phase angle bins stating from 0 degreens < phase angle < 7.5 degrees to 120 degrees < phase angle < 150 degrees for the first and second basis functions as seen in Eq. 18 of Muinonen et al. 2010.
    #The second pair of rows contains the basis function cubic spline at the location of the 6 phase angle bins in the first rows for the first and second basis functions.
    # The third pair of rows contains a cubic spline of the derivative of the basis function spline at the location of the 6 phase angle bins in the first rows for the first and second basis functions.
    basis_function_splines_1_2 = \
        [[[  1.30899694e-01,   1.30899694e-01],
        [  5.23598776e-01,   5.23598776e-01],
        [  1.04719755e+00,   1.04719755e+00],
        [  1.57079633e+00,   1.57079633e+00],
        [  2.09439510e+00,   2.09439510e+00],
        [  2.61799388e+00,   2.61799388e+00]],

       [[  7.50000000e-01,   9.25000000e-01],
        [  3.34860160e-01,   6.28841690e-01],
        [  1.34105600e-01,   3.17554950e-01],
        [  5.11047560e-02,   1.27163670e-01],
        [  2.14656870e-02,   2.23739030e-02],
        [  3.63969890e-03,   1.65056890e-04]],

       [[ -1.90985932e+00,  -5.72957795e-01],
        [ -5.54634324e-01,  -7.67053671e-01],
        [ -2.44045986e-01,  -4.56657890e-01],
        [ -9.49804381e-02,  -2.80718090e-01],
        [ -2.14114236e-02,  -1.11732569e-01],
        [ -9.13286120e-02,  -8.65731380e-08]]]

    # Cubic splines of the third basis function in the HG12 system
    # The third basis function is only used with phase angles of less than 30 degrees.
    #The first row contains phase angle increments starting from 0 degrees < phase angle < 0.3 degrees to 20 < phase angle < 30 degrees.
    #The second row contains a cublic spline of the third basis function evaluted at the phase angle bin locations for phase angle < 30.
    #The third row contains a cublic spline of derivative of the basis function evaluated at the phase angle bin locations for phase angle < 30.
    basis_function_splines_3 = \
        [[  0.00000000e+00,   5.23598776e-03,   1.74532925e-02,
          3.49065850e-02,   6.98131701e-02,   1.39626340e-01,
          2.09439510e-01,   3.49065850e-01,   5.23598776e-01],
       [  1.00000000e+00,   8.33811850e-01,   5.77354240e-01,
          4.21447720e-01,   2.31742300e-01,   1.03481780e-01,
          6.17334730e-02,   1.61070060e-02,   0.00000000e+00],
       [ -1.06300970e-01,  -4.11804395e+01,  -1.03669147e+01,
         -7.57846149e+00,  -3.69609502e+00,  -7.86056516e-01,
         -4.65270118e-01,  -2.04595448e-01,   0.00000000e+00]]


    if system == 'HG12':
        basis_function = [0.0,0.0,0.0]#initiate basis function array
        #defining G1 and G2 phase slope parameters
        if G_or_G_12 < 0.2:#see Eq 28 in Muinonen et al. 2010
            G1 = (0.7527*G_or_G_12) + 0.06164
            G2 = (-0.9612*G_or_G_12) + 0.6270
        if G_or_G_12 >= 0.2:
            G1 = (0.9529*G_or_G_12) + 0.02162
            G2 = (-0.6125*G_or_G_12) + 0.5572

        if phase_angle_radians < 0.1308996938995747: #for phase angles < 7.5 degrees, basis functions 1 and 2
            for i in range(0,2):
                basis_function[i] = ((phase_angle_radians - basis_function_splines_1_2[i][0][1]) * basis_function_splines_1_2[2][0][i]) + basis_function_splines_1_2[1][0][i]

        if phase_angle_radians >=  0.1308996938995747:#for phase angles >= 7.5 degrees
            for i in range(0,2):#loop over basis functions 1 and 2
                for j in range(0, len(basis_function_splines_1_2[0])-1): #loop over phase-angle bins
                    if basis_function_splines_1_2[0][j][i] <= phase_angle_radians and phase_angle_radians <= basis_function_splines_1_2[0][j+1][i]:
                        phase_angle_lo = basis_function_splines_1_2[0][j][i]
                        phase_angle_hi = basis_function_splines_1_2[0][j+1][i]
                        basis_function_lo = basis_function_splines_1_2[1][j][i]
                        basis_function_hi = basis_function_splines_1_2[1][j+1][i]
                        phase_angle_derivative_lo = basis_function_splines_1_2[2][j][i]
                        phase_angle_derivative_hi = basis_function_splines_1_2[2][j+1][i]
                        Delta_phase_angle = phase_angle_hi - phase_angle_lo
                        Delta_basis_function = basis_function_hi - basis_function_lo
                        a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
                        b = (-1.0*phase_angle_derivative_hi*Delta_phase_angle) + Delta_basis_function
                        t = (phase_angle_radians - phase_angle_lo) / Delta_phase_angle
                        basis_function[i] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))

        if phase_angle_radians > 2.6179938779914944:# for phase angles > 150.0 degrees
            #2
            phase_angle_lo = basis_function_splines_1_2[0][5][0]
            phase_angle_hi = phase_angle_radians
            basis_function_lo = basis_function_splines_1_2[1][4][0]
            basis_function_hi = basis_function_splines_1_2[1][5][0]
            phase_angle_derivative_lo = basis_function_splines_1_2[2][4][0]
            phase_angle_derivative_hi = basis_function_splines_1_2[2][5][0]
            Delta_phase_angle = phase_angle_hi - phase_angle_lo
            Delta_basis_function = basis_function_hi - basis_function_lo
            a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
            b = (-1.0*phase_angle_derivative_hi*Delta_phase_angle) + Delta_basis_function
            t = (phase_angle_radians - phase_angle_lo) / Delta_phase_angle
            basis_function[0] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))
            #2
            phase_angle_lo = basis_function_splines_1_2[0][5][1]
            phase_angle_hi = phase_angle_radians
            basis_function_lo = basis_function_splines_1_2[1][4][1]
            basis_function_hi = basis_function_splines_1_2[1][5][1]
            phase_angle_derivative_lo = basis_function_splines_1_2[2][4][1]
            phase_angle_derivative_hi = basis_function_splines_1_2[2][5][1]
            Delta_phase_angle = phase_angle_hi - phase_angle_lo
            Delta_basis_function = basis_function_hi - basis_function_lo
            a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
            b = (-1.0*phase_angle_derivative_hi*Delta_phase_angle) + Delta_basis_function
            t = (phase_angle_radians - phase_angle_lo) / Delta_phase_angle
            basis_function[1] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))

        #basis function 3
        if phase_angle_radians > 0.52359877559829882:#for phase angles > 30 degrees
            basis_function[2] = 0.0
        if phase_angle_radians <= 0.52359877559829882:#for phase angles <= 30 degrees
            for i in range(0, len(basis_function_splines_3[0]) -1):#loop through the phase angle bins for phase angle <= 30 degrees
                if basis_function_splines_3[0][i] <= phase_angle_radians and phase_angle_radians <= basis_function_splines_3[0][i+1]:
                    phase_angle_lo = basis_function_splines_3[0][i]
                    phase_angle_hi = basis_function_splines_3[0][i+1]
                    basis_function_lo = basis_function_splines_3[1][i]
                    basis_function_hi = basis_function_splines_3[1][i+1]
                    phase_angle_derivative_lo = basis_function_splines_3[2][i]
                    phase_angle_derivative_hi = basis_function_splines_3[2][i+1]
                    Delta_phase_angle = phase_angle_hi - phase_angle_lo
                    Delta_basis_function = basis_function_hi - basis_function_lo
                    a = (phase_angle_derivative_lo*Delta_phase_angle) - Delta_basis_function
                    b = ((-1.0 *phase_angle_derivative_hi)*Delta_phase_angle) + Delta_basis_function
                    t = (phase_angle_radians - phase_angle_lo) / Delta_phase_angle
                    basis_function[2] = ((1-t)*basis_function_lo) + (t*basis_function_hi) + (t*(1-t)*((a*(1-t))+(b*t)))

        V_minus_H_mag = -2.5 * math.log10( (G1 * basis_function[0]) + (G2 * basis_function[1]) + ((1.0 - G1 - G2) * basis_function[2])) + 5.0*math.log10(sun_to_asteroid_distance_au*observer_to_asteroid_distance_au)

    if system == 'HG':
        #basis function variables
        #variables for basis function 1
        A_1 = 3.332
        B_1 = 0.631
        C_1 = 0.986
        #variables for basis function 2
        A_2 = 1.862
        B_2 = 1.218
        C_2 = 0.238

        basis_function_1 =  (math.exp(-90.56*math.tan(0.5*phase_angle_radians)**2) * (1.0 - (C_1*math.sin(phase_angle_radians)) / ( 0.119 + (1.341 * math.sin(phase_angle_radians)) - (0.754*(math.sin(phase_angle_radians)**2)) ) ) ) + ( (1 - math.exp(-90.56*math.tan(0.5*phase_angle_radians)**2)) * ( math.exp(-1.0 * A_1 *math.tan(0.5*phase_angle_radians)**B_1) ) )
        basis_function_2 =  (math.exp(-90.56*math.tan(0.5*phase_angle_radians)**2) * (1.0 - (C_2*math.sin(phase_angle_radians)) / ( 0.119 + (1.341 * math.sin(phase_angle_radians)) - (0.754*(math.sin(phase_angle_radians)**2)) ) ) ) + ( (1 - math.exp(-90.56*math.tan(0.5*phase_angle_radians)**2)) * ( math.exp(-1.0 * A_2 *math.tan(0.5*phase_angle_radians)**B_2) ) )

        V_minus_H_mag = -2.5*math.log10( ((1.0 - G_or_G_12)*basis_function_1) + (G_or_G_12*basis_function_2) ) + 5.0*math.log10(sun_to_asteroid_distance_au*observer_to_asteroid_distance_au)

    return V_minus_H_mag

def relative_tracking_rate_for_streak_length_arcsec_p_minute(streak_length_arcsec,exposure_time_s):
    return (streak_length_arcsec / (exposure_time_s/minutes_to_seconds))

def sign_to_plus_minus(number):
    sign = np.sign(number)
    if sign >= 0:
        sign_str = "+"
    if sign < 0:
        sign_str = "-"
    return sign_str

def split_mpc_utc_hours(mpc_hour_minute_ut):
    index_0 = np.ones(len(mpc_hour_minute_ut))*0
    index_0 = index_0.astype('int')
    index_1 = np.ones(len(mpc_hour_minute_ut))*2
    index_1 = index_1.astype('int')
    index_2 = np.ones(len(mpc_hour_minute_ut))*2
    index_2 = index_2.astype('int')
    index_3 = np.ones(len(mpc_hour_minute_ut))*4
    index_3 = index_3.astype('int')
    hour = np.array(map(split_string,mpc_hour_minute_ut,index_0,index_1))
    hour = hour.astype('float')
    minute = np.array(map(split_string,mpc_hour_minute_ut,index_2,index_3))
    minute = minute.astype('float')
    return hour, minute

def streak_length_arcsec_given_SNRs_SNRt(seeing_FWHM_arcsec, SNR_stationary, SNR_trailed):
    streak_length_arcsec = seeing_FWHM_arcsec * ((SNR_stationary/SNR_trailed) - 1.0)
    return streak_length_arcsec

def time_coverage_by_seeing_limited_FWHM_arc_sec_and_streak_length_arcsec(seeing_FWHM_arcsec, streak_length_arcsec, exposure_time_s):
    num_seeing_FWHM_in_streak_length_arcsec = streak_length_arcsec/seeing_FWHM_arcsec
    time_coverage_per_seeing_FWHM_streak_sec = exposure_time_s / num_seeing_FWHM_in_streak_length_arcsec
    return time_coverage_per_seeing_FWHM_streak_sec

def time_of_perhelion_epoch_period_to_mean_anomaly_at_epoch_degrees(time_perihelion_mjd, epoch_mjd, period_seconds):
    #Mean Anomaly         = n * (t - T)  =  (t - T) * 360_deg / P

    period_days = period_seconds * seconds_to_days
    M_radians = (epoch_mjd - time_perihelion_mjd) * ((2 * np.pi)/period_days)
    M_degrees = np.degrees(M_radians)
    if M_degrees <0:
        M_degrees += 360.
    return M_degrees

def total_track_length(streak_length_arcsec, exposure_time_s, readout_time_s, num_exposures, num_filters, cfht_filter_change_time_s):
    tracking_rate_arcsec_p_s = streak_length_arcsec / exposure_time_s
    total_track_time_s = (num_filters * ((num_exposures * exposure_time_s) + (num_exposures * readout_time_s))) + ((num_filters-1) * cfht_filter_change_time_s) - readout_time_s #since we do not count the last readout time, nor do we count the first filter for the purposes of a filter change time
    total_track_length_arcsec = tracking_rate_arcsec_p_s * total_track_time_s
    return total_track_time_s, total_track_length_arcsec

def visir_exp_time_s(SNR,flux_Jy):
    return (SNR/(visir_constant * flux_Jy))**2

