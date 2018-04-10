#Bryce Bolin
'''

Generates commands for tracking of solar system objects with APO

e.g.

tcc track 136.72666, 17.43467, -0.00000146, 0.00000046 Fk5=2000.0

#astropy convert utc to mjd and back 
t = Time('2010-01-01 00:00:00', scale='utc')
t.mjd
t = Time(55197.0, format='mjd', scale='utc')
t.iso

get coords from JPL Horizons


The BIG, IMPORTANT difference is that JPL gives dRA and dDec in arcsec/hour, whereas our software wants DEG?/SEC!!!!!!!!!!!!!!!!. That means YOU HAVE TO DIVIDE THE SPEED BY 12960000

make out put to have timesteps
2008-Apr-24 04:00:00 154.3043103 12.5870917 -2.51 0.58 tcc track 154.3043103, 12.5870917, -0.000000194, 0.000000047 Fk5=2000.0 /Rotation=0.0 /Rottype=Object

dRA has cos dec multiplied
dRA*cosD

ipython -i -- apo_observing_functions.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 306 -stetss 44113.0 44114 0.002

'''



parser = argparse.ArgumentParser()
parser.add_argument("-mpcorbf", "--mpc_orb_file", help="location of lightcurve shape models, e.g., /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT", nargs='*')
parser.add_argument("-ooloc", "--oorb_location", help="location of oorb, e.g., /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb", nargs='*')
parser.add_argument("-an", "--asteroid_name", help="numbered name of asteriod", nargs='*')
parser.add_argument("-stetss", "--start_time_end_time_step_size", help="start_time,end_time in MJD and step size in days 57303 57335 0.0025", nargs='*')
args = parser.parse_args()

mpc_orb_file = str(args.mpc_orb_file[0])
oorb_location = str(args.oorb_location[0])
asteroid_name = str(args.asteroid_name[0])
start_time_mjd, end_time_mjd, step_size_days = string_seperated_to_array_spaces(args.start_time_end_time_step_size,'float')

orbit_line =  grep_asteroid_from_MPCORBDAT_to_KEP_DES_format(asteroid_name,mpc_orb_file)
id_generator_orb = id_generator()
orbit_file_name = '''orbit_''' + id_generator_orb + '''.des'''
echo_orbit_line_to_des = '''echo "''' + orbit_line + '''"''' + ''' > ''' + orbit_file_name
os.system(echo_orbit_line_to_des)
#os.system('rm *' + id_generator_orb + '*')

JD_light_time_corrected_m_astro_hel_x_au_astro_hel_y_au_astro_hel_z_au_astro_toppo_x_au_astro_toppo_y_au_astro_toppo_z = compute_oorb_astroidcentric_helio_and_toppo_vectors_with_JD(oorb_location,orbit_file_name, start_time_mjd, end_time_mjd, step_size_days, orbit_file_name)

run_lightcurve_code(JD_light_time_corrected_m_astro_hel_x_au_astro_hel_y_au_astro_hel_z_au_astro_toppo_x_au_astro_toppo_y_au_astro_toppo_z, asteroid_name, shape_model_directory, lightcurve_code, start_time_mjd, end_time_mjd)
os.system('rm *' + id_generator_orb + '*')