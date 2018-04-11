#Bryce Bolin

import sys
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/APO_observing')
from apo_observing_functions import *

'''
mjd for 2018 04 10
58218.0
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 4231 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 4231
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 2322 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 2322
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 933 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 933
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 1700 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 1700
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 2776 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 2776
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 2536 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 2536
python tracking_command_solar_system_APO.py -mpcorbf /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT -ooloc /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb -an 2792 -stetss 58219.0 58220.0 0.002 -rot -67.5 > 2792
'''

parser = argparse.ArgumentParser()
parser.add_argument("-mpcorbf", "--mpc_orb_file", help="location of lightcurve shape models, e.g., /Users/bolin/Thermal/asteroid_lists/MPCORB.DAT", nargs='*')
parser.add_argument("-ooloc", "--oorb_location", help="location of oorb, e.g., /Users/bolin/NEO/OpenOrb/oorb-master/main/oorb", nargs='*')
parser.add_argument("-an", "--asteroid_name", help="numbered name of asteriod", nargs='*')
parser.add_argument("-stetss", "--start_time_end_time_step_size", help="start_time,end_time in MJD and step size in days 57303 57335 0.0025", nargs='*')
parser.add_argument("-rot", "--rotation_deg", help="rotation of the slit, should be 90-parallactic angle", nargs='*')
args = parser.parse_args()

mpc_orb_file = str(args.mpc_orb_file[0])
oorb_location = str(args.oorb_location[0])
asteroid_name = str(args.asteroid_name[0])
start_time_mjd, end_time_mjd, step_size_days = string_seperated_to_array_spaces(args.start_time_end_time_step_size,'float')
rot_angle_str = str(args.rotation_deg[0])

orbit_line =  grep_asteroid_from_MPCORBDAT_to_KEP_DES_format(asteroid_name,mpc_orb_file)
id_generator_orb = id_generator()
orbit_file_name = '''orbit_''' + id_generator_orb + '''.des'''
echo_orbit_line_to_des = '''echo "''' + orbit_line + '''"''' + ''' > ''' + orbit_file_name
os.system(echo_orbit_line_to_des)

MJD_RA_deg_dec_deg_dRA_dt_cos_dec_deg_per_day_ddec_dt_deg_per_day = compute_oorb_sky_coords(oorb_location,orbit_file_name, start_time_mjd, end_time_mjd, step_size_days, orbit_file_name)
os.system('rm *' + id_generator_orb + '*')


MJD = MJD_RA_deg_dec_deg_dRA_dt_cos_dec_deg_per_day_ddec_dt_deg_per_day[:,0]
RA_deg = MJD_RA_deg_dec_deg_dRA_dt_cos_dec_deg_per_day_ddec_dt_deg_per_day[:,1]
dec_deg = MJD_RA_deg_dec_deg_dRA_dt_cos_dec_deg_per_day_ddec_dt_deg_per_day[:,2]
dRA_dt_cos_dec_deg_sec = MJD_RA_deg_dec_deg_dRA_dt_cos_dec_deg_per_day_ddec_dt_deg_per_day[:,3]/(24.*3600.)
ddec_dt_deg_per_deg_sec = MJD_RA_deg_dec_deg_dRA_dt_cos_dec_deg_per_day_ddec_dt_deg_per_day[:,4]/(24.*3600.)

t = Time(MJD, format='mjd', scale='utc')
time_stamp = t.iso

for i in range(0, len(MJD)):
    print ' %5s %2s %5s %7.5f, %7.5f, %7.8f, %7.8f %7s %7s '%(time_stamp[i],'tcc','track',RA_deg[i], dec_deg[i], dRA_dt_cos_dec_deg_sec[i], ddec_dt_deg_per_deg_sec[i], 'Fk5=2000.0', '/Rottype=Object /Rotation='+rot_angle_str)
