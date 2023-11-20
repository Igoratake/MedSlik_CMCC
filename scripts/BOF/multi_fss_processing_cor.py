import glob
import os
import re
import os
from datetime import  *
import shutil
import pandas as pd
from BayesianOptimizationLib_.bayes_opt.path_ import Path
import numpy as np
import argparse

dt = 12

parser = argparse.ArgumentParser()
parser.add_argument('--outdataframe', '-o', help='Out dataframe name')
parser.add_argument('--parentscript', '-p', help='Parent script name')
args = parser.parse_args()

''' Punto in cui agire per restringere le osservazioni a pochi giorni '''
list_of_obs = glob.glob(Path.SAT_FOLDER + Path.DAYS_GROUP)

p = open(Path.PARTICLE_FILE)
lines = p.readlines()
n_particles = lines[0].split('\n')[0]
p.close()
	
f = open(Path.MEDSLIKII_OUT_DIR + '/config1.txt')
lines = f.readlines()
f.close()

# get sim_length
ss = lines[3]
rr = re.search('=(.*)\n',ss)
sim_length = float(lines[1][11:15])

# get simulation date
# day
ss=lines[5]
rr =re.search('=(.*)\n',ss)
#print(rr.group(1))
dd=(rr.group(1))
# month
ss=lines[6]
rr =re.search('=(.*)\n',ss)
mm=(rr.group(1))
# year
ss=lines[7]
rr =re.search('=(.*)\n',ss)
yy=(rr.group(1))
# hour
ss=lines[8]
rr =re.search('=(.*)\n',ss)
hh=(rr.group(1))

sim_date = date(2000+int(yy),int(mm),int(dd)).toordinal()
sim_date = sim_date + float(hh)/24.
print(sim_date)

days = float(sim_length)/24
print(days)

fss_df = pd.DataFrame(columns=['folder', 'fss'])

obs = 0
for slick_folder in list_of_obs:
	null, slick_id=os.path.split(slick_folder)
	slick_date=date(int(slick_id[0:4]),int(slick_id[4:6]),int(slick_id[6:8])).toordinal()
	slick_date = slick_date + float(slick_id[9:11])/24.

	if (slick_date-sim_date) < days:
     
		if(args.parentscript == 'MDKII_single_run'):
			out_folder = Path.OUT_FOLDER_SIN_RUN
			detection_folder = Path.DETECTION_FOLDER_SIN_RUN + str(obs)
		else:
			out_folder = Path.OUT_FOLDER
			detection_folder = Path.DETECTION_FOLDER + str(obs)

# per RS
		if not os.path.exists(out_folder):
			os.mkdir(out_folder)
			print(out_folder)

		if os.path.exists(detection_folder):
			shutil.rmtree(detection_folder)
			os.mkdir(detection_folder)
		else:
			os.mkdir(detection_folder)
			print(detection_folder)

		pysteps_string = 'python ' + Path.SINGLE_DET_FSS
		args_string = '-f ' + Path.MEDSLIKII_OUT_DIR + ' -o ' + slick_folder + ' -u ' + detection_folder + ' 2> /dev/null'
		os.system(pysteps_string + args_string)

		files = glob.glob(detection_folder + '/fss_syria_*.txt')

		for file in files:
      
			f = pd.read_csv(file, sep = ' ', header = None)
   
			row = pd.DataFrame([{"folder" : out_folder, "fss" : round(f.iloc[0,1], 4)}])
			fss_df = pd.concat([fss_df, row])

		obs += 1

fss_df.to_csv(out_folder + '/' + args.outdataframe)
