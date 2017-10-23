# Diamond I16 Analysis Script
# [Experiment Name]
# [Script Description]
# 
# By [User]
# 23/10/2017


"--------------Import Stuff-------------"
import sys,os
import numpy as np
import matplotlib.pyplot as plt # Plotting
from mpl_toolkits.mplot3d import Axes3D # 3D plotting

# Load Py16progs
sys.path.insert(0,'/dls_sw/i16/software/python/Py16/Py16progs.py') # location of Py16progs
import Py16progs as dp


"-----Set Directories and parameters----"
# Current Directory of this file - if you need it
cf=os.path.dirname(__file__)

# Directory to load data from ***CHANGE THIS***
dp.filedir = '/dls/i16/data/2015/cm12169-2' 

# Directory to save files to ***CHANGE THIS***
dp.savedir='/home/i16user/Desktop' 

# Update default save location for exported plots
plt.rcParams["savefig.directory"] = dp.savedir

# Experiment Parameters (feel free to ignore or remove unless you want to change them)
dp.exp_ring_current = 300.0 # Average ring current for normalisation
dp.exp_monitor = 800.0 # Average monitor current for normalisation
dp.normby = 'rc' # ring current ('rc'), monitor ('ic1') or none ('none')
dp.pil_centre = [110, 242] # Find the current value in /dls_sw/i16/software/gda/config/scripts/localStation.py (search for "ci=")
dp.peakregion=[7,153,186,332] # 'nroi_peak' will search for peaks within this area of the detector [min_y,min_x,max_y,max_x]
dp.exp_title = '' # exp_title will appear in plot titles

"---------------------------------------"
"------------Analysis Stuff-------------"
"---------------------------------------"
### Single Scan ###
d = dp.readscan(123456)
eta_array = d.eta
initial_eta_value = d.metadata.eta

x,y,dy,varx,vary,ttl,d = dp.getdata(scn,vary='roi1_sum')

### Multiple Scans ###
scans = range(512404,512664,5) + [512665,512666]

# Multi-scan plots
help(dp.plotscan) # See the function documentation!
dp.plotscan(scans)
dp.plotscans3D(scans)
dp.plotscansSURF(scans)

# Automatic Peak Fitting & Integration
fit,err = dp.fit_scans(scans,vary='roi1_sum',depvar='Ta',peaktest=1,fit_type='pVoight',plot=['all','int'],saveFIT=None,savePLOT=None)

# Load fitted data:
#fit,err = dp.load_fits(scans,depvar='Ta',fit_type='pVoight')

### Manual analysis of data ###
temp,value = [],[]
for scn in scans:
	x,y,dy,varx,vary,ttl,d = dp.getdata(scn,vary='roi1_sum') # read normalised data
	# x = array of scanned values (e.g. eta)
	# y = array of recoreded values, normalised by time, attenuation + ring current (e.g. 'roi1_sum')
	# dy = errors on y
	# varx/y = string label for x/y
	# ttl = automatically generated scan title
	# d = dataholder as from d = dp.readscan(scn)
	temp += [d.metadata.Ta]
	value += [sum(y)]

# Basic plotting with Matplotlib
plt.figure()
plt.plot(temp,value,'r-o')
dp.labels('Title','temp','value')
