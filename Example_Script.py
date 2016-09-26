# Diamond I16 Analysis Script
# [Experiment Name]
# [Script Description]
# 
# By [User]
# 18/08/2016

import sys,os
import matplotlib.pyplot as plt # Plotting
from mpl_toolkits.mplot3d import Axes3D # 3D plotting

# Load Py16progs
sys.path.insert(0,'/dls_sw/i16/software/python/Py16/Py16progs.py') # location of Py16progs
import Py16progs as dp

# Current Directory
cf=os.path.dirname(__file__)

# Directory to load data from
dp.filedir = '/dls/i16/data/2015/cm12169-2' 

# Directory to save files to
dp.savedir='/home/i16user/Desktop' 

# Experiment Parameters
dp.exp_ring_current = 300.0 # Average ring current for normalisation
dp.exp_monitor = 800.0 # Average monitor current for normalisation
dp.normby = 'rc' # ring current ('rc'), monitor ('ic1') or none ('none')
dp.pil_centre = [110, 242] 
dp.peakregion=[7,153,186,332] # Search for peaks within this area of the detector [min_y,min_x,max_y,max_x]
dp.exp_title = ''

# Scan numbers
scans = range(512404,512664,5) + [512665,512666]

# Multi-scan plots
help(dp.plotscan) # See the function documentation!
dp.plotscan(scans)
dp.plotscans3D(scans)
dp.plotscansSURF(scans)

# Automatic Peak Fitting & Integration
fit,err = dp.fit_scans(scans,vary='roi1_sum',depvar='Ta',peaktest=60,fit_type='pVoight',saveFIT=None,save=None)

# Load fitted data:
#fit,err = dp.load_fits(scans,depvar='Ta',fit_type='pVoight')

# Manual analysis of data
temp,value = [],[]
for scn in scans:
	x,y,dy,varx,vary,ttl,d = dp.getdata(scn,vary='roi1_sum') # read normalised data
	temp += [d.metadata.Ta]
	value += [sum(y)]

# Basic plotting with Matplotlib
plt.figure()
plt.plot(temp,value)
dp.labels('Title','temp','value')
