# -*- coding: utf-8 -*-
"""
Module: I16 Data analysis programs "Py16Progs.py"

By Dan Porter, PhD
Diamond
2016

Usage: 
******In Script********
include the following lines at the top of your script

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/Py16') # location of Py16Progs
import Py16Progs as p16
p16.filedir = '/dls/i16/data/2015/mt0000-1' # your experiment folder
p16.savedir = '/dls/i16/data/2015/mt0000-1/processing' # where to save output

# Use functions:
p16.checkexp()
num = p16.latest()
d= p16.readscan(num)
p16.plotscan(num)
p16.plotpil(num)

******In Console*******
Run the file in the current console:
    >> cd /dls_sw/i16/software/python/Py16/
    >> ipython -i --matplotlib tk
In the python console:
    >> import Py16progs as p16
Change filedir: 
    >> p16.filedir = '/dls/i16/data/2015/cm12169-2/CsFeSe' # your experiment folder

# Use functions:
    >> d= p16.readscan(12345) # generate data structure for scan number 12345
    >> d.eta # returns scan data for eta psudo-device
    >> d.metadata.Energy # returns metadata for scan (at start of scan) 
    >> p16.checkexp() # experiment dates + scans numbers
    >> num = p16.latest() # Latest scan number

**********************
Some Useful Functions:
    d = readscan(num) 
    x,y,dy,varx,vary,ttl,d = getdata(num/d,'varx','vary',save=None)
    x,y,z,varx,vary,varz,ttl = joindata([nums],'varx','vary','varz')
    vol = getvol(num)
    ROI_sum,ROI_maxval,ROI_bkg = pilroi(num,ROIcen,ROIsize,findpeak,peakregion)
    A = getmeta([nums],'Field')
    
    num = latest()
    checkexp()
    checkscan(num1=None,num2=None)
    checklog(time=None,mins=2,cmd=False)
    prend(start,end)
    polflip(sigsig,sigpi,fit='Gauss',output=False)
    metaprint(d1,d2=None)
    savedata(num/d,'varx','vary')
    
    plotscan(num=None,vary='',fit=None,save=None)
    plotpil(num,cax=None,imnum=None)
    plotscans3D(scans,depvar='Ta',vary='',save=None)
    
    fit,err = fit_scans(scans,depvar='Ta',vary='',save=None)
    
    out,err = peakfit(x,y,dy,type='dVoight')
    
    HHH,KKK,LLL = pixel2hkl(num)
    XXX,YYY,ZZZ = pixel2xyz(num)
    TTH,INT = pixel2tth(num)
    
    ROIcen_ij,ROIcen_frame = pilpeak(Y,test = 1,disp=False)
    ispeak(Y,test = 1,disp=False)
    peak = peakfind(Y,cutoff=0.01)*
    
    labels(title,xlabel,ylabel,zlabel)
    vals = frange(start,stop=None,step=1)
    str = stfm(val,err)
    

Version 4.8.7
Last updated: 31/07/22

Version History:
07/02/16 0.9    Program created from DansI16progs.py V3.0
18/02/16 1.0    Functions renamed, readscan updated for errors, simplifications to outputs
03/03/16 1.1    Tested on Beamline, removed errors, cleaned up
08/03/16 1.2    Refinements after a week of beamline use
15/04/16 1.3    New masks method in peak_fit
05/05/16 1.4    Fixed bug that did integer normalisation of APD/counttime
10/05/16 1.5    Added print buffer function
19/05/16 1.6    Fixed 3D plotting bugs, new outputs on fit_scans, catch psi='unavailable' bug (psi=0)
12/07/16 1.7    Ungraded peakfit with variable backgrouns, estimate input and new default values.
10/08/16 1.8    Added logplot and diffplot options to plotscan
08/09/16 1.9    Generalised getvol for any detector by pre-loading the first image
21/09/16 2.0    Removed dnp.io.load. New title includes folder and new functions for lists of scan numbers
07/10/16 2.1    Ordered keys in dataloader, some minor fixes
17/10/16 2.2    Improved ROI functionality - _sfm to remove first frame
14/12/16 2.3    Improved multiscan titles, funtions to write + load previous experiment directories for fast access
20/12/16 2.4    Added checkscans and auto_varx/vary functions, added dont_normalise for none-normalisable values
02/02/17 2.5    Added checkhkl and checkwl. Bugs fixed when using DifCalc (psi+azih not recorded)
22/02/17 2.6    Added dead pixel mask for pilatus for experiments in Feb 2017
13/06/17 2.7    Added getmeta() function
24/07/17 2.8    Added d.metadata.hkl_str, fixing incorrect position bug. Various other fixes.
01/08/17 2.9    Added check for large pilatus arrays.
02/10/17 2.9    Added various misc functions, other bugs fixed
06/10/17 3.0    Added pixel2hkl functions, incluing plots, other bug fixes
17/10/17 3.1    Added choice of temperature sensor, d.metadata.temperature
23/10/17 3.2    Added new plotting features, better multi-dimensional fits
17/11/17 3.3    Added normalisation to pixel2tth, added default bpm behaviour
06/12/17 3.4    Added new plots for pixel2tth
23/01/18 3.5    Added d.counter to scans, made pixel2hkl work for hkl scans
23/02/18 3.6    Added bin_pixel_hkl_cut, new fit methods
09/03/18 3.6    Added, scanpol, updated to correct for python3.6 test errors
18/03/18 3.7    Added exp_parameters_save and _load
01/05/18 3.8    Added sort to joinscans, updated for new PA + pilatus3, added scanimage
21/05/18 3.8    Changed automatic titles to include psi angle and reference
01/08/18 3.9    Removed psutil, added getRAM, various updates and fixes, added plotmeta and pilmaxval, corrected joinscan save
19/10/18 3.9    Corrected type input of getvol.
26/11/18 4.0    Output of checkscan, checklog now str
14/12/18 4.1    Update to simpfit, giving better estimates of peak position, changed default tmpdir
21/02/19 4.2    Update to pcolorplot, removed some bugs
18/03/19 4.2    Update to plotscans, correcting normalisation for multi-vary
15/04/19 4.3    Update to readscan, additional checks on metadata added
15/05/19 4.4    Updated labels function, added newplot, multiplot, sliderplot, sliderplot2D
02/10/19 4.5    Added plotqbpm, added defaults for phase plate scans
23/10/19 4.6    Fixed exec compatibility, now python3 compatible, readnexus added using nexusformat or h5py
29/11/19 4.7    Changed fit_scans save name to include vary, added nexus_rsremap and nexus_plot_rsremap
21/02/20 4.7    Added findscans
29/02/20 4.7    Added output to plotqbpm
27/05/20 4.7    Added licence
01/07/20 4.7.5	Updated plotqbpm to correct error
11/02/21 4.8    Added colormap options
23/04/21 4.8    Added detector to findscans
24/05/21 4.8.1  Added file input to readscan
29/09/21 4.8.2  Added fig_dpi parameter
01/10/21 4.8.3  Corrected error from os.path.isfile(d)
21/01/22 4.8.4  Corrected error on dat files with "=" in ubMeta
26/04/22 4.8.5  Corrected error in CheckScan for python3
03/05/22 4.8.6  Corrected error for merlinroi1 in getdata
30/07/22 4.8.7  Corrected error in polflip plotting, remove plt.show from plotscan

###FEEDBACK### Please submit your bug reports, feature requests or queries to: dan.porter@diamond.ac.uk

@author: Dan Porter
I16, Diamond Light Source
2016

-----------------------------------------------------------------------------
   Copyright 2020 Diamond Light Source Ltd.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Dr Daniel G Porter, dan.porter@diamond.ac.uk
www.diamond.ac.uk
Diamond Light Source, Chilton, Didcot, Oxon, OX11 0DE, U.K.
"""

#Ideas for the future:
# - Pilatus peak finding (2D fitting?) as separate function to dpilroi
# - dataloader contain functions (more pythony...)
# - Refactor fitting routines in new Py16fit module
# - Create dead_pixel files for detectors, load based on which detector used
# - Better naming convention of fit files, including HKL
# - save parameters per experiment


from __future__ import print_function
import sys,os
import glob # find files
import re # regular expressions
import datetime # Dates and times
import time # For timing
import tempfile # Find system temp directory
import numpy as np
import h5py # read hdf5 files such as nexus (.nxs)
#import scisoftpy as dnp # Make sure this is in your python path
import matplotlib.pyplot as plt # Plotting
import matplotlib.ticker as mtick # formatting of tick labels
from matplotlib.colors import LogNorm # logarithmic colormaps
from mpl_toolkits.mplot3d import Axes3D # 3D plotting
from scipy.optimize import curve_fit # Peak fitting
#from scipy import misc # read pilatus images
from imageio import imread # read Tiff images
from scipy.signal import convolve
from itertools import product
from collections import OrderedDict

# read hdf5 files such as nexus (.nxs)
try:
    from nexusformat.nexus import nxload
except ImportError:
    print('nexusformat not available, nexus files will be read using h5py')
    from h5py import File as nxload 


###########################################################################
#############################PARAMETERS####################################
###########################################################################
"-----------------------Default Experiment Directory----------------------"
# Variable filedir is called from the namespace and can be changed at any 
# time,even once the module has been imported, if you change it in the current namespace

filedir = '/dls/i16/data/2021' 
savedir = '/home/i16user/Desktop'

#tmpdir = tempfile.gettempdir() # I16 user accounts don't have access
tmpdir = os.path.expanduser('~')

"-----------------------------Data file format----------------------------"
datfile_format = '%i.dat'
nxsfile_format = '%i.nxs'

"-----------------------Error Estimation Parameters-----------------------"
error_func = lambda x: np.sqrt(np.abs(x)+0.1) # Define how the error on each intensity is estimated
#error_func = rolling_fun
# error_func = lambda x: 0*x + 1 # Switch errors off

"-------------------------Normalisation Parameters------------------------"
exp_ring_current = 300.0 # Standard ring current for current experiment for normalisation
exp_monitor = 800.0 # Standard ic1monitor value for current experiment for normalisation
normby = 'rc' # Incident beam normalisation option: 'rc','ic1' or 'none'
dont_normalise = ['Ta','Tb','Tc','Td','ic1monitor','rc','Cmes','Vmes','Rmes']
detectors = ['APD','sum','maxval'] # error bars are only calculated for counting detectors
default_sensor = 'Ta' # d.metadata.temperature will read this sensor as the standard

"----------------------------Pilatus Parameters---------------------------"
pil_centre = [106,238] # Centre of pilatus images [y,x], find the current value in /dls_sw/i16/software/gda/config/scripts/localStation.py (search for "ci=")
hot_pixel = 2**20-100 # Pixels on the pilatus greater than this are set to the intensity chosen by dead_pixel_func
peakregion=[7,153,186,332] # Search for peaks within this area of the detector [min_y,min_x,max_y,max_x]
pilpara=[119.536,1904.17,44.4698,0.106948,-0.738038,412.19,-0.175,-0.175] # pilatus position parameters for pixel2hkl
dead_pixel_func = np.median # Define how to choose the replaced intensity for hot/broken pixels 
pil_max_size = 1e8 # Maximum volume size that will be loaded. 1e8~1679*1475*41 -> 41 points of pil2M, ~ 800MB
# pilatus_dead_pixels = np.array([[101, 244],
#        [102, 242],
#        [102, 243],
#        [102, 244],
#        [103, 242],
#        [103, 243],
#        [103, 244],
#        [104, 241],
#        [104, 242],
#        [104, 243],
#        [104, 244],
#        [105, 240],
#        [105, 242],
#        [105, 243],
#        [105, 244],
#        [106, 240],
#        [107, 242],
#        [107, 243]]) # Dead pixels at centre of pilatus ***Feb 2017***
pilatus_dead_pixels = np.zeros(shape=(0,2),dtype=int) # Blank

"----------------------------Plotting Parameters--------------------------"
plot_colors = ['b','g','m','c','y','k','r'] # change the default colours of plotscan 
exp_title = '' # Experiment title
default_colormap = 'viridis'  # maplotlib colormap
fig_dpi = 80
#plt.xkcd() # cool plots!


"-------------------------------Data Defaults-----------------------------"
# These parameters will always be defined in the metadata
# Structure meta_defaults['name']=[['list of possible names'],default value if not available]
meta_defaults = {}
meta_defaults['Energy'] = [['Energy','Energy2','en','pgm_energy'],0]
meta_defaults['stokes'] = [['stokes','stoke'],0]
meta_defaults['azih'] = [['azih'],0]
meta_defaults['azik'] = [['azik'],0]
meta_defaults['azil'] = [['azil'],0]
meta_defaults['psi'] = [['psi'],666]
meta_defaults['Ta'] = [['Ta'],300.0]
meta_defaults['h'] = [['h'],0]
meta_defaults['k'] = [['k'],0]
meta_defaults['l'] = [['l'],0]

# If these 
scan_defaults = {}
scan_defaults['t'] = [['t', 'count_time']]
scan_defaults['energy'] = [['energy','energy2']]
scan_defaults['rc'] = [['rc','energy2']]
scan_defaults['TimeSec'] = [['TimeSec','TimeFromEpoch']]


###########################################################################
##############################FUNCTIONS####################################
###########################################################################

"----------------------------Loading Functions----------------------------"


def read_dat_file(filename):
    """
    Reads #####.dat files from instrument, returns class instance containing all data
    Input: 
      filename = string filename of data file
    Output:
      d = class instance with parameters associated to scanned values in the data file, plus:
         d.metadata - class containing all metadata from datafile
         d.keys() - returns all parameter names
         d.values() - returns all parameter values
         d.items() - returns parameter (name,value) tuples
    """
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    
    # Read metadata
    meta = OrderedDict()
    lineno = 0
    for ln in lines:
        lineno += 1
        if '&END' in ln: break
        ln = ln.strip(' ,\n')
        neq = ln.count('=')
        if neq == 1:
            'e.g. cmd = "scan x 1 10 1"'
            inlines = [ln]
        elif neq > 1 and '{' not in ln:
            'e.g. SRSRUN=571664,SRSDAT=201624,SRSTIM=183757'
            ' but not ubMeta={"name": "crystal=big", ...}'
            inlines = ln.split(',')
        else:
            'e.g. <MetaDataAtStart>'
            continue
        
        for inln in inlines:
            vals = inln.split('=')
            if len(vals) != 2: continue
            try:
                meta[vals[0]] = eval( vals[1] )
            except:
                meta[vals[0]] = vals[1]
    
    
    # Read Main data
    # previous loop ended at &END, now starting on list of names
    names = lines[lineno].split()
    # Load 2D arrays of scanned values
    vals = np.loadtxt(lines[lineno+1:],ndmin=2)
    # Assign arrays to a dictionary
    main = OrderedDict()
    for name,value in zip(names,vals.T):
        main[name] = value
    
    # Convert to class instance
    d = dict2obj(main, order=names)
    d.metadata = dict2obj(meta)
    return d

def readscan(num):
    """
        Loads Diamond #.dat data file
            d=readscan(#)
            d=readscan(0) - gives current run data
            d=readscan(-n) - gives previous run data
        
        d is a dataholder storing all scan data. For example:
            d.eta = array of eta values in an eta scan
            d.APD = array of intensities when using the APD
            d.roi2_sum = array of pilatus region sums when using this option
        d.metadata gives ancillary data. For example:
            d.metadata.cmd = scan command
            d.metadata.Ta = current temperature
            d.metadata.en = current energy
        use d.keys() or d.metadata.keys() to see all values
    """
    
    if type(num) is dict2obj:
        return num
    elif type(num) is str and os.path.isfile(num):
        file = num
    else:
        if os.path.isdir(filedir) == False: 
            print( "I can't find the directory: {}".format(filedir) )
            return None

        if num < 1: 
            if latest() is None: return None
            num = latest()+num
        
        file = os.path.join(filedir, datfile_format %num)
        #print(file)
    try:
        d = read_dat_file(file)
        #d = dnp.io.load(file,warn=False) # from SciSoftPi
    except:
        print( "Scan {} doesn't exist or can't be read".format(num) )
        return None
    
    " Shorten the scan command if it is very long to help with titles etc."
    try:
        cmd = d.metadata.cmd
    except:
        print( "Scan {} doesn't contain any metadata".format(num) )
        return None
    if len(cmd) > 80:
        " Shorten long commands"
        " This little functions finds long strings and rounds them to three dp"
        def caller(match):
            return '{:1.4g}'.format(round(float(match.group(0)),2))
        d.metadata.cmd_short=re.sub('[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',caller,cmd) # find numbers in command and round them to 3dp
    else:
        d.metadata.cmd_short = cmd
    
    " Correct re-assigned values"
    keys = list(d.keys())
    metakeys = list(d.metadata.keys())
    " For some reason, parameter names change occsaionaly, add correction here"
    # Correct d.scannables
    if 'count_time' in keys: d.t = d.count_time
    if 'energy2' in keys: d.energy = d.energy2
    if 'rc' not in keys: d.rc = exp_ring_current*np.ones(len(d[keys[0]]))
    if 'TimeSec' not in keys: d.TimeSec = np.arange(0,len(d[keys[0]]))
    # Correct d.metadata
    if 'en' in metakeys: d.metadata.Energy = d.metadata.en; metakeys+=['Energy']
    if 'pgm_energy' in metakeys: d.metadata.Energy = d.metadata.pgm_energy; metakeys+=['Energy']
    if 'stokes' in metakeys: d.metadata.stoke = d.metadata.stokes
    # Missing metaddata
    if 'SRSRUN' not in metakeys: d.metadata.SRSRUN = num
    if 'a' not in metakeys: d.metadata.a = 1.
    if 'b' not in metakeys: d.metadata.b = 1.
    if 'c' not in metakeys: d.metadata.c = 1.
    if 'alpha1' not in metakeys: d.metadata.alpha1 = 90.0
    if 'alpha2' not in metakeys: d.metadata.alpha2 = 90.0
    if 'alpha3' not in metakeys: d.metadata.alpha3 = 90.0
    if 'thp' not in metakeys: d.metadata.thp = 0.0
    if 'delta_axis_offset' not in metakeys: d.metadata.delta_axis_offset = 0.0
    if 'Energy' not in metakeys: d.metadata.Energy = np.NaN
    if 's5xgap' not in metakeys: d.metadata.s5xgap=-1;d.metadata.s5ygap=-1
    if 's6xgap' not in metakeys: d.metadata.s6xgap=-1;d.metadata.s6ygap=-1
    if 'm4pitch' not in metakeys: d.metadata.m4pitch=0.0
    if 'Atten' not in metakeys: d.metadata.Atten=-1;d.metadata.Transmission=1.0
    if 'gam' not in metakeys: d.metadata.gam=0
    if 'azih' not in metakeys: d.metadata.azih=0;d.metadata.azik=0;d.metadata.azil=0
    if 'psi' not in metakeys: d.metadata.psi = 666
    if 'h' not in metakeys: d.metadata.h = 0.0
    if 'k' not in metakeys: d.metadata.k = 0.0
    if 'l' not in metakeys: d.metadata.l = 0.0
    if 'Ta' not in metakeys: d.metadata.Ta = 300.0
    
    " Add filename"
    d.metadata.filename = file
    " Add a scan position counter"
    d.counter = np.arange(0,len(d.TimeSec))
    " Add the scan time value"
    d.ScanTime = d.TimeSec - d.TimeSec[0]
    " Add a hkl string"
    d.metadata.hkl_str = scanhkl(d)
    " Add a temperature string"
    d.metadata.temperature = scantemp(d,default_sensor)
    " Add energy string"
    d.metadata.energy_str = scanenergy(d)
    " Add minimirror string"
    d.metadata.minimirrors = scanminimirrors(d)
    " Correct psi values"
    if type(d.metadata.psi) is str: d.metadata.psi = 666
    if d.metadata.psi < -1000: d.metadata.psi = 0.0
    if type(d.metadata.azih) is str: d.metadata.azih=0;d.metadata.azik=0;d.metadata.azil=0

    # Add additional parameters
    try:
        " Add azimuth string"
        d.metadata.azimuth = scanazimuth(d)
        " Add title string"
        d.metadata.title = scantitle(d)
    except:
        print('The scantitle or azimuth couldn''t be generated')
        pass
    " Array items"
    d.scannables = [x for x in keys if type(d[x]) == np.ndarray ] # only keys linked to arrays
    " Update keys"
    d.update(d.__dict__)
    d.metadata.update(d.metadata.__dict__)
    return d

def getdata(num=None, varx='', vary='', norm=True, abscor=None):
    """
    Get useful values from scans
         x,y,dy,varx,vary,ttl,d = getdata(#/d,'varx','vary')
    
        vary can also call re-intergrate pilatus regions of interest:
            vary = 'nroipeak' - search for peak and create roi around it
            vary = 'nroipeakbkg' - subtract average background defined outside roi
            vary = 'nroi[110,242,75,67]' - creates roi [ceni,cenj,widi,widj]
            vary = 'nroi[31,31]' - creates centred roi with size [widi,widj]
            vary = 'nroi' - creates roi2 like roi
            vary = 'nroibkg' - as above but with background subtraction
    """
    
    "---Handle inputs---"
    if num is None:
        num = latest()
    
    " Multiple nums given"
    if type(num)==list:
        nums = num[1:]
        num = num[0]
    else:
        nums = []
    
    "---Load data---"
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    if d is None:
        return [],[],[],'','','',d
    
    "---Get metadata---"
    keys = [x for x in d.keys() if type(d[x]) == np.ndarray ] # only keys linked to arrays
    m = d.metadata
    cmd = m.cmd # Scan command
    
    """
    # Is this a 2D scan?
    nfloats = 0
    for s in cmd.split():
        try:
            float(s)
            nfloats+=1
        except ValueError:
            pass
    
    if nfloats > 6:
        print( "This may be a 2D Scan - I'm afraid this isn't implemented yet!" )
    """
    
    ttl = scantitle(d)
    
    "---Determine scan variables from scan command---"
    if varx not in keys:
        try:
            x = np.asarray(getattr(m, varx)).reshape(-1) # allow single value x
        except:
            varx = auto_varx(d)
            x = getattr(d, varx)
    else:
        " y values from data file"
        x = getattr(d, varx)
    
    
    "***Get y values***"
    if vary.lower().startswith('nroi'):
        " y values from automatic peak search in pilatus"
        # e.g. 'nroi'           > Defaults to ROI2 in centre of pilatus
        # e.g. 'nroi[11,11]'    > in centre with size [11,11]
        # e.g. 'nroi[109,241,11,11] > centred at [109,241] with size [11,11]
        # e.g. 'nroibkg         > As 'nroi', with background subtraction
        # e.g. 'nroipeak'       > As 'nroi', with peak search 
        # e.g. 'nroipeak[11,11]'> at peak centre, with size [11,11]
        # e.g. 'nroipeakbkg[31,31]' > at peak centre with background subtraction
        # e.g. 'nroi_sfm' > subtract frame 0 from all images before ROI is generated
        
        vol = getvol(m.SRSRUN) # Load pilatus images as 3D volume
        
        " Subtract first frame"
        if 'sfm' in vary.lower():
            bkgcut = vol[:,:,0].copy()
            for n in range(len(x)):
                vol[:,:,n] = vol[:,:,n] - bkgcut
        
        " Find requested position in string"
        roival = re.findall('\d+',vary) # get numbers in vary
        if len(roival) == 2: roival = [str(pil_centre[0]),str(pil_centre[1]),roival[0],roival[1]] # Only two values given for size
        if len(roival) < 4: roival = [str(pil_centre[0]),str(pil_centre[1]),'75','67'] # No values given - default to I16 ROI2
        ROIcen = [int(roival[0]),int(roival[1])]
        ROIsize = [int(roival[2]),int(roival[3])]
        
        " Search for peak in pilatus"
        if 'peak' in vary.lower():
            ROIcen,frame = pilpeak(vol,disp=True)
        
        y,ROI_maxval,ROI_bkg = pilroi(vol,ROIcen=ROIcen,ROIsize=ROIsize)
        
        if 'maxval' in vary.lower():
            y = ROI_maxval
        
        " Handle single value x"
        if len(x) < len(y) and len(x)==1:
            y = np.asarray(np.mean(y)).reshape(-1)
        
        " Calculate errors"
        dy = error_func(y)
        
        " Background"
        if 'bkg' in vary.lower():
            y = y-ROI_bkg
            dy = np.sqrt( dy**2 + error_func(ROI_bkg)**2 )
        
    else:
        
        if vary not in keys:
            try:
                y = getattr(m,vary)*np.ones(len(x))
            except:
                vary = auto_vary(d)
                y = getattr(d,vary)
        else:
            " y values from data file"
            y = getattr(d,vary)
        
        " Handle single value x"
        if len(x) < len(y) and len(x)==1:
            y = np.asarray(np.mean(y)).reshape(-1)
        
        " Only apply errors to counting detectors"
        if 'roi' in vary or vary in detectors:
            dy = error_func(y)
        else:
            dy = 0*y
        
    
    "---Get count time---"
    # Count time can have several names
    if 't' in keys:
        cnt = 1.0*d.t
    elif 'counttime' in keys:
        cnt = 1.0*d.counttime
    elif 'count_time' in keys:
        cnt = 1.0*d.count_time
    else:
        cnt=1.0
    
    "---Get incident beam normalisation---"
    if normby.lower() in ['ic1','ic1monitor']:
        inorm = d.ic1monitor/exp_monitor
        normtxt = '('+normby+'/'+str(exp_monitor)+')'
    elif normby.lower() in ['rc','ring current']:
        inorm = d.rc/exp_ring_current
        normtxt = '('+normby+'/'+str(exp_ring_current)+')'
    else:
        inorm = np.ones(len(d.rc))
        normtxt = ''
    
    " Normalise the data by time, transmission and incident beam"
    if norm and vary not in dont_normalise:
        y = y/d.metadata.Transmission/cnt/inorm
        dy = dy/d.metadata.Transmission/cnt/inorm
        vary += '/Tran./t/'+normtxt
    
    "---Apply absorption correction---"
    if abscor is not None and abscor is not False:
        # Absorption correction based on a flat plate perpendicular to the scattering plane
        # Requires abscor = absorption coefficient of material
        if abscor == True:
            abscor = 1.0
        A = scanvolcor(d,abscor)
        y = y/A
        dy = dy/A
        vary += '/A'
    return x,y,dy,varx,vary,ttl,d

def getmeta(nums=None, field='Energy'):
    """
    Get metadata from multiple scans and return array of values
    Usage: A = getmeta( [nums], 'field')
        [nums] = list of scan numbers
        field = name of metadata field
        A = array of values from each scan with same length as nums
            Any scans that do not contain field are returned 0
    """
    
    nums = np.asarray(nums).reshape([-1])
    
    # Prepare output array
    metavals = np.zeros(len(nums))
    for n,num in enumerate(nums):
        d = readscan(num)
        if d is None: continue
        
        # Get values from dat file
        if field in d.keys():
            val = np.mean(getattr(d,field))
        elif field in d.metadata.keys():
            val = getattr(d.metadata,field)
        else:
            print( 'ERROR: Scan ',num,' contains no ',field,' data' )
            continue
        
        metavals[n] = val
    if len(metavals) == 1:
        metavals = metavals[0]
    return metavals

def getmetas(nums=[0], fields=['Energy']):
    """
    Get several metadata values from multiple scans and return list of values
    Usage: A = getmetas( [nums], 'field')
        [nums] = list of scan numbers / or nums=12345
        fields = name or list of metadata field/ fields
        A = list of values from each scan with same length as nums
            Any scans that do not contain field are returned NaN
    """
    
    " Multiple nums given as input"
    try:
        nums = nums[0:]
    except TypeError:
        nums=[nums]
    
    " single field given"
    if type(fields) is str:
        fields = [fields]
    
    "Loop over scans"
    metavals = []
    for n,num in enumerate(nums):
        " scan object given as input"
        try:
            d = readscan(num)
        except TypeError:
            d = num
        if d is None: continue
        
        fieldvals = []
        for field in fields:
            # Get values from dat file
            if field.lower() == 'hkl':
                val = scanhkl(d)
            elif field in d.keys():
                val = np.mean(getattr(d,field))
            elif field in d.metadata.keys():
                val = getattr(d.metadata,field)
            else:
                val = np.NaN
                print( 'ERROR: Scan ',num,' contains no ',field,' data' )
            
            fieldvals += [val]
        
        metavals += [fieldvals]
    return metavals

def joindata(nums=None, varx='', vary='Energy', varz='', norm=True, sort=False, abscor=None, save=False):
    """
    Get useful values from scans, join multiple scans together, output joined arrays
     x,y,z,varx,vary,varz,ttl = joindata([nums],'varx','vary','varz')
     
      x = nxm array of scanned values, where n is the length of the scan and m is the number of scans
      y = nxm array of values that change with each scan
      z = nxm array of measured values
      varx,vary,varz = names of scan values
      ttl = title generated for these scans
      
      Options:
          norm    True/False    Normalisation options as in getdata
          sort    True/False    Sort the data by y
          abscor  True/False    Absorption correction as in getdata
          save    True/False    Save the resulting arrays in a numpy file
    
      If the output is saved as a numpy file, it can be loaded using the command:
          x,y,z,dz = np.load('path/outputfile.npy')
    """
    
    " Single 2D scan"
    if type(nums) == int or len(nums)==1:
        d = readscan(nums)
        keys = d.keys()
        if varx not in keys: varx = keys[0]
        if vary not in keys: vary = keys[1]
        x,z,dz,out_varx,out_varz,ttl,d = getdata(nums,varx=varx,vary=varz,norm=norm,abscor=abscor)
        y,z,dz,out_vary,out_varz,ttl,d = getdata(nums,varx=vary,vary=varz,norm=norm,abscor=abscor)
        
        " Determine steps in 1st dimension"
        cmd = d.metadata.cmd.split()
        ln=len(frange(float(cmd[2]),float(cmd[3]),float(cmd[4])))
        n_steps = len(x)/ln
        # Reshape
        storex = x.reshape([-1,n_steps])
        storey = y.reshape([-1,n_steps])
        storez = z.reshape([-1,n_steps])
        storedz = dz.reshape([-1,n_steps])
        skip = []
        nums = [nums]
    else:
        " Load first scan to assign data sizes"
        num = nums[0]
        x,z,dz,out_varx,out_varz,ttl,d = getdata(num,varx=varx,vary=varz,norm=norm,abscor=abscor)
        scn_points = len(x)
        n_scans = len(nums)
        ini_cmd = d.metadata.cmd_short
        
        " Assign storage arrays"
        storex = np.zeros([scn_points,n_scans])
        storey = np.zeros([scn_points,n_scans])
        storez = np.zeros([scn_points,n_scans])
        storedz = np.zeros([scn_points,n_scans])
        skip = []
    
        " Loop over each run and get data"
        for n,num in enumerate(nums):
            # Get x and z values using getdata
            x,z,dz,varx2,varz2,ttl2,d = getdata(num,varx,varz,norm,abscor)
            
            # Get y values from dat file
            if vary in d.keys():
                y = getattr(d,vary)
            elif vary in d.metadata.keys():
                yval = getattr(d.metadata,vary)
                if vary == 'psi' and yval < 0:
                    yval = yval+360.0
                y = np.ones(scn_points)*yval
            else:
                print( 'ERROR: Scan ',num,' contains no ',vary,' data' )
                return
            
            # Check correct length of scan
            if len(x) != scn_points:
                print( 'Error: Scan ',num,' has the wrong number of data points. Skipping...' )
                skip += [n]
                continue
            # Check correct type of scan
            if varx2 != out_varx:
                print( 'Error: Scan ',num,' is the wrong type of scan. Skipping...' )
                skip += [n]
                continue
            
            storex[:,n] = x # Scanned variable 
            storey[:,n] = y # dependent variable
            storez[:,n] = z # measured variable
            storedz[:,n] = dz # error on measured variable
        
        " Remove skipped columns"
        storex = np.delete(storex,skip,1)
        storey = np.delete(storey,skip,1)
        storez = np.delete(storez,skip,1)
        storedz = np.delete(storedz,skip,1)
        
        " Sort y"
        if sort:
            sort_index = np.argsort(storey[0,:])
            storex = storex[:,sort_index]
            storey = storey[:,sort_index]
            storez = storez[:,sort_index]
            storedz = storedz[:,sort_index]
    
    "Generate title"
    """
    vals=ttl.split()
    ettl = vals[0]
    rn = vals[1]
    hkl = vals[2]
    energy = vals[3]
    temp = vals[4]
    if len(vals) > 5: 
        pol = vals[5]
    else: 
        pol = ''
    fmt = '{ttl} {xvar}-{yvar} #{rn1:1.0f}-{rn2:1.0f} {hkl} {en} {temp} {pol}'
    if out_varx in ['Energy','en'] or vary in ['Energy','en']:
        energy=''
    if out_varx in ['Ta','Tb'] or vary in ['Ta','Tb']:
        temp=''
    out_ttl = fmt.format(ttl=ettl,rn1=nums[0],rn2=nums[-1],
                         xvar=out_varx,yvar=vary,
                         hkl=hkl,en=energy,temp=temp,pol=pol).strip()
    """
    out_ttl = scantitle(nums)
    
    " Save a .npy file"
    " Data can be loaded with x,y,z,dz = np.load(file.npy)"
    if save not in [None, False, '']:
        # save data to file
        if type(save) is str:
            savefile = os.path.join(savedir, '{}.npy'.format(save))
            head = '{}\n{}\n{}, {}, {}, {}'.format(out_ttl,ini_cmd,out_varx,vary,out_varz,'error_'+out_varz)
        else:
            savefile = os.path.join(savedir, '{} dep {}-{}.npy'.format(vary,nums[0],nums[-1]))
            head = '{}\n{}\n{}, {}, {}, {}'.format(out_ttl,ini_cmd,out_varx,vary,out_varz,'error_'+out_varz)
        np.save(savefile,(storex,storey,storez,storedz))#,header=head)
        print( 'Scans saved as {}'.format(savefile) )
        print( 'Load with: x,y,z,dz = np.load(\'{}\')'.format(savefile))
    
    return storex,storey,storez,out_varx,vary,out_varz,out_ttl

def getvol(num, ROIcen=None, ROIsize=None):
    """
    Load Pilatus images into a single volume: 
        vol=getvol(#/d)
        # is the scan number or dataloader from readscan
        vol is a [n x m x o] sized numpy array, where n and m are the size of the 
        pilatus detector (195 x 487) and o is the length of the scan.
    
    The volume from a region of interest can be obtained using:
        vol=getvol(#/d,ROIcen=[110,242],ROIsize=[31,31])
        
    """
    
    " Load data"
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    " Check 100k (small) or 2m (large) detector"
    try:
        pilname = [s for s in d.metadata.keys() if "_path_template" in s][0]
        pilpath = getattr(d.metadata,pilname)
    except IndexError:
        print( 'Not a pilatus file!' )
        return
    
    " Load first image to get the detector size"
    tif=pilpath % d.path[0]
    file = os.path.join(filedir,tif)
    file=file.replace('/',os.path.sep)
    #im=misc.imread(file,flatten=True) # flatten required to read zylar 16 bit files
    im = imread(file) # imageio
    pil_size = im.shape
    
    
    """
    if 'pilatus100k_path_template' in d.metadata.keys():
        pilpath = d.metadata.pilatus100k_path_template
        pil_size = [195,487]
    elif 'pilatus2m_path_template' in d.metadata.keys():
        pilpath = d.metadata.pilatus2m_path_template
        pil_size = [1679,1475]
    elif 'medipix_path_template' in d.metadata.keys():
        pilpath = d.metadata.medipix_path_template
        pil_size = [256,256]
    elif 'cam1_path_template' in d.metadata.keys():
        pilpath = d.metadata.cam1_path_template
        pil_size = [946,1292]
    else:
        print( 'Not a pilatus file!' )
        return
    """
    
    "Check volume size - don't load large volumes as it is too slow"
    numimag = len(d.path)
    if pil_size[0]*pil_size[1]*numimag > pil_max_size:
        print('The array size is too large! Returning only the first image')
        numimag = 1
    
    "Load each image into a single volume"
    vol = np.zeros([pil_size[0],pil_size[1],numimag]) # [195,487,~31]
    for n in range(numimag):
        " Prepare file name"
        tif=pilpath % d.path[n]
        file = os.path.join(filedir,tif)
        file=file.replace('/',os.path.sep)
        
        " Load image "
        #t=dnp.io.load(file,warn=False)
        #vol[:,:,n] = t.image0 #"image0" varies for some reason
        #im=misc.imread(file,flatten=True) # this is more reliable than dnp.io.load
        im = imread(file) # imageio
        
        " Flip image"
        #im = np.flipud(im)
        
        " Convert to double"
        im = 1.0*im
        " Assign dead pixels to median intensity, or zero for bad image (all -1)"
        im[im<0] = dead_pixel_func(im)
        im[im<0] = 0 
        " Assign hot pixels to median intensity, or zero for bad image (all hot)"
        im[im>hot_pixel] = dead_pixel_func(im)
        im[im>hot_pixel] = 0
        " Mask dead pixels ***Feb 2017***"
        im[pilatus_dead_pixels[:,0],pilatus_dead_pixels[:,1]] = 0
        " Add to 3D array"
        vol[:,:,n] = im
    
    if ROIsize is not None:
        if ROIcen is None:
            ROIcen = pil_centre
        " Only load the volume within the region of interest"
        idxi = [ROIcen[0]-ROIsize[0]//2,ROIcen[0]+ROIsize[0]//2+1]
        idxj = [ROIcen[1]-ROIsize[1]//2,ROIcen[1]+ROIsize[1]//2+1]
        #print( 'ROI = [{0},{1},{2},{3}]'.format(idxi[0],idxj[0],idxi[1],idxj[1]) )
        
        " Check the box is within the detector"
        if idxi[0] < 0: idxi[1] = idxi[1] - idxi[0]; idxi[0] = 0
        if idxi[1] > pil_size[0]: idxi[0] = idxi[0] - (idxi[1]-pil_size[0]); idxi[1] = pil_size[0]
        if idxj[0] < 0: idxj[1] = idxj[1] - idxj[0]; idxj[0] = 0
        if idxj[1] >  pil_size[1]: idxj[0] = idxj[0] - (idxj[1]-pil_size[1]); idxj[1] = pil_size[1]
        #print( 'new ROI = [{0},{1},{2},{3}]'.format(idxi[0],idxj[0],idxi[1],idxj[1]) )
        
        vol = vol[idxi[0]:idxi[1],idxj[0]:idxj[1],:]
    
    return vol

def pilroi(vol, ROIcen=None, ROIsize=[31,31], disp=False):
    """
    Define new ROI in Pilatus Detector
     ROI_sum,ROI_maxval,ROI_bkg = pilroi(vol,ROIcen,ROIsize)
    """
    
    
    "---Get basic data---"
    Nframe = vol.shape[2]
    pil_size = vol.shape[:2]
    
    if ROIcen is None:
        ROIcen = pil_centre
        
    pilsize = vol.shape[:2]
    
    "---Define ROI box---"
    idxi = [ROIcen[0]-ROIsize[0]//2, ROIcen[0]+ROIsize[0]//2+1] # [min, max] short axis
    idxj = [ROIcen[1]-ROIsize[1]//2, ROIcen[1]+ROIsize[1]//2+1] # [min, max] long axis
    if disp: print( 'ROI = [{0},{1},{2},{3}]'.format(idxi[0],idxj[0],idxi[1],idxj[1]) )
    
    " Check the box is within the detector"
    " Move the region, keeping ROIsize constant"
    if idxi[0] < 0: idxi[1] = idxi[1] - idxi[0]; idxi[0] = 0
    if idxi[1] > pil_size[0]: idxi[0] = idxi[0] - (idxi[1]-pil_size[0]); idxi[1] = pil_size[0]
    if idxj[0] < 0: idxj[1] = idxj[1] - idxj[0]; idxj[0] = 0
    if idxj[1] >  pil_size[1]: idxj[0] = idxj[0] - (idxj[1]-pil_size[1]); idxj[1] = pil_size[1]
    if disp: print( 'new ROI = [{0},{1},{2},{3}]'.format(idxi[0],idxj[0],idxi[1],idxj[1]) )
    
    "---Background ROI---"
    " Determine background of image and subtract from the region by pixel"
    " The background region is twice the area of the required region"
    idxbi = [ROIcen[0]-ROIsize[0],ROIcen[0]+ROIsize[0]]
    idxbj = [ROIcen[1]-ROIsize[1],ROIcen[1]+ROIsize[1]]
    
    " Check the box is within the detector"
    if idxbi[0] < 0: idxbi[1] = idxbi[1] - idxbi[0]; idxbi[0] = 0
    if idxbi[1] > pil_size[0]: idxbi[0] = idxbi[0] - (idxbi[1]-pil_size[0]); idxbi[1] = pil_size[0]
    if idxbj[0] < 0: idxbj[1] = idxbj[1] - idxbj[0]; idxbj[0] = 0
    if idxbj[1] > pil_size[1]: idxbj[0] = idxbj[0] - (idxbj[1]-pil_size[1]); idxbj[1] = pil_size[1]
    if disp: print( 'Background ROI = [{0},{1},{2},{3}]'.format(idxbi[0],idxbj[0],idxbi[1],idxbj[1]) )
    
    " Create background mask of the background ROI and remove the requested ROI"
    MASK = np.zeros(pil_size,dtype=int)
    MASK[idxbi[0]:idxbi[1],idxbj[0]:idxbj[1]] = 1 # Add background area
    MASK[idxi[0]:idxi[1],idxj[0]:idxj[1]] = 0 # remove peak ROI
    BKGidx = np.nonzero(MASK)
    BKGnorm = ROIsize[0]*ROIsize[1] # number of pixels in ROI
    
    "---Determine the ROI values for each frame---"
    ROI_sum = np.zeros(Nframe)
    ROI_maxval = np.zeros(Nframe)
    ROI_bkg = np.zeros(Nframe)
    for n in range(Nframe):
        ROI = vol[idxi[0]:idxi[1],idxj[0]:idxj[1],n]
        ROI_sum[n] = ROI.sum()
        ROI_maxval[n] = ROI.max()
        
        " Subtract the ROI from the background area"
        BKG = vol[BKGidx[0],BKGidx[1],n]
        " Instead of taking the mean, use the median"
        " this stops rogue peaks creating a false background"
        pixel_bkg = np.median(BKG)
        " However, if the pilatus gain is high, most pixels will be zero, so take the mean"
        if pixel_bkg == 0: pixel_bkg = np.mean(BKG)
        " Multiply the background by the number of pixels"
        ROI_bkg[n] = pixel_bkg*BKGnorm
    
    return ROI_sum,ROI_maxval,ROI_bkg

def pixel2hkl(num, detdim=[195,487], UB=None):
    """
    Generate hkl coordinates of detector pixel positions
      HHH,KKK,LLL = pixel2hkl(#/d, detdim,UB)
      
      detdim = detector dimesions def=[195,487] (Pilatus 100K)
      UB = orientation matrix, def=None (uses UB stored in metadata)
    
    Uses the pilatus parameters local:"pilpara" to define detector distance, centre etc.
        pilpara=[119.536,1904.17,44.4698,0.106948,-0.738038,412.19,-0.175,-0.175]
    
    Code copied and adapted from Alessandro's readpil.hklmatrix()
    
    Only works on delta, gam, eta, chi, phi scans - not energy or hkl scans!
    """
    
    "---Handle inputs---"
    if num is None:
        num = latest()
    
    "---Load data---"
    try:
        d = readscan(num)
    except TypeError:
        d = num
        num = d.metadata.SRSRUN
    
    varx = auto_varx(d)
    numimag = len(d.path)
    m = d.metadata
    scan_arrays = [x for x in d.keys() if type(d[x]) == np.ndarray ] # only keys linked to arrays
    pi=np.pi
    t1=time.clock()
    
    if varx in ['h','k','l','hkl']:
        # for hkl scans, convert motor angles to eulerian
        # Calculate eulerian angles
        # Taken from: /dls_sw/i16/software/gda/config/scripts/EulerianKconversion.py
        Kalpha=np.radians(50) # alpha angle on Kappa stage
        kap = np.radians(d.kap)
        alpha = -np.degrees(np.arctan( np.cos(Kalpha)*np.tan(kap/2) ))
        chi = -2*np.degrees(np.arcsin( np.sin(Kalpha)*np.sin(kap/2) ))
        theta = d.kth - alpha
        phi = d.kphi - alpha
        
        if 'delta' not in scan_arrays:
            d.delta = d.kdelta
        if 'gam' not in scan_arrays:
            d.gam = d.kgam
        if 'mu' not in scan_arrays:
            d.mu = d.kmu
        if 'eta' not in scan_arrays:
            d.eta = theta
        if 'chi' not in scan_arrays:
            d.chi = chi
        if 'phi' not in scan_arrays:
            d.phi = phi
        d.update(d.__dict__)
        scan_arrays = [x for x in d.keys() if type(d[x]) == np.ndarray ] # only keys linked to arrays
    elif varx not in ['delta','gam','eta','mu','chi','phi']:
        print("Doesnt work on this type of scan!")
    
    ########
    # Everything below here is from Alessandro's readpil.hklmatrix()
    ########
    
    # Rotation functions
    def R_x_r(alpha):
       "right rotation about x"
       alpha = np.float(alpha)
       ALPHA=np.matrix([[1.0, 0.0, 0.0],[0.0, np.cos(alpha),-np.sin(alpha)], [0.0, np.sin(alpha), np.cos(alpha)]])
       return ALPHA
    
    def R_y_r(alpha):
       "right rotation about y"
       alpha = np.float(alpha)
       ALPHA=np.matrix([[np.cos(alpha),0.0, np.sin(alpha)],[0.0, 1.0, 0.0], [-np.sin(alpha), 0.0, np.cos(alpha)]])
       return ALPHA
    
    def R_z_l(alpha):
       "left rotation about z"
       alpha = np.float(alpha)
       ALPHA=np.matrix([[np.cos(alpha),np.sin(alpha), 0.0],[-np.sin(alpha), np.cos(alpha), 0.0], [0.0, 0.0, 1.0]])
       return ALPHA
    
    # detector  and energy initialization
    #indipendent from the scan type
    temp=np.ones((numimag,1))
    H=np.ones([numimag, detdim[0], detdim[1]])
    K=np.ones([numimag, detdim[0], detdim[1]])
    L=np.ones([numimag, detdim[0], detdim[1]])
    
    wl=12.39842/m.Energy
    
    ni=pilpara[0] 
    nj=pilpara[1] 
    delp=pilpara[2] 
    gamp=pilpara[3] 
    pSi=pilpara[4] 
    ndist=pilpara[5]
    x_pixel=pilpara[6] 
    y_pixel=pilpara[7]
    
    x=np.matrix([1,0,0]).transpose()
    y=np.matrix([0,1,0]).transpose()
    z=np.matrix([0,0,1]).transpose()
    
    rotgammaprime=R_x_r(gamp*pi/180)
    rotdeltaprime=R_z_l(delp*pi/180)
    
    
    n0=rotgammaprime*rotdeltaprime*y*ndist
    n0=n0.transpose()
    
    #UB initialization 
    if UB is None:
        UB=np.matrix([[m.UB11,m.UB12,m.UB13],[m.UB21,m.UB22,m.UB23],[m.UB31,m.UB32,m.UB33]])
        invUB=UB.I
    else:
        invUB=np.linalg.inv(UB)
    
    
    Ex=np.cross(x.transpose(),n0)/np.linalg.linalg.norm(np.cross(x.transpose(),n0))*x_pixel;
    Ey=np.cross(n0,Ex)/np.linalg.linalg.norm(np.cross(n0,Ex))*y_pixel;  
    
    Dx0=Ex*np.cos(pSi*pi/180)-Ey*np.sin(pSi*pi/180);
    Dy0=Ex*np.sin(pSi*pi/180)+Ey*np.cos(pSi*pi/180);
    
    #initialization HHH,KKK,LLL matrices
    
    HHH=np.arange(0,detdim[0]*detdim[1]*numimag,1,float).reshape(detdim[0],detdim[1],numimag)
    KKK=np.arange(0,detdim[0]*detdim[1]*numimag,1,float).reshape(detdim[0],detdim[1],numimag)
    LLL=np.arange(0,detdim[0]*detdim[1]*numimag,1,float).reshape(detdim[0],detdim[1],numimag)
    
    #inizialization of arrays in the single image
    #~ PiPj=np.arange(0,detdim[0]*detdim[1],1,float).reshape(detdim[0],detdim[1])
    Pjg,Pig=np.meshgrid(range(0,detdim[1],1),range(0,detdim[0],1))
    
    normP=np.arange(0,detdim[0]*detdim[1],1,float).reshape(detdim[0],detdim[1])
    nval=np.arange(0,detdim[0]*detdim[1],1,float).reshape(detdim[0],detdim[1])
    values=np.arange(0,detdim[0]*detdim[1]*3,1,float).reshape(detdim[0],detdim[1],3)
    khat=np.arange(0,detdim[0]*detdim[1]*3,1,float).reshape(detdim[0],detdim[1],3)
    k=np.arange(0,detdim[0]*detdim[1]*3,1,float).reshape(detdim[0],detdim[1],3)    
    
    #######################################
    ##Rotation Matrices - pre-calculation##
    #######################################
    floatsdet=False
    floatschiphi=False
    
    rotdelta=R_z_l((m.delta-abs(m.delta_axis_offset))*pi/180)
    rotgamma=R_x_r(m.gam*pi/180)
    roteta=R_z_l(m.eta*pi/180)
    rotmu=R_x_r(m.mu*pi/180)
    rotchi=R_y_r(m.chi*pi/180)
    rotphi=R_z_l(m.phi*pi/180)
    rotchiphi=rotchi*rotphi
    
    # Not chi or phi scan
    if 'chi' not in scan_arrays and 'phi' not in scan_arrays:
        floatschiphi=True
    
    # Not delta or gam scan - pre-calculate matrices for speed
    if 'delta' not in scan_arrays and 'gam' not in scan_arrays:
        # Calculate the final wavevector, k
        floatsdet=True
        n=rotgamma*rotdelta*n0.transpose();
        Dx=rotgamma*rotdelta*Dx0.transpose();
        Dy=rotgamma*rotdelta*Dy0.transpose();
        
        nt=n.transpose()
        Dxt=Dx.transpose()
        Dyt=Dy.transpose()
        
        temp1=((Pig-ni).reshape(detdim[0],detdim[1]).repeat(3)).reshape(detdim[0],detdim[1],3)*np.array(Dxt)
        temp2=(Pjg-nj).repeat(3).reshape(detdim[0],detdim[1],3)*np.array(Dyt)+np.array(nt)
        values=temp1+temp2
        
        nval=np.sqrt(values[:,:,0]**2+values[:,:,1]**2+values[:,:,2]**2)
        
        for indice in [0,1,2]:
            values[:,:,indice]=values[:,:,indice]/nval;
        
        Phat_m_y=-values
        Phat_m_y[:,:,1]=Phat_m_y[:,:,1]+1
        normP=np.sqrt(Phat_m_y[:,:,0]**2+Phat_m_y[:,:,1]**2+Phat_m_y[:,:,2]**2)
        
        for indice in [0,1,2]:
            khat[:,:,indice]=Phat_m_y[:,:,indice]/normP
        
        for indice in [0,1,2]:
            k[:,:,indice]=-khat[:,:,1]*2/wl*khat[:,:,indice]
    
    ############################
    ##Rotation Matrices - Scan##
    ############################
    for index_r in range(numimag):
        # delta or gam move:
        if 'gam' in scan_arrays: 
            rotgamma=R_x_r(d.gam[index_r]*pi/180)
        if 'delta' in scan_arrays: 
            rotdelta=R_z_l((d.delta[index_r]-abs(d.delta_axis_offset[index_r]))*pi/180) 
        if floatsdet is False:
            # This is identical to the part above, except rotdelta or rotgamma changes on each iteration
            # Calculate the final wavevector, k
            n=rotgamma*rotdelta*n0.transpose();
            Dx=rotgamma*rotdelta*Dx0.transpose();
            Dy=rotgamma*rotdelta*Dy0.transpose();
            
            nt=n.transpose()
            Dxt=Dx.transpose()
            Dyt=Dy.transpose()
            
            temp1=((Pig-ni).reshape(detdim[0],detdim[1]).repeat(3)).reshape(detdim[0],detdim[1],3)*np.array(Dxt)
            temp2=(Pjg-nj).repeat(3).reshape(detdim[0],detdim[1],3)*np.array(Dyt)+np.array(nt)
            values=temp1+temp2
            
            nval=np.sqrt(values[:,:,0]**2+values[:,:,1]**2+values[:,:,2]**2)
            
            for indice in [0,1,2]:
                values[:,:,indice]=values[:,:,indice]/nval
            
            Phat_m_y=-values
            Phat_m_y[:,:,1]=Phat_m_y[:,:,1]+1;
            normP=np.sqrt(Phat_m_y[:,:,0]**2+Phat_m_y[:,:,1]**2+Phat_m_y[:,:,2]**2)
            
            for indice in [0,1,2]:
                khat[:,:,indice]=Phat_m_y[:,:,indice]/normP
            
            for indice in [0,1,2]:
                k[:,:,indice]=-khat[:,:,1]*2/wl*khat[:,:,indice]
        
        # change the rotation matrices for the current scan
        if 'eta' in scan_arrays:
            roteta=R_z_l(d.eta[index_r]*pi/180) 
        if 'mu' in scan_arrays: 
            rotmu=R_x_r(d.mu[index_r]*pi/180)
        if 'phi' in scan_arrays: 
            rotphi=R_z_l(d.phi[index_r]*pi/180)
        if 'chi' in scan_arrays: 
            rotchi=R_y_r(d.chi[index_r]*pi/180)
        if floatschiphi is not True:
            rotchiphi=rotchi*rotphi
        
        # Generate the general rotation matrix, then calculate the hkl mesh
        Z=rotmu*(roteta*(rotchiphi))
        invUBinvZ=invUB*Z.I;
        HHH[:,:,index_r]=invUBinvZ[0,0]*k[:,:,0] + invUBinvZ[0,1]*k[:,:,1] + invUBinvZ[0,2]*k[:,:,2]
        KKK[:,:,index_r]=invUBinvZ[1,0]*k[:,:,0] + invUBinvZ[1,1]*k[:,:,1] + invUBinvZ[1,2]*k[:,:,2];
        LLL[:,:,index_r]=invUBinvZ[2,0]*k[:,:,0] + invUBinvZ[2,1]*k[:,:,1] + invUBinvZ[2,2]*k[:,:,2];
    
    t2=time.clock()
    print("#{} Time spent generating hkl positions={}".format(num,t2-t1))
    return HHH,KKK,LLL

def pixel2xyz(num, detdim=[195,487]):
    """
    Generate cartesian coordinates of detector pixel positions, in A-1
      XXX,YYY,ZZZ = pixel2xyz(#/d, detdim)
      
      detdim = detector dimesions def=[195,487] (Pilatus 100K)
    
    Uses the pilatus parameters local:"pilpara" to define detector distance, centre etc.
        pilpara=[119.536,1904.17,44.4698,0.106948,-0.738038,412.19,-0.175,-0.175]
    
    Code copied and adapted from Alessandro's readpil.hklmatrix()
    
    Only works on delta, gam, eta, chi, phi scans - not energy or hkl scans!
    """
    
    UB = np.eye(3)/(2*np.pi) # B matrix if a/b/c=2pi
    XXX,YYY,ZZZ = pixel2hkl(num,detdim,UB)
    return XXX,YYY,ZZZ

def pixel2tth(num, detdim=[195,487], centre_only=False, norm=True):
    """
    Generate two-theta coordinates of detector pixel positions
      TTH,INT = pixel2tth(#/d, detdim)
      
      detdim = detector dimesions def=[195,487] (Pilatus 100K)
    
    Uses the pilatus parameters local:"pilpara" to define detector distance, centre etc.
        pilpara=[119.536,1904.17,44.4698,0.106948,-0.738038,412.19,-0.175,-0.175]
    
    To bin the data, use: 
        bin_tth,bint_int = bindata(TTH,INT,0.01) 
    
    Only works on delta, gam, eta, chi, phi scans - not energy or hkl scans!
    """
    
    XXX,YYY,ZZZ = pixel2xyz(num,detdim)
    t1=time.clock()
    Qmag = np.sqrt(XXX**2 + YYY**2 + ZZZ**2)
    
    TTH = q2tth(Qmag,getmeta(num,'Energy'))
    vol = getvol(num)
    norm_val,norm_txt = normalise(num)
    
    # Normalise vs time etc
    if norm:
        vol = vol/norm_val
    
    if centre_only:
        cen = pil_centre[0]
        TTH = TTH[cen-10:cen+10,:,:]
        vol = vol[cen-10:cen+10,:,:]
    
    TTH = TTH.flatten()
    vol = vol.flatten()
    sortindex = np.argsort(TTH)
    TTH = TTH[sortindex]
    vol = vol[sortindex]
    
    t2=time.clock()
    print("time spent sorting tth values={}".format(t2-t1))
    
    return TTH,vol

def pixel2tth2(num, detdim=[195,487], norm=True, pixel_centre=[104, 205], frame_centre=0, pixel_width=10, frame_width=4):
    """
    Generate two-theta coordinates of detector pixel positions
      TTH,INT = pixel2tth2(#/d, detdim, norm, pixel_centre, frame_centre, pixel_width, frame_width)
      
      detdim = detector dimesions def=[195,487] (Pilatus 100K)
      pixel_centre = [i,j] coordinates of the peak on the detector
      frame_centre = frame number of peak maximum
      pixel_width = number of pixels to integrate over in the short direction on the detector
      frame_width = number of frames to integrate over, around the frame_centre
    
    This function works differently to pixel2tth as the tth values are taken from a single frame, with
    intensity values integrated around a central point and plotted against the long axis of the detector.
    
    E.G.
        vol = getvol(12345)
        [i,j], frame = pilpeak(vol)
        tth, int = pixel2tth2(12345, pixel_centre=[i,j], frame_centre=frame)
    
    Shortcuts:
        pixel_centre = 'peak' >> takes [i,j], frame from pilpeak
        pixel_centre = 'centre' >> takes [i,j] from pil_centre attribute, frame as half way throught scan
    
    Uses the pilatus parameters local:"pilpara" to define detector distance, centre etc.
        pilpara=[119.536,1904.17,44.4698,0.106948,-0.738038,412.19,-0.175,-0.175]
    Only works on delta, gam, eta, chi, phi scans - not energy or hkl scans!
    """
    
    XXX,YYY,ZZZ = pixel2xyz(num,detdim)
    t1=time.clock()
    Qmag = np.sqrt(XXX**2 + YYY**2 + ZZZ**2)
    TTH = q2tth(Qmag,getmeta(num,'Energy'))
    vol = getvol(num)
    norm_val,norm_txt = normalise(num)
    
    # Normalise vs time etc
    if norm:
        vol = vol/norm_val
    
    if pixel_centre == 'peak':
        pixel_centre, frame_centre = pilpeak(vol,disp=True)
    elif pixel_centre == 'centre':
        pixel_centre = pil_centre
        frame = vol.shape[2]//2
    
    maxdif_i = np.max(TTH[:,pixel_centre[1],frame_centre])-np.min(TTH[:,pixel_centre[1],frame_centre])
    maxdif_j = np.max(TTH[pixel_centre[0],:,frame_centre])-np.min(TTH[pixel_centre[0],:,frame_centre])
    if maxdif_i > maxdif_j:
        # Horizontal geometry (larger change in TTH along short axis, i)
        tth = TTH[:,pixel_centre[1],frame_centre]
        voltth = vol[:,pixel_centre[1]-pixel_width//2:pixel_centre[1]+pixel_width//2,frame_centre-frame_width//2:frame_centre+frame_width//2]
        I = voltth.sum(axis=1).sum(axis=1)
    else:
        # Vertical geometry (larger change in TTH along long axis, j)
        tth = TTH[pixel_centre[0],:,frame_centre]
        voltth = vol[pixel_centre[0]-pixel_width//2:pixel_centre[0]+pixel_width//2,:,frame_centre-frame_width//2:frame_centre+frame_width//2]
        I = voltth.sum(axis=0).sum(axis=1)
    t2=time.clock()
    print("time spent sorting tth values={}".format(t2-t1))
    
    return tth,I

def pixel2chi(num, detdim=[195,487], norm=True, pixel_centre=[104, 205], frame_centre=0, pixel_width=10, frame_width=4):
    """
    Generate two-theta coordinates of detector pixel positions
      TTH,INT = pixel2tth2(#/d, detdim, norm, pixel_centre, frame_centre, pixel_width, frame_width)
      
      detdim = detector dimesions def=[195,487] (Pilatus 100K)
      pixel_centre = [i,j] coordinates of the peak on the detector
      frame_centre = frame number of peak maximum
      pixel_width = number of pixels to integrate over in the short direction on the detector
      frame_width = number of frames to integrate over, around the frame_centre
    
    This function works differently to pixel2tth as the tth values are taken from a single frame, with
    intensity values integrated around a central point and plotted against the long axis of the detector.
    
    E.G.
        vol = getvol(12345)
        [i,j], frame = pilpeak(vol)
        tth, int = pixel2tth2(12345, pixel_centre=[i,j], frame_centre=frame)
    
    Shortcuts:
        pixel_centre = 'peak' >> takes [i,j], frame from pilpeak
        pixel_centre = 'centre' >> takes [i,j] from pil_centre attribute, frame as half way throught scan
    
    Uses the pilatus parameters local:"pilpara" to define detector distance, centre etc.
        pilpara=[119.536,1904.17,44.4698,0.106948,-0.738038,412.19,-0.175,-0.175]
    Only works on delta, gam, eta, chi, phi scans - not energy or hkl scans!
    """
    
    XXX,YYY,ZZZ = pixel2xyz(num,detdim)
    t1=time.clock()
    CHI=90+np.rad2deg(np.arcsin(YYY/ZZZ))
    vol = getvol(num)
    norm_val,norm_txt = normalise(num)
    
    # Normalise vs time etc
    if norm:
        vol = vol/norm_val
    
    if pixel_centre == 'peak':
        pixel_centre, frame_centre = pilpeak(vol,disp=True)
    elif pixel_centre == 'centre':
        pixel_centre = pil_centre
        frame = vol.shape[2]//2
    
    chi = CHI[:,pixel_centre[1],frame_centre]
    volchi = vol[:,pixel_centre[1]-pixel_width//2:pixel_centre[1]+pixel_width//2,frame_centre-frame_width//2:frame_centre+frame_width//2]
    I = volchi.sum(axis=1).sum(axis=1)
    t2=time.clock()
    print("time spent sorting chi values={}".format(t2-t1))
    
    return chi,I

def readnexus(num, nexusformat=False):
    """
    Read nexus file (.nxs)
     n = readnexus(scan_no)

    The Nexus format is a hdf5 type file with many more options for storing data
    Nexus files are read with the h5py package, however the nexusformat package
    provides additional functionality:

    If using nexusformat:
        from nexusformat.nexus import nxload
        n = nxload('12345.nxs')
        n = readnexus(12345, nexusformat=True)
        eta = n.entry1.measurement.eta
        sum = n['entry1/measurement/sum']
        print(n.tree)
        n.entry1.measurement.plot()
    If using h5py:
        import h5py
        n = h5py.File('12345.nxs')
        n = readnexus(12345)
        eta = n['entry1/measurement/eta'][:]
        sum = n['entry1/measurement/sum'][:]
    """

    if os.path.isdir(filedir) == False: 
        print( "I can't find the directory: {}".format(filedir) )
        return None
    
    if num < 1: 
        if latest() is None: return None
        num = latest()+num
    
    file = os.path.join(filedir, nxsfile_format %num)
    #print(file)
    try:
        if nexusformat:
            n = nxload(file)
        else:
            n = h5py.File(file)
    except:
        print( "Scan {} doesn't exist or can't be read".format(num) )
        return None
    return n

"------------------------Experiment Check Functions-----------------------"


def latest():
    "Get the latest run during and experiment"
    
    if os.path.isdir(filedir) == False: 
        print( "I can't find the directory: {}".format(filedir) )
        return None
    
    # Get all data files in folder
    ls=glob.glob(filedir+os.path.sep+'*.dat')
    ls = np.sort(ls)
    
    if len(ls) < 1:
        print( "No files in directory: {}".format(filedir) )
        return
    
    #newest = ls[-1] # file with highest number
    newest = max(ls, key=os.path.getctime) # file created last
    #num = np.int(os.path.split(newest)[-1][:-4])
    num = np.abs(int(os.path.split(newest)[-1][-10:-4])) # handles i10-#####.dat and ######.dat
    return num

def checkexp():
    "Check experiment folder data"
    
    # Get all data files in folder
    ls=glob.glob(filedir+os.path.sep+'*.dat')
    ls = np.sort(ls)
    
    if len(ls) < 1:
        print( "I can't find the directory: {}".format(filedir) )
        return
    
    # First + last run numbers
    fn = np.int(re.findall('\d+',os.path.split(ls[0])[1])[0])
    ln = np.int(re.findall('\d+',os.path.split(ls[-1])[1])[0])
    
    # Load data
    fd = readscan(fn)
    ld = readscan(ln)
    
    # Get first and last file creation dates (not reliable)
    #ft = datetime.datetime.fromtimestamp(os.path.getctime(ls[0]))
    #lt = datetime.datetime.fromtimestamp(os.path.getctime(ls[-1]))
    
    # Use SRSDAT and SRSTIM - because date isn't in every file
    # Note that SRSDAT stores months and days without 0's, 
    # so 11 Feb (2015112) and 1 Dec (2015112) look the same! Not sure how to fix this!
    ft = datetime.datetime.strptime(str(fd.metadata.SRSDAT)+' '+str(fd.metadata.SRSTIM),'%Y%m%d %H%M%S')
    lt = datetime.datetime.strptime(str(ld.metadata.SRSDAT)+' '+str(ld.metadata.SRSTIM),'%Y%m%d %H%M%S')
    
    # Format times
    ft = ft.strftime('%Y-%m-%d %H:%M:%S,%f')
    lt = lt.strftime('%Y-%m-%d %H:%M:%S,%f')
    
    #if 'date' in fd.metadata.keys(): ft = fd.metadata.date
    #if 'date' in ld.metadata.keys(): lt = ld.metadata.date
    
    # Print data to screen
    print( 'Experiment: ',filedir )
    print( ' First scan: #{0}    {1}'.format(fn,ft) )
    print( '  Last scan: #{0}    {1}'.format(ln,lt) )
    print( '  No. scans: {}'.format(len(ls)) )
    print( 'Experiment ring current: ',exp_ring_current )
    print( 'Experiment monitor: ',exp_monitor )
    print( 'Normalisation option: ',normby )
    print( 'Pilatus Centre: ',pil_centre )
    return

def checkscan(num=None, showval=None):
    """
    Get run number information,  returns a string
        checknum(#)             Display lots of information about a single scan
        checknum([run1,run2,...]) List all the scans in the list or array
    
    if # <= 0, # will be given as latest()-#
    
    checknum(...,showval='phi'): also display infomation from a variable, such as 'phi'
    
    """
    
    if os.path.isdir(filedir) == False: 
        print( "I can't find the directory: {}".format(filedir) )
        return ''
    
    if num is None:
        num = latest()
    
    # turn into array
    num = np.asarray(num).reshape(-1)
    
    "-----------------Multi run------------------"
    if len(num)>1:
        outstr = checkscans(num,showval=showval)
        return outstr
    
    "----------------Single run------------------"
    d = readscan(num[0])
    if d is None:
        print( "File does not exist!" )
        return 'File does not exist!'
    m = d.metadata
    ks = d.keys()
    
    # Print information
    outstr = ''
    outstr += '-----------Run %s-----------\n' % m.SRSRUN
    outstr += '  File Dir: {}\n'.format(filedir)
    outstr += '   Command: {}\n'.format(m.cmd)
    outstr += '   Npoints: {}\n'.format(len(d.TimeSec))
    outstr += '       HKL: ({0},{1},{2})\n'.format(m.h,m.k,m.l)
    outstr += '    Energy: {0} keV\n'.format(m.Energy)
    outstr += '      Temp: {0} K\n\n'.format(m.temperature)
    
    # Check for vertical or horizontal geopmetry
    if m.gam > 0.1:
        outstr += '    Horizontal Geometry\n'
        outstr += '        Mu: {0}\n'.format(m.mu)
        outstr += '       Chi: {0}\n'.format(m.chi)
        outstr += '     Gamma: {0}\n'.format(m.gam)
    else:
        outstr += '      Vertical Geometry'
        outstr += '       Eta: {0}\n'.format(m.eta)
        outstr += '       Chi: {0}\n'.format(m.chi)
        outstr += '     Delta: {0}\n'.format(m.delta)
    outstr += 'Psi({0},{1},{2}): {3}\n\n'.format(m.azih,m.azik,m.azil,m.psi)

    outstr += '       Sx: {0} mm\n'.format(m.sx)
    outstr += '       Sy: {0} mm\n'.format(m.sy)
    outstr += '       Sz: {0} mm\n\n'.format(m.sz)
    
    outstr += '  Sample Slits: {}\n'.format(scanss(num[0]))
    outstr += 'Detector Slits: {}\n\n'.format(scands(num[0]))
    
    # Minimirrors
    if m.m4x > 0.1: 
        mm = 'in' 
    else: 
        mm = 'out'
    outstr += 'Minimirrors: {} ({:4.2f} deg)\n'.format(mm,m.m4pitch)
    
    if 'sum' in ks:
        # Pilatus
        outstr += '\n'
        outstr += 'Detector: Pilatus (do={})\n'.format(m.delta_axis_offset)
        outstr += '   Count: {0:1.3g}s\n'.format(np.mean(d.t))
        outstr += ' Max val: {0:5.3g}\n'.format(max(d.maxval))
        
    if 'APD' in ks:
        # APD
        outstr += '\n'
        outstr += 'Detector: APD\n'
        if m.gam > 0.:
            # Horizontal
            if m.stoke < 45.:
                outstr += '   Pol: {0} (pi-pi)\n'.format(m.stoke)
            else:
                outstr += '   Pol: {0} (pi-sigma)\n'.format(m.stoke)
        else:
            # Vertical
            if m.stoke < 45.:
                outstr += '   Pol: {0} (sigma-sigma)\n'.format(m.stoke)
            else:
                outstr += '   Pol: {0} (sigma-pi)\n'.format(m.stoke)
        outstr += '   Count: {0:1.3g}s\n'.format(np.mean(d.counttime))
        outstr += ' Max val: {0:5.3g}\n'.format(max(d.APD))
    
    if 'FF' in ks:
        outstr += '\n'
        outstr += 'Detector: Vortex\n'
        outstr += '   Count: {0:1.3g}s\n'.format(np.mean(d.count_time))
        outstr += '   ROIs (maxval)\n:'
        ROIs = [n for n in d.keys() if 'Element' in n]
        for roi in ROIs:
            outstr += '    {} ({})\n'.format(roi,max(getattr(d,roi)))
        
    # Attenuation
    outstr += '\n    Atten: {0} ({1}%)\n\n'.format(m.Atten,m.Transmission*100)
    
    # additional info
    if showval is not None:
        if showval in ks:
            val = np.mean(getattr(d,showval))
        elif showval in m.keys():
            val = getattr(m,showval)
        else:
            val = 'No Data'
        outstr += '{} = {}\n\n'.format(showval, val)

    # Timing
    time = d.TimeSec[-1]-d.TimeSec[0]
    hours = np.int(np.floor(time/3600))
    mins = np.int(np.floor(np.remainder(time,3600)/60))
    secs = np.remainder(np.remainder(time,3600),60)
    outstr += 'Ran on {}\n'.format(m.date)
    outstr += 'Time taken: {} hours, {} mins, {} seconds\n'.format(hours,mins,secs)
    return outstr

def checkscans(num1=None,num2=None,showval=None,find_scans=None):
    """
    Get brief info about each scan in a list, returns a string
    
    checkscans([587892,587893,587894])
    
    checkscans(...,showval='phi') >> also display infomation from a variable, such as 'phi'
    
    checkscans(...,find_scans='eta') >> only display 'eta' scans
    """
    
    if num1 is None:
        num1 = range(-10,1)
    if num2 is not None:
        num1 = range(num1,num2+1)
    nums = np.asarray(num1,dtype=int).reshape(-1)
    
    showvals = np.asarray(showval).reshape(-1)
    
    # Print brief info
    outstr = ''
    #fmt = '{nums} | {date} | {mode:4s} {energy:5.3g}keV {temp:5.3g}K ({h:1.2g},{k:1.2g},{l:1.2g}) | {cmd}'
    fmt = '{nums} | {date} | {mode:4s} {ss} {ds} {energy:5.3g}keV {temp:6s} {hkl:17s} {showval} | {cmd}\n'
    #showval_dict = {'show':'','equal':'','val':''}
    for n in range(len(nums)):
        d = readscan(nums[n])
        if d is None:
            continue
        
        m = d.metadata
        ks = d.keys()
        varx = auto_varx(d)
        
        # Only display scans of type find_scans
        if find_scans is not None and find_scans != varx:
            continue
        
        sampsl = '{0:4.2g}x{1:<4.2g}'.format(m.s5xgap,m.s5ygap)
        detsl = '{0:4.2g}x{1:<4.2g}'.format(m.s6xgap,m.s6ygap)
        
        h = round(m.h*10)/10.0 + 0.0 # + 0.0 to remove -0 
        k = round(m.k*10)/10.0 + 0.0
        l = round(m.l*10)/10.0 + 0.0
        hkl = '({h:1.2g},{k:1.2g},{l:1.2g})'.format(h=h,k=k,l=l)
        if m.gam>0.1:
            mode= 'H'
            pol = 'p'
        else:
            mode='V'
            pol='s'
        if 'sum' in ks: mode=mode+'Pil'
        if 'FF' in ks: mode=mode+'Vtx'
        if 'xmapc' in ks: mode=mode+'xma'
        if 'APD' in ks:
            if m.stoke < 45.: 
                if m.gam>0.0:
                    pol = pol+'p'
                else:
                    pol = pol+'s'
            else:
                if m.gam>0.0:
                    pol = pol+'s'
                else:
                    pol = pol+'p'
            mode = mode+pol
        
        showvaltxt = []
        for showval in showvals:
            if showval is None: continue
            
            if showval == 'hkl':
                val = scanhkl(m.SRSRUN)
            elif showval in ks:
                val = np.mean(getattr(d,showval))
            elif showval in m.keys():
                val = getattr(m,showval)
            else:
                val = 'No Data'
            showvaltxt += ['{} = {}'.format(showval,val)]
        showvaltxt = ', '.join(showvaltxt)
        
        outstr += fmt.format(nums=m.SRSRUN,date=m.date,mode=mode,ss=sampsl,ds=detsl,energy=m.Energy,temp=m.temperature,hkl=hkl,showval=showvaltxt,cmd=m.cmd_short)
    return outstr

def checklog(time=None,mins=2,cmd=False,find=None):
    """
    Look at experiment log file, returning specific part as str

    returns log file from time-mins to time

    checklog(time,mins,commands,findstr)
      time = None - Uses current time (default)
      time = scan number (uses last modified time of this scan), 0 gives latest scan
      time = [hour,(min),(day),(month),(year)], () values set to default, defaults are time of last scan
      time = '2015-07-07 13:50:08,000' 
      mins = log file will display from time-mins to time
      mins = 'all' - display whole log from first scan time
      commands = True/False - only display command strings
      find = only display lines that include findstr
    """

    filename = os.path.join(filedir, 'gdaterminal.log')
    
    # time is a datetime object
    if time is None:
        #date = readscan(0).metadata.date
        #time = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
        #time = datetime.datetime.now()
        time = datetime.datetime.fromtimestamp(os.stat(filename).st_mtime) # last modified time

    
    if type(time) is str:
        time=datetime.datetime.strptime(time,'%Y-%m-%d %H:%M:%S,%f')
    
    if type(time) is not type(datetime.datetime.now()):
        if type(time) is int:
            # input is a scan number (0 for latest)
            #date = readscan(time).metadata.date
            #time = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
            file = scanfile(time)
            time = datetime.datetime.fromtimestamp(os.stat(file).st_mtime) # last modified time
            # Add 10s for good measure
            time = time + datetime.timedelta(0,10) # days, seconds
        else:
            date = readscan(0).metadata.date
            now = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
            if len(time) < 1: time+=[now.hour]
            if len(time) < 2: time+=[now.minute]
            if len(time) < 3: time+=[now.day]
            if len(time) < 4: time+=[now.month]
            if len(time) < 5: time+=[now.year]
            
            if time[4] < 1000: time[4]+=2000 # 15 > 2015
            
            hour = time[0]
            min = time[1]
            day = time[2]
            month = time[3]
            year = time[4]
            time = datetime.datetime(year,month,day,hour,min,0)
        print(time)
    
    if mins == 'all':
        # Get all data files in folder
        ls=glob.glob(filedir+os.path.sep+'*.dat')
        ls = np.sort(ls)
        fn = np.int(re.findall('\d+',os.path.split(ls[0])[1])[0])
        # Load first run
        fd = readscan(fn)
        ft = datetime.datetime.strptime(str(fd.metadata.SRSDAT)+str(fd.metadata.SRSTIM),'%Y%m%d%H%M%S')
        
        # Ge number of minutes
        diff = time - ft
        mins = np.ceil(diff.total_seconds()/60.0)
        print(mins)
    
    mintime = time - datetime.timedelta(minutes=mins)
    
    outstr = ''
    with open(filename) as file:
        for line in file:
            if cmd and '>>>' not in line:
                continue
            
            if find is not None and find not in line:
                continue
            
            spt = line.split()
            tim=datetime.datetime.strptime(spt[0]+spt[1],'%Y-%m-%d%H:%M:%S,%f')
            
            if tim >= mintime and tim <= time:
                outstr += line
    
    return outstr

def auto_varx(num):
    "Determine the scanned variable"
    
    "---Load data---"
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    
    "---Get metadata---"
    keys = [x for x in d.keys() if type(d[x]) == np.ndarray ] # only keys linked to arrays
    m = d.metadata
    cmd = m.cmd # Scan command
    
    varx = cmd.split()[1]
    if varx == 'hkl':
        if cmd.split()[0] == 'scan':
            hvar,kvar,lvar=cmd.split()[8:11]
        elif cmd.split()[0] == 'scancn':
            hvar,kvar,lvar=cmd.split()[2:5]
        else:
            print( 'Error: Wrong type of scan' )
        hvar = float(re.sub("[^0-9.]", "",hvar))
        kvar = float(re.sub("[^0-9.]", "",kvar))
        lvar = float(re.sub("[^0-9.]", "",lvar))
        if hvar > 0.0:
            varx = 'h'
        elif kvar > 0.0:
            varx = 'k'
        else:
            varx = 'l'
    elif varx == 'sr2':
        if 'azimuthal' in keys: # DiffCalc sr2 scan
            varx = 'azimuthal'
        else:
            varx = 'phi'
    elif varx == 'th2th':
        varx = 'delta'
        
    if varx not in keys:
        varx = keys[0]
    if '_energy' in varx:
        varx = varx.replace('_energy', '_offset')
    return varx

def auto_vary(num):
    "Determine the default scannable"
    
    "---Load data---"
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    
    "---Get metadata---"
    keys = d.scannables
    m = d.metadata
    cmd = m.cmd # Scan command
    
    vary = cmd.split()[-1]
    if vary[0:3] == 'roi':
        vary = vary + '_sum'
    elif 'pil100k' in cmd:
        vary = 'sum'
    elif 'pil2m' in cmd:
        vary = 'sum'
    elif 'APD' in keys:
        vary = 'APD'
    elif 'xmapc' in keys:
        vary = 'xmroi2'
    elif 'bpm' in cmd:
        vary = 'sum'
    elif 'QBPM6' in cmd:
        vary = 'C1'
    
    if vary not in keys:
        vary = keys[-1]
    return vary

def scantitle(num):
    "Generate a formatted title for the current scan"
    
    "---Handle inputs---"
    if num is None:
        num = latest()
    
    " Multiple nums given"
    num = np.asarray(num).reshape(-1)
    if len(num) > 1:
        nums = num[1:]
        num = num[0]
    else:
        nums = []
        num = num[0]
    
    "---Load data---"
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    if d is None:
        return 'No Scan'
    
    "---Get metadata---"
    #keys = [x for x in d.keys() if type(d[x]) == np.ndarray ] # only keys linked to arrays
    m = d.metadata
    cmd = m.cmd # Scan command
    
    " Use scan command to determine variables of scan and title"
    HKL = m.hkl_str
    Energy = scanenergy(d)
    Temp = scantemp(d)
    psi = scanazimuth(d,'short')
    pol = scanpol(d)
    
    "---exp_title---"
    if exp_title is '':
        etitle = os.path.basename(filedir)
    else:
        etitle = exp_title
    
    "---Generate title---"
    if len(nums) < 1:
        #ttl = '{} #{} {} {} {}{}'.format(etitle,m.SRSRUN,HKL,Energy,Temp,pol).strip()
        ttl='{} #{} {} {} {} {}{}'.format(etitle,m.SRSRUN,Energy,Temp,HKL,psi,pol).strip()
        return ttl
    else:
        "Multiple run title"
        scan_range = numbers2string([num]+list(nums))
        
        " Use last scan to determine what has changed between scans"
        d2 = readscan(nums[-1])
        last_HKL = d2.metadata.hkl_str
        last_Energy = scanenergy(d)
        last_Temp = scantemp(d)
        last_psi = scanazimuth(d,'short')
        last_pol = scanpol(d)
        
        if Temp != last_Temp:
            ttl = '{} #{} {} {} {}'.format(etitle,scan_range,HKL,Energy,pol).strip()
        elif Energy != last_Energy:
            ttl = '{} #{} {} {} {}'.format(etitle,scan_range,HKL,Temp,pol).strip()
        else:
            ttl = '{} #{} {} {} {}{}'.format(etitle,scan_range,HKL,Energy,Temp,pol).strip()
        return ttl

def scanlatt(num):
    "Displays the lattice parameters of the chosen scan"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    print('---Lattice Parameters #{}---'.format(m.SRSRUN))
    print('    a = {:7.4f},     b = {:7.4f},      c = {:7.4f}'.format(m.a,m.b,m.c))
    print('alpha = {:7.2f},  beta = {:7.2f},  gamma = {:7.2f}'.format(m.alpha1,m.alpha2,m.alpha3))

def scanhkl(num):
    "Returns the average HKL of the chosen scan as a formatted string"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    if 'h' in d.keys():
        h = np.mean(d.h)
    else:
        h = m.h
    
    if 'k' in d.keys():
        k = np.mean(d.k)
    else:
        k = m.k
    
    if 'l' in d.keys():
        l = np.mean(d.l)
    else:
        l = m.l
    
    return '({0:1.3g},{1:1.3g},{2:1.3g})'.format(round(h,2)+0.0,round(k,2)+0.0,round(l,2)+0.0)

def scanazimuth(num,style=None):
    "Returns the azimuthal angle and zero reference as a formatted string"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    psi = m.psi
    h,k,l = m.azih, m.azik, m.azil
    if style in ['short','simple','s']:
        return 'psi%1.0f%1.0f%1.0f=%1.1f'%(m.azih,m.azik,m.azil,m.psi)
    return '{:7.2f} ({:1.3g},{:1.3g},{:1.3g})'.format(psi,h,k,l)

def scanpol(num,latex=False):
    "Returns the polarisation of the current scan"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    pol = ''
    if m.delta_axis_offset < 1 and np.abs(m.thp) > 10:
        if m.gam > 0.01:
            if m.stoke < 45.:
                pol = ' p-p'
                ltx = '$\pi-\pi$'
            else:
                pol = ' p-s'
                ltx = '$\pi-\sigma$'
        else:
            if m.stoke < 45.:
                pol = ' s-s'
                ltx = '$\sigma-\sigma$'
            else:
                pol = ' s-p'
                ltx = '$\sigma-\pi$'
    if latex:
        return ltx
    return pol

def scantemp(num, sensor=None, return_value=False):
    "Returns the average temperature of the chosen scan as a formatted string"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    if sensor is None:
        sensor = default_sensor
    
    m = d.metadata
    if sensor in d.keys():
        T = np.mean(getattr(d,sensor))
    elif sensor in d.metadata.keys():
        T = getattr(m,sensor)
    else:
        T = 300
    
    # temperature given as 0 if lakeshore not connected
    if T < 0.1:
        T = 300

    if return_value:
        return T    
    return '{:1.3g}K'.format(T)

def scanss(num):
    "Returns the sample slit gaps as a formatted string"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    return '{0:4.2f}x{1:4.2f} mm'.format(m.s5xgap,m.s5ygap)

def scands(num):
    "Returns the MiniMirror Pitch"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    # analyser slits s6 changed name to s7 in April 2018
    if 's7xgap' in m.keys():
        return '{0:4.2f}x{1:4.2f} mm'.format(m.s7xgap,m.s7ygap)
    else:
        return '{0:4.2f}x{1:4.2f} mm'.format(m.s6xgap,m.s6ygap)

def scanminimirrors(num):
    "Returns the minimirror position a formatted string"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    if m.m4pitch > 0.02:
        mmin = 'in'
    else:
        mmin = 'out'
    return '{} ({:4.2f} deg)'.format(mmin, m.m4pitch)

def scanenergy(num):
    "Returns the average energy of the chosen scan as a formatted string"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    if 'Energy' in d.keys():
        E = np.mean(d.Energy)
    else:
        E = m.Energy
    
    return '{:1.5g}keV'.format(E)

def scaneuler(num, mean=False):
    """
    Returns the euler angles eta chi phi mu delta gamma
    """
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    # Calculate eulerian angles
    # Taken from: /dls_sw/i16/software/gda/config/scripts/EulerianKconversion.py
    Kalpha=np.radians(50) # alpha angle on Kappa stage
    kap = np.radians(d.kap)
    alpha = -np.degrees(np.arctan( np.cos(Kalpha)*np.tan(kap/2) ))
    
    chi = -2*np.degrees(np.arcsin( np.sin(Kalpha)*np.sin(kap/2) ))
    eta = d.kth - alpha
    phi = d.kphi - alpha
    delta = d.kdelta
    gam = d.kgam
    mu = d.kmu
    if mean:
        eta = np.mean(eta)
        chi = np.mean(chi)
        phi = np.mean(phi)
        delta = np.mean(delta)
        gam = np.mean(gam)
        mu = np.mean(mu)
    return eta, chi, phi, mu, delta, gam

def scanwl(num):
    "Returns the initial wavelength of the given scan in Angstroms"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    m = d.metadata
    
    e = 1.6021733E-19;  # C  electron charge
    h = 6.62606868E-34; # Js  Plank consant
    c = 299792458;      # m/s   Speed of light
    A = 1e-10;          # m Angstrom      
    
    # Energy(keV)-Energy(eV):
    E = m.Energy*1000*e
    
    # SI: E = hc/lambda
    lam = h*c/E
    wl = lam/A
    return wl

def scantime(num):
    "Returns the time the scan file was last modified, as a datetime object"

    #date = readscan(num).metadata.date
    #time = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
    file = scanfile(num)
    time = datetime.datetime.fromtimestamp(os.stat(file).st_mtime) # last modified time
    # Add 1s for good measure
    time = time + datetime.timedelta(0,1) # days, seconds
    return time

def scanfile(num):
    "Returns the full file name of scan #num"

    try:
        return num.metadata.filename
    except:
        pass
    
    if os.path.isdir(filedir) == False: 
        print( "I can't find the directory: {}".format(filedir) )
        return None
    
    try:
        num = num.metadata.SRSRUN
    except:
        pass
    
    if num < 1: 
        if latest() is None: return None
        num = latest()+num
    
    file = os.path.join(filedir, '%i.dat' %num)
    return file

def scanimagefile(num, idx=0):
    " Return the filename of an image file in the scan"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    try:
        pilname = [s for s in d.metadata.keys() if "_path_template" in s][0]
        pilpath = getattr(d.metadata,pilname)
    except IndexError:
        print( 'Not a pilatus file!' )
        return
    
    " Load first image to get the detector size"
    tif=pilpath % d.path[idx]
    file = os.path.join(filedir,tif)
    file=file.replace('/',os.path.sep)
    return file

def scanimagefiles(num, setvarx=None):
    "Return inputs for image_gui: image_file_list, names_list, ttl"
    
    x, y, dy, varx, vary, ttl, d = getdata(num, varx=setvarx)
    file_list = [scanimagefile(d, n) for n in range(len(x))]
    name_list = ['%s = %10.5g' % (varx, v) for v in x]
    return file_list, name_list, ttl

def scanimage(num, idx=0):
    "Returns a single image array from scan #num"
    
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    file = scanimagefile(num, idx)
    #im=misc.imread(file, flatten=True) # flatten required to read zylar 16 bit files
    im = imread(file) # imageio
    return im

def nexustree(tree, name=[]):
    """
    Returns string of the complete nexus tree
        n = readnexus(12345)
        ss = nexustree(n)
        print(ss)
    """
    outstr = ''
    try:
        for branch in tree.keys():
            branchname = name + [branch]
            try:
                val = tree[branch][...]
                if np.prod(val.shape) > 1:
                    pval = val.shape
                elif np.prod(val.shape) == 1:
                    pval = val.reshape(-1)[0]
                else:
                    pval = val
                pname = '/'.join(branchname)
                outstr += '{} = {}\n'.format(pname, pval)
            except AttributeError:
                outstr += nexustree(tree[branch], branchname)
    except AttributeError:
        return outstr
    return outstr

def nexussearch(find, tree, name=[], Whole=False, Case=False):
    """
    Search nexus tree for named value
        addresses, values = nexussearch('Energy', Whole=False, Case=False)

    Returns list of nexus addresses and the associated value
    """
    outnames = []
    outvals = []
    try:
        for branch in tree.keys():
            branchname = name + [branch]
            if (
                find == "*" or 
                Whole and Case and find == branch or 
                Whole and not Case and find.lower() == branch.lower() or 
                not Whole and Case and find in branch or 
                not Whole and not Case and find.lower() in branch.lower()
                ):
                try:
                    val = tree[branch][...]
                    if np.prod(val.shape) > 1:
                        pval = val.shape
                    else:
                        pval = val.reshape(-1)[0]
                    pname = '/'.join(branchname)
                    outnames += [pname]
                    outvals += [pval]
                    continue
                except AttributeError:
                    pass
            pname, pval = nexussearch(find, tree[branch], branchname)
            outnames += pname
            outvals += pval
    except AttributeError:
        return outnames, outvals
    return outnames, outvals

def prend(start=0,end=None):
    "Calculate the end time of a run, return str"

    outstr = ''
    if end is None and ( start == 0 or start == latest() ):
        # End of current scan is required
        st = readscan(start)
        m = st.metadata
        t1=datetime.datetime(m.Year,m.Month,m.Day,m.Hours,m.Minutes,m.Seconds,0)
        t2=datetime.datetime.now()
        
        # Number of points completed
        Ncomplete = len(st.TimeSec)
        
        # Number of points in scan
        cmd = m.cmd.split()
        scantype = cmd[0]
        varx = cmd[1]
        if scantype == 'scancn':
            # e.g. scancn x 0.1 10 t 1
            scanlen = int(cmd[3])
        elif scantype == 'scan':
            # e.g. scan x 1 10 0.1 t 1
            scan_start = float(cmd[2])
            scan_end = float(cmd[3])
            scan_step = float(cmd[4])
            scanlen = len(np.arange(scan_start,scan_end,scan_step))
        
        nrem = scanlen-Ncomplete
        
        tdif = t2-t1
        
        # Predict end time
        tperrun = tdif / Ncomplete
        trem = nrem*tperrun
        tend = t2 + trem
        
        outstr += 'Scan number: #{}\n'.format(latest())
        outstr += 'Scan Started: {}\n'.format(t1)
        outstr += 'Points complete: {} / {}\n'.format(Ncomplete,scanlen)
        outstr += 'Time per point: {}\n'.format(tperrun)
        outstr += 'Still to go: {} ({})\n'.format(nrem,trem)
        outstr += 'Scan will end: {}\n'.format(tend)
        return outstr
    
    st = readscan(start)
    nd = readscan(end)
    
    m = st.metadata
    t1=datetime.datetime(m.Year,m.Month,m.Day,m.Hours,m.Minutes,m.Seconds,0)
    
    """ If run has already finished """
    if nd is not None:
        m = nd.metadata
        t2=datetime.datetime(m.Year,m.Month,m.Day,m.Hours,m.Minutes,m.Seconds,0)
        
        tdif = t2-t1
        
        outstr += 'Run started: {}\n'.format(t1)
        outstr += 'Run ended: {}\n'.format(t2)
        outstr += 'Run took: {}\n'.format(tdif)
        return outstr
    
    """ If run is still going """
    # Get current scan
    cur = latest()
    cr = readscan(cur)
    
    m = cr.metadata
    t2=datetime.datetime(m.Year,m.Month,m.Day,m.Hours,m.Minutes,m.Seconds,0)
    
    tdif = t2-t1
    ndif = cur-start
    nrem = end-cur
    
    # Predict end time
    tperrun = tdif / ndif
    trem = nrem*tperrun
    tend = t2 + trem
    
    outstr += 'Run Started: {}\n'.format(t1)
    outstr += 'Latest scan: #{}\n'.format(cur)
    outstr += 'Time per scan: {}\n'.format(tperrun)
    outstr += 'Runs completed: {} ({})\n'.format(ndif,tdif)
    outstr += 'Still to go: {} ({})\n'.format(nrem,trem)
    outstr += 'Run will end: {}\n'.format(tend)
    return outstr

def findfile(num,topdir=None):
    "Find the experiment directory of a scan number"
    
    if topdir == None:
        topdir = filedir
        
    file = '{}.dat'.format(num)
    
    for root, subFolders, files in os.walk(topdir):
        # Remove pilatus directories for speed
        for n in subFolders[:]:
            if '-files' in n: subFolders.remove(n) # ####-pilatus100k-files
            if '-data' in n: subFolders.remove(n) # snapped-data
        #print( subFolders )
        #print( files )
        #a=raw_input('press enter')
        #if a == 'exit': break
        
        if file in files:
            print( 'Scan {} was in directory: {}'.format(num,root) )
            return root

def polflip(sigsig, sigpi, fit='Gauss', output=False, plot=False):
    "Calculate flipping ratio ect."
    # Flipping ratio - measured on the straight through beam
    x1,y1,dy1,varx1,vary1,ttl1,d1 = getdata(sigsig) # scan 1
    x2,y2,dy2,varx2,vary2,ttl2,d2 = getdata(sigpi) # scan 2 
    
    # Find polarisation of scan1
    if 's-s' in ttl1:
        nprx = x1
        npry = y1
        npr = 'ss'
        nprt = ttl1
        nprN = sigsig
    elif 'p-p' in ttl1: 
        nprx = x1
        npry = y1
        npr = 'pp'
        nprt = ttl1
        nprN = sigsig
    elif 's-p' in ttl1:
        prx = x1
        pry = y1
        pr = 'sp'
        prt = ttl1
        prN = sigsig
    elif 'p-s' in ttl1:
        prx = x1
        pry = y1
        pr = 'ps'
        prt = ttl1
        prN = sigsig
    else:
        print( 'Run {} is not polarized'.format(sigsig) )
        return
    
    # Find polarisation of scan2
    if 's-s' in ttl2:
        nprx = x2
        npry = y2
        npr = 'ss'
        nprt = ttl2
        nprN = sigpi
    elif 'p-p' in ttl2: 
        nprx = x2
        npry = y2
        npr = 'pp'
        nprt = ttl2
        nprN = sigpi
    elif 's-p' in ttl2:
        prx = x2
        pry = y2
        pr = 'sp'
        prt = ttl2
        prN = sigpi
    elif 'p-s' in ttl2:
        prx = x2
        pry = y2
        pr = 'ps'
        prt = ttl2
        prN = sigpi
    else:
        print( 'Run {} is not polarized'.format(sigpi) )
        return
    
    
    nprfit,nprerr = peakfit(nprx,npry,type=fit)
    prfit,prerr = peakfit(prx,pry,type=fit)
    
    npr_ara = nprfit['Area']
    pr_ara = prfit['Area']
    #npr_amp,npr_cen,npr_wid,npr_bkg,npr_ara,npr_damp,npr_dcen,npr_dwid,npr_dbkg,npr_dara,npr_yfit = lorzfit(nprx,npry)
     #   pr_amp,pr_cen,pr_wid,pr_bkg,pr_ara,pr_damp,pr_dcen,pr_dwid,pr_dbkg,pr_dara,pr_yfit = lorzfit(prx,pry)

    FR = pr_ara / npr_ara
    CT = (pr_ara-npr_ara)/(pr_ara+npr_ara)
    
    if output: return FR,CT
    
    print( '--------------Polarisation Analysis-----------------' )
    print( 'Fit type: {}'.format(fit) )
    print( '{}: {}'.format(npr,nprt) )
    print( '{}: {}'.format(pr,prt) )
    print( '      {0} / {1} = {2:6.3f}  ({3:3.1f}%)'.format(pr,npr,FR,FR*100) )
    print( '{0}-{1} / {0}+{1} = {2:6.3f}'.format(pr,npr,CT) )
    print()
    
    if plot:
        plotscan([nprN,prN],fit=fit,fits=True)
        plt.legend([r'#{}: $\sigma-\sigma$'.format(nprN),r'#{}: $\sigma-\pi$'.format(prN)])
        # plt.show()
        # saveplot('POLFLIP '+saveable(ttl1))

def polenergy(sigsig, sigpi, background=None, vary='', bkg_scale=None, flipping_ratio=None, low_points=5, save=False):
    "Create Plot of energy-polarisation scans and calculate the subtraction"
    
    " Get the signal data - measured at a resonant feature in ss and sp"
    x1,y1,dy1,varx1,vary1,ttl1,d1 = getdata(sigsig,vary=vary) # scan 1
    x2,y2,dy2,varx2,vary2,ttl2,d2 = getdata(sigpi,vary=vary) # scan 2
    " Get the background data - measured away from the Bragg peak in sp"
    if background is None:
        xb,yb,dyb=1*x1,0*x1,0*x1
    else:
        xb,yb,dyb,varxb,varyb,ttlb,db = getdata(background,vary=vary) # Background  
    
    " Get metadata"
    m = d2.metadata
    cmd = m.cmd_short # Scan command
    hkl = m.hkl_str
    T = m.temperature
    sampsl = scanss(d2)
    detsl = scands(d2)
    atten1 = '{0:1.0f}'.format(m.Atten)
    psival = '{:1.3g}'.format(m.psi)
    azir = r'$\angle$ ({:1.0f}{:1.0f}{:1.0f})'.format(m.azih,m.azik,m.azil)
    
    ttl = '{} {} $\Psi$={} {}\nss={} ds={} atten={}'.format(hkl,T,psival,azir,sampsl,detsl,atten1)
    
    " Use the title to confirm the polarisation, so that pi-pi and pi-sigma are also caught"
    " npr = non-polarisation rotation = ss, pp"
    " pr = polarisation rotation = sp, ps"
    
    " Find polarisation of scan1"
    if 's-s' in ttl1:
        nprx = x1
        npry = y1
        npr = 'ss'
        nprt = ttl1
        nprN = sigsig
        nprlab = '$\sigma\sigma$'
    elif 'p-p' in ttl1: 
        nprx = x1
        npry = y1
        npr = 'pp'
        nprt = ttl1
        nprN = sigsig
        nprlab = '$\pi\pi$'
    elif 's-p' in ttl1:
        prx = x1
        pry = y1
        pr = 'sp'
        prt = ttl1
        prN = sigsig
        prlab = '$\sigma\pi'
    elif 'p-s' in ttl1:
        prx = x1
        pry = y1
        pr = 'ps'
        prt = ttl1
        prN = sigsig
        prlab = '$\pi\sigma$'
    else:
        print( 'Run {} is not polarized'.format(sigsig) )
        return
    
    # Find polarisation of scan2
    if 's-s' in ttl2:
        nprx = x2
        npry = y2
        npr = 'ss'
        nprt = ttl2
        nprN = sigpi
        nprlab = '$\sigma\sigma$'
    elif 'p-p' in ttl2: 
        nprx = x2
        npry = y2
        npr = 'pp'
        nprt = ttl2
        nprN = sigpi
        nprlab = '$\pi\pi$'
    elif 's-p' in ttl2:
        prx = x2
        pry = y2
        pr = 'sp'
        prt = ttl2
        prN = sigpi
        prlab = '$\sigma\pi$'
    elif 'p-s' in ttl2:
        prx = x2
        pry = y2
        pr = 'ps'
        prt = ttl2
        prN = sigpi
        prlab = '$\pi\sigma$'
    else:
        print( 'Run {} is not polarized'.format(sigpi) )
        return
    
    " npr = non-polarisation rotation = ss, pp"
    " pr = polarisation rotation = sp, ps"
    
    if len(pry) != len(yb):
        yb = np.interp(prx,xb,yb)
        xb = prx
        
    if len(pry) != len(npry):
        npry = np.interp(prx,nprx,npry)
        nprx= prx
    
    " Determine the background scale"
    if bkg_scale is None:
        bkg_scale = 1
    
    " Determine the Flipping ratio from the low energy data"
    if flipping_ratio is None:
        flipping_ratio = sum(npry[:low_points])/sum(pry[:low_points]) 
    
    " Calculate the subtractions"
    DIF = pry - yb*bkg_scale - npry/flipping_ratio
    
    " Labels"
    diflab = '{}-{}/{:1.3g}-bkg/{:1.3g}'.format(prlab,nprlab,flipping_ratio,bkg_scale)
    nprlab = '{}/{:1.3g} #{}'.format(nprlab,flipping_ratio,nprN)
    prlab = '{} #{}'.format(prlab,prN)
    bkglab = 'BKG/{:1.3g} #{}'.format(bkg_scale,background)
    
    " Create Plot"
    fig = plt.figure(figsize=[10,8], dpi=fig_dpi)
    
    plt.plot(nprx,npry/flipping_ratio,'-ob',linewidth=2,label=nprlab)
    plt.plot(prx,pry,'-og',linewidth=2,label=prlab)
    plt.plot(xb,yb*bkg_scale,'-ok',linewidth=2,label=bkglab)
    plt.plot(prx,DIF,'-or',linewidth=2,label=diflab)
    
    plt.legend(loc=0, fontsize=16)
    plt.ylabel(vary2, fontsize=18)
    plt.xlabel('Energy (keV)', fontsize=18)
    #plttl = ttl2+'\n'+cmd+'\npsi={}, ss ={}, ds ={}, atten = {}'.format(sampsl,detsl,atten1)
    plt.title(ttl, fontsize=14)
    plt.show()
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('PolEng_{}_{}_{}.png'.format(hkl,T,psival))

def scanabscor(num=0,u=1,eta_offset=0.0,chi_offset=0.0):
    """
    Calculate absorption correction
     A = scanabscor(num,u)
     Icorrected = Iexp/A
     See abscor for more details
    """
    
    # Get data
    try: 
        d = readscan(num)
    except:
        d = num
    
    # Calculate eulerian angles
    eta, chi, phi, mu, delta, gam = scaneuler(d)
    A = [abscor(eta[n],chi[n],delta[n],u=u) for n in range(len(eta))]
    return A

def scanvolcor(num=0,u=1,eta_offset=0.0,chi_offset=0.0):
    """
    Calculate absorption correction
     A = scanvolcor(num,u)
     Icorrected = Iexp/A
     See abscor for more details
    """
    
    # Get data
    try: 
        d = readscan(num)
    except:
        d = num
    
    u = np.asarray(u).reshape(-1)
    
    # Calculate eulerian angles
    eta, chi, phi, mu, delta, gam = scaneuler(d)
    
    if len(u) == 1:
        u = np.repeat(u, len(eta))
    
    A = [volcor(eta[n],chi[n],delta[n],u=u[n]) for n in range(len(eta))]
    return A

def metaprint(d1,d2=None):
    """ 
    Print metadata of a scan:
        metadata(d)
    OR compare metadata of two scans:
        metadata(d1,d2)
    """
    
    try: 
        d1 = readscan(d1)
    except:
        'Entered data'
    
    if d2 is None:
        for k in d1.metadata.keys():
            if len(str(d1.metadata[k])) > 1000:
                showstr = '...to long...'
            else:
                showstr = str(d1.metadata[k])
            print( '{:>20} : {:<20}'.format(k,showstr) )
    else:
        try:
            d2= readscan(d2)
        except:
            'Entered data'
        
        # Compare meta data scans
        print( 'Key                  : #{:<10}: #{:<10}'.format(d1.metadata.SRSRUN,d2.metadata.SRSRUN) )
        for k in d1.metadata.keys():
            try:
                m1 = d1.metadata[k]
                m2 = d2.metadata[k]
                diff = ''
                if m1 != m2: 
                    diff = '***'
                    #print( m1,' does not equal ',m2 )
                if len(str(m1)) > 1000:
                    showstr1 = '...to long...'
                else:
                    showstr1 = str(m1)
                if len(str(m2)) > 1000:
                    showstr2 = '...to long...'
                else:
                    showstr2 = str(m2)
                print( '{:>20} : {:10} : {:<10} {}'.format(k,showstr1,showstr2,diff) )
            except:
                print( '{} does not exist in #{}'.format(k,d2.metadata.SRSRUN) )

def checkpeaks(num,test=1,vary=''):
    "Check multiple scans for peaks"
    
    # Turn num into a list
    try:
        N = len(num)
    except:
        num=[num]
    
    for n in num:
        x,y,dy = getdata(n,vary=vary)[:3]
        rat = ispeak(y,dy,return_rat=True)
        print( '{} {:8.2f} {}'.format(n,rat,rat>test) )

def savescan(num=None,varx='',vary='',norm=True,abscor=None):
    "Save scan as .dat file"
    
    # load data
    x,y,dy,varx,vary,ttl,d = getdata(num,varx,vary,norm,abscor)
    
    # save data to file
    savefile = os.path.join(savedir, '{}_{}.dat'.format(num,saveable(vary)))
    head = '{}\n{}\n{}, {}, {}'.format(ttl,d.metadata.cmd_short,varx,vary,'error_'+vary)
    np.savetxt(savefile,np.transpose([x,y,dy]),header=head)
    print( 'Scan #{} has been saved as {}'.format(num,savefile) )

def loadscan(num,vary=None):
    "Load a scan.dat file saved with savescan from the analysis folder"
    
    if vary is None:
        # Get all data files in folder
        ls=glob.glob(savedir+os.path.sep+'*.dat')
        ls = np.sort(ls)
        for file in ls:
            if str(num) in file:
                filename = file
                break
    else:
        filename = os.path.join(savedir, '{}_{}.dat'.format(num,saveable(vary)))
    try:
        print('Reading ',filename)
    except:
        print('There is no file for scan #',num)
    
    # Get Header Data
    with open(filename,'r') as ff:
        # Line 1 = title
        ttl = ff.readline().strip('# ').strip('\n')
        # Line 2 = scan command
        cmd = ff.readline().strip('# ')
        # Line 3 = variable names
        varx,vary,dvary = ff.readline().strip('# ').split(', ')
        
    # Get x, y, dy data
    data = np.loadtxt(filename)
    if data.shape[0] == 3:
        # old files - saved as 3 rows
        x,y,dy = data
    else:
        # new files - saved as 3 columns
        x,y,dy = data.T
    
    return x,y,dy,varx,vary,ttl

def create_analysis_file(scans,depvar='Ta',vary='',varx='',fit_type = 'pVoight',bkg_type='flat',peaktest=1,
                  norm=True,abscor=None,plot='all',show_fits=True,mask_cmd=None,estvals=None,xrange=None,sortdep=True,
                  Nloop=10, Binit=1e-5, Tinc=2, change_factor=0.5, converge_max=100, min_change=0.01,
                  savePLOT=False,saveFIT=False):
    """
    Creates a new python analysis script in the anaysis folder
    The script will contain all the python code to call the Py16progs module, load the scan data and
    run an automated fit routine.
    Inputs: - See fit_scans for a description of the inputs
    """
    
    # Turn depvar into a list (encase of multiple depvar values)
    if type(depvar) is str:
        depvar = [depvar]
    Ndep = len(depvar)
    
    ttl = saveable(exp_title)
    filename = savedir + '/I16_Analysis_{}.py'.format('_'.join([ttl,depvar[0]]))
    
    # Add str quotes for printing, but not if list
    if type(vary) is str:
        vary = '\''+vary+'\''
    if type(plot) is str:
        plot = '\''+plot+'\''
    
    # Check file doesn't already exist, if it does, create new file name
    n = 1
    while os.path.isfile(filename):
        n += 1
        filename = '_'.join([ttl, depvar[0], str(n)])
        filename = 'I16_Analysis_{}.py'.format(filename)
        filename = os.path.join(savedir, filename)
        #filename = savedir + '/I16_Analysis_{}_{}_{}.py'.format(ttl,depvar[0],n)
    
    with open(filename,'w') as f:
        # Write comments at top
        f.write('# Diamond I16 Analysis Script\n')
        f.write('# {}\n'.format(exp_title))
        f.write('# Analysis of {} dependence\n'.format(' '.join(depvar)))
        f.write('# \n')
        f.write('# By User\n')
        f.write('# {}\n\n'.format(datetime.datetime.strftime(datetime.datetime.now(),'%d/%m/%Y')))
        
        # Import stuff
        f.write('import sys,os\n')
        f.write('import numpy as np\n')
        #f.write('import scisoftpy as dnp # Make sure this is in your python path\n')
        f.write('import matplotlib.pyplot as plt # Plotting\n')
        f.write('from mpl_toolkits.mplot3d import Axes3D # 3D plotting\n\n')
        
        
        f.write('# Load Py16progs\n')
        f.write('sys.path.insert(0,\'{}\') # location of Py16progs\n'.format(os.path.dirname(__file__)))
        f.write('import Py16progs as dp\n\n')
        
        # Data and save directories
        f.write('# Current Directory\n')
        f.write('cf = os.path.dirname(__file__)\n\n')
        f.write('# Directory to load data from\n')
        f.write('dp.filedir = \'{}\' \n\n'.format(filedir))
        f.write('# Directory to save files to\n')
        f.write('dp.savedir=\'{}\' \n\n'.format(savedir))
        f.write('# Update default save location for exported plots\n')
        f.write('plt.rcParams["savefig.directory"] = dp.savedir \n\n')
        
        # Normalisation Options and Pilatus Centre
        f.write('# Experiment Parameters\n')
        f.write('dp.exp_ring_current = {}\n'.format(exp_ring_current))
        f.write('dp.exp_monitor = {}\n'.format(exp_monitor))
        f.write('dp.normby = \'{}\'\n'.format(normby))
        f.write('dp.pil_centre = {}\n'.format(pil_centre))
        f.write('dp.peakregion = {}\n'.format(peakregion))
        f.write('dp.exp_title = \'{}\'\n\n'.format(exp_title))
        
        # Scan numbers
        f.write('# Scan numbers\n')
        try:
            if np.max(np.diff(np.diff(scans))) == 0:
                f.write('scans = range({},{},{})\n\n'.format(scans[0],scans[-1]+1,scans[1]-scans[0]))
            else:
                f.write('scans = {} \n\n'.format(str(list(scans))))
        except:
            f.write('scans = {} \n\n'.format(str(list(scans))))
        
        if mask_cmd is not None:
            try:
                mask_cmd = str( mask_cmd.tolist() )
            except:
                mask_cmd = str( list(mask_cmd) )
        
        # Fit
        f.write('# Fitting\n')
        f.write('mask_cmd = {}\n'.format(mask_cmd))
        f.write('estvals = {}\n'.format(str(estvals)))
        f.write('fitopt = dict(depvar={},vary={},varx=\'{}\',fit_type = \'{}\',bkg_type=\'{}\',peaktest={},\n'.format(depvar,vary,varx,fit_type,bkg_type,peaktest))
        f.write('              norm={},abscor={},plot={},show_fits={},mask_cmd=mask_cmd,estvals=estvals,xrange={},sortdep={},\n'.format(norm,abscor,plot,show_fits,xrange,sortdep))
        f.write('              Nloop={}, Binit={}, Tinc={}, change_factor={}, converge_max={}, min_change={},\n'.format(Nloop,Binit,Tinc,change_factor,converge_max,min_change))
        f.write('              savePLOT={},saveFIT={})\n\n'.format(savePLOT,saveFIT))
        f.write('fit,err = dp.fit_scans(scans,**fitopt)\n\n')
        # Load
        f.write('# Load fitted data:\n')
        f.write('#fit,err = dp.load_fits([{},{}], depvar={}, vary={}, fit_type=\'{}\')\n\n'.format(scans[0],scans[-1],depvar,vary,fit_type))
    
    print( 'New Analysis file written to ',filename )

def example_script(scanno=None):
    """
    Returns example script as string
    """
    if scanno is None:
        scanno = latest()
    time = datetime.datetime.now()

    outstr = '# Example Python Scipt for I16_Data_Viewer\n'
    outstr += '# {}\n\n'.format(time)
    outstr += '#pp.filedir = ''{}''\n'.format(filedir)
    outstr += '#pp.savedir = ''{}''\n\n'.format(savedir)

    outstr += 'd = pp.readscan({}) # read scan file (e.g. d.eta, d.metadata.psi,...)\n\n'.format(scanno)

    outstr += 'd = pp.plotscan({}) # automatic plot, see help(pp.plotscan)\n\n'.format(scanno)

    outstr += 'x,y,dy, varx, vary, ttl, d = pp.getdata({}, varx=\'Auto\', vary=\'Auto\', norm=True) # automatic axes and normalisation\n\n'.format(scanno)

    outstr += 'pp.newplot(x, y, label={})\npp.labels(ttl, varx, vary, legend=True)\n\n'.format(scanno)

    outstr += 'x2,y2,dy2, varx2, vary2, ttl2, d2 = pp.getdata({}, varx=\'Auto\', vary=\'Auto\', norm=True) \n'.format(scanno-1)    
    outstr += 'pp.multiplot([x,x2], [y,y2], labels=[{},{}])\npp.labels(ttl, varx, vary, legend=True)\n\n'.format(scanno,scanno-1)
    return outstr

def get_all_scannos():
    "Returns the scan numbers available in filedir"
    
    # Get all data files in folder
    ls=glob.glob(filedir+os.path.sep+'*.dat')
    ls = np.sort(ls)
    
    # Convert to int
    scannos = [np.int(os.path.split(file)[-1][:-4]) for file in ls]
    return scannos

def findscans(scannos=None, hkl=None, energy=None, temperature=None, scan_type=None, detector=None, stokes=None, psi=None, endtime=None, hours_before=None):
    """
    Find scans with certain properties within a given set of scan numbers
      scans = findscans(scannos, hkl=[0,0,2])
      - returns scans with metadata recorded at [0,0,2]

     - endtime: datetime.datetime object
     - hours_before - finds scans between endtime-hours_before and endtime
    """

    hkl_tol = 0.05
    energy_tol = 0.001
    temp_tol = 0.5
    stokes_tol = 0.1
    psi_tol = 0.1

    if scannos is None:
        scannos = get_all_scannos()

    if hkl:
        hkl = np.asarray(hkl)

    if endtime is None:
        endtime = datetime.datetime.now() + datetime.timedelta(seconds=5)
    if hours_before:
        hours_before = datetime.timedelta(hours=hours_before)
        starttime = endtime - hours_before

    outscans = []
    for scan in scannos:
        d = readscan(scan)
        if d is None: continue
        m = d.metadata

        findall = []
        if hkl is not None:
            scanhkl = np.array([m.h, m.k, m.l])
            diff = np.sqrt(np.sum( (scanhkl - hkl)**2 ))
            if diff < hkl_tol:
                findall += [True]
            else:
                findall += [False]

        if energy:
            if np.abs(m.Energy - energy) < energy_tol:
                findall += [True]
            else:
                findall += [False]

        if temperature:
            T = scantemp(d, return_value=True)
            if np.abs(T - temperature) < temp_tol:
                findall += [True]
            else:
                findall += [False]

        if scan_type:
            varx = auto_varx(d)
            if scan_type == varx:
                findall += [True]
            else:
                findall += [False]

        if detector:
            cmd = m.cmd
            if detector in cmd:
                findall += [True]
            else:
                findall += [False]

        if stokes:
            if np.abs(m.stokes - stokes) < stokes_tol:
                findall += [True]
            else:
                findall += [False]

        if psi:
            if np.abs(m.psi - psi) < psi_tol:
                findall += [True]
            else:
                findall += [False]

        if hours_before is not None:
            scandate = scantime(scan)
            if starttime < scandate < endtime:
                findall += [True]
            else:
                findall += [False]

        if np.all(findall):
            outscans += [scan]
    return outscans



"----------------------Previous Experiment Functions----------------------"


def exp_list_add():
    "Creates/appends files in the system Temp directory to store previously used experiment folders"
    
    exp_dir = '{}\n'.format(filedir)
    
    if filedir in exp_list_get():
        return
    
    with open(os.path.join(tmpdir,'Py16_experiment_directories.txt'),'a') as f:
        f.write(exp_dir)
    print('Added {} to experiment list'.format(filedir))

def sav_list_add():
    "Creates/appends files in the system Temp directory to store previously used save folder"
    
    sav_dir = '{}\n'.format(savedir)
    
    if savedir in sav_list_get():
        return
    
    with open(os.path.join(tmpdir,'Py16_analysis_directories.txt'),'a') as f:
        f.write(sav_dir)
    print('Added {} to analysis list'.format(savedir))

def exp_list_get():
    "Loads a list of directories of previous experiments"
    
    if not os.path.isfile(os.path.join(tmpdir,'Py16_experiment_directories.txt')):
        return []
    
    with open(os.path.join(tmpdir,'Py16_experiment_directories.txt'),'r') as f:
        exp_list = f.readlines()
        exp_list = [x.strip() for x in exp_list] # remove \n
    return exp_list

def sav_list_get():
    "Loads a list of analysis directories used previously"
    
    if not os.path.isfile(os.path.join(tmpdir,'Py16_analysis_directories.txt')):
        return []
    
    with open(os.path.join(tmpdir,'Py16_analysis_directories.txt'),'r') as f:
        sav_list = f.readlines()
        sav_list = [x.strip() for x in sav_list] # remove \n
    return sav_list

def recall_last_exp():
    "Loads the last saved experiment & analysis directories from the temp directory and sets filedir and savedir"
    
    exp_list = exp_list_get()
    sav_list = sav_list_get()
    if len(exp_list) > 0:
        global filedir
        print('Recalling previous experiment folder: {}'.format(exp_list[-1]))
        filedir = exp_list[-1]
    if len(sav_list) > 0:
        global savedir
        print('Recalling previous analysis folder: {}'.format(sav_list[-1]))
        savedir = sav_list[-1]

def clear_prev_exp():
    "Clear the previous experiment + analysis directories files"
    
    with open(os.path.join(tmpdir,'Py16_experiment_directories.txt'),'w') as f:
        f.write('')
    with open(os.path.join(tmpdir,'Py16_analysis_directories.txt'),'w') as f:
        f.write('')

def exp_parameters_save(saveto=None):
    """
    Save experimental parameters in savedir directory
    """
    
    if saveto is None:
        saveto = savedir
    
    if os.path.isdir(saveto) is False:
        print('This directory doesn\'t exist')
        return
    
    savefile = os.path.join(saveto, 'Py16_parameters.txt')
    txt = ''
    txt += 'exp_title %s\n' %exp_title
    txt += 'filedir %s\n' %filedir
    txt += 'savedir %s\n' %savedir
    txt += 'datfile_format %s\n' %datfile_format
    txt += 'exp_ring_current %6.2f\n' %exp_ring_current
    txt += 'exp_monitor %6.2f\n' %exp_monitor
    txt += 'normby %s\n' %normby
    txt += 'dont_normalise %s\n' %(' '.join(dont_normalise))
    txt += 'detectors %s\n' %(' '.join(detectors))
    txt += 'default_sensor %s\n' %default_sensor
    txt += 'pil_centre %3.0f %3.0f\n' %(pil_centre[0],pil_centre[1])
    txt += 'hot_pixel %i\n' %hot_pixel
    txt += 'peakregion %i %i %i %i\n' %(peakregion[0],peakregion[1],peakregion[2],peakregion[3])
    txt += 'pilpara %f %f %f %f %f %f %f %f\n' %(pilpara[0],pilpara[1],pilpara[2],pilpara[3],pilpara[4],pilpara[5],pilpara[6],pilpara[7])
    txt += 'pil_max_size %f\n' %pil_max_size
    txt += 'plot_colors %s\n' %(' '.join(plot_colors))
    
    with open(savefile,'w') as file:
        file.write(txt)
    
    print('Py16 parameters written to %s' %(savefile))

def exp_parameters_load(loadfrom=None):
    """
    Load Py16 parameters from file
    """
    
    global filedir,savedir,datfile_format,exp_ring_current,exp_monitor,normby,dont_normalise
    global detectors,default_sensor,pil_centre,hot_pixel,peakregion,pilpara,pil_max_size,plot_colors
    
    if loadfrom is None:
        loadfrom = savedir
    
    savefile = os.path.join(loadfrom, 'Py16_parameters.txt')
    
    if not os.path.isfile(savefile):
        print('No Parameter file in this directory!')
        return
    
    print('Loading Py16 parameters from %s' %(savefile))
    
    with open(savefile) as file:
        lines = file.readlines()
    
    params = {}
    for line in lines:
        split_line = line.split()
        params[split_line[0]] = split_line[1:]
    
    exp_title = ' '.join(params['exp_title'])
    #filedir = ' '.join(params['filedir'])
    #savedir = ' '.join(params['savedir'])
    datfile_format = ' '.join(params['datfile_format'])
    exp_ring_current = float(params['exp_ring_current'][0])
    exp_monitor = float(params['exp_monitor'][0])
    normby = params['normby'][0]
    dont_normalise = params['dont_normalise']
    detectors = params['detectors']
    default_sensor = params['default_sensor'][0]
    pil_centre = [int(x) for x in params['pil_centre']]
    hot_pixel = int(params['hot_pixel'][0])
    peakregion=[int(x) for x in params['peakregion']]
    pilpara=[float(x) for x in params['pilpara']]
    #dead_pixel_func = np.median # Define how to choose the replaced intensity for hot/broken pixels 
    pil_max_size = float(params['pil_max_size'][0])
    #pilatus_dead_pixels = np.zeros(shape=(0,2),dtype=int) # Blank
    plot_colors = params['plot_colors']
    

"----------------------------Analysis Functions---------------------------"


def fit_scans(scans,depvar='Ta',vary='',varx='',fit_type = 'pVoight',bkg_type='flat',peaktest=1,
                  norm=True,abscor=None,plot='all',show_fits=True,mask_cmd=None,estvals=None,xrange=None,sortdep=True,
                  Nloop=10, Binit=1e-5, Tinc=2, change_factor=0.5, converge_max=100, min_change=0.01,
                  savePLOT=False,saveFIT=False):
    """ 
     Automated routine to fit peaks and plot results within runs of multiple scans dependent  
     on another variable, such as temperature, energy or azimuthal depdences.
     
         fit_scans(scans,depvar='Ta',vary='',fit_type = 'pVoight',plot='all')
     OR
         val,err = fit_scans(scans,depvar='Ta',vary='',fit_type = 'pVoight',plot=None)
     
     INPUTS:
             scans : list/array of scan numbers to use
           depvar : ('Ta')    : Independent variable, e.g. 'Ta','Energy','psi'
             vary : ('')      : Dependent variable e.g. 'APD', 'roi1_sum', 'peakroi' or '' for automatic
             varx : ('')      : Independent axis variable e.g. 'eta','chi','delta', or '' for automatic
         fit_type : ('pVoight): Type of fit to perform e.g. 'Simple','pVoight','Lorentz','Gauss'
         bkg_type : ('flat')  : Type of background to use e.g. 'flat','slope','step'
         peaktest : (10)      : Parameter to determine whether a peak exists, see help(ispeak) 
             norm : (True)    : Data normalised by time, transmission + ring current
           abscor : (None)    : Absorption coefficient for material, None switches off correction
             plot : ('all')   : Create plot of fits at end. 'all' plots all paramters, 'int' plot integrated area, None doesn't plot
        show_fits : (False)   : Generate plots for all scans and show fits (True/False)
         mask_cmd : (None)    : Define regions to remove, e.g. [ ['x<-1','x>1','np.abs(x-34)>0.1'], [x<-1.2','x>1.2'] ]
          estvals : (None)    : Define initial estimates of peak parameters, e.g. [height,cen,wid,bkg,slope]
           xrange : (None)    : Give the plots a specific xrange e.g. [0,300]
          sortdep : (True)    : Sort the final result by the depdent variable to give better plots (True/False)
            Nloop : (100)     : Number of RMC iterations of fit to perform, 0 to turn off RMC fitting
            Binit : (1e-5)    : Initial inverse-temperature for RMC
             Tinc : (2)       : RMC temperature increase factor
    change_factor : (0.5)     : RMC - width of normal distribution about which the estimate parameters are randomly varied
     converge_max : (100)     : RMC is converged when converge_max is reached
       min_change : (0.01)    : RMC - converge value in increased by one when the relative change is less than this
         savePLOT : (None)    : Save a jpg of the plot in the default directory (True/'Name.jpg'/None)
          saveFIT : (None)    : Save a text file of the fits in the default directory (True/'Name.txt'/None)
        
     OUTPUTS:
              val : Fitted values:
                      val['Scan Number']
                    val[depvar]
                    val['Peak Height']
                    val['Peak Centre']
                    val['FWHM']
                    val['Lorz frac']
                    val['Background']
                    val['Area']
                    val['CHI2 per dof']
              err : Errors on fitted values
                      err['Scan Number']
                    err['Peak Height']
                    err['Peak Centre']
                    err['FWHM']
                    err['Lorz frac']
                    err['Background']
                    err['Area']
                    out['CHI2 per dof']
    
    NOTES:
         > for psi scans, angles will automatically be placed between 0 and 360 degrees, unless xrange is set
         > if data has been saved with saveFIT=True, fitted values can be regained using:
             val,err = load_fits(scans,depvar='Ta',fit_type='pVoight')
    
    PLOTS:
        "plot" keyword accepts the following as a single string or list of strings:
            > 'all' - a single figure with subplots of integrated intensity, width, centre, background etc.
            > 'int' - a figure with a single axis showing the integrated intensity
            > 'wid' - a figure with a single axis showing the fitted FWHM
            > 'cen' - a figure with a single axis showing the fitted centre
            > 'surface' - *2D fits only!* plots integrated intensities vs x and y data. scans should vary each time in y and every n times in x
    
    MASKS:
        > Masks define regions of each scan to remove
        > input "mask_cmd" takes string arguments that will be evaluated 
        > These string arguments should be boolean statements about the variables x, y or dy
        > Any TRUE values will be removed
        > mask_cmd requires a list of lists as an input, where the each scan has its own list of string arguments
        > If only 1 item is given, it will be used in all scans
        > E.G.
        >    mask_cmd = [['x>1']] # will remove values in all scans where x>1
        >    mask_cmd = [[np.abs(x-1.4)<0.1]] # will remove values in all scans within 0.1 of 1.4
        >    mask_cmd = [['x>1','x<-1']] # will remove values in all scans where x>1 and where x<-1
        >    mask_cmd = [['x>1'],['x>2'],['x>1'],...] # Remove different regions in different scans, MUST be same length as scans
    
    ESTIMATES:
        > Estimate the initial parameters of each peak
        > input "estvals" takes list of lists, where each list contains initial parameters for each scan
        > E.G. For gaussian + slope:
        >    estvals = [height,cen,wid,bkg,slope] # same for all scans
        >    estvals = [[height1,cen1,wid1,bkg1,slope1],[height2,cen2,wid2,bkg2,slope2],...] # Different for each scan
    
    FITTING PROCEDURE:
        >
    """

    scans = np.asarray(scans).astype(int)
    
    # Turn depvar into a list
    if type(depvar) is str:
        depvar = [depvar]
    Ndep = len(depvar)
    
    "Masks should be a 2D list of strings"
    "Each string should define an evaluatable statement that defines a region of x"
    "Each defined region will be REMOVED"
    "e.g. [ ['x<-1','x>1','np.abs(x-34)>0.1'], [x<-1.2','x>1.2'] ]"
    if mask_cmd is not None:
        if type(mask_cmd) is str:
            mask_cmd = [[mask_cmd]] # str->2D list
        if type(mask_cmd[0]) is str:
            mask_cmd = [mask_cmd] # 1D list->2D list
        if len(mask_cmd) == 1:
            mask_cmd *= len(scans) # 2D list with single element -> 2D list same length as scans
    
    "estvals should be a 2D list of peak values estimates"
    "There should be the same number of rows as scans"
    "Each row should contain initial estimates for each of the fit parameters"
    "E.G. if gause+slope: [height,cen,wid,bkg,slope]"
    if estvals is not None:
        if len(estvals) < len(scans):
            estvals = [estvals]
        if len(estvals) == 1:
            estvals *= len(scans)
    else:
        estvals = [None]*len(scans)
    
    # Define output dictionary names
    dict_names = ['Scan Number'] + depvar + ['Peak Height','Peak Centre','FWHM','Lorz frac','Background','Area','CHI2 per dof']
    
    # Pre-allocate variables
    valstore = np.zeros([len(scans),8+Ndep])
    errstore = np.zeros([len(scans),8+Ndep])
    x_exp,y_exp = [],[]
    x_fit,y_fit = [],[]
    
    "-----Loading-----"
    for n,run in enumerate(scans):
        d = readscan(run)
        if d is None: print( 'File for run #{} does not exist!'.format(run) ); return
        x,y,dy,labvarx,labvary,ttl = getdata(d,varx=varx,vary=vary,norm=norm,abscor=abscor)[:6]
        
        if mask_cmd is not None: x,y,dy = maskvals(x,y,dy,mask_cmd[n])
            
        valstore[n,0] = run
        errstore[n,0] = run
        
        depstr = ''
        for v in range(Ndep):
            if depvar[v] in d.keys():
                dep_value = np.mean(getattr(d,depvar[v]))
            elif depvar[v] in d.metadata.keys():
                dep_value = getattr(d.metadata,depvar[v])
            else:
                dep_value = n
            valstore[n,v+1] = dep_value
            errstore[n,v+1] = dep_value
            depstr += '{} = {:5.1f} '.format(depvar[v],dep_value)
        
        "-----Run Fit-----"
        out,err = peakfit(x,y,dy,fit_type,bkg_type,peaktest,estvals[n],
                          Nloop,Binit,Tinc,change_factor,converge_max,min_change)
        
        "-----Output-----"
        if 'Peak Height' in out.keys():
            valstore[n,Ndep+1] = out['Peak Height']
            valstore[n,Ndep+2] = out['Peak Centre']
            valstore[n,Ndep+3] = out['FWHM']
            valstore[n,Ndep+4] = out['Lorz frac']
            valstore[n,Ndep+5] = out['Background']
            valstore[n,Ndep+6] = out['Area']
            valstore[n,Ndep+7] = out['CHI2 per dof']
            
            errstore[n,Ndep+1] = err['Peak Height']
            errstore[n,Ndep+2] = err['Peak Centre']
            errstore[n,Ndep+3] = err['FWHM']
            errstore[n,Ndep+4] = err['Lorz frac']
            errstore[n,Ndep+5] = err['Background']
            errstore[n,Ndep+6] = err['Area']
            errstore[n,Ndep+7] = out['CHI2 per dof']
        else:
            print( 'Peak height not found... something went wrong?' )
        
        x_exp += [x]
        y_exp += [y]
        x_fit += [out['x']]
        y_fit += [out['y']]
        peak_rat = ispeak(y,dy,test = peaktest,return_rat=True)
        print( '{0:3.0f}/{1} {2} {3}: Peak={4:7.3g}  Amp={5:<8.0f}  Cen={6:<6.2f}  Wid={7: <5.2g}  Frac={8:<5.2g}  Bkg={9:<8.2g}  Int={10:<8.2g}    CHI{11:8.2g}'.format(n,len(scans)-1,run,depstr,peak_rat,*valstore[n,Ndep+1:]) )
    
    # Data range
    if xrange is None and 'Ta' in depvar:
        Tadep = depvar.index('Ta')+1
        xrange = [min(valstore[:,Tadep])-1,max(valstore[:,Tadep]+1)]
    elif xrange is None and 'psi' in depvar:
        xrange = [0,360]
    elif xrange is None:
        xrange = [min(valstore[:,1]),max(valstore[:,1])]
    
    # PSI range
    if 'psi' in depvar:
        psidep = depvar.index('psi')+1
        valstore[ abs(valstore[:,psidep]) > 1000 ,psidep] = 180.0 # correct psi calculation error
        valstore[ valstore[:,psidep] < xrange[0] ,psidep] = valstore[ valstore[:,psidep] < xrange[0] ,psidep] + 360
        valstore[ valstore[:,psidep] > xrange[0]+360 ,psidep] = valstore[ valstore[:,psidep] > xrange[0]+360 ,psidep] - 360
    
    # Title
    ttl = scantitle(scans)
    fttl = '{} {}\n{} vs {}\n{}'.format(ttl, fit_type, labvarx, labvary, d.metadata.cmd_short)
    
    "---Plot individual fits---"
    if show_fits == True:
        # Plot each scan in an axis of a 16*16 figure
        Nplots = 25
        for n in range(int(np.ceil(len(scans)/float(Nplots)))):
            fig, axs = plt.subplots(5,5,figsize=[24,14],dpi=80)
            fig.subplots_adjust(hspace = 0.35,wspace=0.32,left=0.07,right=0.97)
            plt.suptitle(fttl,fontsize=14)
            #fig.subplots_adjust(wspace = 0.25)
            axs = axs.flatten()
            pltruns = scans[n*Nplots:(n+1)*Nplots]
            for axn,rn in enumerate(pltruns):
                calno = axn+n*Nplots
                axs[axn].plot(x_exp[calno],y_exp[calno],'b-+',linewidth=2)
                axs[axn].plot(x_fit[calno],y_fit[calno],'r-',linewidth=2)
                axs[axn].set_title('#{}: {}={}'.format(rn,depvar[0],valstore[calno,1]))
            "---Save Multiplot figure---"
            if savePLOT not in [None, False, '']:
                if type(savePLOT) is str:
                    saveplot('{} FITS {}'.format(savePLOT,n))
                else:
                    saveplot('{0} ScansFIT {1:1.0f}-{2:1.0f} {3} {4} FITS {5}'.format(' '.join(depvar),scans[0],scans[-1],vary,fit_type,n))
                
    # Sort by depvar
    unsorted_valstore = 1.0*valstore
    unsorted_errstore = 1.0*errstore
    if sortdep:
        idx = np.argsort(valstore[:,1])
        valstore = valstore[idx,:]
        errstore = errstore[idx,:]
    
    # Save valstore & errstore values in text files
    if saveFIT not in [None, False, '']:
        header = '{}\n{} vs {}\n{}'.format(fttl, labvarx, labvary, ','.join(dict_names))
        if type(saveFIT) is str:
            savefile = os.path.join(savedir, '{}.dat'.format(saveFIT))
            esavefile = os.path.join(savedir, '{}_errors.dat'.format(saveFIT))
            np.savetxt(savefile,unsorted_valstore,header=header)
            np.savetxt(esavefile,unsorted_errstore,header=header)
            print( 'Saved as {}'.format(savefile) )
        else:
            savefile = os.path.join(savedir, '{0} ScansFIT {1:1.0f}-{2:1.0f} {3} {4}.dat'.format(' '.join(depvar),scans[0],scans[-1],vary,fit_type))
            esavefile = os.path.join(savedir, '{0} ScansFIT {1:1.0f}-{2:1.0f} {3} {4}_errors.dat'.format(' '.join(depvar),scans[0],scans[-1],vary,fit_type))
            np.savetxt(savefile,unsorted_valstore,header=header)
            np.savetxt(esavefile,unsorted_errstore,header=header)
            print( 'Saved as {}'.format(savefile) )
            print( 'Reload this scan with:\n val,err = load_fits([{},{}], depvar={}, vary=\'{}\', fit_type=\'{}\')'.format(scans[0], scans[-1], depvar, vary, fit_type) )
        
    
    "------Plotting------"
    try:
        plot.__iter__; # test if array
    except AttributeError:
        plot = [plot]
    for nplot in plot:
        if nplot == 'all':
            # 2D Plots of area, width and position
            fig = plt.figure(figsize=[18,12], dpi=fig_dpi)
            ax1 = fig.add_subplot(231) # Area
            #plt.plot(valstore[:,0],valstore[:,6],'-o',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+6],errstore[:,Ndep+6],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Integrated Area',fontsize=18)
            plt.xlim(xrange)
            plt.ylim([0,max(valstore[:,Ndep+6])*1.1])
             
            ax2 = fig.add_subplot(232) # Width
            #plt.plot(valstore[:,0],valstore[:,3],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+3],errstore[:,Ndep+3],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Fitted Width',fontsize=18)
            plt.xlim(xrange)
            #plt.ylim([0,plt.ylim()[1]])
             
            ax3=fig.add_subplot(233) # Centre
            #plt.plot(valstore[:,0],valstore[:,2],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+2],errstore[:,Ndep+2],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Fitted Centre',fontsize=18)
            plt.xlim(xrange)
            
            fig.add_subplot(234) # Area
            #plt.plot(valstore[:,0],valstore[:,5],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+5],errstore[:,Ndep+5],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Fitted Background',fontsize=18)
            plt.xlim(xrange)
             
            ax2=fig.add_subplot(235) # Lorentz fraction
            #plt.plot(valstore[:,0],valstore[:,4],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+4],errstore[:,Ndep+4],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Lorentz fraction',fontsize=18)
            plt.xlim(xrange)
            #plt.ylim([0,plt.ylim()[1]])
             
            ax3=fig.add_subplot(236) # CHI^2
            plt.plot(valstore[:,1],valstore[:,Ndep+7],'-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('CHI^2 per dof',fontsize=18)
            plt.xlim(xrange)
            
            plt.suptitle(fttl,fontsize=14)
            fig.subplots_adjust(wspace = 0.25,hspace=0.3)
            #plt.tight_layout()
            # Stop offset in x tick labels
            #x_formatter = ticker.ScalarFormatter(useOffset=False)
            #ax3.xaxis.set_major_formatter(x_formatter)
        elif nplot == 'int':
            # 2D Plots of area
            fig = plt.figure(figsize=[8,8], dpi=fig_dpi)
            #plt.plot(valstore[:,0],valstore[:,6],'-o',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+6],errstore[:,Ndep+6],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Integrated Sum('+labvary+')',fontsize=18)
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            plt.ylim([0,max(valstore[:,Ndep+6])*1.1])
            fig.subplots_adjust(left=0.15)
        elif nplot == 'cen':
            # 2D Plots of centre
            fig = plt.figure(figsize=[8,8], dpi=fig_dpi)
            #plt.plot(valstore[:,0],valstore[:,2],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+2],errstore[:,Ndep+2],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel(labvarx+' Centre',fontsize=18)
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        elif nplot == 'wid':
            # 2D Plots of width
            fig = plt.figure(figsize=[8,8], dpi=fig_dpi)
            #plt.plot(valstore[:,0],valstore[:,3],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+3],errstore[:,Ndep+3],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel(labvarx+' Width',fontsize=18)
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        elif nplot == 'surface':
            # 2D surface plot of multi-dimension scans
            sx = unsorted_valstore[:,Ndep-1]
            sy = unsorted_valstore[:,Ndep]
            inte = unsorted_valstore[:,Ndep+6]
            
            # Determine the repeat length of the scans
            delta_x = np.abs(np.diff(sx))
            ch_idx_x = np.where(delta_x > delta_x.max()*0.9) # find biggest changes
            ch_delta_x = np.diff(ch_idx_x)
            rep_len_x = np.round(np.mean(ch_delta_x))
            delta_y = np.abs(np.diff(sy))
            ch_idx_y = np.where(delta_y > delta_y.max()*0.9) # find biggest changes
            ch_delta_y = np.diff(ch_idx_y)
            rep_len_y = np.round(np.mean(ch_delta_y))
            print('Scans in {} are repeating every {} iterations'.format(depvar[0], rep_len_x))
            print('Scans in {} are repeating every {} iterations'.format(depvar[1], rep_len_y))
            rep_len = int(max(rep_len_x, rep_len_y))
            
            # Reshape into square arrays
            # If this is problematic, look at scipy.interpolate.griddata
            sx_squareA = sx[:rep_len*(len(sx)//rep_len)].reshape(-1,rep_len)
            sy_squareA = sy[:rep_len*(len(sy)//rep_len)].reshape(-1,rep_len)
            int_squareA = inte[:rep_len*(len(inte)//rep_len)].reshape(-1,rep_len)
            
            plt.figure(figsize=[12,10], dpi=fig_dpi)
            plt.pcolormesh(sx_squareA, sy_squareA, int_squareA)
            plt.axis('image')
            cb = plt.colorbar()
            plt.xlabel(dict_names[Ndep-1],fontsize=22)
            plt.ylabel(dict_names[Ndep],fontsize=22)
            plt.title(fttl,fontsize=18)
            cb.set_label('Integrated Area',fontsize=22)
        
        "---Save Figure---"
        if savePLOT not in [None, False, '']:
            if type(savePLOT) is str:
                saveplot('{} {}'.format(savePLOT,nplot))
            else:
                saveplot('{0} ScansFIT {1} {2:1.0f}-{3:1.0f} {4} {5}'.format(' '.join(depvar),nplot,scans[0],scans[-1],vary,fit_type))
    plt.show()
    
    " Prepare output dicts"
    out_values = {}
    out_errors = {}
    for n,name in enumerate(dict_names):
        out_values[name] = valstore[:,n]
        out_errors[name] = errstore[:,n]
    
    return out_values,out_errors

def load_fits(scans=[0], depvar='Ta', plot=None, vary='sum', fit_type = 'pVoight', file=None, disp=False, save=False,**fitopt):
    """ 
    New
     Load previously fitted data from fit_scans(), assuming it completed succesfully and was saved.
      > The inputs are used to determine the automatic filename
      > Both the values and error files are read
      > Plots can be re-generated
      > A typical filename is: 'Ta ScansFIT 00001-00010 pVoight.dat'
     
         val,err = load_fits(scans,depvar='Ta',fit_type = 'pVoight',plot=None)
     
     INPUTS:
             scans : [N,M]     : list of scan numbers used in fit_scans(), only the first and final values are required
           depvar : 'Ta'      : Independent variable, e.g. 'Ta','Energy','psi'
             vary : 'roi2_sum': Scanned variable e.g. 'APD', 'sum', 'nroi[31,31]'
         fit_type : 'pVoight  : Type of fit to perform e.g. 'Simple','pVoight','Lorentz','Gauss'
             plot : None      : Plot the loaded fit data, available: [None,'all','int','cen','wid']
             file : None      : If given, this filename will be used to open the data files
             disp : False     : If True, prints the data on the console
             save : False     : If true, saves images of the plots
       ALTERNATIVE:
             scans : 'filename': Name of .dat file where fit data is stored (no other parameters are required) 
     OUTPUTS:
              val : Fitted values
                    val['Scan Number']
                    val[depvar]
                    val['Peak Height']
                    val['Peak Centre']
                    val['FWHM']
                    val['Lorz frac']
                    val['Background']
                    val['Area']
                    val['CHI2 per dof']
                    
              err : Errors on fitted values
                    err['Scan Number']
                    err['Peak Height']
                    err['Peak Centre']
                    err['FWHM']
                    err['Lorz frac']
                    err['Background']
                    err['Area']
                    err['CHI2 per dof']
    
    E.G. val,err = load_fits([572138,572368],depvar=['sx', 'sy'],fit_type='Gauss')
    E.G. val,err = load_fits('sx sy ScansFIT 572138-572368 Gauss')
    
    """
    
    # Turn depvar into a list
    if type(depvar) is str:
        depvar = [depvar]
    Ndep = len(depvar)
    
    if type(scans) is str:
        # scans == filename
        if scans[-4:] == '.dat':
            scans = scans[:-4]
        file = os.path.join(savedir,scans+'.dat')
        efile = os.path.join(savedir,scans+'_errors.dat')
    
    if file is None:
        scans = np.asarray(scans).astype(int)
        fname = '{0} ScansFIT {1:1.0f}-{2:1.0f} {3} {4}.dat'.format(' '.join(depvar), scans[0], scans[-1], vary, fit_type)
        ename = '{0} ScansFIT {1:1.0f}-{2:1.0f} {3} {4}_errors.dat'.format(' '.join(depvar), scans[0], scans[-1], vary, fit_type)
        file = os.path.join(savedir, fname)
        efile = os.path.join(savedir, ename)
        print(file)
        print(efile)
    
    valstore = np.loadtxt(file)
    errstore = np.loadtxt(efile)
    with open(file,'r') as ff:
        first_line = ff.readline().strip()
        # if title given in header, get next line
        while 'Scan' not in first_line:
            first_line = ff.readline().strip()
    
    first_line = first_line.strip('# ')
    names = first_line.split(',')
    "-----Printing-----"
    if disp:
        fmt = '{:3.0f}/{:<3.0f} '
        for name in names:
            fmt += name+'={:<8.2g} '
        for n in range(len(valstore)):
            #print( '{0:3.0f}/{1} {2} {3}: Amp={5:<8.0f}  Cen={6:<6.2f}  Wid={7: <5.2g}  Frac={8:<5.2g}  Bkg={9:<8.2g}  Int={10:<8.2g}    CHI{11:8.2g}'.format(n,len(scans)-1,*valstore[n,Ndep+1:]) )
            print( fmt.format(n,len(scans)-1,*valstore[n,:]) )
    
    "------Plotting------"
    try:
        plot.__iter__; # test if array
    except AttributeError:
        plot = [plot]
    fttl = '#{} {} {}'.format(numbers2string(scans),vary,fit_type)
    xrange = [valstore[:,1].min(),valstore[:,1].max()]
    depvar = names[1:]
    for nplot in plot:
        if nplot is None:
            break
        if nplot == 'all':
            # 2D Plots of area, width and position
            fig = plt.figure(figsize=[18,12], dpi=fig_dpi)
            ax1 = fig.add_subplot(231) # Area
            #plt.plot(valstore[:,0],valstore[:,6],'-o',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+6],errstore[:,Ndep+6],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Integrated Area',fontsize=18)
            plt.xlim(xrange)
            plt.ylim([0,max(valstore[:,Ndep+6])*1.1])
             
            ax2 = fig.add_subplot(232) # Width
            #plt.plot(valstore[:,0],valstore[:,3],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+3],errstore[:,Ndep+3],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Fitted Width',fontsize=18)
            plt.xlim(xrange)
            #plt.ylim([0,plt.ylim()[1]])
             
            ax3=fig.add_subplot(233) # Centre
            #plt.plot(valstore[:,0],valstore[:,2],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+2],errstore[:,Ndep+2],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Fitted Centre',fontsize=18)
            plt.xlim(xrange)
            
            fig.add_subplot(234) # Area
            #plt.plot(valstore[:,0],valstore[:,5],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+5],errstore[:,Ndep+5],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Fitted Background',fontsize=18)
            plt.xlim(xrange)
             
            ax2=fig.add_subplot(235) # Lorentz fraction
            #plt.plot(valstore[:,0],valstore[:,4],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+4],errstore[:,Ndep+4],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Lorentz fraction',fontsize=18)
            plt.xlim(xrange)
            #plt.ylim([0,plt.ylim()[1]])
             
            ax3=fig.add_subplot(236) # CHI^2
            plt.plot(valstore[:,1],valstore[:,Ndep+7],'-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('CHI^2 per dof',fontsize=18)
            plt.xlim(xrange)
            
            plt.suptitle(fttl,fontsize=14)
            fig.subplots_adjust(wspace = 0.25,hspace=0.3)
            #plt.tight_layout()
            # Stop offset in x tick labels
            #x_formatter = ticker.ScalarFormatter(useOffset=False)
            #ax3.xaxis.set_major_formatter(x_formatter)
        elif nplot == 'int':
            # 2D Plots of area
            fig = plt.figure(figsize=[8,8], dpi=fig_dpi)
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,6],'-o',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+6],errstore[:,Ndep+6],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Integrated Sum',fontsize=18)
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            plt.ylim([0,max(valstore[:,Ndep+6])*1.1])
            fig.subplots_adjust(left=0.15)
        elif nplot == 'cen':
            # 2D Plots of centre
            fig = plt.figure(figsize=[8,8], dpi=fig_dpi)
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,2],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+2],errstore[:,Ndep+2],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Centre',fontsize=18)
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        elif nplot == 'wid':
            # 2D Plots of width
            fig = plt.figure(figsize=[8,8], dpi=fig_dpi)
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,3],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+3],errstore[:,Ndep+3],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Width',fontsize=18)
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        elif nplot == 'surface':
            # 2D surface plot of multi-dimension scans
            sx = valstore[:,Ndep-1]
            sy = valstore[:,Ndep]
            inte = valstore[:,Ndep+6]
            
             # Determine the repeat length of the scans
            delta_x = np.abs(np.diff(sx))
            ch_idx_x = np.where(delta_x > delta_x.max()*0.9) # find biggest changes
            ch_delta_x = np.diff(ch_idx_x)
            rep_len_x = np.round(np.mean(ch_delta_x))
            delta_y = np.abs(np.diff(sy))
            ch_idx_y = np.where(delta_y > delta_y.max()*0.9) # find biggest changes
            ch_delta_y = np.diff(ch_idx_y)
            rep_len_y = np.round(np.mean(ch_delta_y))
            print('Scans in {} are repeating every {} iterations'.format(depvar[0], rep_len_x))
            print('Scans in {} are repeating every {} iterations'.format(depvar[1], rep_len_y))
            rep_len = int(max(rep_len_x, rep_len_y))
            
            # Reshape into square arrays
            # If this is problematic, look at scipy.interpolate.griddata
            sx_squareA = sx[:rep_len*(len(sx)//rep_len)].reshape(-1,rep_len)
            sy_squareA = sy[:rep_len*(len(sy)//rep_len)].reshape(-1,rep_len)
            int_squareA = inte[:rep_len*(len(inte)//rep_len)].reshape(-1,rep_len)
            
            plt.figure(figsize=[12,10], dpi=fig_dpi)
            plt.pcolormesh(sx_squareA, sy_squareA, int_squareA)
            plt.axis('image')
            cb = plt.colorbar()
            plt.xlabel(dict_names[Ndep-1],fontsize=22)
            plt.ylabel(dict_names[Ndep],fontsize=22)
            plt.title(fttl,fontsize=18)
            cb.set_label('Integrated Area',fontsize=22)
        
        "---Save Figure---"
        if save not in [None, False, '']:
            if type(save) is str:
                saveplot('{} {}'.format(save,nplot))
            else:
                saveplot('{0} ScansFIT {1} {2:1.0f}-{3:1.0f} {4} {5}'.format(' '.join(depvar),nplot,scans[0],scans[-1],vary,fit_type))
    plt.show()
    
    " Prepare output dicts"
    out_values = {}
    out_errors = {}
    for n,name in enumerate(names):
        out_values[name] = valstore[:,n]
        out_errors[name] = errstore[:,n]
    
    return out_values,out_errors

def pil_peaks(scans,depvar='Ta',ROIsize=[31,31],cax=None,save=False):
    "Show pilatus peak search results"
    
    # Load data
    x,y,z,varx,vary,varz,ttl = joindata(scans,vary=depvar)
    
    # Plot each scan in an axis of a 16*16 figure
    Nplots = 25
    ROIcenvals=np.zeros([len(scans),2])
    for n in range(int(np.ceil(len(scans)/float(Nplots)))):
        fig, axs = plt.subplots(5,5,figsize=[24,14],dpi=80)
        fig.subplots_adjust(hspace = 0.35,wspace=0.32,left=0.07,right=0.97)
        plt.suptitle(ttl,fontsize=14)
        #fig.subplots_adjust(wspace = 0.25)
        axs = axs.flatten()
        pltruns = scans[n*Nplots:(n+1)*Nplots]
        for axn,rn in enumerate(pltruns):
            calno = axn+n*Nplots
            
            vol=getvol(rn)
            ROIcen,frame = pilpeak(vol,disp=True)
            ROIcenvals[calno] = ROIcen
            
            im = axs[axn].imshow(vol[:,:,frame])
            
            " Default colour thresholds"
            if cax is None:
                md = np.median(vol[:,:,frame])
                mx = np.max(vol[:,:,frame])
                cmax = md + 10**(0.7*np.log10(mx-md))
                if cmax <= 0: cmax = 1
                cax = [0,cmax]
            im.set_clim(cax)
            
            pil_size = vol.shape[:2]
            idxi = np.array([ROIcen[0]-ROIsize[0]//2,ROIcen[0]+ROIsize[0]//2+1])
            idxj = np.array([ROIcen[1]-ROIsize[1]//2,ROIcen[1]+ROIsize[1]//2+1])
            print( 'ROI = [{0},{1},{2},{3}]'.format(idxi[0],idxj[0],idxi[1],idxj[1]) )
            axs[axn].plot(idxj[[0,1,1,0,0]],idxi[[0,0,1,1,0]],'k-',linewidth=2)
            axs[axn].plot([pil_centre[1],pil_centre[1]],[0,pil_size[0]],'k:',linewidth=2)
            axs[axn].plot([0,pil_size[1]],[pil_centre[0],pil_centre[0]],'k:',linewidth=2)
            axs[axn].set_aspect('equal')
            axs[axn].autoscale(tight=True)
            axs[axn].set_title('#{}: {}={}\nROIcen=[{},{}]'.format(rn,depvar,y[0,calno],ROIcen[0],ROIcen[1]))
            
        
        "---Save Multiplot figure---"
        if save not in [None, False, '']:
            if type(save) is str:
                saveplot('{} {}'.format(save,n))
            else:
                saveplot('{} PilatusPeaks {:1.0f}-{:1.0f} {}'.format(depvar,scans[0],scans[-1],n))
        
    plt.figure(figsize=[10,6], dpi=fig_dpi)
    plt.plot(ROIcenvals[:,1],ROIcenvals[:,0],'b-o',linewidth=2)
    plt.plot(ROIcenvals[0,1],ROIcenvals[0,0],'g+',markeredgewidth=3,markersize=22,label='Start')
    plt.plot(ROIcenvals[-1,1],ROIcenvals[-1,0],'r+',markeredgewidth=3,markersize=22,label='Finish')
    plt.plot([pil_centre[1],pil_centre[1]],[0,pil_size[0]],'k:',linewidth=2)
    plt.plot([0,pil_size[1]],[pil_centre[0],pil_centre[0]],'k:',linewidth=2)
    plt.axis([1,487,1,195])
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal')
    plt.legend(loc=0)
    labels(ttl,'Pilatus i','Pilatus j')
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot('{} Pilatus i-j'.format(save))
        else:
            saveplot('{} PilatusPeaks {:1.0f}-{:1.0f} i-j'.format(depvar,scans[0],scans[-1]))


"---------------------------Plotting Functions----------------------------"


def plotscan(num=None,vary='',varx='',fit=None,norm=True,sum=False,subtract=False,fits=False,
             labels=None,logplot=False,diffplot=False,normalise=False,save=False):
    """
    Default plot of I16 data, plotscan(#), or plotscan(#,save=True)
    
    Plotting with fits: plotscan(#,fit='Gauss') or plotscan(#,fit='Lorentz') or plotscan(#,fit='pVoight')
    
    Plotting multiple regions of interest: plotscan(#,vary=['roi2_sum','roi1_sum'])
    
    Subtract regions of interest: plotscan(#,vary=['droi[110,240,51,51','droi[110,340,51,51]],subtract=True) = droi1-droi2
    
    Multi-run plot: plotscan([#1,#2,#3,...])
    
    Multi-run plot with fits: plotscan([#1,#2,#3,...],fit='Gauss',fits=True)
    
    Sum of several scans: plotscan([#1,#2],sum=True)
    
    """
    
    "---Handle inputs---"
    if num is None:
        num = latest()
    
    " Multiple nums given"
    try:
        nums = num[1:]
        num = num[0]
    except:
        nums=[]
    
    " Multiple vary's given"
    if type(vary) is str or vary is None:
        varys = []
    else:
        varys = vary[1:]
        vary = vary[0]
    
    "---Load data---"
    try:
        d = readscan(num)
        if d is None: return
    except:
        d = num
        num = d.metadata.SRSRUN
    
    " Get metadata"
    m = d.metadata
    cmd = m.cmd_short # Scan command
    sampsl = '{0:4.2g}x{1:<4.2g}'.format(m.s5xgap,m.s5ygap)
    detsl = '{0:4.2g}x{1:<4.2g}'.format(m.s6xgap,m.s6ygap)
    atten1 = '{0:1.0f}'.format(m.Atten)
    if labels is None:
        lbl = str(m.SRSRUN)
        if len(varys)>0: lbl=vary
    else:
        # make list
        labels = np.asarray(labels).reshape(-1)
        
        lbl = '{:6.0f}'.format(m.SRSRUN)
        for label in labels:
            if label in m.keys():
                lbl += ', {} = {}'.format(label,getattr(m,label))
            elif label in d.keys():
                lbl += ', {} = {}'.format(label,np.mean(getattr(d,label)))
            elif label == 'hkl':
                lbl += ', {}'.format(scanhkl(d))
            elif label == 'T':
                lbl += ', {}'.format(scantemp(d))
    
    " Get data" 
    x,y,dy,varx,varynew,ttl = getdata(d,vary=vary,varx=varx,norm=norm)[:6]
    
    "---Sum scans---"
    if len(nums) > 0 and sum==True:
        lgd = '{}'.format(m.SRSRUN)
        dy = dy**2 # add errors in quadrature
        for rn in nums:
            # Get data 
            x2,y2,dy2,varx2,vary2,ttl2,d2 = getdata(rn,vary=vary,varx=varx,norm=norm)
            lgd += '+\n{}'.format(rn)
            y += y2
            dy += dy2**2
        dy = np.sqrt(dy)
    
    "---Subtract vary's---"
    if len(varys) > 0 and subtract==True:
        lgd = vary
        dy = dy**2 # add errors in quadrature
        for nvary in varys:
            # Get data 
            x2,y2,dy2,varx2,vary2,ttl2,d2 = getdata(d,vary=nvary,varx=varx,norm=norm)
            len(y2)
            lgd += '-\n{}'.format(nvary)
            varynew += '-{}'.format(nvary)
            y -= y2
            dy += dy2**2
        dy = np.sqrt(dy)
    
    "---Differentiate---"
    if diffplot:
        y = np.abs(np.gradient(y))
        varynew = 'd {} / d {}'.format(varynew,varx)
    
    "---Normalise---"
    if normalise:
        ynorm = np.sqrt(np.sum( y**2 ))/len(y)
        y = y/ynorm
        dy = dy/ynorm
        varynew = 'Normalised Intensity'
    
    "---Create plot---"
    fig = plt.figure(figsize=[10,8], dpi=fig_dpi)
    
    plt.errorbar(x,y,dy,fmt='-o',c=plot_colors[0],linewidth=2,label=lbl)
    
    plt.xlabel(varx, fontsize=18)
    plt.ylabel(varynew, fontsize=18)
    plttl = ttl+'\n'+cmd+'\nss ={}, ds ={}, atten = {}'.format(sampsl,detsl,atten1)
    plt.title(plttl, fontsize=14)
    
    if logplot: plt.gca().set_yscale(u'log')
    if len(nums) > 0 and sum==True: plt.legend([lgd],loc='best')
    #if len(varys) > 0 and subtract==True: plt.legend([lgd],loc='best')
    fig.subplots_adjust(left=0.15)
    fig.subplots_adjust(top=0.85)
    
    # Change formats in x & y axes so they are nicer
    plt.gca().get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%8.3g'))
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    
    "---Fit Data---"
    if fit != None:
        out,err = peakfit(x,y,type=fit,Nloop=100,peaktest=-100,interpolate=True,disp=True)
        
        xfit,yfit=out['x'],out['y']
        amp,damp = out['Peak Height'],err['Peak Height']
        cen,dcen = out['Peak Centre'],err['Peak Centre']
        wid,dwid = out['FWHM'],err['FWHM']
        bkg,dbkg = out['Background'],err['Background']
        ara,dara = out['Area'],err['Area']
        
        " add to plot"
        ax = plt.gca()
        fitcol = plot_colors[-1]
        txtcol = 'k'
        if fits: 
            fitcol = plot_colors[0]
            txtcol = plot_colors[0]
        ax.plot(xfit,yfit,':',linewidth=3,c=fitcol,label='Fit')
        
        " Print results on plot"
        t0 = '   ---- {} Fit ----    '.format(fit)
        t1 = ' Amp = {}'.format(stfm(amp,damp))
        t2 = ' Cen = {}'.format(stfm(cen,dcen))
        t3 = ' Wid = {}'.format(stfm(wid,dwid))
        t4 = ' Bkg = {}'.format(stfm(bkg,dbkg))
        t5 = ' Int = {}'.format(stfm(ara,dara))
        
        axs = ax.axis()
        axw = axs[1]-axs[0]
        axh = axs[3]-axs[2]
        if cen > axs[0] + axw/2.:
            tx = axs[0] + axw*0.05
        else:
            tx = axs[0] + axw*0.65
        ty = min(yfit) + 0.9*(max(yfit)-min(yfit))
        txt = t0+'\n'+t1+'\n'+t2+'\n'+t3+'\n'+t4+'\n'+t5
        plt.text(tx,ty,txt,color=txtcol,
         horizontalalignment = 'left',
         verticalalignment   = 'top',
         multialignment      = 'left')
    
    "---Plot multiple scans---"
    if len(nums) > 0 and sum==False:
        for n,rn in enumerate(nums):
            " Get data" 
            x2,y2,dy2,varx2,vary2,ttl2,dn = getdata(rn,vary=vary,varx=varx,norm=norm)
            mn = dn.metadata
            
            " Generate label"
            if labels is None:
                lbl = str(mn.SRSRUN)
            else:
                lbl = '{:6.0f}'.format(mn.SRSRUN)
                for label in labels:
                    if label in m.keys():
                        lbl += ', {} = {}'.format(label,getattr(mn,label))
                    elif label in d.keys():
                        lbl += ', {} = {}'.format(label,np.mean(getattr(dn,label)))
                    elif label == 'hkl':
                        lbl += ', {}'.format(scanhkl(dn))
                    elif label == 'T':
                        lbl += ', {}'.format(scantemp(dn))
            
            " Subtract ROIs"
            if len(varys) > 0 and subtract==True:
                dy2 = dy2**2 # add errors in quadrature
                for nvary in varys:
                    # Get data 
                    x3,y3,dy3 = getdata(dn,vary=nvary,varx=varx,norm=norm)[:3]
                    y2 -= y3
                    dy2 += dy3**2
                dy2 = np.sqrt(dy2)
            
            "---Differentiate---"
            if diffplot:
                y2 = np.abs(np.gradient(y2))
            
            " Normalise"
            if normalise:
                y2norm = np.sqrt(np.sum( y2**2 ))/len(y2)
                y2 = y2/y2norm
                dy2 = dy2/y2norm
            
            " Plot data"
            plotcol = plot_colors[(n+1)-(n+1)//len(plot_colors)*len(plot_colors)]
            plt.errorbar(x2,y2,dy2,fmt='-o',c=plotcol,linewidth=2,label=lbl)
            
            " Fit data"
            if fits == True:
                out,err = peakfit(x2,y2,type=fit,peaktest=-100,interpolate=True,disp=True)
        
                xfit,yfit=out['x'],out['y']
                amp,damp = out['Peak Height'],err['Peak Height']
                cen,dcen = out['Peak Centre'],err['Peak Centre']
                wid,dwid = out['FWHM'],err['FWHM']
                bkg,dbkg = out['Background'],err['Background']
                ara,dara = out['Area'],err['Area']
                
                ax = plt.gca()
                ax.plot(xfit,yfit,':',linewidth=2,c=plotcol,label='Fit')
                
                # Print results to plot
                t0 = '   ---- {} Fit ----    '.format(fit)
                t1 = ' Amp = {}'.format(stfm(amp,damp))
                t2 = ' Cen = {}'.format(stfm(cen,dcen))
                t3 = ' Wid = {}'.format(stfm(wid,dwid))
                t4 = ' Bkg = {}'.format(stfm(bkg,dbkg))
                t5 = ' Int = {}'.format(stfm(ara,dara))
                
                axs = ax.axis()
                axw = axs[1]-axs[0]
                axh = axs[3]-axs[2]
                if cen > axs[0] + axw/2.:
                    tx = axs[0] + axw*0.05
                else:
                    tx = axs[0] + axw*0.65
                #ty = axs[2] + axh*0.4
                ty = min(yfit) + 0.9*(max(yfit)-min(yfit))
                txt = t0+'\n'+t1+'\n'+t2+'\n'+t3+'\n'+t4+'\n'+t5
                plt.text(tx,ty,txt,color=plotcol,
                 horizontalalignment = 'left',
                 verticalalignment   = 'top',
                 multialignment      = 'left')
            
        plt.legend(loc='best', frameon=False)
        #plt.legend(loc='upper left')
        #plt.legend(loc='upper right')
        
    "---Plot multiple vary values---"
    if len(varys) > 0 and subtract==False:
        for n,nvary in enumerate(varys):
            " Get data" 
            x2,y2,dy2,varx2,vary2,ttl2,dn = getdata(num,vary=nvary,varx=varx,norm=norm)
            mn = dn.metadata
            
            " Generate label"
            if labels is None:
                lbl = '{}'.format(nvary)
            else:
                lbl = '{}, {} = {}'.format(str(mn.SRSRUN),labels,getattr(mn,labels))
            
            " Plot data"
            plotcol = plot_colors[(n+1)-((n+1)//len(plot_colors))*len(plot_colors)]
            plt.errorbar(x2,y2,dy2,fmt='-o',c=plotcol,linewidth=2,label=lbl)
            
            " Fit data"
            if fits == True:
                out,err = peakfit(x2,y2,type=fit,interpolate=True,disp=True)
        
                xfit,yfit=out['x'],out['y']
                amp,damp = out['Peak Height'],err['Peak Height']
                cen,dcen = out['Peak Centre'],err['Peak Centre']
                wid,dwid = out['FWHM'],err['FWHM']
                bkg,dbkg = out['Background'],err['Background']
                ara,dara = out['Area'],err['Area']
                
                ax = plt.gca()
                ax.plot(xfit,yfit,':',c=plotcol,label='Fit')
                
                # Print results to plot
                t0 = '   ---- {} Fit ----    '.format(fit)
                t1 = ' Amp = {}'.format(stfm(amp,damp))
                t2 = ' Cen = {}'.format(stfm(cen,dcen))
                t3 = ' Wid = {}'.format(stfm(wid,dwid))
                t4 = ' Bkg = {}'.format(stfm(bkg,dbkg))
                t5 = ' Int = {}'.format(stfm(ara,dara))
                
                axs = ax.axis()
                axw = axs[1]-axs[0]
                axh = axs[3]-axs[2]
                if cen > axs[0] + axw/2.:
                    tx = axs[0] + axw*0.05
                else:
                    tx = axs[0] + axw*0.65
                #ty = axs[2] + axh*0.4
                ty = min(yfit) + 0.9*(max(yfit)-min(yfit))
                txt = t0+'\n'+t1+'\n'+t2+'\n'+t3+'\n'+t4+'\n'+t5
                plt.text(tx,ty,txt,color=plotcol,
                 horizontalalignment = 'left',
                 verticalalignment   = 'top',
                 multialignment      = 'left')
        
        plt.legend(loc='best', frameon=False)
        #plt.legend(loc='upper left')
        #plt.legend(loc='upper right')
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot(ttl)
    # plt.show()

def plotpil(num,cax=None,varx='',imnum=None,bkg_img=None,ROIcen=None,ROIsize=[75,67],show_ROIbkg=False,show_peakregion=False,log_colors=False,save=False, colormap=None):
    """
    Pilatus image viewer, plotpil(#)
    Displays pilatus image with a slider to looks throught the whole scan
    Options:
      num             : Scan number (0 for latest scan)
      cax             : [min,max] colour map cut offs, or None for automatic (Default=None)
      varx            : scannable to use as slider axis, or None for automatic (Default=None)
      imnum           : Image number to start at - 0-scanlength, or None for middle (Default=None)
      bkg_img         : If not None, image number bkg_img will be subtracted from all the other images in the scan (Default=None)
      ROIcen          : Centre of region of interest [-y,x], or None for pp.pil_centre (Default=None)
      ROIsize         : Size of region of interest to plot [ywid,xwid] (Default=[75,67]-ROI2)
      show_ROIbkg     : If True, plots the background region of the ROI (Default=False)
      show_peakregion : If True, plots the pp.peakregion area (Default=False)
      save            : If True, saves an image of the resulting figure to pp.savedir (Default=False)
    
    E.G. plotpil(0,[0,1E3],imnum=0,bkg_img=0,ROIcen=[110,242],ROIsize=[21,10],show_peakregion=True)
    """
    
    " Load data file"
    d = readscan(num)
    if d is None: print( 'File for run #{} does not exist!'.format(num) ); return
    
    " Get data"
    x,y,dy,varx,varynew,ttl = getdata(d,varx=varx)[:6]
    
    " Load pilatus images as volume"
    vol = getvol(num)
    
    " Get scan information"
    Nframe = len(d.path)
    cmd = d.metadata.cmd_short
    
    " Subtract one frame from all the others for background subtraction"
    if bkg_img is not None:
        bkgcut = vol[:,:,bkg_img].copy()
        for n in range(len(x)):
            vol[:,:,n] = vol[:,:,n] - bkgcut
    
    " Default initial frame"
    if imnum is None:
        imnum = int(Nframe//2)
    
    if colormap is None:
        colormap = default_colormap
    " Colormap name"
    if type(colormap) == str:
        colormap = plt.get_cmap(colormap)
    
    " Default colour thresholds"
    if cax is None:
        if log_colors:
            cax = [1,vol.max()]
        else:
            md = np.median(vol[:,:,imnum])
            mx = np.max(vol[:,:,imnum])
            cmax = md + 10**(0.7*np.log10(mx-md))
            if cmax <= 0: cmax = 1
            cax = [0,cmax]
        print( 'caxis set at [{0:1.3g},{1:1.3g}]'.format(cax[0],cax[1]) )
    
    " Create figure & plot 1st image"
    fig = plt.figure(figsize=[10,6], dpi=fig_dpi)
    ax = fig.add_subplot(111)
    if log_colors:
        p = plt.imshow(vol[:,:,imnum]+1, cmap=colormap, norm=LogNorm())
    else:
        p = plt.imshow(vol[:,:,imnum])
    p.set_clim(cax)
    ttl = ttl+'\n'+cmd
    ax.set_title(ttl, fontsize=20)
    
    " Add roi lines to plot"
    if ROIcen is None:
        ROIcen = pil_centre
    
    pil_size = vol.shape[:2]
    idxi = np.array([ROIcen[0]-ROIsize[0]//2,ROIcen[0]+ROIsize[0]//2+1])
    idxj = np.array([ROIcen[1]-ROIsize[1]//2,ROIcen[1]+ROIsize[1]//2+1])
    print( 'ROI = [{0},{1},{2},{3}]'.format(idxi[0],idxj[0],idxi[1],idxj[1]) )
    ax.plot(idxj[[0,1,1,0,0]],idxi[[0,0,1,1,0]],'k-',linewidth=2)
    ax.plot([pil_centre[1],pil_centre[1]],[0,pil_size[0]],'k:',linewidth=2)
    ax.plot([0,pil_size[1]],[pil_centre[0],pil_centre[0]],'k:',linewidth=2)
    
    if show_peakregion:
        pry1,prx1,pry2,prx2 = peakregion
        print('peakregion = [{},{},{},{}]'.format(*peakregion))
        ax.plot([prx1,prx2,prx2,prx1,prx1],[pry1,pry1,pry2,pry2,pry1],'y-',linewidth=2)
    
    if show_ROIbkg:
        "---Background ROI---" "FROM fun pilroi()"
        " Determine background of image and subtract from the region by pixel"
        " The background region is twice the area of the required region"
        idxbi = [ROIcen[0]-ROIsize[0],ROIcen[0]+ROIsize[0]]
        idxbj = [ROIcen[1]-ROIsize[1],ROIcen[1]+ROIsize[1]]
        
        " Check the box is within the detector"
        if idxbi[0] < 0: idxbi[1] = idxbi[1] - idxbi[0]; idxbi[0] = 0
        if idxbi[1] > pil_size[0]: idxbi[0] = idxbi[0] - (idxbi[1]-pil_size[0]); idxbi[1] = pil_size[0]
        if idxbj[0] < 0: idxbj[1] = idxbj[1] - idxbj[0]; idxbj[0] = 0
        if idxbj[1] > pil_size[1]: idxbj[0] = idxbj[0] - (idxbj[1]-pil_size[1]); idxbj[1] = pil_size[1]
        print( 'Background ROI = [{0},{1},{2},{3}]'.format(idxbi[0],idxbj[0],idxbi[1],idxbj[1]))
        idxbi = np.asarray(idxbi); idxbj = np.asarray(idxbj);
        ax.plot(idxbj[[0,1,1,0,0]],idxbi[[0,0,1,1,0]],'r-',linewidth=2)
    
    ax.set_aspect('equal')
    ax.autoscale(tight=True)
    
    " Create slider on plot"
    axsldr = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    
    if pil_size[0] > 195:
        "2M pilatus"
        ax.set_position([0.25,0.16,0.5,0.7])
        axsldr.set_position([0.15,0.06,0.65,0.03])
    
    sldr = plt.Slider(axsldr, varx, 1, vol.shape[2], \
            valinit=imnum, valfmt = '%0.0f')
    txt = plt.xlabel('{0} = {1} [{2:1.0f}]'.format(varx,x[imnum],imnum),\
                   fontsize=18 )
    
    " Slider update function"
    def update(val):
        "Update function for pilatus image"
        imgno = int(round(sldr.val))
        p.set_data(vol[:,:,imgno-1]+log_colors)
        txt.set_text('{0} = {1}'.format(varx,x[imgno-1]))
        #p.set_clim(cax)
        plt.draw()
        #fig.canvas.draw()
    sldr.on_changed(update)
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot(ttl+'_'+str(imnum))

def plotmeta(scans, fields='Energy', use_time=False, plot_against=None):
    """
    Plot a single value of metadata against run number or time
    """
    
    # X-axis
    if use_time:
        x = getmeta(scans,'TimeSec')
        x = x - x[0]
        xlab = 'Time [s]'
    elif plot_against is not None:
        x = getmeta(scans,plot_against)
        xlab = plot_against
    else:
        x = scans
        xlab = 'Scan number'
    
    # create figure
    fig = plt.figure(figsize=[10,8], dpi=fig_dpi)
    
    # multiple fields
    fields = np.asarray(fields).reshape(-1)
    for field in fields:
        # y-axis
        y = getmeta(scans, field)
        plt.plot(x,y,'-o',linewidth=3,label=field)
    
    if len(fields) > 1:
        plt.legend(loc=0, frameon=False, fontsize=18)
    plt.xlabel(xlab, fontsize=18)
    plt.ylabel(fields[0], fontsize=18)
    plt.suptitle(scantitle(scans), fontsize=20)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ticklabel_format(useOffset=False)
    plt.ticklabel_format(style='sci',scilimits=(-3,3))
    
    # set border space outside axes
    plt.subplots_adjust(left=0.15,bottom=0.12)

def plotdwn(num,save=None):
    "Default plot of I16 data, plotscan(#), or plotscan(#,save=1)"
    # Load data
    d = readscan(num)
    # Get principle variable data
    cmd = d.metadata.cmd_short # Scan command
    # Get variables
    varx = cmd.split()[1]
    vary = cmd.split()[-1]
    if vary[0:3] == 'roi':
        vary = vary + '_sum'
    x = getattr(d,varx)
    y = getattr(d,vary)
    
    #Create plot
    #dnp.plot.plot(x,y,title='#{0}: {1}'.format(num,d.metadata.cmd),name = '#{0}'.format(num))
    return

def plotscans(scans=[],depvar=None,vary='',varx='',fit=None,norm=True,logplot=False,save=False):
    """
    Plot multiple scans of I16 data, plotscans([#1,#2,...])
    
    Plotting with fits: plotscans([#1,#2,...],fit='Gauss'), fit can also be 'Lorentz' or 'pVoight'
    
    Good for plotting many scans as a continuous colormap is used to vary the line colours.
    
    """
    
    "---Create plot---"
    fig = plt.figure(figsize=[10,8], dpi=fig_dpi)
    
    " Define Colourmap"
    cm_plot = plt.get_cmap('rainbow')
    cm_fit = plt.get_cmap('spring') 
    depvals = getmeta(scans,depvar)
    colrange = (depvals-depvals.min())/(depvals.max()-depvals.min())
    
    "---Load data---"
    for n,num in enumerate(scans):
        try:
            d = readscan(num)
            if d is None: print( 'File for run #{} does not exist!'.format(num) ); return 
        except:
            d = num
            num = d.metadata.SRSRUN
        
        " Get metadata"
        m = d.metadata
        cmd = m.cmd_short # Scan command
        sampsl = '{0:4.2g}x{1:<4.2g}'.format(m.s5xgap,m.s5ygap)
        detsl = '{0:4.2g}x{1:<4.2g}'.format(m.s6xgap,m.s6ygap)
        atten1 = '{0:1.0f}'.format(m.Atten)
        if depvar not in m.keys():
            lbl = str(m.SRSRUN)
        else:
            lbl = '{}, {} = {}'.format(str(m.SRSRUN),depvar,getattr(m,depvar))
        
        " Get data" 
        x,y,dy,varxnew,varynew,ttl = getdata(d,vary=vary,varx=varx,norm=norm)[:6]
        
        " Plot data"
        plt.errorbar(x,y,dy,fmt='-o',c=cm_plot(colrange[n]),linewidth=2,label=lbl)
        
        "---Fit Data---"
        if fit != None:
            out,err = peakfit(x,y,type=fit,Nloop=100,peaktest=-100,interpolate=True,disp=True)
            
            xfit,yfit=out['x'],out['y']
            amp,damp = out['Peak Height'],err['Peak Height']
            cen,dcen = out['Peak Centre'],err['Peak Centre']
            wid,dwid = out['FWHM'],err['FWHM']
            bkg,dbkg = out['Background'],err['Background']
            ara,dara = out['Area'],err['Area']
            
            " add to plot"
            ax = plt.gca()
            #txtcol = 'k'
            ax.plot(xfit,yfit,':',linewidth=2,c=cm_fit(colrange[n]),label='#{} Fit'.format(str(m.SRSRUN)))
    
    plt.xlabel(varxnew, fontsize=18)
    plt.ylabel(varynew, fontsize=18)
    #plttl = ttl+'\n'+cmd+'\nss ={}, ds ={}, atten = {}'.format(sampsl,detsl,atten1)
    #plttl = '#{} - #{}'.format(scans[0],scans[-1])
    plttl = scantitle(scans)
    plt.suptitle(plttl, fontsize=16)
    
    if logplot: plt.gca().set_yscale(u'log')
    #if len(varys) > 0 and subtract==True: plt.legend([lgd],loc='best')
    fig.subplots_adjust(left=0.15)
    fig.subplots_adjust(top=0.85)
    
    # Change formats in x & y axes so they are nicer
    plt.gca().get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%8.3g'))
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    
    if n < 6:
        # Legend
        plt.legend(loc='best')
        #plt.legend(loc='upper left')
        #plt.legend(loc='upper right')
    else:
        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cm_plot)
        sm.set_array(depvals)
        cbar = plt.colorbar(sm)
        cbar.set_label(depvar,fontsize=24)
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot(ttl)
    plt.show()

def plotscans3D(scans,depvar='Ta',vary='',varx='',norm=True,logplot=False,save=False):
    """
    Plot 3D eta scans of energy or temp dependence
        scans = list of scan numbers
        depvar = device that varies between scans
        varx = device that was scanned, e.g. 'eta'
        vary = device that was recorded, e.g. 'roi2_sum'
        norm = True/False, normalisation option
        sort = False/True, if true, sorts the scans by depvar
        logplot = False/True, if true, logs the vary data
        save = False/True, if true, saves a png image of the figure to the savedir directory
    """
    
    "---Create Figure---"
    fig = plt.figure(figsize=[14,12], dpi=fig_dpi)
    ax = fig.add_subplot(111,projection='3d')
    
    "-----Load & Plot-----"
    if type(scans) == int:
        "2D scan"
        "e.g. scan sx -1.0 1.0 0.05 sy -1.0 1.0 0.05 BeamOK pil100k 0.1 roi1 roi2"
        d = readscan(scans)
        cmds = d.metadata.cmd.split()
        if depvar not in d.keys():
            depvar = cmds[1]
        if varx not in d.keys():
            varx = cmds[5]
        x,y,dy,labvarx,labvary,ttl,d = getdata(scans,varx=varx,vary=vary,norm=norm)
        z,y2,dy2,labvarx2,labvary2,ttl2,d2 = getdata(scans,varx=depvar,vary=vary,norm=norm)
        ax.plot(z,x,y)
        scans = [scans]
    else:
        for n,run in enumerate(scans):
            x,y,dy,labvarx,labvary,ttl,d = getdata(scans[n],vary=vary,varx=varx,norm=norm)
            if logplot==True: 
                y = np.log10(y)
                labvary='log10('+labvary+')'
            if depvar in d.keys():
                z = getattr(d,depvar)
            elif depvar in d.metadata.keys():
                z = np.ones(x.shape)*getattr(d.metadata,depvar)
            else:
                continue
            
            ax.plot(z,x,y)
    
    # Axis labels
    ax.set_xlabel(depvar, fontsize=18)
    ax.set_ylabel(labvarx, fontsize=18)
    ax.set_zlabel(labvary, fontsize=18)
    #ax.set_zscale('log')
    
    ttl = scantitle(scans)
    ax.set_title(ttl,fontsize=18)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} Scans {1:1.0f}-{2:1.0f} 3D'.format(depvar,scans[0],scans[-1]))
    plt.show()

def plotscans2D(scans,depvar='Ta',vary='',varx='',norm=True,sort=False,logplot=False,save=False):
    """
    Plot pcolor of multiple scans
        scans = list of scan numbers
        depvar = device that varies between scans
        varx = device that was scanned, e.g. 'eta'
        vary = device that was recorded, e.g. 'roi2_sum'
        norm = True/False, normalisation option
        sort = False/True, if true, sorts the scans by depvar
        logplot = False/True, if true, logs the vary data
        save = False/True, if true, saves a png image of the figure to the savedir directory
    """
    
    "-----Loading-----"
    x,y,z,varx,vary,varz,ttl = joindata(scans,varx,depvar,vary,norm,sort)
    
    # Create Plot
    fig = plt.figure(figsize=[14,12], dpi=fig_dpi)
    ax = plt.subplot(1,1,1)
    if logplot:
        plt.pcolor(x,y,np.log10(z))
        varz='log10('+varz+')'
    else:
        plt.pcolor(x,y,z)
    #plt.axis('image')
    plt.axis('tight')
    cb = plt.colorbar()
    # Axis labels
    ax.set_xlabel(varx, fontsize=18)
    ax.set_ylabel(vary, fontsize=18)
    cb.set_label(varz, fontsize=18)
    plt.suptitle(ttl,fontsize=20)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} Scans {1:1.0f}-{2:1.0f} 2D'.format(depvar,scans[0],scans[-1]))
    plt.show()
    
def plotscansSURF(scans,depvar='Ta',vary='',varx='',norm=True,sort=False,logplot=False,save=False):
    """
    Plot 3D surface built from data from multiple scans
        scans = list of scan numbers
        depvar = device that varies between scans
        varx = device that was scanned, e.g. 'eta'
        vary = device that was recorded, e.g. 'roi2_sum'
        norm = True/False, normalisation option
        sort = False/True, if true, sorts the scans by depvar
        logplot = False/True, if true, logs the vary data
        save = False/True, if true, saves a png image of the figure to the savedir directory
    """
    
    "-----Loading-----"
    x,y,z,varx,vary,varz,ttl = joindata(scans,varx,depvar,vary,norm,sort)
    
    # Create Plot
    fig = plt.figure(figsize=[14,12], dpi=fig_dpi)
    ax = plt.subplot(1,1,1,projection='3d')
    if logplot:
        surf = ax.plot_surface(y,x,np.log10(z),rstride=1,cstride=1, cmap=plt.cm.jet)
        varz='log10('+varz+')'
    else:
        surf = ax.plot_surface(y,x,z,rstride=1,cstride=1, cmap=plt.cm.jet)
    #cmap=plt.cm.ScalarMappable(cmap=plt.cm.jet)
    #cmap.set_array(slice)
    plt.axis('tight')
    ax.view_init(30, -40) # elev, azim
    #cb = plt.colorbar(surf)
    
    # Axis labels
    ax.set_xlabel(vary, fontsize=18)
    ax.set_ylabel(varx, fontsize=18)
    ax.set_zlabel(varz, fontsize=18)
    plt.suptitle(ttl,fontsize=20)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} Scans {1:1.0f}-{2:1.0f} SURF'.format(depvar,scans[0],scans[-1]))
    plt.show()

def plotpilSURF(num,varx='',ROIcen=None,wid=10,save=False):
    """ 
    Plot the horizontal (chi direction) pixels of a scan's area detector vs the scan direction
        num = scan number
        varx = scan device, '' = automatic
        ROIcen = detector centre [i,j], if None, pilpeak will be used to determine the peak centre
        wid = width along the detector vertical (delta) pixels to sum
        save = if True, saves a png of the figure to the savedir directory
    """
    
    x,y,dy,varx,vary,ttl,d = getdata(num,varx=varx)
    vol = getvol(num)
    if ROIcen is None:
        ROIcen,frame = pilpeak(vol,disp=True)
    
    deltavals = range(ROIcen[1]-wid//2,ROIcen[1]+wid//2)
    slice = np.sum(vol[:,deltavals,:],1)
    X,Y = np.meshgrid(x,range(vol.shape[0]))
    
    
    fig = plt.figure(figsize=[14,12], dpi=fig_dpi)
    ax = plt.subplot(1,1,1,projection='3d')
    surf = ax.plot_surface(X,Y,slice,rstride=1,cstride=1, cmap=plt.cm.jet)
    #cmap=plt.cm.ScalarMappable(cmap=plt.cm.jet)
    #cmap.set_array(slice)
    plt.axis('tight')
    cb = plt.colorbar(surf)
    
    # Axis labels
    ax.set_xlabel(varx, fontsize=18)
    ax.set_ylabel('Pilatus Horizontal pixel', fontsize=18)
    ax.set_zlabel(vary, fontsize=18)
    cb.set_label(vary, fontsize=18)
    plt.suptitle(ttl,fontsize=20)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} PILSURF'.format(num))
    plt.show()

def plotpilhist(num, frame=None, bins=100, save=False):
    """
    Generate histogram of pixel values in area detector
    """

    x,y,dy,varx,vary,ttl,d = getdata(num)
    vol = getvol(num)
    
    if frame is None:
        vol = vol.reshape(-1)
        ttl += '\n {:s} = {:1.4g} - {:1.4g}'.format(varx, x[0], x[-1])
    else:
        vol = vol[:,:,frame].reshape(-1)
        ttl += '\n {:s} = {:1.4g}'.format(varx, x[frame])

    plt.figure(figsize=[12,10], dpi=fig_dpi)
    plt.hist(np.log10(vol+1), bins)
    labels(ttl, 'Pixel Intensity', 'No. Pixels')
    plt.yscale('log')
    
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} PILHist'.format(num))
    plt.show()


def plotpilhkl_cuts(num,hkl_centre=None,image_centre=None,sum_tolarance=[0.05,0.05,0.05],max_points=301,fit_type='Gauss'):
    """
    Generate h,k,l cuts through a pilatus image scan
        plotpilhkl_cuts(num, [h,k,l], None, [0.01,0.01,0.01],301
    Inputs:
      num = scan number
      image_centre = [i,j,k] or None* - centre of the cuts on the detector e.g. [104,205,31]
      hkl_centre = [h,k,l] or None* - centre of the cuts in hkl, if None the peak position will be used
      sum_tolarance = [dh,dk,dl] - distance around hkl position to sum, in reciprocal lattice units
      max_points = 100 - maximum number of points in each cut (reduces calculation time)
    
    If image_centre is None, the centre of the detector and scan will be used.
    
    At a choosen central (h,k,l), 3 cuts are generated:
            Cut 1: (H, k-dk:k+dk, l-dl:l+dl)
            Cut 2: (h-dh:h+dh, K, l-dl:l+dl)
            Cut 3: (h-dh:h+dh, k-dk:k+dk, L)
    The step dh,dk,dl are determined by the maximum step size between pixels
    At each point H/K/L along each cut, pixels matching the following are summed:
            Pixel sum at H:
                abs(h_pixels-H) < H_step/2
                abs(k_pixels-k) < dk
                abs(l_pixels-l) < dl
            Where H_step is defined by the number of cut_points
    cut_points defines the number of points along axis H/K/L, generated from the minimum to maximum h,k,l
    
    A single figure with 3 subplots is created.
    """
    
    [h_list,h_scan],[k_list,k_scan],[l_list,l_scan] = bin_pixel_hkl_cut(num,hkl_centre,image_centre,sum_tolarance,max_points)
    h_centre = np.mean(h_list)
    k_centre = np.mean(k_list)
    l_centre = np.mean(l_list)
    
    # Fitting
    if fit_type is not None:
        # Fit h
        print('Binned hkl fitting #{}'.format(num))
        print('h fit:')
        fit_h,err_h = peakfit(h_list,h_scan,type=fit_type,disp=True)
        print('\nk fit:')
        fit_k,err_k = peakfit(k_list,k_scan,type=fit_type,disp=True)
        print('\nl fit:')
        fit_l,err_l = peakfit(l_list,l_scan,type=fit_type,disp=True)
        
    
    """
     # Load hkl matrices
    hhh,kkk,lll = pixel2hkl(num)
    print('Max hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.max(),kkk.max(),lll.max()))
    print('Min hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.min(),kkk.min(),lll.min()))
    
    # Load pilatus volume
    vol = getvol(num)
    
    # Load data
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    # h,k,l min tol
    min_tol_h = np.max(np.abs(np.diff(hhh)))
    min_tol_k = np.max(np.abs(np.diff(kkk)))
    min_tol_l = np.max(np.abs(np.diff(lll)))
    print('Maximum hkl step per pixel = [%1.3g, %1.3g, %1.3g]' % (min_tol_h,min_tol_k,min_tol_l))
    
    # Centre of cuts
    if hkl_centre is None:
        if image_centre is None:
            i,j = pil_centre
            k = len(x)//2
        else:
            i,j,k = image_centre
        h_centre = hhh[i,j,k]
        k_centre = kkk[i,j,k]
        l_centre = lll[i,j,k]
        print('Peak Centre: (%1.3g,%1.3g,%1.3g)'%(h_centre,k_centre,l_centre))
    else:
        h_centre,k_centre,l_centre = hkl_centre
    
    # Pixels close to centre
    h_cen_idx = np.abs(hhh-h_centre) < sum_tolarance[0]
    k_cen_idx = np.abs(kkk-k_centre) < sum_tolarance[1]
    l_cen_idx = np.abs(lll-l_centre) < sum_tolarance[2]
    
    kl_cen_idx = k_cen_idx*l_cen_idx
    hl_cen_idx = h_cen_idx*l_cen_idx
    hk_cen_idx = h_cen_idx*k_cen_idx
    
    # cut ranges
    h_list = np.linspace(hhh[kl_cen_idx].min(),hhh[kl_cen_idx].max(),cut_points[0])
    k_list = np.linspace(kkk[hl_cen_idx].min(),kkk[hl_cen_idx].max(),cut_points[1])
    l_list = np.linspace(lll[hk_cen_idx].min(),lll[hk_cen_idx].max(),cut_points[2])
    h_step = h_list[1]-h_list[0]
    k_step = k_list[1]-k_list[0]
    l_step = l_list[1]-l_list[0]
    
    h_scan = np.zeros(len(h_list))
    k_scan = np.zeros(len(k_list))
    l_scan = np.zeros(len(l_list))
    print('Binning h axis at (H,%1.3g,%1.3g) in %1.0f steps, summing <%1.0f pixels per step' % (k_centre,l_centre,len(h_list),np.sum(kl_cen_idx)))
    for n in range(len(h_list)):
        hval = h_list[n]
        hidx = np.abs(hhh-hval) < h_step/2
        small_vol = vol[ hidx*kl_cen_idx ]
        h_scan[n] = np.sum(small_vol)
    print('Binning k axis at (%1.3g,K,%1.3g) in %1.0f steps, summing <%1.0f pixels per step' % (h_centre,l_centre,len(k_list),np.sum(hl_cen_idx)))
    for n in range(len(k_list)):
        kval = k_list[n]
        kidx = np.abs(kkk-kval) < k_step/2
        small_vol = vol[ kidx*hl_cen_idx ]
        k_scan[n] = np.sum(small_vol)
    print('Binning l axis at (%1.3g,%1.3g,L) in %1.0f steps, summing <%1.0f pixels per step' % (h_centre,k_centre,len(l_list),np.sum(hk_cen_idx)))
    for n in range(len(l_list)):
        lval = l_list[n]
        lidx = np.abs(lll-lval) < l_step/2
        small_vol = vol[ lidx*hk_cen_idx ]
        l_scan[n] = np.sum(small_vol)
    """
    plt.figure(figsize=[8,12], dpi=fig_dpi)
    plt.subplot(311)
    plt.errorbar(h_list,h_scan,np.sqrt(h_scan+1),fmt='-o',lw=2,ms=12)
    labels(None,'(h,0,0)','Intensity')
    if fit_type is not None:
        plt.plot(fit_h['x'],fit_h['y'],'r-')
        plt.text(0.55,0.5, fit_h['Results'],transform=plt.gca().transAxes,fontsize=12)
    plt.subplot(312)
    plt.errorbar(k_list,k_scan,np.sqrt(k_scan+1),fmt='-o',lw=2,ms=12)
    labels(None,'(0,k,0)','Intensity')
    if fit_type is not None:
        plt.plot(fit_k['x'],fit_k['y'],'r-')
        plt.text(0.55,0.5, fit_k['Results'],transform=plt.gca().transAxes,fontsize=12)
    plt.subplot(313)
    plt.errorbar(l_list,l_scan,np.sqrt(l_scan+1),fmt='-o',lw=2,ms=12)
    labels(None,'(0,0,l)','Intensity')
    if fit_type is not None:
        plt.plot(fit_l['x'],fit_l['y'],'r-')
        plt.text(0.55,0.5, fit_l['Results'],transform=plt.gca().transAxes,fontsize=12)
    ttl = scantitle(num)
    add_title = '\n(%6.3f,%6.3f,%6.3f) +/- (%1.1g,%1.1g,%1.1g)' %(h_centre,k_centre,l_centre,sum_tolarance[0],sum_tolarance[1],sum_tolarance[2])
    plt.suptitle(ttl+add_title,fontsize=20,fontweight='bold')
    plt.subplots_adjust(left=0.2,hspace=0.25,top=0.91)

def plotpilhkl_surf3d(num,cax=None,log_colors=False,initial_image=None):
    """
    Plot pilatus frames as 3D surface in reciprocal space
       plotpilhkl_surf3d(num)
    
    Options:
        cax = [min,max] or None* | Intensity cut offs
        log_colors = True/False* | Log the intensities
        initial_image = int or None | Number of initial image to display, None = half way
    
    pixel2hkl used to generate reciprocal space coordinates, using UB matrix stored
    in metadata. Detector position calibration defined by locals: pilpara
    
    The images are rebinned to allow fast plotting.
    
    A slider is used to move through the images. If the slider is unresponsive, try running the program again.
    """
    
    # Load hkl matrices
    hhh,kkk,lll = pixel2hkl(num)
    
    # Load pilatus volume
    vol = getvol(num)
    vol = rebin(vol,[2,2,1])
    
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    cm = plt.get_cmap('jet')
    nframes = vol.shape[2]
    if initial_image is None:
        initial_image = nframes//2
    slice_h = hhh[::2,::2,initial_image]
    slice_k = kkk[::2,::2,initial_image]
    slice_l = lll[::2,::2,initial_image]
    #slice_h = hhh[:,:,initial_image]
    #slice_k = kkk[:,:,initial_image]
    #slice_l = lll[:,:,initial_image]
    #slice_v = vol[:,:,initial_image]
    
    if cax is None:
        if log_colors:
            cax = [1,vol.max()]
        else:
            md = np.median(vol[:,:,initial_image])
            mx = np.max(vol[:,:,initial_image])
            cmax = md + 10**(0.7*np.log10(mx-md))
            if cmax <= 0: cmax = 1
            cax = [0,cmax]
        print( 'caxis set at [{0:1.3g},{1:1.3g}]'.format(cax[0],cax[1]) )
    
    # Normalise
    norm = plt.Normalize(cax[0],cax[1])
    colours = []
    for imgno in range(nframes):
        colours += [cm(norm(vol[:,:,imgno]))]
    slice_c = colours[imgno]
    
    # Create Plot
    fig = plt.figure(figsize=[14,12], dpi=fig_dpi)
    ax = plt.subplot(1,1,1,projection='3d')
    surf = ax.plot_surface(slice_h,slice_k,slice_l,rstride=2,cstride=2, facecolors=slice_c)
    plt.axis('tight')
    #ax.view_init(30, 60) # elev, azim
    ax.set_position([0.1,0.2,0.8,0.7])
    
    ttl = scantitle(num)
    labels(ttl,'h','k','l',size='normal')
    plt.show()
    
    " Create slider on plot"
    axsldr = plt.axes([0.15, 0.15, 0.65, 0.03], axisbg='lightgoldenrodyellow')
    sldr = plt.Slider(axsldr, varx, 1, vol.shape[2], \
            valinit=initial_image, valfmt = '%0.0f')
    txt = plt.xlabel('{0} = {1} [{2:1.0f}]'.format(varx,x[initial_image],initial_image),\
                   fontsize=18 )
    
    " Slider update function"
    def update(val):
        "Update function for pilatus image"
        imgno = int(round(sldr.val))
        slice_h = hhh[::2,::2,imgno]
        slice_k = kkk[::2,::2,imgno]
        slice_l = lll[::2,::2,imgno]
        slice_c = colours[imgno]
        ax.collections[0].remove() # ax.collections is location of surface object
        #surf.remove()
        #ax.clear()
        ax.plot_surface(slice_h,slice_k,slice_l,rstride=2,cstride=2, facecolors=slice_c)
        txt.set_text('{0} = {1}'.format(varx,x[imgno-1]))
        plt.show()
        #plt.draw()
        #fig.canvas.draw()
    sldr.on_changed(update)

def plotpilxyz_surf3d(num,cax=None,log_colors=False,initial_image=None):
    """
    Plot pilatus frames as 3D surface in reciprocal space
       plotpilxyz_surf3d(num)
    
    Options:
        cax = [min,max] or None* | Intensity cut offs
        log_colors = True/False* | Log the intensities
        initial_image = int or None | Number of initial image to display, None = half way
    
    pixel2xyz used to generate reciprocal space coordinates, using UB matrix stored
    in metadata. Detector position calibration defined by locals: pilpara
    
    The images are rebinned to allow fast plotting.
    
    A slider is used to move through the images. If the slider is unresponsive, try running the program again.
    """
    
    # Load hkl matrices
    xxx,yyy,zzz = pixel2xyz(num)
    
    # Load pilatus volume
    vol = getvol(num)
    vol = rebin(vol,[2,2,1])
    
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    cm = plt.get_cmap('jet')
    nframes = vol.shape[2]
    if initial_image is None:
        initial_image = nframes//2
    slice_x = xxx[::2,::2,initial_image]
    slice_y = yyy[::2,::2,initial_image]
    slice_z = zzz[::2,::2,initial_image]
    
    if cax is None:
        if log_colors:
            cax = [1,vol.max()]
        else:
            md = np.median(vol[:,:,initial_image])
            mx = np.max(vol[:,:,initial_image])
            cmax = md + 10**(0.7*np.log10(mx-md))
            if cmax <= 0: cmax = 1
            cax = [0,cmax]
        print( 'caxis set at [{0:1.3g},{1:1.3g}]'.format(cax[0],cax[1]) )
    
    # Normalise
    norm = plt.Normalize(cax[0],cax[1])
    colours = []
    for imgno in range(nframes):
        colours += [cm(norm(vol[:,:,imgno]))]
    slice_c = colours[imgno]
    
    # Create Plot
    fig = plt.figure(figsize=[14,12], dpi=fig_dpi)
    ax = plt.subplot(1,1,1,projection='3d')
    surf = ax.plot_surface(slice_x,slice_y,slice_z,rstride=2,cstride=2, facecolors=slice_c)
    plt.axis('tight')
    #ax.view_init(30, 60) # elev, azim
    ax.set_position([0.1,0.2,0.8,0.7])
    
    ttl = scantitle(num)
    labels(ttl,'Q$_x$ [$\AA^{-1}$]','Q$_y$ [$\AA^{-1}$]','Q$_z$ [$\AA^{-1}$]',size='normal')
    plt.show()
    
    " Create slider on plot"
    axsldr = plt.axes([0.15, 0.15, 0.65, 0.03], axisbg='lightgoldenrodyellow')
    sldr = plt.Slider(axsldr, varx, 1, vol.shape[2], \
            valinit=initial_image, valfmt = '%0.0f')
    txt = plt.xlabel('{0} = {1} [{2:1.0f}]'.format(varx,x[initial_image],initial_image),\
                   fontsize=18 )
    
    " Slider update function"
    def update(val):
        "Update function for pilatus image"
        imgno = int(round(sldr.val))
        slice_x = xxx[::2,::2,imgno]
        slice_y = yyy[::2,::2,imgno]
        slice_z = zzz[::2,::2,imgno]
        slice_c = colours[imgno]
        print(ax.collections)
        ax.collections[0].remove() # ax.collections is location of surface object
        ax.plot_surface(slice_x,slice_y,slice_z,rstride=2,cstride=2, facecolors=slice_c)
        txt.set_text('{0} = {1}'.format(varx,x[imgno-1]))
        plt.show()
    sldr.on_changed(update)

def plotpiltth(num=None,binsep=0.1,centre_only=False):
    """
    Plot binned two-theta representation of pilatus frames
       plotpiltth(num,binsep,centre_only)
        num = scan number or list of scan numbers for muliplt overlaid plots
        binsep = width of bins
        centre_only = if True, only bins pixels around the centre of the pilatus.
    
    pixel2tth used to generate reciprocal space coordinates, using UB matrix stored
    in metadata. Detector position calibration defined by locals: pilpara
    """
    
    "---Handle inputs---"
    if num is None or num == 0:
        num = latest()
    
    " Multiple nums given"
    nums = np.asarray(num,dtype=int).reshape(-1)
    
    plt.figure(figsize=[10,8], dpi=fig_dpi)
    for num in nums:
        # Load hkl matrices
        tth,ival = pixel2tth(num,centre_only=centre_only)
        tth,ival = bindata(tth,ival,binsep)
        plt.plot(tth,ival,'-o',lw=2,ms=12)
    ttl = scantitle(nums)
    labels(ttl,'Two-Theta [Deg]','Pilatus Intensity',size='normal')
    plt.show()


def plotqbpm(scannos, xvar=None, normaliseby='left', plot=True, show_difference=False, save=False, omit_left=0, omit_right=0):
    """
    Plots 4 currents from QBPM scans
        e.g. scancn ppth1 0.001 w .5 qbpm6

        neg, pos, centre, offset = plotqbpm(123456)
      scannos = scan to plot, multiple scans are appended together
      normaliseby = min, left*, right, mean, none - method of normalisation
      plot = True/ False - create plot
      show_difference=False/ True - display diff on plot
      save = True/ False
      omit_left = 0 number of values to omit on left
      omit_right = 0 number of values to omit of right

      neg = negative offset
      pos = positive offset
      centre = average point between neg and pos
      offset = offset value from cenre to pos
    """
    
    scannos = np.asarray(scannos).reshape(-1)
    if xvar is None:
        xvar = auto_varx(scannos[0])
    xval = np.empty(0)
    mon = np.empty(0)
    C1=np.empty(0)
    C2=np.empty(0)
    C3=np.empty(0)
    C4=np.empty(0)
    for scan in scannos:
        d = readscan(scan)
        norm = d.ic1monitor/(d.rc/exp_ring_current)
        ln = len(norm) - omit_right
        norm = norm[omit_left:ln]
        
        xval = np.append(xval, getattr(d, xvar)[omit_left:ln])
        mon = np.append(mon, norm)
        C1 = np.append(C1, d.C1[omit_left:ln]/norm)
        C2 = np.append(C2, d.C2[omit_left:ln]/norm)
        C3 = np.append(C3, d.C3[omit_left:ln]/norm)
        C4 = np.append(C4, d.C4[omit_left:ln]/norm)
    
    if normaliseby.lower() in ['none']:
        minC1 = 0
        minC2 = 0
        minC3 = 0
        minC4 = 0
        maxC1 = 1
        maxC2 = 1
        maxC3 = 1
        maxC4 = 1
    elif normaliseby.lower() in ['mean']:
        minC1 = C1.min()
        minC2 = np.mean(np.append(C2[:5], C2[-5:]))
        minC3 = C3.min()
        minC4 = np.mean(np.append(C4[:5], C4[-5:]))
        maxC1 = np.mean(np.append(C1[:5], C1[-5:]))
        maxC2 = C2.max()
        maxC3 = np.mean(np.append(C3[:5], C3[-5:]))
    elif normaliseby.lower() in ['right', 'r']:
        minC1 = C1.min()
        minC2 = np.mean(C2[-5:])
        minC3 = C3.min()
        minC4 = np.mean(C4[-5:])
        maxC1 = np.mean(C1[-5:])
        maxC2 = C2.max()
        maxC3 = np.mean(C3[-5:])
        maxC4 = C4.max()
    elif normaliseby.lower() in ['left', 'l', 'start']:
        minC1 = C1.min()
        minC2 = np.mean(C2[:5])
        minC3 = C3.min()
        minC4 = np.mean(C4[:5])
        maxC1 = np.mean(C1[:5])
        maxC2 = C2.max()
        maxC3 = np.mean(C3[:5])
        maxC4 = C4.max()
    else:
        maxC4 = C4.max()
        minC1 = C1.min()
        minC2 = C2.min()
        minC3 = C3.min()
        minC4 = C4.min()
        maxC1 = C1.max()
        maxC2 = C2.max()
        maxC3 = C3.max()
        maxC4 = C4.max()

    C1 = (C1-minC1)/np.max(C1-minC1)
    C2 = (C2-minC2)/np.max(C2-minC2)
    C3 = (C3-minC3)/np.max(C3-minC3)
    C4 = (C4-minC4)/np.max(C4-minC4)
    
    #C1 = C1/maxC1
    #C2 = (C2/minC2) -1
    #C3 = C3/maxC3
    #C4 = (C4/minC4) -1
    
    # interpolate
    ival = np.linspace(np.min(xval), np.max(xval), 100*len(xval))
    iC1 = np.interp(ival, xval, C1)
    iC2 = np.interp(ival, xval, C2)
    iC3 = np.interp(ival, xval, C3)
    iC4 = np.interp(ival, xval, C4)
    
    diff = np.abs((iC1+iC3)/2 - (iC2+iC4)/2)

    if plot:
        fig, ax1 = plt.subplots(figsize=[12,10])
        plt.plot(xval, C1, 'b-', lw=2, label='C1')
        plt.plot(xval, C2, 'r-', lw=2, label='C2')
        plt.plot(xval, C3, 'c-', lw=2, label='C3')
        plt.plot(xval, C4, 'm-', lw=2, label='C4')
        if show_difference:
            plt.plot(ival, diff, 'k-', lw=0.5, label='|C1+C3-C2-C4|')
        cmd = d.metadata.cmd_short
        ttl = scantitle(scannos) + '\n%s\nppchi = %5.1f, ppz1 = %5.2f, ppz2 = %5.2f'%(cmd, d.metadata.ppchi,d.metadata.ppz1, d.metadata.ppz2)
        labels(ttl, xvar, 'QBPM6', legend=True)
        
        ax2 = ax1.twinx()
        plt.plot(xval, mon, 'g:', lw=2, label='ic1monitor')
        labels(None,None,'ic1monitor')
        ax2.tick_params(axis='y', labelcolor='g')
        ax2.set_ylabel('ic1monitor',color='g')
    
    """
    # Find intercepts
    midpoint = np.mean(ival)
    midpoint = xval[np.argmin(mon)]
    step = np.mean(np.diff(xval))
    # Difference
    d1 = np.abs(iC1-iC2)
    d2 = np.abs(iC1-iC4)
    d3 = np.abs(iC3-iC2)
    d4 = np.abs(iC3-iC4)
    
    ivaln = ival[ival<midpoint-step]
    #neg1 = ivaln[np.argmin(d1[ival<midpoint])]
    #neg2 = ivaln[np.argmin(d2[ival<midpoint])]
    #neg3 = ivaln[np.argmin(d3[ival<midpoint])]
    #neg4 = ivaln[np.argmin(d4[ival<midpoint])]
    #neg = np.mean([neg1,neg2,neg3,neg4])
    #neg = ivaln[np.argmin(diff[ival<midpoint-step])]

    ivalp = ival[ival>midpoint+step]
    #pos1 = ivalp[np.argmin(d1[ival>midpoint])]
    #pos2 = ivalp[np.argmin(d2[ival>midpoint])]
    #pos3 = ivalp[np.argmin(d3[ival>midpoint])]
    #pos4 = ivalp[np.argmin(d4[ival>midpoint])]
    #pos = np.mean([pos1,pos2,pos3,pos4])
    pos = ivalp[np.argmin(diff[ival>midpoint+step])]
    """
    # new method 29/2/2020
    # find smallest differences furthest appart
    npoints = 0
    percentile = 0
    while npoints < 2:
        percentile += 1
        minthresh = np.percentile(diff, percentile)
        minxvals = ival[ diff < minthresh ]
        npoints = len(minxvals)
    neg = minxvals[0]
    pos = minxvals[-1]

    avmid = (pos+neg)/2
    negoff = neg-avmid
    posoff = pos-avmid
    midpoint = xval[np.argmin(mon)]
    print(ttl)
    print('Estimated Midpoint: %s=%7.3f'%(xvar,midpoint))
    print('   Lower intercept: %s=%8.4f (%+6.4f)'%(xvar,neg,negoff))
    print('   Upper intercept: %s=%8.4f (%+6.4f)'%(xvar,pos,posoff))
    print('   Actual midpoint: %s=%8.4f'%(xvar,avmid))
    
    if plot:
        plt.sca(ax1)
        plt.axvline(neg, c='k', lw=0.5)
        plt.axvline(pos, c='k', lw=0.5)
    
        "---Save---"
        if save not in [None, False, '']:
            if type(save) is str:
                saveplot(save)
            else:
                saveplot('{0} QBPM6'.format(num))
        plt.show()
    return neg, pos, avmid, posoff
        


"-----------------------Peak Fitting Functions----------------------------"


def FWHM(x,y,interpolate=False):
    """
    Calculate a simple FWHM from a peak
    Finds the peak centre using the maximum values
    finds the width by finding the positions to the left and right of the centre
    at half the height of the maximum value.
    """
    
    if interpolate:
        interx = np.linspace(x[0],x[-1],len(x)*100)
        intery = np.interp(interx,x,y) 
        x, y = interx, intery
    
    mx = max(y)
    ln = len(y)
    
    # Peak position
    pkpos = y.argmax()
    
    # Split into two parts - before and after the peak
    hfxx1 = x[:pkpos+1]
    hfxx2 = x[pkpos:]
    
    # Find the half-max positions
    hfmx1 = abs(y[:pkpos+1]-mx//2)
    hfmx2 = abs(y[pkpos:]-mx//2)
    
    hfpos1 = hfxx1[hfmx1.argmin()]
    hfpos2 = hfxx2[hfmx2.argmin()]
    
    # Return FWHM
    return abs(hfpos2-hfpos1)

def centre(x,y):
    """
    Calcualte centre of a peak using the weighted average of the largest 20% of points
    """
    srt = np.argsort(y)
    cen = np.average( x[ srt[ -len(x)//5: ] ] ,weights=y[ srt[ -len(x)//5: ] ]**2)
    return cen

def straightline(x,grad=1.0,inter=0.0):
    "Staigh line"
    return grad*x + inter

def linefit(x,y,dy=None,disp=False):
    """
    Fit a line to data, y = mx + c
    
    fit,err = linefit(x,y,dy,disp=False)
    x,y = arrays of data to fit
    dy = error bars on each y value (or leave as None)
    disp = True/False display results
    
    fit/err = dicts of results with entries:
        'Gradient'    - the gradient of the line (m)
        'Intercept'   - the intercept (c)
        'x'           - x values, same as x
        'y'           - the fitted y values for each x
    
    Note: a matching line can be generated with:
        y = straightline(x,grad=fit['Gradient'],inter=fit['Intercept'])
    """
    
    # Set dy to 1 if not given
    if dy is None: dy=np.ones(len(y))
    
    # Remove zeros from x - causes errors in covariance matrix
    xold = x
    offset = 0.
    if any(np.abs(x)<0.001):
        print( 'Zero detected - adding 0.001 to x values' )
        offset = 0.001
        x = x + offset
    if any(np.isnan(dy)):
        print( 'Ignoring errors due to NaNs' )
        dy=np.ones(len(y))
    
    # Handle zero intensities
    y[ y<0.01 ] = 0.01
    dy[ dy<0.01 ] = 0.01
    
    # Starting parameters
    grad = 0.0
    inter = np.mean(y)
    
    try:
        vals, covmat = curve_fit(straightline,x,y,[grad,inter],sigma=dy)
    except RuntimeError:
        vals = [0,0]
        covmat = np.diag([np.nan,np.nan])
    
    # Values
    grad = vals[0]
    inter = vals[1]
    # Errors
    perr = np.sqrt(np.diag(covmat))
    dgrad = perr[0]
    dinter = perr[1]
    
    # Calculate fit
    yfit = straightline(xold,grad,inter)
    
    # Calculate CHI^2
    chi = np.sum( (y-yfit)**2 / dy)
    dof = len(y) - 2 # Number of degrees of freedom (Nobs-Npar)
    chinfp = chi/dof
    
    fit,err = {},{}
    fit['Gradient'] = grad
    fit['Intercept'] = inter
    fit['x'] = xold
    fit['y'] = yfit
    err['Gradient'] = dgrad
    err['Intercept'] = dinter
    
    # Print Results
    if disp:
        print( ' ------Line Fit:----- ' )
        print( '  Gradient = {0:10.3G} +/- {1:10.3G}'.format(grad,dgrad) )
        print( ' Intercept = {0:10.3G} +/- {1:10.3G}'.format(inter,dinter) )
        print( '     CHI^2 = {0:10.3G}'.format(chi) )
        print( '  CHI^2 per free par = {0:10.3G}'.format(chinfp) )
    return fit,err

def simpplt(x,height=1,cen=0,FWHM=0.5,bkg=0):
    "Plot an Illustration of simpfit"
    
    minpos = cen-FWHM
    maxpos = cen+FWHM
    y = np.ones(len(x))*bkg
    y[len(x)//5:-len(x)//5] += height/2
    y[np.logical_and(x>minpos, x<maxpos)] += height/2
    return y

def simpfit(x,y,disp=None):
    "Simple peak parameters"
    
    # Starting parameters
    wid = FWHM(x,y,interpolate=True)
    
    #bkgrgn = np.concatenate( (y[:len(x)//5],y[-len(x)//5:]) ) # background method 1 - wrong if peak is off centre
    #bkgrgn = np.percentile(y,range(0,20)) # background method 2 - average lowest 5th of data
    #bkg = np.mean(bkgrgn) 
    h,bin = np.histogram(y,10)
    bincen = (bin[1:] + bin[:-1]) / 2.0
    bkg = bincen[np.argmax(h)]
    amp = max(y) - bkg
    #if amp > 5*bkg:
    #    cen = x[y.argmax()]
    #else:
    #    cen = x[len(x)//2]
    # Alternative centre method 9/2/16
    #srt = np.argsort(y)
    #cen = np.average( x[ srt[ -len(x)//5: ] ] ,weights=y[ srt[ -len(x)//5: ] ])
    # Adapted 3/12/18, square weights to correct for long scans of sharp peaks.
    cen = centre(x,y)

    # Errors
    damp = np.sqrt(amp)
    dwid = abs(x[1]-x[0])
    #dbkg = np.sqrt(np.sum(bkgrgn**2))//len(bkgrgn)
    dbkg = np.sqrt(bkg)
    dcen = dwid
    
    # Integrated area
    scanwid = abs(x[-1]-x[0])
    ara = np.sum(y-bkg)*scanwid/len(x)
    dara = np.sqrt(np.sum(y))*scanwid/len(x)
    
    # Print Results
    if disp is not None:
        print( ' ------Simple Fit:----- ' )
        print( ' Amplitude = {0:10.3G} +/- {1:10.3G}'.format(amp,damp) )
        print( '    Centre = {0:10.3G} +/- {1:10.3G}'.format(cen,dcen) )
        print( '      FWHM = {0:10.3G} +/- {1:10.3G}'.format(wid,dwid) )
        print( 'Background = {0:10.3G} +/- {1:10.3G}'.format(bkg,dbkg) )
        print( '      Area = {0:10.3G} +/- {1:10.3G}'.format(ara,dara) )
    
    return amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara

def gauss(x,height=1,cen=0,FWHM=0.5,bkg=0):
    "Define Gaussian"
    "From http://fityk.nieto.pl/model.html"
    return height*np.exp(-np.log(2)*((x-cen)/(FWHM/2))**2) + bkg

def lorentz(x,height=1,cen=0,FWHM=0.5,bkg=0):
    "Define Lorentzian"
    "From http://fityk.nieto.pl/model.html"
    return height/(1 + ((x-cen)/(FWHM/2))**2 ) + bkg

def pvoight(x,height=1,cen=0,FWHM=0.5,LorFrac=0.5,bkg=0):
    "Define pseudo-Voight"
    "From http://fityk.nieto.pl/model.html"
    HWHM = FWHM/2.0
    ln2 = 0.69314718055994529
    pos = x-cen
    L = LorFrac/( 1 + (pos/HWHM)**2 )
    G = (1-LorFrac)*np.exp( -ln2*(pos/HWHM)**2 )
    return height*(G + L) + bkg

def create_peak_fun(text_fn, params):
    """
    Create a function from a string, return the function
     func = create_peak_fun(text_fn, params)

    text_fn = str function acting on variable 'x'
    params = list of variables in function other than 'x'

    e.g.
      func = create_peak_fun('x**2+y*z', ['y','z'])
    Returns func, which definition:
    def func(x,y,z):
        return x**2+y*z
    """
    inputs = ','.join(params)
    funcstr = 'def func(x,{}):\n    return {}'.format(inputs,text_fn)
    
    fitlocals = {}
    exec(funcstr,globals(),fitlocals) # python >2.7.9
    func = fitlocals['func']
    return func

def peakfit(x,y,dy=None,type='pVoight',bkg_type='flat',peaktest=1,estvals=None,
            Nloop=10,Binit=1e-5,Tinc=2,change_factor=0.5,converge_max = 100,
            min_change=0.01,interpolate=False,debug=False,plot=False,disp=False):
    """ General Peak Fitting function to fit a profile to a peak in y = f(x)
    Allows several possible profiles to be used and can try to find the best estimates for 
    fitting parameters using an RMC-based least-squares routine.
    
    out,err = peakfit(x,y)
    out,err = peakfit(x,y,dy=None,type='pVoight',**fitoptions)
    
    **fitoptions:
    Basic parameters:
        x = array of the dependent variable, e.g. eta, mu, phi
        y = array of the independent variable, e.g. maxval,roi1_sum, APD
        dy = errors on y (default = None)
        type = function type. Allowed: 'pVoight' (default), 'Gauss', 'Lorentz', 'Simple'*
        bkg_type = background type. Allowed: 'flat' (default), 'slope', 'step'
    RMC options:
        Nloop = Number of iterations per temperature, default = 0 (RMC off)**
        Binit = Initial values of 1/kbT used for RMC, default = 1e-3 (lower = Higher temp)
        Tinc = After Nloop steps, the temperature is increased by factor Tinc, default = 2
        change_factor = Each parameter is multiplied by a normal distribution around 1 with width change_factor. (Default = 0.5)
        converge_max = Iterations will end when convergece reaches converge_max. (Default = 100)
        min_change = Relative variation of parameters must be below min_change to increase convergence value. (Default = 0.01)
    Output options:
        interpolate = True: The output fit will have interpolated (much finer) values in x and y. (Default = False)
        debug = True: Output of each iteration is displayed. (Default = False)
        disp = True: The final fitted parameters will be displayed in the command line. (Dafault = False)
    
    Output:
        out = dict with fitted parameters
        err = dict with errors on fitted paramters
        
        out.keys() = ['Peak Height','Peak Centre','FWHM','Lorz frac','Background','Area','CHI**2','CHI2 per dof','x','y']
    
    * selecting type='simple' will not fit the data, just provide a very simple estimation.
    ** Nloop must be set > 0 for the RMC routine to be used, for Nloop=0 or converge_max=0, a simple gradient decend method from a simple estimation is used.
    
    Notes on the RMC routine:
     - see the code
    
    """
    
    # Set dy to 1 if not given
    if dy is None: dy=np.ones(len(y))
    
    # Remove zeros from x - causes errors in covariance matrix
    xold = 1.0*x
    offset = 0.
    if any(np.abs(x)<0.0001):
        print( 'Zero detected - adding 0.0001 to x values' )
        offset = 0.0001
        x = x + offset
    if any(np.isnan(dy)):
        print( 'Ignoring errors due to NaNs' )
        dy=np.ones(len(y))
    
    # Handle zero intensities
    y[ y<0.01 ] = y[ y<0.01 ]+0.01
    dy[ dy<0.01 ] = dy[ dy<0.01 ]+0.01
    
    # Estimate starting parameters
    amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(x,y)
    frac = 0.5
    
    
    '-----------------------------------------------------------'
    '---------------------CHOOSE FUNCTIONS----------------------'
    '-----------------------------------------------------------'
    
    # Define background function
    if bkg_type.lower() in ['slope','sloping','grad']:
        bkgfunc = '+slope*(x-cen)'
        inpvals = ['slope']
        defestvals = [ (y[-1]-y[0])/(xold[-1]-xold[0]) ]
        valnames = ['Background Slope']
        minvals = [-np.inf]
        maxvals = [np.inf]
    elif bkg_type.lower() in ['step']:
        bkgfunc = '+np.append(bkg-step*np.ones(np.floor(len(x)/2.0)),bkg+step*np.ones(np.ceil(len(x)/2.0)))'
        inpvals = ['step']
        defestvals = [ (y[-1]-y[0])/2.0 ]
        valnames = ['Background Step']
        minvals = [-np.inf]
        maxvals = [np.inf]
    else:
        bkgfunc = ''
        inpvals = []
        defestvals = []
        valnames = []
        minvals = []
        maxvals = []
        
    
    # Define starting parameters for choosen function
    if type.lower() in ['gauss','gaussian','g']:
        txtfunc = 'height*np.exp(-np.log(2)*((x-cen)/(FWHM/2))**2)+bkg'+bkgfunc
        inpvals = ['height','cen','FWHM','bkg']+inpvals
        fitfunc = create_peak_fun(txtfunc,inpvals)
        defestvals = [amp,cen,wid,bkg]+defestvals
        valnames= ['Peak Height','Peak Centre','FWHM','Background']+valnames
        minvals = [np.std(y),min(x),abs(x[1]-x[0]),-np.inf]+minvals
        maxvals = [5*amp,max(x),2*(max(x)-min(x)),np.inf]+maxvals
    elif type.lower() in ['lorz','lorentz','lorentzian','l']:
        txtfunc = 'height/(1 + ((x-cen)/(FWHM/2))**2 )+bkg'+bkgfunc
        inpvals = ['height','cen','FWHM','bkg']+inpvals
        fitfunc = create_peak_fun(txtfunc,inpvals)
        defestvals = [amp,cen,wid,bkg]+defestvals
        valnames= ['Peak Height','Peak Centre','FWHM','Background']+valnames
        minvals = [0,min(x),abs(x[1]-x[0]),-np.inf]+minvals
        maxvals = [5*amp,max(x),2*(max(x)-min(x)),np.inf]+maxvals
    elif type.lower() in ['simp','simple','basic','s','sum','total','max','maxval','maximum']:
        fitfunc = simpplt
        defestvals = [amp,cen,wid,bkg]
        valnames= ['Peak Height','Peak Centre','FWHM','Background']
        minvals = [0,min(x),abs(x[1]-x[0]),-np.inf]
        maxvals = [5*amp,max(x),2*(max(x)-min(x)),np.inf]
    else:
        txtfunc = 'height*( LorFrac/( 1.0 + (2.0*(x-cen)/FWHM)**2 ) + (1.0-LorFrac)*np.exp( -np.log(2)*(2.*(x-cen)/FWHM)**2 ) )+bkg'+bkgfunc
        inpvals = ['height','cen','FWHM','LorFrac','bkg']+inpvals
        fitfunc = create_peak_fun(txtfunc,inpvals)
        defestvals = [amp,cen,wid,frac,bkg]+defestvals
        valnames= ['Peak Height','Peak Centre','FWHM','Lorz frac','Background']+valnames
        minvals = [0,min(x),abs(x[1]-x[0]),-0.5,-np.inf]+minvals
        maxvals = [5*amp,max(x),2*(max(x)-min(x)),2,np.inf]+maxvals
    
    if estvals is None:
        estvals = defestvals[:]
    
    '-----------------------------------------------------------'
    '-------------------------FIT DATA--------------------------'
    '-----------------------------------------------------------'
    # Fitting not reuqired
    if type.lower() in ['simp','simple','basic','s']:
        amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(xold,y)
        if ara < 0: ara=0
        fitvals = [amp,cen,wid,bkg]
        errvals = [damp,dcen,dwid,dbkg]
        chi=0
    
    elif type.lower() in ['sum','total']:
        amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(xold,y)
        ara = y.sum()
        dara = np.sqrt(ara)
        fitvals = [amp,cen,wid,bkg]
        errvals = [damp,dcen,dwid,dbkg]
        chi=0
    
    elif type.lower() in ['max','maxval','maximum']:
        amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(xold,y)
        ara = y.max()
        dara = np.sqrt(ara)
        cen = xold[y.argmax()]
        fitvals = [amp,cen,wid,bkg]
        errvals = [damp,dcen,dwid,dbkg]
        chi=0
    
    # Perform fitting
    else:
        # Check if a peak exists to fit
        peak_rat = ispeak(y,dy,test=peaktest,disp=False,return_rat=True)
        if debug: print( 'Peak ratio: {:1.2g} ({:1.2g})'.format(peak_rat,peaktest) )
        if peak_rat < peaktest:
            if debug: print( 'No peak here (rat={:1.2g}). Fitting background instead!'.format(peak_rat) )
            amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(xold,y)
            type = 'Background'
            fitfunc = straightline
            valnames = ['Slope','Background']
            estvals = [0,bkg]
            minvals = [-np.inf,-np.inf]
            maxvals = [np.inf,np.inf]
        
        # Perform fitting
        # Initial Fit (but don't update the estimators yet)
        try:
            fitvals, covmat = curve_fit(fitfunc,x,y,estvals,sigma=dy,absolute_sigma=True)
        except RuntimeError:
            if debug: print( 'Initial fit failed!' )
            fitvals = 1*estvals
            covmat = np.nan*np.eye(len(estvals))
        yfit = fitfunc(xold,*fitvals) # New curve
        chi = np.sum( (y-yfit)**2 / dy) # Calculate CHI^2
        if debug: print( 'Initial Fit CHI**2 = ',chi )
        
        # Check errors are reasonable
        errvals = np.sqrt(np.diag(covmat))
        if any(np.isnan(errvals)):
            chi = np.inf
        
        if debug: print( 'Estimates: ',estvals )
        if debug: print( 'Initial Fit: ',list(fitvals),'CHI**2 = ',chi )
        
        # Check new values are reasonable
        for n,val in enumerate(fitvals):
            if val < minvals[n] or val > maxvals[n]:
                if debug: print( 'Initial value out of range: {} = {} ({}:{})'.format(valnames[n],val,minvals[n],maxvals[n]) )
                chi = np.inf # will not accept change if fitvalues fall out of range
        
        "----------------RMC-------------------"
        changes = np.zeros(len(estvals))
        converge = 0
        Ntemp = 0
        while converge < converge_max:
            beta = Binit*Tinc**Ntemp
            if debug: print( 'New Temperature: ',Ntemp,beta )
            Ntemp += 1
            if Ntemp > Nloop:
                break
            for MCloop in range(Nloop):
                ini_estvals = 1*estvals # 1*estvals copies the array rather than links to it!
                if debug: print( Ntemp,MCloop,'Current estimates: ',list(ini_estvals) )
                # Loop over each estimator and randomly vary it
                for estn in range(len(estvals)):
                    inc_factor = np.random.normal(1,change_factor)
                    est_new = 1*estvals
                    est_new[estn] = est_new[estn]*inc_factor
                    if debug: print( '\tNew {} = {}'.format(valnames[estn],est_new[estn]) )
                    try:
                        fitvals, covmat = curve_fit(fitfunc,x,y,est_new,sigma=dy,absolute_sigma=True)
                    except RuntimeError:
                        if debug: print( beta,MCloop,estn,'Fit failed.' )
                        continue
                    yfit = fitfunc(xold,*fitvals) # New curve
                    chi_new = np.sum( (y-yfit)**2 / dy) # Calculate CHI^2
                    
                    # Check errors are reasonable
                    errvals = np.sqrt(np.diag(covmat))
                    if any(np.isnan(errvals)) or any(np.isinf(errvals)):
                        chi_new = np.inf
                    
                    # Check new values are reasonable
                    for n,val in enumerate(fitvals):
                        #if debug: print( beta,MCloop,estn,'CheckVal: ',n,val,minvals[n],maxvals[n] )
                        if val < minvals[n] or val > maxvals[n]:
                            if debug: print( '\t\tValue out of range: {} = {} ({}:{})'.format(valnames[n],val,minvals[n],maxvals[n]) )
                            chi_new = np.inf # will not accept change if fitvalues fall out of range
                    if debug: print( '\tFits: {}'.format(list(fitvals)) )
                    if debug: print( '\tErrors: {}'.format(list(errvals)) )
                    if debug: print( '\tCHI**2: {}'.format(chi_new) )
                    
                    # Metropolis Algorithm
                    if chi_new < chi or np.exp(beta*(chi-chi_new)) > np.random.rand():
                        if debug: print( '\tFits Kept!' )
                        estvals = 1*fitvals # = 1*est_new
                        chi = 1*chi_new
                        changes[estn] += 1
                
                # Track changes
                chvals = np.divide(np.abs(np.subtract(estvals,ini_estvals)),ini_estvals)
                if np.any(chvals > min_change):
                    converge = 0
                else:
                    converge += 1
                
                if debug: print( beta,MCloop,chi,'Changes: ',changes,chvals,converge )
                
                # break the loop if the solution has converged
                if converge >= converge_max:
                    if debug: print( 'Fit converged in {} temps!'.format(Ntemp-1) )
                    break
        
        # After the loop, perform a final check
        try:
            fitvals, covmat = curve_fit(fitfunc,x,y,estvals,sigma=dy,absolute_sigma=True)
        except RuntimeError:
            fitvals = 1*estvals
            fitvals[0] = 0.0
            covmat = np.nan*np.eye(len(estvals))
        
        errvals = np.sqrt(np.diag(covmat))
    
    # Check fit has worked
    if any(np.isnan(errvals)) or chi == np.inf:
        print( 'Fit didnt work: use summation instead' )
        amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(xold,y)
        if ara < 0: ara=0
        type = 'Simple'
        fitfunc = simpplt
        valnames = ['Peak Height','Peak Centre','FWHM','Background']
        fitvals = [amp,cen,wid,bkg]
        errvals = [damp,dcen,dwid,dbkg]
    
    # create output dict
    output = dict(zip(valnames,fitvals))
    outerr = dict(zip(valnames,errvals))
    
    '-----------------------------------------------------------'
    '-------------Calulate area (profile dependent)-------------'
    '-----------------------------------------------------------'
    if type.lower() in ['gauss','gaussian','g']:
        output['Lorz frac'] = 0.0
        outerr['Lorz frac'] = 0.0
        output['Peak Height'] = abs(output['Peak Height'])
        output['FWHM'] = abs(output['FWHM'])
        amp,damp = abs(fitvals[0]),abs(errvals[0])
        wid,dwid = abs(fitvals[2]),abs(errvals[2])
        sig = wid/(2*np.sqrt(2*np.log(2))) # Gaussian sigma
        dsig = dwid/((2*np.sqrt(2*np.log(2))))
        ara = np.abs(amp*sig*np.sqrt(2*np.pi))
        dara = ara*np.sqrt( (damp/amp)**2 + (dsig/sig)**2 )
    elif type.lower() in ['lorz','lorentz','lorentzian','l']:
        output['Lorz frac'] = 1.0
        outerr['Lorz frac'] = 0.0
        output['Peak Height'] = abs(output['Peak Height'])
        output['FWHM'] = abs(output['FWHM'])
        amp,damp = abs(fitvals[0]),abs(errvals[0])
        wid,dwid = abs(fitvals[2]),abs(errvals[2])
        ara = np.pi*amp*wid/2
        dara = ara*np.sqrt( (damp/amp)**2 + (dwid/wid)**2 )
    elif type.lower() in ['simp','simple','basic','s','max','maxval','maximum','sum','total']:
        output['Lorz frac'] = -1.0
        outerr['Lorz frac'] = 0.0
    elif type.lower() in ['voight','pvoight','pseudovoight','v']:
        output['Peak Height'] = abs(output['Peak Height'])
        output['FWHM'] = abs(output['FWHM'])
        amp,damp = abs(fitvals[0]),abs(errvals[0])
        wid,dwid = abs(fitvals[2]),abs(errvals[2])
        frac,dfrac = fitvals[3],errvals[3]
        
        # Calculated Voight area = Gaussian + Voight
        sig = wid/(2*np.sqrt(2*np.log(2))) # Gaussian sigma
        dsig = dwid/((2*np.sqrt(2*np.log(2))))
        Gara = np.abs(amp*sig*np.sqrt(2*np.pi))
        Lara = np.pi*amp*wid/2
        ara = frac*Lara + (1-frac)*Gara
    
        # Error on area
        dGara = Gara*np.sqrt( (damp/amp)**2 + (dsig/sig)**2 )
        dLara = Lara*np.sqrt( (damp/amp)**2 + (dwid/wid)**2 )
        dVara1= (1-frac)*Gara*np.sqrt( (dfrac/(1-frac))**2 + (dGara/Gara)**2 )
        dVara2= frac*Lara*np.sqrt( (dfrac/frac)**2 + (dLara/Lara)**2 )
        dara = np.sqrt( dVara1**2 + dVara2**2 )
    elif type.lower() in ['background']:
        output['Lorz frac'] = np.nan
        outerr['Lorz frac'] = 0.0
        output['Peak Height'] = np.nan
        outerr['Peak Height'] = 0.0
        output['FWHM'] = np.nan
        outerr['FWHM'] = 0.0
        output['Peak Centre'] = np.nan
        outerr['Peak Centre'] = 0.0
        output['Background'] = fitfunc(xold[len(xold)//2],*fitvals)
        outerr['Background'] = np.std(y)
        ara = 0.0
        ara = 0.0
    output['Area'] = ara
    outerr['Area'] = dara
    
    '-----------------------------------------------------'
    '----------------------Extra data---------------------'
    '-----------------------------------------------------'
    # Calculate fit
    if interpolate:
        xfit = np.linspace(min(xold),max(xold),50*len(xold))
    else:
        xfit = xold
    yfit = fitfunc(xfit,*fitvals)
    output['x'] = xfit
    output['y'] = yfit
    
    # Calculate CHI^2
    ycomp = fitfunc(xold,*fitvals)
    chi = np.sum( (y-ycomp)**2 / dy**2)
    dof = len(y) - len(fitvals) # Number of degrees of freedom (Nobs-Npar)
    chinfp = chi/dof
    output['CHI**2'] = chi
    output['CHI2 per dof'] = chinfp
    
    # Results String
    res_str = ' ------{} Fit:----- \n'.format(type)
    for estn in range(len(fitvals)):
        res_str +='{0:12s} = {1:20s}\n'.format(valnames[estn],stfm(fitvals[estn],errvals[estn]))
    res_str +='        Area = {0:20s}\n'.format(stfm(ara,dara))
    output['Results'] = res_str
    
    # Print Results
    if disp:
        res_str += '       CHI^2 = {0:10.8G}\n'.format(chi)
        res_str += 'CHI^2 per free par = {0:10.3G}\n'.format(chinfp) 
        print(res_str)
    
    # Plot Results
    if plot:
        plt.figure(figsize=[12,10], dpi=fig_dpi)
        plt.errorbar(x,y,dy,fmt='b-o',lw=2,label='Data')
        plt.plot(xfit,yfit,'r-',lw=2,label='Fit')
        plt.legend(loc=0, frameon=False, fontsize=18)
        plt.text(0.7,0.6, res_str,transform=plt.gca().transAxes,fontsize=16,fontname='Times New Roman')
        plt.show()
    return output,outerr

def fittest(x,y,dy=None,tryall=False,disp=False):
    "Attempt multiple fit types and return the best"
    
    # Set dy to 1 if not given
    if dy is None: dy=np.ones(len(y))
    
    # Remove zeros from x - causes errors in covariance matrix
    xold = x
    offset = 0.
    if any(np.abs(x)<0.001):
        print( 'Zero detected - adding 0.001 to x values' )
        offset = 0.001
        x = x + offset
    if any(np.isnan(dy)):
        print( 'Ignoring errors due to NaNs' )
        dy=np.ones(len(y))
    
    # Handle zero intensities
    y[ y<0.01 ] = y[ y<0.01 ]+0.01
    dy[ dy<0.01 ] = dy[ dy<0.01 ]+0.01
    
    # Estimate starting parameters
    amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara = simpfit(x,y)
    
    # Define functions and starting paramters
    fitname = ['Line','Gaussian','Lorentzian','pVoight']
    fitfuncs = [straightline,gauss,lorentz,pvoight]
    estvals = [[0,bkg],
               [amp,cen,wid,bkg],
               [amp,cen,wid,bkg],
               [amp,cen,wid,0.5,bkg]]
    valnames = [['Gradient','Intercept'],
                ['Peak Height','Peak Centre','FWHM','Background'],
                ['Peak Height','Peak Centre','FWHM','Background'],
                ['Peak Height','Peak Centre','FWHM','Lorz frac','Background']]
    trialvals=[[np.arange(-2,2,0.5),np.linspace(min(y),max(y),5)],
               [np.linspace(min(y),max(y),5),np.linspace(min(x),max(x),5),np.linspace(x[1]-x[0],x[-1]-x[0],5),np.linspace(min(y),max(y),5)],
               [np.linspace(min(y),max(y),5),np.linspace(min(x),max(x),5),np.linspace(x[1]-x[0],x[-1]-x[0],5),np.linspace(min(y),max(y),5)],
               [np.linspace(min(y),max(y),5),np.linspace(min(x),max(x),5),np.linspace(x[1]-x[0],x[-1]-x[0],5),np.linspace(0,1,5),np.linspace(min(y),max(y),5)]]
    minvals = [[-np.inf,0],
               [0,min(x),x[1]-x[0],-np.inf],
               [0,min(x),x[1]-x[0],-np.inf],
               [0,min(x),x[1]-x[0],-np.inf,-np.inf]]
    maxvals = [[np.inf,np.inf],
               [np.inf,max(x),5*(x[-1]-x[0]),np.inf],
               [np.inf,max(x),5*(x[-1]-x[0]),np.inf],
               [np.inf,max(x),5*(x[-1]-x[0]),np.inf,np.inf]]
    chival = np.zeros(len(fitfuncs))
    
    # Loop through each function and determine the best fit 
    for n in range(len(fitfuncs)):
        if tryall:
            # Attemp every compbination of trial values VERY SLOW!!! 
            trials = product(*trialvals[n]) # from itertools
            chival[n] = np.inf
            fitvals = np.inf*np.ones(len(valnames[n]))
            errvals = np.inf*np.ones(len(valnames[n]))
            for tt in trials:
                # attempt the fit
                try:
                    tfitvals, tcovmat = curve_fit(fitfuncs[n],x,y,tt,sigma=dy)
                except RuntimeError:
                    continue
                
                # Check fit is within the allowed range
                for v in range(len(valnames[n])):
                    if tfitvals[v] < minvals[n][v] or tfitvals[v] > maxvals[n][v]:
                        continue
                
                yfit = fitfuncs[n](xold,*tfitvals) # New curve
                newchi = np.sum( (y-yfit)**2 / dy) # Calculate CHI^2
                if newchi < chival[n]:
                    fitvals = tfitvals
                    errvals = np.sqrt(np.diag(tcovmat))
                    chival[n] = 1*newchi
        else:
            # attempt the fit
            try:
                fitvals, covmat = curve_fit(fitfuncs[n],x,y,estvals[n],sigma=dy)
            except RuntimeError:
                fitvals = 1*estvals[n]
                covmat = np.nan*np.eye(len(estvals[n]))
            
            errvals = np.sqrt(np.diag(covmat))
            
            yfit = fitfuncs[n](xold,*fitvals) # New curve
            chival[n] = np.sum( (y-yfit)**2 / dy) # Calculate CHI^2
            
            # Check fit is within the allowed range
            for v in range(len(valnames[n])):
                if fitvals[v] < minvals[n][v] or fitvals[v] > maxvals[n][v]:
                    valnames[n][v] += '*'
                    chival[n] = np.inf
        
        if disp:
            print( '----{}: {}----'.format(n,fitname[n]) )
            for v in range(len(valnames[n])):
                print( '{0:10s} = {1:10.3G} +/- {2:10.3G}'.format(valnames[n][v],fitvals[v],errvals[v]) )
            #print( '      Area = {0:10.3G} +/- {1:10.3G}'.format(ara,dara) )
            print( '     CHI^2 = {0:10.8G}'.format(chival[n]) )
    
    # Find the minimum chi
    minval = np.argmin(chival)
    return fitname[minval]

def gauss2D(XY,height=1,cen_x=0,cen_y=0,FWHM_x=.5,FWHM_y=.5,bkg=0):
    "Define 2D Gaussian"
    X,Y = XY
    G= height*np.exp(-np.log(2)*( ((X-cen_x)/(FWHM_x/2.0))**2 + ((Y-cen_y)/(FWHM_y/2.0))**2 )) + bkg
    return G.ravel()

def gaussfit2D(cut):
    "Fit a peak in a pilatus image ***Unfinished***"
    
    X,Y = np.meshgrid(np.arange(487.0),np.arange(195.0))
    iniguess=[max(cut.ravel()),487//2,195//2,20,20,median(cut)]
    popt,pcov = opt.curve_fit(gauss2D,(X,Y),cut.ravel(),p0=iniguess)
    G = gauss2D((X,Y),*popt).reshape(cut.shape)
    return G

def orderpar(x,Tc=100,beta=0.5,amp=1):
    "Generate an order paramter"
    #op = amp*np.real(np.power(np.complex(Tc-x),beta))
    op = amp*np.power(Tc-x,beta)
    op[np.isnan(op)] = 0.0
    return op

def orderparfit(x,y,dy=None,Tc=None,disp=None):
    "Fit an order parameter to a temperature dependence y = f(T)"
    
    # Set dy to 1 if not given
    if dy is None: dy=np.ones(len(y))
    
    # Remove zeros from x - causes errors in covariance matrix
    xold = x
    offset = 0.
    if any(np.abs(x)<0.001):
        print( 'Zero detected - adding 0.001 to x values' )
        offset = 0.001
        x = x + offset
    if any(np.isnan(dy)):
        print( 'Ignoring errors due to NaNs' )
        dy=np.ones(len(y))
    
    # Handle zero intensities
    y[ y<0.01 ] = 0.01
    dy[ dy<0.01 ] = 0.01
    
    # Starting parameters
    if Tc is None:
        Tc = x[len(x)//2]
    beta = 0.5
    amp = np.mean(y[:len(y)//10])
    print( Tc,beta,amp )
    
    try:
        vals, covmat = curve_fit(orderpar,x,y,[Tc,beta,amp],sigma=dy)
    except RuntimeError:
        vals = [0,beta,amp]
        covmat = np.diag([np.nan,np.nan,np.nan])
    # Values
    Tc = vals[0]-offset
    beta = vals[1]
    amp = vals[2]
    # Errors
    perr = np.sqrt(np.diag(covmat))
    dTc = perr[0]
    dbeta = perr[1]
    damp = perr[2]
    
    # Calculate fit
    yfit = orderpar(xold,Tc,beta,amp)
    
    # Calculate CHI^2
    chi = np.sum( (y-yfit)**2 / dy)
    dof = len(y) - 4 # Number of degrees of freedom (Nobs-Npar)
    chinfp = chi/dof
    
    # Check fit has worked
    if Tc <= 0 or any(np.isnan([dTc,dbeta,damp])):
        print( 'Fit didn''t work: oh dear' )
        return
    
    # Print Results
    if disp:
        print( ' ------Order Parameter Fit:----- ' )
        print( '        Tc = {0:10.3G} +/- {1:10.3G}'.format(Tc,dTc) )
        print( '      Beta = {0:10.3G} +/- {1:10.3G}'.format(beta,dbeta) )
        print( '       Amp = {0:10.3G} +/- {1:10.3G}'.format(amp,damp) )
        print( '     CHI^2 = {0:10.3G}'.format(chi) )
        print( '  CHI^2 per free par = {0:10.3G}'.format(chinfp) )
    return Tc,beta,amp,dTc,dbeta,damp,yfit


"----------------------------Smoothing Functions--------------------------"


def rebin(vol,step=[4,4,1]):
    """
    V = rebin(volume,[4,4,1])
    Rebins the Pilatus array, taking the aveage of multiple pixels
    """
    
    sh = vol.shape
    
    # remove trailing elements
    vol2=vol[:sh[0]-sh[0]%step[0],:sh[1]-sh[1]%step[1],:sh[2]-sh[2]%step[2]]
    
    # sum pixels
    ct = 0.0
    volsum = np.zeros(list(np.array(vol2.shape)//np.array(step)))
    for n in range(step[0]):
        for m in range(step[1]):
            for o in range(step[2]):
                ct+=1.0
                volsum += vol2[n::step[0],m::step[1],o::step[2]]
    
    # Take the average
    return volsum / ct

def conv_gauss(y, window_size=31, widval=0.2):
    " smooth data using covolution with gaussian"
    
    half_window = (window_size -1) // 2
    
    x = np.arange(-half_window,half_window+1)
    sig = widval*window_size
    amp = 1.0/(sig*np.sqrt(2*np.pi)) # normalised to area of 1
    cen = 0
    bkg = 0
    G = amp*np.exp(-(x-cen)**2/(2*sig**2)) + bkg
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    #return convolve(G,y,mode= 'valid')
    return convolve(y,G,mode= 'valid')


"--------------------------Peak Finding Functions-------------------------"


def pilpeak(vol,disp=False):
    """
    Find peak in pilatus detector, within peakregion
     [pk_i,pk_j],frame = pilpeak(vol)
     
    Currently uses a very basic search that rebins the volume into courser bins,
    averaging the intensities then finding the position of the max bin.
    """
    
    step = [2,2,2] # rebining step size
    bvol = rebin(vol,step) # rebin the volume and average the values (removing some noise)
    pki=peakregion[0:4:2]
    pkj=peakregion[1:4:2]
    volpeak = bvol[pki[0]//step[0]:pki[1]//step[1],pkj[0]//step[0]:pkj[1]//step[1],:]
    pos = np.unravel_index(volpeak.argmax(),volpeak.shape) # max pixel within the region
    pos = np.multiply(pos,step) # recover unbinned index
    ROIcen = [pos[0]+peakregion[0],pos[1]+peakregion[1]]
    if disp:
        print( 'Peak found at [{0},{1}], frame #{2}, in region [{3},{4},{5},{6}]'.format(ROIcen[0],ROIcen[1],pos[2],peakregion[0],peakregion[1],peakregion[2],peakregion[3]) )
    return ROIcen,pos[2]

def pilmaxval(vol,step=[4,4,1]):
    """
    Rebin the area detector images and return the array of maxvals of the rebinned image
      maxval = pilmaxval( vol, step )
      vol = [nxmxs] volume from getvol
      step = [dN,dM,dS] number of steps to average over in each direction
      maxval = [dS] array of maxvals
    """
    
    bvol = rebin(vol,step) # rebin the volume and average the values (removing some noise)
    maxvals = bvol.max(axis=0).max(axis=0)
    return maxvals

def pilpeaks(vol,rebin_step=[4,4,2],peakheight=100,disp=False):
    """
    Find peak in pilatus detector, within peakregion
     [pk_i,pk_j],frame = pilpeak(vol)
     
    Currently uses a very basic search that rebins the volume into courser bins,
    averaging the intensities then finding the position of the max bin.
    """
    
    background = np.median(vol)
    bvol = rebin(vol,rebin_step) # rebin the volume and average the values (removing some noise)
    pki=peakregion[0:4:2]
    pkj=peakregion[1:4:2]
    volpeak = bvol[pki[0]//rebin_step[0]:pki[1]//rebin_step[1],pkj[0]//rebin_step[0]:pkj[1]//rebin_step[1],:]
    
    volpeak_bkg = volpeak-background # background subtraction
    idx = volpeak_bkg > peakheight # find peaks
    pos_volpeak = np.array(np.where(idx)).T # gen array positions of peaks in the volpeak array
    peak_inten=volpeak_bkg[pos_volpeak[:,0],pos_volpeak[:,1],pos_volpeak[:,2]] # get intensity values
    
    # Remove neighboring pixels
    # Take the largest pposition
    peak_sep = np.sum(np.abs(np.diff(pos_volpeak,axis=0)),axis=1) # compare peak positions
    peak_sep = np.insert(peak_sep,0,100) # add first value
    
    
    # Recover original intensities
    pos_vol = np.multiply(pos_volpeak,rebin_step) + np.array(rebin_step)//2 # gen array positions of peaks in the vol array, recoving unbinned index
    ROIcen = [pos_vol[:,0]+peakregion[0],pos_vol[:,1]+peakregion[1]]
    if disp:
        for n in range(pos_vol.shape[0]):
            print( 'Peak found at [{:3.0f},{:3.0f}], frame #{:3.0f}, with intensity {:1.0f}'.format(ROIcen[0][n],ROIcen[1][n],pos_vol[n,2],peak_inten[n]) )
    return ROIcen,pos_vol[:,2],peak_inten

def peakfind(Y,cutoff=0.01):
    "Finds peaks in the given 1D spectrum"
    " Doesn't really work at the moment and needs testing"
    
    # Code from "deinonychusaur" http://stackoverflow.com/questions/24656367/find-peaks-location-in-a-spectrum-numpy
    # Requires scipy.signal.convolve
    
    
    
    #Obtaining derivative
    kernel = [1, 0, -1]
    dY = convolve(Y, kernel, 'valid')
    #dY = np.diff(Y)
    
    #Checking for sign-flipping
    S = np.sign(dY)
    ddS = convolve(S, kernel, 'valid')
    
    #These candidates are basically all negative slope positions
    #Add one since using 'valid' shrinks the arrays
    candidates = np.where(dY > 0)[0] + (len(kernel) - 1)
    
    #Here they are filtered on actually being the final such position in a run of
    #negative slopes
    peaks = sorted(set(candidates).intersection(np.where(ddS == -2)[0] + 1))
    
    #If you need a simple filter on peak size you could use:
    peaks = np.array(peaks)[Y[peaks] > cutoff]
    return peaks

def ispeak(Y,dY=None,test = 1,disp=False,return_rat=False):
    "Determines whether a peak exists in the given dataset"
    
    if dY is None:
        dY = error_func(Y)
    
    "From Blessing, J. Appl. Cryst. (1997). 30, 421-426"
    "EQU: (1) + (6)"
    " Background estimation added by me"
    s = np.mean(Y)
    bkg = np.min(Y)
    wi = 1/dY**2
    signal = np.sum( wi*(Y-bkg) )/np.sum( wi )
    err = np.sqrt( len(Y) / np.sum( wi ))
    
    #s = np.sum(Y)/len(Y)
    #h,bin = np.histogram(Y,10)
    #bkg = bin[np.argmax(h)]
    #signal = np.sum(Y-bkg)/len(Y)
    #srt = np.sort(Y)
    #err = 3*np.mean(np.abs(np.diff(srt[:len(Y)//2])))
    
    #sig=np.average(Y,weights=np.sqrt(np.abs(Y)))
    #err=np.average(np.sqrt(np.abs(Y)),weights=np.sqrt(np.abs(Y)))
    
    #sig = np.mean(Y)
    #err = np.mean(np.sqrt(np.abs(Y)))
    
    rat = signal / err
    # A peak exists if the signal/background ratio is greater than about 15
    if disp:
        print( 'avg: ',s )
        print( 'bkg: ',bkg )
        print( 'signal: ',signal )
        print( 'error: ',err )
        print( 'rat: ',rat )
    if return_rat:
        return rat
    return rat > test


"----------------------------Misc Functions-------------------------------"

def normalise(num):
    """
    Automatically determine normalisation factor for given scan
        Returns norm_factor, norm_txt
        
        In = Io/norm_factor
        dIn = dIo/norm_factor
        vary += norm_txt
    """
    
    "---Load data---"
    try:
        d = readscan(num)
    except TypeError:
        d = num
    
    "---Get metadata---"
    keys = d.scannables
    m = d.metadata
    cmd = m.cmd # Scan command
    
    "---Get count time---"
    # Count time can have several names
    if 't' in keys:
        cnt = 1.0*d.t
    elif 'counttime' in keys:
        cnt = 1.0*d.counttime
    elif 'count_time' in keys:
        cnt = 1.0*d.count_time
    else:
        cnt=1.0
    
    "---Get incident beam normalisation---"
    if normby.lower() in ['ic1','ic1monitor']:
        inorm = d.ic1monitor/exp_monitor
        normtxt = '('+normby+'/'+str(exp_monitor)+')'
    elif normby.lower() in ['rc','ring current']:
        inorm = d.rc/exp_ring_current
        normtxt = '('+normby+'/'+str(exp_ring_current)+')'
    else:
        inorm = np.ones(len(d.rc))
        normtxt = ''
    
    full_norm = d.metadata.Transmission*cnt*inorm
    normtxt = '/Tran./t/'+normtxt
    return full_norm,normtxt

def abscor(eta=0,chi=90,delta=0,mu=0,gamma=0,u=1.0,disp=False,plot=False):
    """
    Calculate absorption correction
     A = abscor(eta,chi,delta,mu,gamma,u, disp, plot)
         eta/chi/delta/mu/gamma = diffractometer angles
         u = absoprtion coefficient
         disp = True/False - display information
         plot = True/False - plot diffractometer/sample orientation
    
    Usage:
        scans = range(584448,584481,1)
        eta = getmeta(scans,'eta')
        delta = getmeta(scans,'delta')
        chi = getmeta(scans,'chi')
        A=np.zeros(len(eta))
        for n in range(len(eta)): 
            A[n] = abscor(eta[n]-4.41, chi[n]+2.43, delta[n],disp=True)
        Icorrected = Iexp/A
    
     From: IT-C Secion 6.3.3 Absorption corrections
     Table 6.3.3.1. (1) Reflection from a crystal slab with negligible transmission
          A = sin(theta - phi) / mu( sin(theta-phi) + sin(theta+phi) )
     
     A = Transmission coefficient. (Absorption correction A*=1/A) 
         "The reduction in the intensity of an x-ray reflection" 
             Iobs = A*Iorig
             Icorrected = Iobs/A
     theta = Bragg angle
     phi = "the crystal planes are inclined at an angle phi to the extended face and the normal in the plane of the incidnend and diffracted beams"
     mu = absorption coefficient
    """
        
    " Convert angles"
    eta = np.deg2rad(eta)
    mu = np.deg2rad(mu)
    chi = np.deg2rad(chi)
    delta = np.deg2rad(delta)
    gamma = np.deg2rad(gamma)
    
    " Determine Wavevector, Q"
    " Beamline coordinates z->Along beam, towards beamstop, y->upwards, x->away from synchrotron"
    tth = delta + gamma # one of these should be zero
    theta = tth/2.0
    ki = np.array([0,0,1]) # incident wavevector
    #kf = np.array([0,np.sin(tth),np.cos(tth)]) # final wavevector, vertical scattering
    #kf = np.array([np.sin(tth),0,np.cos(tth)]) # final wavevector, horizontal scattering
    #kf = np.array([np.sin(gamma),np.sin(delta),np.cos(gamma)*np.cos(delta)])
    kf = rot3D(ki,0,np.rad2deg(gamma),np.rad2deg(-delta))
    
    Q = kf - ki # wavevector transfer
    mQ = np.sqrt( np.sum( Q**2 ) )
    
    " Determine surface normal, N"
    #N = np.array([np.cos(eta)*np.cos(chi) , np.cos(eta)*np.sin(chi) , np.sin(-eta)]) # vertical scattering, sample surface perp to [0,1,0] at eta 0, chi 90
    #N = np.array([np.cos( mu)*np.sin(chi) ,-np.cos(chi) , np.sin(-mu)*np.sin( chi)]) # Horizontal scattering, sample surface perp to [1,0,0] at mu 0, chi 90
    def normal(eta,mu,chi):
        #N = np.array([np.cos(eta)*np.cos(chi) , np.cos(eta)*np.sin(chi) , np.sin(-eta)])
        #N = np.array([np.cos(eta)*np.sin(chi)+np.cos(mu)*np.cos(chi),
        #              np.cos(eta)*np.sin(chi),
        #              np.sin(-eta)*np.sin(chi)+np.sin(-mu)*np.cos(chi)])
        N = np.array([np.sin(mu)*np.sin(eta)*np.sin(chi) + np.cos(mu)*np.cos(chi),
                  np.cos(eta)*np.sin(chi),
                 -np.cos(mu)*np.sin(eta)*np.sin(chi) - np.sin(mu)*np.cos(chi)])
        #N = rot3D(np.array([1,0,0]),np.rad2deg(chi),np.rad2deg(mu),np.rad2deg(-eta))
        return N
    N = normal(eta,mu,chi)
    mN = np.sqrt( np.sum( N**2 ) )
    
    " Calculate angle between Q and N"
    phi = np.arccos( np.dot(Q,N) /(mQ*mN) )
    
    " Calcualte the absorption correction"
    A = np.sin(theta-phi) / (u*( np.sin(theta-phi) + np.sin(theta+phi) )) # IT-C Table 6.3.3.1. A = Transmission coefficient. "The reduction in the intensity of an x-ray reflection" Iobs = A*Iorig, Icorrected = Iobs/A
    #print( np.round(np.rad2deg(eta),2),np.round(np.rad2deg(chi),2),np.round(np.rad2deg(delta),2),np.round(Q,3),np.round(N,3),' phi = ',round(np.rad2deg(phi),1),' A = ',np.round(A,3) )
    
    if disp:
        " Convert angles"
        deta = np.rad2deg(eta)
        dmu = np.rad2deg(mu)
        dchi = np.rad2deg(chi)
        ddelta = np.rad2deg(delta)
        dgamma = np.rad2deg(gamma)
        dphi = np.rad2deg(phi)
        sQ = str(np.round(Q,2))
        sN = str(np.round(N,2))
        print( 'eta: {:5.2f} mu: {:5.2f} chi: {:5.2f} delta: {:5.2f} gamma: {:5.2f}  Q={:16s}  N={:16s}  phi = {:5.2f}  A = {:5.2g}'.format(deta,dmu,dchi,ddelta,dgamma,sQ,sN,dphi,A) )
        
    
    if plot:
        " Create Figure"
        fig = plt.figure(figsize=[12,14], dpi=fig_dpi)
        ax = fig.add_subplot(211,projection='3d')
        " Plot axes lines"
        ax.plot3D([-1,1],[0,0],[0,0],'k-',linewidth=1)
        ax.plot3D([0,0],[-1,1],[0,0],'k-',linewidth=1)
        ax.plot3D([0,0],[0,0],[-1,1],'k-',linewidth=1)
        " swap axes z -> y, y -> z"
        " Plot ki, kf & Q"
        ax.plot3D([0,0],[-ki[2],0],[0,0],'r-',linewidth=2,label='$k_i$')
        ax.plot3D([0,kf[0]],[0,kf[2]],[0,kf[1]],'r-o',linewidth=2,label='$k_f$')
        ax.plot3D([0,Q[0]],[0,Q[2]],[0,Q[1]],'r:o',linewidth=3,label='Q')
        " Plot surface & Normal"
        pN, = ax.plot([0,N[0]],[0,N[2]],[0,N[1]],'b-',linewidth=3,label='N')
        
        " Labels"
        ax.invert_yaxis()
        ax.set_xlabel('x',fontsize=16)
        ax.set_ylabel('z',fontsize=16)
        ax.set_zlabel('y',fontsize=16)
        
        ddelta = np.rad2deg(delta)
        dgamma = np.rad2deg(gamma)
        sQ = str(np.round(Q,2))
        ttl = 'delta: {:5.2f} gamma: {:5.2f}\nQ={:16s}'.format(ddelta,dgamma,sQ)
        ax.set_title(ttl)
        plt.legend(loc=0)
        
        " Create slider on plot"
        #plt.subplot(212)
        #plt.axes('off')
        axsldr1 = plt.axes([0.15, 0.15, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        axsldr2 = plt.axes([0.15, 0.25, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        axsldr3 = plt.axes([0.15, 0.35, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        
        sldr1 = plt.Slider(axsldr1, 'eta', 0, 120,valinit=np.rad2deg(eta), valfmt = '%0.0f')
        sldr2 = plt.Slider(axsldr2, 'mu', 0, 120,valinit=np.rad2deg(mu), valfmt = '%0.0f')
        sldr3 = plt.Slider(axsldr3, 'chi', -10, 100,valinit=np.rad2deg(chi), valfmt = '%0.0f')
        txt = plt.title('N = [{:1.2g},{:1.2g},{:1.2g}]\nphi = {:1.2g}, A = {:1.2g}'.format(N[0],N[1],N[2],np.rad2deg(phi),A),fontsize=18 )
        
        " Slider update function"
        def update(val):
            "Update function for pilatus image"
            eta = np.deg2rad(sldr1.val)
            mu  = np.deg2rad(sldr2.val)
            chi = np.deg2rad(sldr3.val)
            N = normal(eta,mu,chi)
            mN = np.sqrt( np.sum( N**2 ) )
            uphi = np.arccos( np.dot(Q,N) /(mQ*mN) )
            uA = np.sin(theta-uphi) / (u*( np.sin(theta-uphi) + np.sin(theta+uphi) )) # IT-C Table 6.3.3.1
            pN.set_data( [0,N[0]],[0,N[2]])
            pN.set_3d_properties([0,N[1]])
            txt.set_text('N = [{:1.2g},{:1.2g},{:1.2g}]\nphi = {:1.2g}, A = {:1.2g}'.format(N[0],N[1],N[2],np.rad2deg(uphi),uA))
            plt.draw()
            #fig.canvas.draw()
        sldr1.on_changed(update)
        sldr2.on_changed(update)
        sldr3.on_changed(update)
    return abs(A)

def volcor(eta=0,chi=90,delta=0,mu=0,gamma=0,u=1.0,disp=False,plot=False):
    """
    Calculate scattering volume correction, or self-absorption
     A = volcor(eta,chi,delta,mu,gamma,u, disp, plot)
         eta/chi/delta/mu/gamma = diffractometer angles
         u = absoprtion coefficient
         disp = True/False - display information
         plot = True/False - plot diffractometer/sample orientation
    
    Usage:
        scans = range(584448,584481,1)
        eta = getmeta(scans,'eta')
        delta = getmeta(scans,'delta')
        chi = getmeta(scans,'chi')
        A=np.zeros(len(eta))
        for n in range(len(eta)): 
            A[n] = abscor(eta[n]-4.41, chi[n]+2.43, delta[n],disp=True)
        Icorrected = Iexp/A
    
     
     
     A = Transmission coefficient. (Absorption correction A*=1/A) 
         "The reduction in the intensity of an x-ray reflection" 
             Iobs = A*Iorig
             Icorrected = Iobs/A
     theta = Bragg angle
     phi = "the crystal planes are inclined at an angle phi to the extended face and the normal in the plane of the incidnend and diffracted beams"
     mu = absorption coefficient
    """
    
    " Determine Wavevector, Q"
    " Beamline coordinates z->Along beam, towards beamstop, y->upwards, x->away from synchrotron"
    tth = delta + gamma # one of these should be zero
    theta = tth/2.0
    ki = np.array([0,0,1]) # incident wavevector
    kf = rot3D(ki,0,gamma,-delta)
    
    Q = kf - ki # wavevector transfer
    Q = np.sqrt( np.sum( Q**2 ) )
    
    N = surface_normal(eta,mu,chi)
    phi = np.rad2deg(np.arccos(np.dot(Q,N)))
    alpha = 90-np.rad2deg(np.arccos(np.dot(-ki,N)))
    beta = 90-np.rad2deg(np.arccos(np.dot(kf,N)))
    A = (1/u)/(1 + (np.sin(np.deg2rad(alpha))/np.sin(np.deg2rad(beta)))) # Beutier et al PRB 2008 HoMn2O5 eq. 4, Iobs = Ical*A
    
    if disp:
        sQ = str(np.round(Q,2))
        sN = str(np.round(N,2))
        print( 'eta: {:5.2f} mu: {:5.2f} chi: {:5.2f} delta: {:5.2f} gamma: {:5.2f}  Q={:16s}  N={:16s}  alpha = {:5.2f} beta = {:5.2f} phi = {:5.2f}  A = {:5.2g}'.format(eta,mu,chi,delta,gamma,sQ,sN,alpha,beta,phi,A) )
    
    if plot:
        " Create Figure"
        fig = plt.figure(figsize=[12,14], dpi=fig_dpi)
        ax = fig.add_subplot(211,projection='3d')
        " Plot axes lines"
        ax.plot3D([-1,1],[0,0],[0,0],'k-',linewidth=1)
        ax.plot3D([0,0],[-1,1],[0,0],'k-',linewidth=1)
        ax.plot3D([0,0],[0,0],[-1,1],'k-',linewidth=1)
        " swap axes z -> y, y -> z"
        " Plot ki, kf & Q"
        ax.plot3D([0,0],[-ki[2],0],[0,0],'r-',linewidth=2,label='$k_i$')
        ax.plot3D([0,kf[0]],[0,kf[2]],[0,kf[1]],'r-o',linewidth=2,label='$k_f$')
        ax.plot3D([0,Q[0]],[0,Q[2]],[0,Q[1]],'r:o',linewidth=3,label='Q')
        " Plot surface & Normal"
        pN, = ax.plot([0,N[0]],[0,N[2]],[0,N[1]],'b-',linewidth=3,label='N')
        
        " Labels"
        ax.invert_yaxis()
        ax.set_xlabel('x',fontsize=16)
        ax.set_ylabel('z',fontsize=16)
        ax.set_zlabel('y',fontsize=16)
        
        sQ = str(np.round(Q,2))
        ttl = 'delta: {:5.2f} gamma: {:5.2f}\nQ={:16s}'.format(delta,gamma,sQ)
        ax.set_title(ttl)
        plt.legend(loc=0)
        
        " Create slider on plot"
        #plt.subplot(212)
        #plt.axes('off')
        axsldr1 = plt.axes([0.15, 0.15, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        axsldr2 = plt.axes([0.15, 0.25, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        axsldr3 = plt.axes([0.15, 0.35, 0.65, 0.03], axisbg='lightgoldenrodyellow')
        
        sldr1 = plt.Slider(axsldr1, 'eta', 0, 120,valinit=eta, valfmt = '%0.0f')
        sldr2 = plt.Slider(axsldr2, 'mu', 0, 120,valinit=mu, valfmt = '%0.0f')
        sldr3 = plt.Slider(axsldr3, 'chi', -10, 100,valinit=chi, valfmt = '%0.0f')
        txt = plt.title('N = [{:1.2g},{:1.2g},{:1.2g}]\nphi = {:1.2g}, A = {:1.2g}'.format(N[0],N[1],N[2],phi,A),fontsize=18 )
        
        " Slider update function"
        def update(val):
            "Update function for pilatus image"
            eta = sldr1.val
            mu  = sldr2.val
            chi = sldr3.val
            N = normal(eta,mu,chi)
            uphi = ang(Q,N[n,:])
            ualpha = 90-ang(N[n,:],-ki)
            ubeta = ang(-kf,N[n,:])-90
            uA = (1/u)/(1 + (np.sin(np.deg2rad(alpha[n]))/np.sin(np.deg2rad(beta[n])))) # Beutier et al PRB 2008 HoMn2O5 eq. 4, Iobs = Ical*A
            pN.set_data( [0,N[0]],[0,N[2]])
            pN.set_3d_properties([0,N[1]])
            txt.set_text('N = [{:1.2g},{:1.2g},{:1.2g}]\nphi = {:1.2g}, A = {:1.2g}'.format(N[0],N[1],N[2],uphi,uA))
            plt.draw()
            #fig.canvas.draw()
        sldr1.on_changed(update)
        sldr2.on_changed(update)
        sldr3.on_changed(update)
    return abs(A)

def surface_normal(eta=0, mu=0, chi=90):
    """
    Returns a vector parallel to the phi axis of the diffractometer, normal to the sample surfaace
    """
    reta = np.deg2rad(eta)
    rmu = np.deg2rad(mu)
    rchi = np.deg2rad(chi)
    N = np.array([np.sin(rmu)*np.sin(reta)*np.sin(rchi) + np.cos(rmu)*np.cos(rchi),
              np.cos(reta)*np.sin(rchi),
             -np.cos(rmu)*np.sin(reta)*np.sin(rchi) - np.sin(rmu)*np.cos(rchi)])
    return N

def surface_angles(eta=0, mu=0, chi=90, delta=0, gamma=0):
    """
    Returns alpha and beta, the angles of the incident and scattered beam to the sample surface
    """
    
    " Determine Wavevector, Q"
    " Beamline coordinates z->Along beam, towards beamstop, y->upwards, x->away from synchrotron"
    tth = delta + gamma # one of these should be zero
    theta = tth/2.0
    ki = np.array([0,0,1]) # incident wavevector
    kf = rot3D(ki,0,gamma,-delta)
    
    Q = kf - ki # wavevector transfer
    
    N = surface_normal(eta,mu,chi)
    alpha = 90-np.rad2deg(np.arccos(np.dot(-ki,N)))
    beta = 90-np.rad2deg(np.arccos(np.dot(kf,N)))
    return alpha, beta

def rot3D(vec,alpha=0.,beta=0.,gamma=0.):
    """Rotate 3D vector
        vec = rot3D(vec,alpha=0.,beta=0.,gamma=0.)
       where alpha = angle from X axis to Y axis (Yaw)
             beta  = angle from Z axis to X axis (Pitch)
             gamma = angle from Y axis to Z axis (Roll)
       angles in degrees
       In a right handed coordinate system.
           Z
          /|\
           |
           |________\Y
           \        /
            \
            _\/X 
    """
    
    # Convert to radians
    alpha = alpha*np.pi/180.
    beta  = beta*np.pi/180.
    gamma = gamma*np.pi/180.
    
    # Define 3D rotation matrix
    Rx = np.array([[1,0,0],[0,np.cos(gamma),-np.sin(gamma)],[0.,np.sin(gamma),np.cos(gamma)]])
    Ry = np.array([[np.cos(beta),0.,np.sin(beta)],[0.,1.,0.],[-np.sin(beta), 0., np.cos(beta)]])
    Rz = np.array([[np.cos(alpha),-np.sin(alpha),0.],[np.sin(alpha),np.cos(alpha),0.],[0.,0.,1.]])
    R = np.dot(np.dot(Rx,Ry),Rz)
    
    # Rotate coordinates
    return np.dot(R,vec.T).T

def newplot(*args, **kwargs):
    """
    Shortcut to creating a simple plot
    E.G.
      x = np.arange(-5,5,0.1)
      y = x**2
      newplot(x,y,'r-',lw=2,label='Line')
    """

    if 'linewidth' and 'lw' not in kwargs.keys():
        kwargs['linewidth'] = 2

    plt.figure(figsize=[12, 12], dpi=fig_dpi)
    plt.plot(*args, **kwargs)

    plt.setp(plt.gca().spines.values(), linewidth=2)
    plt.xticks(fontsize=25, fontname='Times New Roman')
    plt.yticks(fontsize=25, fontname='Times New Roman')
    plt.ticklabel_format(useOffset=False)
    plt.ticklabel_format(style='sci', scilimits=(-3, 3))


def multiplot(xvals, yvals=None, datarange=None, cmap='jet', labels=None, marker=None):
    """
    Shortcut to creating a simple multiplot with either colorbar or legend
    E.G.
      x = np.arange(-5,5,0.1)
      ys = [x**2, 1+x**2, 2+x**2, 3+x**2, 4+x**2]
      datarange = [0,1,2,3,4]
      multiplot(x, ys, datarange, cmap='winter')
    OR:
      x = np.arange(-5,5,0.1)
      ys = [x**2, 1+x**2, 2+x**2, 3+x**2, 4+x**2]
      labels = ['x*x','2+x*x','3+x*x','4+x*x']
      multiplot(x, ys, labels=labels)
    """

    if yvals is None:
        yvals = xvals
        xvals = []
    yvals = np.asarray(yvals)
    xvals = np.asarray(xvals)

    if datarange is None:
        datarange = range(len(yvals))
    datarange = np.asarray(datarange,dtype=np.float)

    cm = plt.get_cmap(cmap)
    colrange = (datarange - datarange.min()) / (datarange.max() - datarange.min())
    
    if marker is None:
        marker = ''
    linearg = '-' + marker

    plt.figure(figsize=[12, 12], dpi=fig_dpi)
    for n in range(len(datarange)):
        col = cm(colrange[n])
        if len(xvals) == 0:
            plt.plot(yvals[n], linearg, lw=2, color=col)
        elif len(xvals.shape) == 1:
            plt.plot(xvals, yvals[n], linearg, lw=2, color=col)
        else:
            plt.plot(xvals[n], yvals[n], linearg, lw=2, color=col)

    plt.setp(plt.gca().spines.values(), linewidth=2)
    plt.xticks(fontsize=25, fontname='Times New Roman')
    plt.yticks(fontsize=25, fontname='Times New Roman')
    plt.ticklabel_format(useOffset=False)
    plt.ticklabel_format(style='sci', scilimits=(-3, 3))

    if labels is None:
        # Add Colorbar
        sm = plt.cm.ScalarMappable(cmap=cm)
        sm.set_array(datarange)
        cbar = plt.colorbar(sm)
        #cbar.set_label('variation [unit]', fontsize=24, fontweight='bold', fontname='Times New Roman')
    else:
        # Add legend
        plt.legend(labels, loc=0, frameon=False, prop={'size':20,'family':'serif'})


def newplot3(*args, **kwargs):
    """
    Shortcut to creating a simple 3D plot
    Automatically tiles 1 dimensional x and y arrays to match 2D z array,
    assuming z.shape = (len(x),len(y))

    E.G.
      newplot3([1,2,3,4],[9,8,7],[[2,4,6],[8,10,12],[14,16,18],[20,22,24]],'-o')
    """

    if 'linewidth' and 'lw' not in kwargs.keys():
        kwargs['linewidth'] = 2

    fig = plt.figure(figsize=[12, 12], dpi=fig_dpi)
    ax = fig.add_subplot(111, projection='3d')

    x = np.asarray(args[0], dtype=np.float)
    y = np.asarray(args[1], dtype=np.float)
    z = np.asarray(args[2], dtype=np.float)

    if z.ndim == 2:
        if x.ndim < 2:
            x = np.tile(x, z.shape[1]).reshape(z.T.shape).T
        if y.ndim < 2:
            y = np.tile(y, z.shape[0]).reshape(z.shape)

        # Plot each array independently
        for n in range(len(z)):
            ax.plot(x[n], y[n], z[n], *args[3:], **kwargs)
    else:
        ax.plot(*args, **kwargs)


def sliderplot(YY, X=None, slidervals=None, *args, **kwargs):
    """
    Shortcut to creating a simple 2D plot with a slider to go through a third dimension
    YY = [nxm]: y axis data (initially plots Y[0,:])
     X = [n] or [nxm]:  x axis data (can be 1D or 2D, either same length or shape as Y)
     slidervals = None or [m]: Values to give in the slider

    E.G.
      sliderplot([1,2,3],[[2,4,6],[8,10,12],[14,16,18],[20,22,24]],slidervals=[3,6,9,12])
    """

    if 'linewidth' and 'lw' not in kwargs.keys():
        kwargs['linewidth'] = 2

    fig = plt.figure(figsize=[12, 12], dpi=fig_dpi)

    X = np.asarray(X, dtype=np.float)
    Y = np.asarray(YY, dtype=np.float)
    if slidervals is None:
        slidervals = range(Y.shape[0])
    slidervals = np.asarray(slidervals, dtype=np.float)

    if X.ndim < 2:
        X = np.tile(X, Y.shape[0]).reshape(Y.shape)

    plotline, = plt.plot(X[0, :], Y[0, :], *args, **kwargs)
    plt.axis([X.min(), X.max(), Y.min(), Y.max()])
    plt.subplots_adjust(bottom=0.2)
    ax = plt.gca()

    " Create slider on plot"
    axsldr = plt.axes([0.15, 0.05, 0.65, 0.03], axisbg='lightgoldenrodyellow')

    sldr = plt.Slider(axsldr, '', 0, len(slidervals) - 1)
    txt = axsldr.set_xlabel('{} [{}]'.format(slidervals[0], 0), fontsize=18)

    plt.sca(ax)

    " Slider update function"

    def update(val):
        "Update function for pilatus image"
        pno = int(np.floor(sldr.val))
        plotline.set_xdata(X[pno, :])
        plotline.set_ydata(Y[pno, :])
        txt.set_text('{} [{}]'.format(slidervals[pno], pno))
        plt.draw()
        plt.gcf().canvas.draw()
        # fig1.canvas.draw()

    sldr.on_changed(update)


def sliderplot2D(ZZZ, XX=None, YY=None, slidervals=None, *args, **kwargs):
    """
    Shortcut to creating an image plot with a slider to go through a third dimension
    ZZZ = [nxmxo]: z axis data
     XX = [nxm] or [n]:  x axis data
     YY = [nxm] or [m]: y axis data
     slidervals = None or [o]: Values to give in the slider

    if XX and/or YY have a single dimension, the 2D values are generated via meshgrid

    E.G.
      sliderplot([1,2,3],[[2,4,6],[8,10,12],[14,16,18],[20,22,24]],slidervals=[3,6,9,12])
    """

    if 'linewidth' and 'lw' not in kwargs.keys():
        kwargs['linewidth'] = 2

    fig = plt.figure(figsize=[12, 12], dpi=fig_dpi)

    ZZZ = np.asarray(ZZZ, dtype=np.float)

    if slidervals is None:
        slidervals = range(ZZZ.shape[2])
    slidervals = np.asarray(slidervals, dtype=np.float)

    if XX is None:
        XX = range(ZZZ.shape[1])
    if YY is None:
        YY = range(ZZZ.shape[0])
    XX = np.asarray(XX, dtype=np.float)
    YY = np.asarray(YY, dtype=np.float)
    if XX.ndim < 2:
        XX, YY = np.meshgrid(XX, YY)

    p = plt.pcolormesh(XX, YY, ZZZ[:, :, 0])
    # p.set_clim(cax)

    plt.subplots_adjust(bottom=0.2)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.autoscale(tight=True)

    " Create slider on plot"
    axsldr = plt.axes([0.15, 0.05, 0.65, 0.03], axisbg='lightgoldenrodyellow')

    sldr = plt.Slider(axsldr, '', 0, len(slidervals) - 1)
    txt = axsldr.set_xlabel('{} [{}]'.format(slidervals[0], 0), fontsize=18)

    plt.sca(ax)

    " Slider update function"

    def update(val):
        "Update function for pilatus image"
        pno = int(np.round(sldr.val))
        p.set_array(ZZZ[:-1, :-1, pno].ravel())
        txt.set_text('{} [{}]'.format(slidervals[pno], pno))
        plt.draw()
        plt.gcf().canvas.draw()
        # fig1.canvas.draw()
    sldr.on_changed(update)

def labels(ttl=None, xvar=None, yvar=None, zvar=None, legend=False, size='Normal', font='Times New Roman'):
    """
    Add formatted labels to current plot, also increases the tick size
    :param ttl: title
    :param xvar: x label
    :param yvar: y label
    :param zvar: z label (3D plots only)
    :param legend: False/ True, adds default legend to plot 
    :param size: 'Normal' or 'Big'
    :param font: str font name, 'Times New Roman'
    :return: None
    """

    if size.lower() in ['big', 'large', 'xxl', 'xl']:
        tik = 30
        tit = 32
        lab = 35
        leg = 25
    else:
        # Normal
        tik = 18
        tit = 20
        lab = 22
        leg = 18

    plt.xticks(fontsize=tik, fontname=font)
    plt.yticks(fontsize=tik, fontname=font)
    plt.setp(plt.gca().spines.values(), linewidth=2)
    if plt.gca().get_yaxis().get_scale() != 'log':
        plt.ticklabel_format(useOffset=False)
        plt.ticklabel_format(style='sci',scilimits=(-3,3))

    if ttl is not None:
        plt.gca().set_title(ttl, fontsize=tit, fontweight='bold', fontname=font)

    if xvar is not None:
        plt.gca().set_xlabel(xvar, fontsize=lab, fontname=font)

    if yvar is not None:
        plt.gca().set_ylabel(yvar, fontsize=lab, fontname=font)

    if zvar is not None:
        # Don't think this works, use ax.set_zaxis
        plt.gca().set_zlabel(zvar, fontsize=lab, fontname=font)

    if legend:
        plt.legend(loc=0, frameon=False, prop={'size':leg,'family':'serif'})

def tth2d(tth,energy_kev):
    "Converts two-theta array in degrees to d-spacing in A"
    
    e = 1.6021733E-19;  # C  electron charge
    h = 6.62606868E-34; # Js  Plank consant
    c = 299792458;      # m/s   Speed of light
    A = 1e-10;          # m Angstrom      
    pi = np.pi          # mmmm tasty Pi
        
    th = tth*pi/360 # theta in radians
    # Calculate |Q|
    Qmag = np.sin(th)*(1000*energy_kev)*e*4*pi/ (h*c*1e10)
    dspace = 2*pi/Qmag
    return dspace

def tth2q(tth,energy_kev):
    "Converts two-theta array in degrees to wavevector transfer, Q in A-1"
    
    e = 1.6021733E-19;  # C  electron charge
    h = 6.62606868E-34; # Js  Plank consant
    c = 299792458;      # m/s   Speed of light
    A = 1e-10;          # m Angstrom      
    pi = np.pi          # mmmm tasty Pi
        
    th = tth*pi/360 # theta in radians
    # Calculate |Q|
    Qmag = np.sin(th)*(1000*energy_kev)*e*4*pi/ (h*c*1e10)
    return Qmag

def d2tth(dspace,energy_kev):
    "Converts array of d-spacings in A to two-theta in degrees"
    
    e = 1.6021733E-19;  # C  electron charge
    h = 6.62606868E-34; # Js  Plank consant
    c = 299792458;      # m/s   Speed of light
    A = 1e-10;          # m Angstrom      
    pi = np.pi          # mmmm tasty Pi
    
    Qmag = 2*pi/dspace
    th = np.arcsin( Qmag*h*c*1e10 / (1000*energy_kev*e*4*pi) )
    tth = 360*th/pi # 2-theta in radians
    return tth

def q2tth(qmag,energy_kev):
    "Converts array of wavevector transfers in A-1 to two-theta in degrees"
    
    e = 1.6021733E-19;  # C  electron charge
    h = 6.62606868E-34; # Js  Plank consant
    c = 299792458;      # m/s   Speed of light
    A = 1e-10;          # m Angstrom      
    pi = np.pi          # mmmm tasty Pi
    
    th = np.arcsin( qmag*h*c*1e10 / (1000*energy_kev*e*4*pi) )
    tth = 360*th/pi # 2-theta in radians
    return tth

def c2th(HKL=[0,0,1],scan=0):
    "Calculate the 2-theta value of a reflection at the given scan's energy"
    
    d = readscan(scan)
    m = d.metadata
    a, b, c = m.a, m.b, m.c
    alpha, beta, gam = m.alpha1, m.alpha2, m.alpha3
    h,k,l = HKL
    sina2 = np.sin(np.deg2rad(alpha))**2
    sinb2 = np.sin(np.deg2rad(beta))**2
    sing2 = np.sin(np.deg2rad(gam))**2
    cosa = np.cos(np.deg2rad(alpha))
    cosb = np.cos(np.deg2rad(beta))
    cosg = np.cos(np.deg2rad(gam))
    sin, cos = np.sin, np.cos
    
    x1 = h**2/(a**2*sina2)
    x2 = 2*k*l*(cosb*cosg-cosa)/(b*c)
    x3 = k**2/(b**2*sinb2)
    x4 = 2*h*l*(cosa*cosg-cosb)/(a*c)
    x5 = l**2/(c**2*sing2)
    x6 = 2*h*k*(cosa*cosb-cosg)/(a*b)
    x7 = 1 - cosa**2 - cosb**2 - cosg**2 + 2*cosa*cosb*cosg
    # 1/d*2
    dd = (x1+x2+x3+x4+x5+x6)/x7
    # d
    d = 1/np.sqrt(dd)
    return d2tth(d,m.Energy)

def bindata(X,Y,binsep=0.01,bin_cen=None,bin_func=np.nanmean):
    """
    Bins the data into bins separated by binsep. Values within the bin are averaged 
      CEN,INT = bindata(X,Y,binsep=0.01)
    or:
      CEN,INT = bindata(X,Y,bin_cen=CEN)
    
    Note that np.nanmean is used to average values within each bin, which can be used
    to only average values with signal in. Other functions can be used by specifying bin_func:
      CEN,INT = bindata(X,Y,0.01,bin_func=np.max)
    
    E.G.
        tth,ival = pixel2tth(num)
        ival[ival<1] = np.nan
        CEN,INT = bindata(tth,ival,0.01)
    """
    
    t0 = time.clock()
    
    # Generate bin positions
    if bin_cen is None:
        mxdata = np.max(X)
        mndata = np.min(X)
        
        bin_edge = np.arange(mndata,mxdata+binsep,binsep)
        bin_cen = bin_edge - binsep/2.0
        bin_cen = bin_cen[1:] # remove unused bin
    else:
        bin_cen = np.asarray(bin_cen)
        bin_sep = bin_cen[1] - bin_cen[0]
        bin_edge = bin_cen + binsep/2.0 
    
    # Histogram the th values
    bin_pos = np.digitize(X,bin_edge) - 1 # digitize indexes values from bin_edege[i-1] to bin_edge[i], hence the - 1
    t1 = time.clock()
    print('Digitise took {} s'.format(t1-t0))
    
    # Loop over binned angles and average the intensities
    bin_int = np.zeros(bin_cen.shape)
    for n in range(len(bin_cen)):
        inbin = bin_pos == n
        bin_int[n] = bin_func(Y[inbin])
    
    t2 = time.clock()
    print('Bin average took {} s'.format(t2-t1))
    
    return bin_cen,bin_int

def bin_pixel_hkl_cut(num,hkl_centre=None,image_centre=None,sum_tolarance=[0.05,0.05,0.05],max_points=301):
    """
    Generate h,k,l cuts through a pilatus image scan
        H,K,L = bin_pixel_hkl_cut(num, [h,k,l], None, [0.01,0.01,0.01],[101,101,101])
    Inputs:
      num = scan number
      image_centre = [i,j,k] or None* - centre of the cuts on the detector e.g. [104,205,31]
      hkl_centre = [h,k,l] or None* - centre of the cuts in hkl, if None the peak position will be used
      sum_tolarance = [dh,dk,dl] - distance around hkl position to sum, in reciprocal lattice units
      max_points = n - maximum number of points in each cut
    Returns:
      H = [h,Ih] - list of values and binned intensities along h
      K = [k,Ik] - list of values and binned intensities along k
      L = [l,Il] - list of values and binned intensities along l
    
    If image_centre is None, the centre of the detector and scan will be used.
    
    At a choosen central (h,k,l), 3 cuts are generated:
            Cut 1: (H, k-dk:k+dk, l-dl:l+dl)
            Cut 2: (h-dh:h+dh, K, l-dl:l+dl)
            Cut 3: (h-dh:h+dh, k-dk:k+dk, L)
    At each point H/K/L along each cut, pixels matching the following are summed:
            Pixel sum at H:
                abs(h_pixels-H) < H_step/2
                abs(k_pixels-k) < dk
                abs(l_pixels-l) < dl
            Where H_step is defined by the number of cut_points
    cut_points defines the number of points along axis H/K/L, generated from the minimum to maximum h,k,l
    
    A single figure with 3 subplots is created.
    """
    
     # Load hkl matrices
    hhh,kkk,lll = pixel2hkl(num)
    print('Max hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.max(),kkk.max(),lll.max()))
    print('Min hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.min(),kkk.min(),lll.min()))
    
    # Load pilatus volume
    vol = getvol(num)
    
    # Load data
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    # h,k,l min tol
    # diff only works along rows, doesn't calcualte difference between rows
    #min_tol_h = np.max(np.abs(np.diff(hhh)))
    #min_tol_k = np.max(np.abs(np.diff(kkk)))
    #min_tol_l = np.max(np.abs(np.diff(lll)))
    min_tol_h = np.max([np.max(np.abs(np.diff(hhh,axis=0))),np.max(np.abs(np.diff(hhh,axis=1))),np.max(np.abs(np.diff(hhh,axis=2)))])
    min_tol_k = np.max([np.max(np.abs(np.diff(kkk,axis=0))),np.max(np.abs(np.diff(kkk,axis=1))),np.max(np.abs(np.diff(kkk,axis=2)))])
    min_tol_l = np.max([np.max(np.abs(np.diff(lll,axis=0))),np.max(np.abs(np.diff(lll,axis=1))),np.max(np.abs(np.diff(lll,axis=2)))])
    print('Maximum hkl step per pixel = [%1.3g, %1.3g, %1.3g]' % (min_tol_h,min_tol_k,min_tol_l))
    
    # Centre of cuts
    if hkl_centre is None:
        if image_centre is None:
            i,j = pil_centre
            k = len(x)//2
        elif image_centre == 'peak':
            [i,j],k = pilpeak(vol, disp=True)
        else:
            i,j,k = image_centre
        h_centre = hhh[i,j,k]
        k_centre = kkk[i,j,k]
        l_centre = lll[i,j,k]
        print('Peak Centre: (%1.3g,%1.3g,%1.3g)'%(h_centre,k_centre,l_centre))
    else:
        h_centre,k_centre,l_centre = hkl_centre
    
    # Pixels close to centre
    h_cen_idx = np.abs(hhh-h_centre) < sum_tolarance[0]
    k_cen_idx = np.abs(kkk-k_centre) < sum_tolarance[1]
    l_cen_idx = np.abs(lll-l_centre) < sum_tolarance[2]
    
    kl_cen_idx = k_cen_idx*l_cen_idx
    hl_cen_idx = h_cen_idx*l_cen_idx
    hk_cen_idx = h_cen_idx*k_cen_idx
    
    # cut ranges
    h_len = (hhh[kl_cen_idx].max()-hhh[kl_cen_idx].min())/min_tol_h
    print('n h points = %1.0f'%h_len)
    if h_len > max_points:
        print('Reducing h range to %i points' %max_points)
        h_list = np.arange(h_centre-min_tol_h*(max_points//2),h_centre+min_tol_h*(max_points//2)+min_tol_h,min_tol_h)
    else: 
        h_list = np.arange(hhh[kl_cen_idx].min(),hhh[kl_cen_idx].max(),min_tol_h)
    k_len = (kkk[hl_cen_idx].max()-kkk[hl_cen_idx].min())/min_tol_k 
    print('n k points = %1.0f'%k_len)
    if k_len > max_points:
        print('Reducing k range to %i points' %max_points)
        k_list = np.arange(k_centre-min_tol_k*(max_points//2),k_centre+min_tol_k*(max_points//2)+min_tol_k,min_tol_k)
    else:
        k_list = np.arange(kkk[hl_cen_idx].min(),kkk[hl_cen_idx].max(),min_tol_k) 
    l_len = (lll[hk_cen_idx].max()-lll[hk_cen_idx].min())/min_tol_l 
    print('n l points = %1.0f'%l_len)
    if l_len > max_points:
        print('Reducing l range to %i points' %max_points)
        l_list = np.arange(l_centre-min_tol_l*(max_points//2),l_centre+min_tol_l*(max_points//2)+min_tol_l,min_tol_l)
    else: 
        l_list = np.arange(lll[hk_cen_idx].min(),lll[hk_cen_idx].max(),min_tol_l)
    h_step = h_list[1]-h_list[0]
    k_step = k_list[1]-k_list[0]
    l_step = l_list[1]-l_list[0]
    
    h_scan = np.zeros(len(h_list))
    k_scan = np.zeros(len(k_list))
    l_scan = np.zeros(len(l_list))
    print('Binning h axis at (H,%1.3g,%1.3g) in %1.0f steps, summing <%1.0f pixels per step' % (k_centre,l_centre,len(h_list),np.sum(kl_cen_idx)))
    for n in range(len(h_list)):
        hval = h_list[n]
        hidx = np.abs(hhh-hval) < h_step/2
        small_vol = vol[ hidx*kl_cen_idx ]
        h_scan[n] = np.sum(small_vol)
    print('Binning k axis at (%1.3g,K,%1.3g) in %1.0f steps, summing <%1.0f pixels per step' % (h_centre,l_centre,len(k_list),np.sum(hl_cen_idx)))
    for n in range(len(k_list)):
        kval = k_list[n]
        kidx = np.abs(kkk-kval) < k_step/2
        small_vol = vol[ kidx*hl_cen_idx ]
        k_scan[n] = np.sum(small_vol)
    print('Binning l axis at (%1.3g,%1.3g,L) in %1.0f steps, summing <%1.0f pixels per step' % (h_centre,k_centre,len(l_list),np.sum(hk_cen_idx)))
    for n in range(len(l_list)):
        lval = l_list[n]
        lidx = np.abs(lll-lval) < l_step/2
        small_vol = vol[ lidx*hk_cen_idx ]
        l_scan[n] = np.sum(small_vol)
    
    return [h_list,h_scan],[k_list,k_scan],[l_list,l_scan]

def bin_pixel_h_cut(num,hkl_centre=None,image_centre=None,sum_tolarance=0.05,max_points=301):
    """
    Generate h cut through a pilatus image scan
        h,Ih = bin_pixel_hkl_cut(num, [h,k,l], None, 0.01,1001)
    Inputs:
      num = scan number
      image_centre = [i,j,k] or None* - centre of the cuts on the detector e.g. [104,205,31]
      hkl_centre = [h,k,l] or None* - centre of the cuts in hkl, if None the peak position will be used
      sum_tolarance = dk - distance around hkl position to sum, in reciprocal lattice units
      max_points = nk - maximum number of points in cut, to reduce calculation time
    Returns:
      h - list of values
      Ih - binned intensities along l
    
    If image_centre is None, the centre of the detector and scan will be used.
    """
    
     # Load hkl matrices
    hhh,kkk,lll = pixel2hkl(num)
    print('Max hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.max(),kkk.max(),lll.max()))
    print('Min hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.min(),kkk.min(),lll.min()))
    
    # Load pilatus volume
    vol = getvol(num)
    
    # Load data
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    # h,k,l min tol
    # diff only works along rows, doesn't calcualte difference between rows
    #min_tol_h = np.max(np.abs(np.diff(hhh)))
    #min_tol_k = np.max(np.abs(np.diff(kkk)))
    #min_tol_l = np.max(np.abs(np.diff(lll)))
    min_tol_h = np.max([np.max(np.abs(np.diff(hhh,axis=0))),np.max(np.abs(np.diff(hhh,axis=1))),np.max(np.abs(np.diff(hhh,axis=2)))])
    min_tol_k = np.max([np.max(np.abs(np.diff(kkk,axis=0))),np.max(np.abs(np.diff(kkk,axis=1))),np.max(np.abs(np.diff(kkk,axis=2)))])
    min_tol_l = np.max([np.max(np.abs(np.diff(lll,axis=0))),np.max(np.abs(np.diff(lll,axis=1))),np.max(np.abs(np.diff(lll,axis=2)))])
    print('Maximum hkl step per pixel = [%1.3g, %1.3g, %1.3g]' % (min_tol_h,min_tol_k,min_tol_l))
    
    # Centre of cuts
    if hkl_centre is None:
        if image_centre is None:
            i,j = pil_centre
            k = len(x)//2
        else:
            i,j,k = image_centre
        h_centre = hhh[i,j,k]
        k_centre = kkk[i,j,k]
        l_centre = lll[i,j,k]
        print('Peak Centre: (%1.3g,%1.3g,%1.3g)'%(h_centre,k_centre,l_centre))
    else:
        h_centre,k_centre,l_centre = hkl_centre
    
    # Pixels close to centre
    h_cen_idx = np.abs(hhh-h_centre) < sum_tolarance
    k_cen_idx = np.abs(kkk-k_centre) < sum_tolarance
    l_cen_idx = np.abs(lll-l_centre) < sum_tolarance
    
    kl_cen_idx = k_cen_idx*l_cen_idx
    
    # cut ranges
    h_list = np.arange(hhh[kl_cen_idx].min(),hhh[kl_cen_idx].max(),min_tol_h) 
    if len(h_list) > max_points:
        print('Reducing range to %i points' %max_points)
        h_list = np.arange(h_centre-min_tol_h*(max_points//2),h_centre+min_tol_h*(max_points//2)+min_tol_h,min_tol_h)
        #h_list = np.linspace(hhh[kl_cen_idx].min(),hhh[kl_cen_idx].max(),max_points) 
    h_step = h_list[1]-h_list[0]
    
    h_scan = np.zeros(len(h_list))
    print('Binning h axis at (H,%1.3g,%1.3g) in %1.0f steps, summing <%1.0f pixels per step' % (k_centre,l_centre,len(h_list),np.sum(kl_cen_idx)))
    for n in range(len(h_list)):
        hval = h_list[n]
        hidx = np.abs(hhh-hval) < h_step/2
        small_vol = vol[ hidx*kl_cen_idx ]
        h_scan[n] = np.sum(small_vol)
    
    return h_list,h_scan

def bin_pixel_k_cut(num,hkl_centre=None,image_centre=None,sum_tolarance=0.05,max_points=301):
    """
    Generate h,k,l cuts through a pilatus image scan
        k,Ik = bin_pixel_hkl_cut(num, [h,k,l], None, 0.01,1001)
    Inputs:
      num = scan number
      image_centre = [i,j,k] or None* - centre of the cuts on the detector e.g. [104,205,31]
      hkl_centre = [h,k,l] or None* - centre of the cuts in hkl, if None the peak position will be used
      sum_tolarance = dk - distance around hkl position to sum, in reciprocal lattice units
      max_points = nk - maximum number of points in each cut, to reduce calculatio time
    Returns:
      k - list of values
      Ik - binned intensities along l
    
    If image_centre is None, the centre of the detector and scan will be used.
    """
    
     # Load hkl matrices
    hhh,kkk,lll = pixel2hkl(num)
    print('Max hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.max(),kkk.max(),lll.max()))
    print('Min hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.min(),kkk.min(),lll.min()))
    
    # Load pilatus volume
    vol = getvol(num)
    
    # Load data
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    # h,k,l min tol
    # diff only works along rows, doesn't calcualte difference between rows
    #min_tol_h = np.max(np.abs(np.diff(hhh)))
    #min_tol_k = np.max(np.abs(np.diff(kkk)))
    #min_tol_l = np.max(np.abs(np.diff(lll)))
    min_tol_h = np.max([np.max(np.abs(np.diff(hhh,axis=0))),np.max(np.abs(np.diff(hhh,axis=1))),np.max(np.abs(np.diff(hhh,axis=2)))])
    min_tol_k = np.max([np.max(np.abs(np.diff(kkk,axis=0))),np.max(np.abs(np.diff(kkk,axis=1))),np.max(np.abs(np.diff(kkk,axis=2)))])
    min_tol_l = np.max([np.max(np.abs(np.diff(lll,axis=0))),np.max(np.abs(np.diff(lll,axis=1))),np.max(np.abs(np.diff(lll,axis=2)))])
    print('Maximum hkl step per pixel = [%1.3g, %1.3g, %1.3g]' % (min_tol_h,min_tol_k,min_tol_l))
    
    # Centre of cuts
    if hkl_centre is None:
        if image_centre is None:
            i,j = pil_centre
            k = len(x)//2
        else:
            i,j,k = image_centre
        h_centre = hhh[i,j,k]
        k_centre = kkk[i,j,k]
        l_centre = lll[i,j,k]
        print('Peak Centre: (%1.3g,%1.3g,%1.3g)'%(h_centre,k_centre,l_centre))
    else:
        h_centre,k_centre,l_centre = hkl_centre
    
    # Pixels close to centre
    h_cen_idx = np.abs(hhh-h_centre) < sum_tolarance
    k_cen_idx = np.abs(kkk-k_centre) < sum_tolarance
    l_cen_idx = np.abs(lll-l_centre) < sum_tolarance
    
    hl_cen_idx = h_cen_idx*l_cen_idx
    
    # cut ranges
    k_list = np.arange(kkk[hl_cen_idx].min(),kkk[hl_cen_idx].max(),min_tol_k)
    if len(k_list) > max_points:
        print('Reducing range to %i points' %max_points)
        k_list = np.arange(k_centre-min_tol_k*(max_points//2),k_centre+min_tol_k*(max_points//2)+min_tol_k,min_tol_k)
        #k_list = np.linspace(kkk[hl_cen_idx].min(),kkk[hl_cen_idx].max(),max_points)
    k_step = k_list[1]-k_list[0]
    
    k_scan = np.zeros(len(k_list))
    print('Binning k axis at (%1.3g,K,%1.3g) in %1.0f steps, summing <%1.0f pixels per step' % (h_centre,l_centre,len(k_list),np.sum(hl_cen_idx)))
    for n in range(len(k_list)):
        kval = k_list[n]
        kidx = np.abs(kkk-kval) < k_step/2
        small_vol = vol[ kidx*hl_cen_idx ]
        k_scan[n] = np.sum(small_vol)
    
    return k_list,k_scan

def bin_pixel_l_cut(num,hkl_centre=None,image_centre=None,sum_tolarance=0.05,max_points=301):
    """
    Generate h,k,l cuts through a pilatus image scan
        l,Il = bin_pixel_hkl_cut(num, [h,k,l], None, 0.01,101)
    Inputs:
      num = scan number
      image_centre = [i,j,k] or None* - centre of the cuts on the detector e.g. [104,205,31]
      hkl_centre = [h,k,l] or None* - centre of the cuts in hkl, if None the peak position will be used
      sum_tolarance = dl - distance around hkl position to sum, in reciprocal lattice units
      max_points = nl - maximum number of points in each cut, to reduce calcualtion time
    Returns:
      l - list of values
      Il - binned intensities along l
    
    If image_centre is None, the centre of the detector and scan will be used.
    """
    
     # Load hkl matrices
    hhh,kkk,lll = pixel2hkl(num)
    print('Max hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.max(),kkk.max(),lll.max()))
    print('Min hkl: (%1.3g,%1.3g,%1.3g)'%(hhh.min(),kkk.min(),lll.min()))
    
    # Load pilatus volume
    vol = getvol(num)
    
    # Load data
    x,y,dy,varx,vary,ttl,d = getdata(num)
    
    # h,k,l min tol
    # diff only works along rows, doesn't calcualte difference between rows
    #min_tol_h = np.max(np.abs(np.diff(hhh)))
    #min_tol_k = np.max(np.abs(np.diff(kkk)))
    #min_tol_l = np.max(np.abs(np.diff(lll)))
    min_tol_h = np.max([np.max(np.abs(np.diff(hhh,axis=0))),np.max(np.abs(np.diff(hhh,axis=1))),np.max(np.abs(np.diff(hhh,axis=2)))])
    min_tol_k = np.max([np.max(np.abs(np.diff(kkk,axis=0))),np.max(np.abs(np.diff(kkk,axis=1))),np.max(np.abs(np.diff(kkk,axis=2)))])
    min_tol_l = np.max([np.max(np.abs(np.diff(lll,axis=0))),np.max(np.abs(np.diff(lll,axis=1))),np.max(np.abs(np.diff(lll,axis=2)))])
    print('Maximum hkl step per pixel = [%1.3g, %1.3g, %1.3g]' % (min_tol_h,min_tol_k,min_tol_l))
    
    # Centre of cuts
    if hkl_centre is None:
        if image_centre is None:
            i,j = pil_centre
            k = len(x)//2
        else:
            i,j,k = image_centre
        h_centre = hhh[i,j,k]
        k_centre = kkk[i,j,k]
        l_centre = lll[i,j,k]
        print('Peak Centre: (%1.3g,%1.3g,%1.3g)'%(h_centre,k_centre,l_centre))
    else:
        h_centre,k_centre,l_centre = hkl_centre
    
    # Pixels close to centre
    h_cen_idx = np.abs(hhh-h_centre) < sum_tolarance
    k_cen_idx = np.abs(kkk-k_centre) < sum_tolarance
    l_cen_idx = np.abs(lll-l_centre) < sum_tolarance
    
    hk_cen_idx = h_cen_idx*k_cen_idx
    
    # cut ranges
    l_list = np.arange(lll[hk_cen_idx].min(),lll[hk_cen_idx].max(),min_tol_l)
    if len(l_list) > max_points:
        print('Reducing range to %i points' %max_points)
        l_list = np.arange(l_centre-min_tol_l*(max_points//2),l_centre+min_tol_l*(max_points//2)+min_tol_l,min_tol_l)
        #l_list = np.linspace(lll[hk_cen_idx].min(),lll[hk_cen_idx].max(),max_points)
    l_step = l_list[1]-l_list[0]
    
    l_scan = np.zeros(len(l_list))
    print('Binning l axis at (%1.3g,%1.3g,L) in %1.0f steps, summing <%1.0f pixels per step' % (h_centre,k_centre,len(l_list),np.sum(hk_cen_idx)))
    for n in range(len(l_list)):
        lval = l_list[n]
        lidx = np.abs(lll-lval) < l_step/2
        small_vol = vol[ lidx*hk_cen_idx ]
        l_scan[n] = np.sum(small_vol)
    
    return l_list,l_scan

def saveplot(name,dpi=None):
    "Saves current figure as a png in the savedir directory"
    
    if type(name) is int:
        name = str(name)
    
    gcf = plt.gcf()
    savefile = os.path.join(savedir, '{}.png'.format(saveable(name)))
    gcf.savefig(savefile,dpi=dpi)
    print( 'Saved Figure {} as {}'.format(gcf.number,savefile) )

def savedata(name,x,y,dy=None,header=None):
    """
    Save data as a formated text file in the analysis directory
      savedata( name, x, y, dy, header )
    
    name = name of file, not including extention (.dat will be added)
    x,y,dy = equal length arrays, dy can be left out
    header = header info can be added if required
    
    Uses: np.savetxt(savefile,(x,y,dy),header=head)
    
    E.G.
        x = [1,2,3,4]
        y = [2.4,56,43,1e6]
        dy = np.sqrt(y)
        savedata('test',x,y,dy,header='lovely data')
    """
    
    savefile = os.path.join(savedir, '{}.dat'.format(saveable(name)))
    if dy is None:
        np.savetxt(savefile,(x,y),header=header)
    else:
        np.savetxt(savefile,(x,y,dy),header=header)
    print( 'Data has been saved as {}'.format(savefile) )

def loaddata(name):
    """
    Load data saved using savedata
      x,y,dy = loaddata(name)
    
    name = str as given to savedata, no extention
    x,y,dy = arrays of data as saved
    
    E.G.
        x = [1,2,3,4]
        y = [2.4,56,43,1e6]
        dy = np.sqrt(y)
        savedata('test',x,y,dy,header='lovely data')
        x,y,dy = loaddata('test')
    """
    
    savefile = os.path.join(savedir, '{}.dat'.format(saveable(name)))
    return np.loadtxt(savefile)

def frange(start,stop=None,step=1):
    "Like np.arange but ends at stop, rather than stop-step"
    " A = frange(0,5,1) = [0.,1.,2.,3.,4.,5.]"
    
    if stop is None:
        stop = start
        start = 0
    
    return list(np.arange(start,stop+step*0.99,step,dtype=float))

def stfm(val,err):
    """
    Create standard form string from value and uncertainty"
     str = stfm(val,err)
     Examples:
          '35.25 (1)' = stfm(35.25,0.01)
          '110 (5)' = stfm(110.25,5)
          '0.0015300 (5)' = stfm(0.00153,0.0000005)
          '1.56(2)E+6' = stfm(1.5632e6,1.53e4)
    
    Notes:
     - Errors less than 0.01% of values will be given as 0
     - The maximum length of string is 13 characters
     - Errors greater then 10x the value will cause the value to be rounded to zero
    """
    
    # Determine the number of significant figures from the error
    if err == 0. or val/float(err) >= 1E5:
        # Zero error - give value to 4 sig. fig.
        out = '{:1.5G}'.format(val)
        if 'E' in out:
            out = '{}(0)E{}'.format(*out.split('E'))
        else:
            out = out + ' (0)'
        return out
    elif np.log10(np.abs(err)) > 0.:
        # Error > 0
        sigfig = np.ceil(np.log10(np.abs(err)))-1
        dec = 0.
    elif np.isnan(err):
        # nan error
        return '{} (-)'.format(val)
    else:
        # error < 0
        sigfig = np.floor(np.log10(np.abs(err))+0.025)
        dec = -sigfig
    
    # Round value and error to the number of significant figures
    rval = round(val/(10.**sigfig))*(10.**sigfig)
    rerr = round(err/(10.**sigfig))*(10.**sigfig)
    # size of value and error
    pw = np.floor(np.log10(np.abs(rval)))
    pwr = np.floor(np.log10(np.abs(rerr)))
    
    max_pw = max(pw,pwr)
    ln = max_pw - sigfig # power difference
    
    if np.log10(np.abs(err)) < 0:
        rerr = err/(10.**sigfig)
    
    # Small numbers - exponential notation
    if max_pw < -3.:
        rval = rval/(10.**max_pw)
        fmt = '{'+'0:1.{:1.0f}f'.format(ln)+'}({1:1.0f})E{2:1.0f}'
        return fmt.format(rval,rerr,max_pw)
    
    # Large numbers - exponential notation
    if max_pw >= 4.:
        rval = rval/(10.**max_pw)
        rerr = rerr/(10.**sigfig)
        fmt = '{'+'0:1.{:1.0f}f'.format(ln)+'}({1:1.0f})E+{2:1.0f}'
        return fmt.format(rval,rerr,max_pw)
    
    fmt = '{'+'0:0.{:1.0f}f'.format(dec+0)+'} ({1:1.0f})'
    return fmt.format(rval,rerr)

def saveable(string):
    " Removes bad characters from a string, so it can be used as a filename"
    # Special - replace # with S for scans
    string = string.replace('#','S')
    # Replace some characters with underscores
    for char in '#%{}\/<>@|':
        string = string.replace(char,'_') 
    # Replace other characters with nothing
    for char in '*$�!':
        string = string.replace(char,'') 
    return string

def rolling_fun(Y,fun=np.std):
    """
    Performs a rolling verions of the given function across and array
      rolling_fun(Y,fun)
       Y = array, e.g. [1,2,3,4]
       fun = array function, e.g. np.mean or np.std
    Returns an array, of the same length as Y, with function fun performed across every 3 items
    """
    
    s=[fun(Y[n-1:n+2]) for n in range(1,len(Y)-1)]
    s=[s[0]]+s+[s[-1]] # add first and last values
    return np.asarray(s)

def findranges(scannos,sep=':'):
    "Convert a list of numbers to a simple string"
    
    scannos = np.sort(scannos).astype(int)
    
    dif = np.diff(scannos)
    
    stt,stp,rng = [scannos[0]],[dif[0]],[1]
    for n in range(1,len(dif)):
        if scannos[n+1] != scannos[n]+dif[n-1]:
            stt += [scannos[n]]
            stp += [dif[n]]
            rng += [1]
        else:
            rng[-1] += 1
        #print(n,stt,stp,rng)
    stt += [scannos[-1]]
    rng += [1]
    
    out = []
    x = 0
    while x < len(stt):
        if rng[x] == 1:
            out += ['{}'.format(stt[x])]
            x += 1
        elif stp[x] == 1:
            out += ['{}{}{}'.format(stt[x],sep,stt[x+1])]
            x += 2
        else:
            out += ['{}{}{}{}{}'.format(stt[x],sep,stp[x],sep,stt[x+1])]
            x += 2
    return ','.join(out)

def mini_string_range(scannos,sep=':'):
    "Convert a list of numbers to a simple string"
    
    scannos = np.sort(scannos).astype(int).astype(str)
    
    n = len(scannos[0])
    while np.all([scannos[0][:-n] == x[:-n] for x in scannos]): 
        n -= 1
    
    if n == len(scannos[0]):
        return '{}-{}'.format(scannos[0],scannos[-1])
    
    inistr = scannos[0][:-(n+1)]
    strc = [i[-(n+1):] for i in scannos]
    liststr = findranges(strc,sep=sep)
    return '{}[{}]'.format(inistr,liststr)

def numbers2string(scannos,sep=':'):
    "Convert a list of numbers to a simple string"
    
    if type(scannos) is str or type(scannos) is int or len(scannos) == 1:
        return str(scannos)
    
    string1 = findranges(scannos,sep)
    if len(string1) < 40:
        return string1
    
    string2 = mini_string_range(scannos,sep)
    if len(string2) < 40:
        return string2
    
    scannos = np.sort(scannos).astype(str)
    string3 = '{}{}{}'.format(scannos[0],sep,scannos[-1])
    return string3 

def maskvals(x,y,dy,mask_cmd):
    """
    Returns arrays masked by requested rules
    Allows you to remove regions of a scan easily.
    e.g. x,y,dy = maskvals(x,y,dy,['x<0'])
    
    Masks should be a list of strings
    Each string should define an evaluatable statement that defines a region of x or y
    Each defined region will be REMOVED
    examples: 
        'x<-1'
        'x>1'
        'np.abs(x-34)>0.1']
        [x<-1.2','x>1.2']
    """
    
    if type(mask_cmd) is str:
        mask_cmd = [mask_cmd]
    
    mask = np.ones(len(x),dtype=bool)
    for m in mask_cmd:
        mask[ eval(m) ] = False
    return x[mask],y[mask],dy[mask]

def getRAM():
    """
    Get available RAM of system, or give 1e8 bits
    """
    
    try:
        import psutil # get RAM available
    except ImportError:
        return 1e8
    return psutil.virtual_memory()[1]*0.5

def nexus_rsremap(filename):
    """
    Load rs_remap file, usually in exp_dir/processed/
    returns hh, kk, ll, volume
    """
    n = nxload(filename)
    hval = np.array(n['/processed/process/reciprocal_space/h-axis'])
    kval = np.array(n['/processed/process/reciprocal_space/k-axis'])
    lval = np.array(n['/processed/process/reciprocal_space/l-axis'])
    volume = np.array(n['/processed/process/reciprocal_space/volume'])
    hh, kk, ll = np.meshgrid(hval, kval, lval)
    return hh, kk, ll, volume

def nexus_plot_rsremap(filename, n=None):
    """
    Example plot of rs_remap
    """

    hh, kk, ll, volume = nexus_rsremap(filename)

    if n is None:
        n = ll.shape[2]//2
    lval = ll[0,0,n]
    ttl = '\n'.join(os.path.split(filename))

    plt.figure(figsize=[12,10], dpi=fig_dpi)
    plt.pcolormesh(hh[:,:,n], kk[:,:,n], volume[:,:,n].T)
    plt.colorbar()
    plt.axis('image');
    labels(ttl, '(h,0,%4.2f)'%lval, '(0,k,%4.2f)'%lval)

class dict2obj(OrderedDict):
    "Convert dictionary object to class instance"
    def __init__(self,dictvals,order=None):
        # Initialise OrderedDict (not sure which of these is correct)
        super(dict2obj, self).__init__()
        #OrderedDict.__init__(self)
        
        if order is None:
            order = dictvals.keys()
        
        for name in order:
            setattr(self,name,dictvals[name])
            self.update({name:dictvals[name]})
