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
sys.path.insert(0,'/dls_sw/i16/software/python/userscripts/python') # location of Py16Progs
import Py16Progs as p16
p16.filedir = '/dls/i16/data/2015/mt0000-1' # your experiment folder
p16.savedir = '' # where to save output

# Use functions:
p16.checkexp()
num = p16.latest()
d= p16.readscan(num)
p16.plotscan(num)
p16.plotpil(num)

******In Console*******
Run the file in the current console

Load the module:
runfile('/dls_sw/i16/software/python/userscripts/I16user/Py16/Py16progs.py')

Change filedir: 
filedir = '/dls/i16/data/2015/cm12169-2/CsFeSe' # your experiment folder

# Use functions:
d= readscan(num)
checkexp()
num = latest()

**********************
Functions:
    d = readscan(num) 
    x,y,dy,varx,vary,ttl,d = getdata(num/d,'varx','vary',save=None)
    x,y,z,varx,vary,varz,ttl = joindata([nums],'varx','vary','varz')
    vol = getvol(num)
    ROI_sum,ROI_maxval,ROI_bkg = pilroi(num,ROIcen,ROIsize,findpeak,peakregion)
    
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
    plotscans3D(runs,depvar='Ta',vary='',save=None)
    
    fit,err = fit_scans(runs,depvar='Ta',vary='',save=None)
    
    out,err = peakfit(x,y,dy,type='dVoight')
    
    ROIcen_ij,ROIcen_frame = pilpeak(Y,test = 1,disp=False)
    ispeak(Y,test = 1,disp=False)
    peak = peakfind(Y,cutoff=0.01)*
    
    labels(title,xlabel,ylabel,zlabel)
    vals = frange(start,stop=None,step=1)
    str = stfm(val,err)
    

Version 2.0
Last updated: 21/09/16

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

###FEEDBACK### Please submit your bug reports, feature requests or queries to: dan.porter@diamond.ac.uk

@author: Dan Porter
I16, Diamond Light Source
2016
"""

"""
Ideas for the future:
 - Pilatus peak finding (2D fitting?) as separate function to dpilroi
 - dataloader contain functions (more pythony...)
 - Refactor fitting routines in new Py16fit module

"""

import os
import glob # find files
import re # regular expressions
import datetime # Dates and times
import numpy as np
#import scisoftpy as dnp # Make sure this is in your python path
import matplotlib.pyplot as plt # Plotting
import matplotlib.ticker as mtick # formatting of tick labels
from mpl_toolkits.mplot3d import Axes3D # 3D plotting
from scipy.optimize import curve_fit # Peak fitting
from scipy import misc # read pilatus images
from scipy.signal import convolve
from itertools import product


###########################################################################
#############################PARAMETERS####################################
###########################################################################
"-----------------------Default Experiment Directory----------------------"
# Variable filedir is called from the namespace and can be changed at any 
# time,even once the module has been imported, if you change it in the current namespace

filedir = '/dls/i16/data/2016' 
savedir = '/home/i16user/Desktop'

"-----------------------Error Estimation Parameters-----------------------"
error_func = lambda x: np.sqrt(np.abs(x)+0.1) # Define how the error on each intensity is estimated
# error_func = lambda x: 0*x + 1 # Switch errors off

"-------------------------Normalisation Parameters------------------------"
exp_ring_current = 300.0 # Standard ring current for current experiment for normalisation
exp_monitor = 800.0 # Standard ic1monitor value for current experiment for normalisation
normby = 'rc' # Incident beam normalisation option: 'rc','ic1' or 'none'

"----------------------------Pilatus Parameters---------------------------"
pil_centre = [110,242] # Centre of pilatus images (for current value see localStation.py
hot_pixel = 2**20-100 # Pixels on the pilatus greater than this are set to the intensity chosen by dead_pixel_func
peakregion=[7,153,186,332] # Search for peaks within this area of the detector [min_y,min_x,max_y,max_x]
dead_pixel_func = np.median # Define how to choose the replaced intensity for hot/broken pixels 

"----------------------------Plotting Parameters--------------------------"
plot_colors = ['b','g','m','c','y','k','r'] # change the default colours of plotscan 
exp_title = '' # Experiment title
#plt.xkcd()

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
    meta = {}
    lineno = 0
    for ln in lines:
        lineno += 1
        if '&END' in ln: break
        ln = ln.strip(' ,\n')
        neq = ln.count('=')
        if neq == 1:
            'e.g. cmd = "scan x 1 10 1"'
            inlines = [ln]
        elif neq > 1:
            'e.g. SRSRUN=571664,SRSDAT=201624,SRSTIM=183757'
            inlines = ln.split(',')
        else:
            'e.g. <MetaDataAtStart>'
            continue
        
        for inln in inlines:
            vals = inln.split('=')
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
    main = {}
    for name,value in zip(names,vals.T):
        main[name] = value
    
    # Convert to class instance
    d = dict2obj(main)
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
    
    if os.path.isdir(filedir) == False: 
        print( "I can't find the directory: {}".format(filedir) )
        return None
    
    if num < 1: 
        if latest() is None: return None
        num = latest()+num
    
    file = os.path.join(filedir, '%i.dat' %num)
    try:
        d = read_dat_file(file)
        #d = dnp.io.load(file,warn=False) # from SciSoftPi
    except:
        print( "Scan {} doesn't exist".format(num) )
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
        d.metadata.cmd=re.sub('[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',caller,cmd) # find numbers in command and round them to 3dp
    
    " Correct re-assigned values"
    " For some reason, parameter names change occsaionaly, add correction here"
    if 'en' in d.metadata.keys(): d.metadata.Energy = d.metadata.en
    if 'count_time' in d.keys(): d.t = d.count_time
    if 'energy2' in d.keys(): d.energy = d.energy2
    if 'rc' not in d.keys(): d.rc = exp_ring_current*np.ones(len(d[d.keys()[0]]))
    if 'TimeSec' not in d.keys(): d.TimeSec = np.arange(0,len(d[d.keys()[0]]))
    
    " Correct psi values"
    if d.metadata.psi < -1000: d.metadata.psi = 0.0
    if d.metadata.psi == 'Unavailable': d.metadata.psi = 0.0
    " Update keys"
    d.update(d.__dict__)
    d.metadata.update(d.metadata.__dict__)
    return d

def getdata(num=None,varx='',vary='',norm=True,abscor=None):
    """
    Get useful values from runs
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
    
    " Use scan command to determine variables of scan and title"
    HKL = '({0:1.3g},{1:1.3g},{2:1.3g})'.format(round(m.h,2)+0.0,round(m.k,2)+0.0,round(m.l,2)+0.0)
    Energy = '{0:1.4g}keV'.format(m.Energy)
    Temp = '{0:1.3g}K'.format(m.Ta)
    if Temp == '0K': Temp = '300K'
    pol = ''
    if m.delta_axis_offset < 1 and m.thp > 10:
        if m.gam > 0.:
            if m.stoke < 45.:
                pol = ' p-p'
            else:
                pol = ' p-s'
        else:
            if m.stoke < 45.:
                pol = ' s-s'
            else:
                pol = ' s-p'
    
    "---exp_title---"
    if exp_title is '':
        etitle = os.path.basename(filedir)
    else:
        etitle = exp_title
    
    "---Generate title---"
    ttl = '{} #{} {} {} {}{}'.format(etitle,m.SRSRUN,HKL,Energy,Temp,pol).strip()
    
    "---Determine scan variables from scan command---"
    if varx not in keys:
        " Determine dependent variable from scan command - should be the second word"
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
            varx = 'phi'
        elif varx == 'th2th':
            varx = 'delta'
            
        if varx not in keys:
            varx = keys[0]
        
    "***Get x values***"
    x = getattr(d,varx)
    
    "***Get y values***"
    if 'nroi' in vary.lower():
        " y values from automatic peak search in pilatus"
        # e.g. 'nroi'           > Defaults to ROI2 in centre of pilatus
        # e.g. 'nroi[11,11]'    > in centre with size [11,11]
        # e.g. 'nroi[109,241,11,11] > centred at [109,241] with size [11,11]
        # e.g. 'nroibkg         > As 'nroi', with background subtraction
        # e.g. 'nroipeak'       > As 'nroi', with peak search 
        # e.g. 'nroipeak[11,11]'> at peak centre, with size [11,11]
        # e.g. 'nroipeakbkg[31,31]' > at peak centre with background subtraction
        
        vol = getvol(m.SRSRUN) # Load pilatus images as 3D volume
        
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
        
        " Calculate errors"
        dy = error_func(y)
        
        " Background"
        if 'bkg' in vary.lower():
            y = y-ROI_bkg
            dy = np.sqrt( dy**2 + error_func(ROI_bkg)**2 )
        
    else:
        if vary not in keys:
            " Determine indepdendent variable from scan command - should be the last word"
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
            
            if vary not in keys:
                vary = keys[-1]
        
        " y values from data file"
        y = getattr(d,vary)
        dy = error_func(y)
        
    
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
    if norm:
        y = y/d.metadata.Transmission/cnt/inorm
        dy = dy/d.metadata.Transmission/cnt/inorm
        vary += '/Tran./t/'+normtxt
    
    "---Apply absorption correction---"
    if abscor is not None:
        # Absorption correction based on a flat plate perpendicular to the scattering plane
        # Requires abscor = absorption coefficient of material
        print( 'Not done yet!' )

    return x,y,dy,varx,vary,ttl,d

def joindata(nums=None,varx='',vary='Energy',varz='',norm=True,abscor=None,save=False):
    """
    Get useful values from runs, join multiple runs together, output joined arrays
     x,y,z,varx,vary,varz,ttl = joindata([nums],'varx','vary','varz')
     
      x = nxm array of scanned values, where n is the length of the scan and m is the number of runs
      y = nxm array of values that change with each scan
      z = nxm array of measured values
      varx,vary,varz = names of scan values
      ttl = title generated for these scans
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
        ini_cmd = d.metadata.cmd
        
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
    
    "Generate title"
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
    fmt = '{ttl}{xvar}-{yvar} #{rn1:1.0f}-{rn2:1.0f} {hkl} {en} {temp} {pol}'
    if out_varx in ['Energy','en'] or vary in ['Energy','en']:
        energy=''
    if out_varx in ['Ta','Tb'] or vary in ['Ta','Tb']:
        temp=''
    out_ttl = fmt.format(ttl=ettl,rn1=nums[0],rn2=nums[-1],
                         xvar=out_varx,yvar=vary,
                         hkl=hkl,en=energy,temp=temp,pol=pol)
    
    " Save a .dat file"
    " Data can be loaded with x,y,z,dz = np.loadtxt(file.dat)"
    if save not in [None, False, '']:
        # save data to file
        if type(save) is str:
            savefile = os.path.join(savedir, '{}.dat'.format(save))
            head = '{}\n{}\n{}, {}, {}, {}'.format(out_ttl,ini_cmd,varx,vary,varz,'error_'+varz)
            np.savetxt(savefile,(x,y,z,dz),header=head)
            print( 'Saved as {}'.format(savefile) )
        else:
            savefile = os.path.join(savedir, '{} dep {}-{}.dat'.format(vary,runs[0],runs[-1]))
            head = '{}\n{}\n{}, {}, {}, {}'.format(out_ttl,ini_cmd,varx,vary,varz,'error_'+varz)
            np.savetxt(savefile,(x,y,z,dz),header=head)
            print( 'Scan #{} has been saved as {}'.format(num,savefile) )
    
    return storex,storey,storez,out_varx,vary,out_varz,out_ttl

def getvol(num,ROIcen=None,ROIsize=None):
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
    if type(num)==int or type(num)==np.int64:
        d = readscan(num)
    else:
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
    im=misc.imread(file)
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
    
    "Load each image into a single volume"
    numimag = len(d.path)
    vol = np.zeros([pil_size[0],pil_size[1],numimag]) # [195,487,~31]
    for n in range(numimag):
        " Prepare file name"
        tif=pilpath % d.path[n]
        file = os.path.join(filedir,tif)
        file=file.replace('/',os.path.sep)
        
        " Load image "
        #t=dnp.io.load(file,warn=False)
        #vol[:,:,n] = t.image0 #"image0" varies for some reason
        im=misc.imread(file) # this is more reliable than dnp.io.load
        
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

def pilroi(vol,ROIcen=None,ROIsize=[31,31],disp=False):
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
        ROI_maxval[n] = error_func( ROI.max() )
        
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
    
    newest = max(ls, key=os.path.getctime)
    num = np.int(os.path.split(newest)[-1][:-4])
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

def checkscan(num1=None,num2=None,showval=None):
    """
    Get run number information, 
        checknum(#)             Display lots of information about a single scan
        checknum(first,last)      List all the scans from first to last with brief info
        checknum([run1,run2,...]) List all the scans in the list or array
    
    if # <= 0, # will be given as latest()-#
    
    checknum(...,showval='phi'): also display infomation from a variable, such as 'phi'
    
    """
    
    if os.path.isdir(filedir) == False: 
        print( "I can't find the directory: {}".format(filedir) )
        return
    
    if num1 is None:
        num1 = latest()
    
    if num2 is None:
        num2 = num1
    
    #if num1 < 1: num1 = latest()+num1
    #if num2 < 1: num2 = latest()+num2
    
    if type(num1) in [int,np.int64] :
        num = range(num1,num2+1)
    else:
        num = num1 # num = list or array
    
    "-----------------Multi run------------------"
    if len(num)>1:
        # Print brief info
        #fmt = '{num} | {date} | {mode:4s} {energy:5.3g}keV {temp:5.3g}K ({h:1.2g},{k:1.2g},{l:1.2g}) | {cmd}'
        fmt = '{num} | {date} | {mode:4s} {ss} {ds} {energy:5.3g}keV {temp:5.3g}K {hkl:17s} {show}{equal}{val} | {cmd}'
        showval_dict = {'show':'','equal':'','val':''}
        for n in range(len(num)):
            d = readscan(num[n])
            if d is None:
                print( "File does not exist" )
                continue
            m = d.metadata
            ks = d.keys()
            
            if m.Ta == 0: m.Ta = 300
            
            sampsl = '{0:4.2g}x{1:<4.2g}'.format(m.s5xgap,m.s5ygap)
            detsl = '{0:4.2g}x{1:<4.2g}'.format(m.s6xgap,m.s6ygap)
            
            h = round(m.h*100)//100 + 0.0 # + 0.0 to remove -0 
            k = round(m.k*100)//100 + 0.0
            l = round(m.l*100)//100 + 0.0
            hkl = '({h:1.2g},{k:1.2g},{l:1.2g})'.format(h=h,k=k,l=l)
            if m.gam>0.0:
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
            
            if showval is not None:
                if showval in ks:
                    val = np.mean(getattr(d,showval))
                elif showval in m.keys():
                    val = getattr(m,showval)
                else:
                    val = 'No Data'
                showval_dict['show'] = showval
                showval_dict['equal'] = '='
                showval_dict['val'] = val
            
            print( fmt.format(num=m.SRSRUN,date=m.date,mode=mode,ss=sampsl,ds=detsl,energy=m.Energy,temp=m.Ta,hkl=hkl,cmd=m.cmd,**showval_dict) )
        return 
    
    "----------------Single run------------------"
    d = readscan(num[0])
    if d is None:
        print( "File does not exist!" )
        return d
    m = d.metadata
    ks = d.keys()
    
    if m.Ta == 0: m.Ta = 300
    
    # Print information
    print( '-----------Run ', m.SRSRUN, '-----------' )
    print( '  File Dir: ',filedir )
    print( '   Command: ',m.cmd )
    print( '   Npoints: ',len(d.TimeSec) )
    print( '       HKL: ({0},{1},{2})'.format(m.h,m.k,m.l) )
    print( '    Energy: {0} keV'.format(m.Energy) )
    print( '      Temp: {0} K'.format(m.Ta) )
    print()
    
    # Check for vertical or horizontal geopmetry
    if m.gam > 0.1:
        print( '    Horizontal Geometry' )
        print( '        Mu: {0}'.format(m.mu) )
        print( '       Chi: {0}'.format(m.chi) )
        print( '     Gamma: {0}'.format(m.gam) )
    else:
        print( '      Vertical Geometry' )
        print( '       Eta: {0}'.format(m.eta) )
        print( '       Chi: {0}'.format(m.chi) )
        print( '     Delta: {0}'.format(m.delta) )
    print( 'Psi({0},{1},{2}): {3}'.format(m.azih,m.azik,m.azil,m.psi) )
    print()
    print( '       Sx: {0} mm'.format(m.sx) )
    print( '       Sy: {0} mm'.format(m.sy) )
    print( '       Sz: {0} mm'.format(m.sz) )
    print()
    print( '  Sample Slits: {0:4.2f}x{1:4.2f} mm'.format(m.s5xgap,m.s5ygap) )
    print( 'Detector Slits: {0:4.2f}x{1:4.2f} mm'.format(m.s6xgap,m.s6ygap) )
    print()
    
    # Minimirrors
    if m.m4x > 0.1: 
        mm = 'in' 
    else: 
        mm = 'out'
    print( 'Minimirrors: {} ({:4.2f} deg)'.format(mm,m.m4pitch) )
    
    if 'sum' in ks:
        # Pilatus
        print()
        print( 'Detector: Pilatus (do={})'.format(m.delta_axis_offset) )
        print( '   Count: {0:1.3g}s'.format(np.mean(d.t)) )
        print( ' Max val: {0:5.3g}'.format(max(d.maxval)) )
        
    if 'APD' in ks:
        # APD
        print()
        print( 'Detector: APD' )
        if m.gam > 0.:
            # Horizontal
            if m.stoke < 45.:
                print( '   Pol: {0} (pi-pi)'.format(m.stoke) )
            else:
                print( '   Pol: {0} (pi-sigma)'.format(m.stoke) )
        else:
            # Vertical
            if m.stoke < 45.:
                print( '   Pol: {0} (sigma-sigma)'.format(m.stoke) )
            else:
                print( '   Pol: {0} (sigma-pi)'.format(m.stoke) )
        print( '   Count: {0:1.3g}s'.format(np.mean(d.counttime)) )
        print( ' Max val: {0:5.3g}'.format(max(d.APD)) )
    
    if 'FF' in ks:
        print()
        print( 'Detector: Vortex' )
        print( '   Count: {0:1.3g}s'.format(np.mean(d.count_time)) )
        print( '   ROIs (maxval):' )
        ROIs = [n for n in d.keys() if 'Element' in n]
        for roi in ROIs:
            print( '    {} ({})'.format(roi,max(getattr(d,roi))) )
        
    # Attenuation
    print()
    print( '    Atten: {0} ({1}%)'.format(m.Atten,m.Transmission*100) )
    print()
    
    # additional info
    if showval is not None:
        if showval in ks:
            val = np.mean(getattr(d,showval))
        elif showval in m.keys():
            val = getattr(m,showval)
        else:
            val = 'No Data'
        print( showval,' = ',val )
        print()
    
    # Timing
    time = d.TimeSec[-1]-d.TimeSec[0]
    hours = np.int(np.floor(time/3600))
    mins = np.int(np.floor(np.remainder(time,3600)/60))
    secs = np.remainder(np.remainder(time,3600),60)
    print( 'Ran on {}'.format(m.date) )
    print( 'Time taken: {} hours, {} mins, {} seconds'.format(hours,mins,secs) )
    return 

def checklog(time=None,mins=2,cmd=False,find=None):
    """Look at experiment log file 
    checklog(time,mins,commands,findstr)
      time = None - Uses current time (default)
      time = scan number (uses start time of this scan), 0 gives latest scan
      time = [hour,(min),(day),(month),(year)], () values set to default, defaults are time of last scan
      time = '2015-07-07 13:50:08,000' 
      mins = log file will display from time-mins to time
      mins = 'all' - display whole log from first scan time
      commands = True/False - only display command strings
      find = only display lines that include findstr
    """
    
    # time is a datetime object
    if time is None:
        #date = readscan(0).metadata.date
        #time = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
        time = datetime.datetime.now()
    
    if type(time) is str:
        time=datetime.datetime.strptime(time,'%Y-%m-%d %H:%M:%S,%f')
    
    if type(time) is not type(datetime.datetime.now()):
        if type(time) is int:
            # input is a scan number (0 for latest)
            date = readscan(time).metadata.date
            time = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
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
    
    filename = os.path.join(filedir, 'gdaterminal.log')
    file = open(filename)
    
    for line in file:
        if cmd and '>>>' not in line:
            continue
        
        if find is not None and find not in line:
            continue
        
        spt = line.split()
        tim=datetime.datetime.strptime(spt[0]+spt[1],'%Y-%m-%d%H:%M:%S,%f')
        
        if tim > mintime and tim < time:
            print( line.rstrip() )
    
    return

def prend(start=0,end=None):
    "Calculate the end time of a run"
    
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
            scanlen = len(np.range(scan_start,scan_end,scan_step))
        
        nrem = scanlen-Ncomplete
        
        tdif = t2-t1
        
        # Predict end time
        tperrun = tdif / Ncomplete
        trem = nrem*tperrun
        tend = t2 + trem
        
        print( 'Scan number: #',latest() )
        print( 'Scan Started: ',t1 )
        print( 'Points complete: ',Ncomplete,'/',scanlen )
        print( 'Time per point: ',tperrun )
        print( 'Still to go: ',nrem,' (',trem,')' )
        print( 'Scan will end: ',tend )
        return
    
    st = readscan(start)
    nd = readscan(end)
    
    m = st.metadata
    t1=datetime.datetime(m.Year,m.Month,m.Day,m.Hours,m.Minutes,m.Seconds,0)
    
    """ If run has already finished """
    if nd is not None:
        m = nd.metadata
        t2=datetime.datetime(m.Year,m.Month,m.Day,m.Hours,m.Minutes,m.Seconds,0)
        
        tdif = t2-t1
        
        print( 'Run started: ',t1 )
        print( 'Run ended: ',t2 )
        print( 'Run took: ',tdif )
        return
    
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
    
    print( 'Run Started: ',t1 )
    print( 'Latest scan: #',cur )
    print( 'Time per scan: ',tperrun )
    print( 'Runs completed: ',ndif,' (',tdif,')' )
    print( 'Still to go: ',nrem,' (',trem,')' )
    print( 'Run will end: ',tend )
    return

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

def polflip(sigsig,sigpi,fit='Gauss',output=False,plot=False):
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
        saveplot('POLFLIP '+saveable(ttl1))

def polenergy(sigsig,sigpi,background=None,vary='',bkg_scale=None,flipping_ratio=None,low_points=5,save=False):
    "Create Plot of energy-polarisation scans and calculate the subtraction"
    
    " Get the signal data - measured at a resonant feature in ss and sp"
    x1,y1,dy1,varx1,vary1,ttl1,d1 = getdata(sigsig,vary=vary) # scan 1
    x2,y2,dy2,varx2,vary2,ttl2,d2 = getdata(sigpi,vary=vary) # scan 2
    " Get the background data - measured away from the Bragg peak in sp"
    xb,yb,dyb,varxb,varyb,ttlb,db = getdata(background,vary=vary) # Background  
    
    " Get metadata"
    m = d2.metadata
    cmd = m.cmd # Scan command
    hkl = '({:3.1f},{:3.1f},{:3.1f})'.format(m.h,m.k,m.l)
    T = '{:1.3g}K'.format(m.Ta)
    sampsl = '{0:4.2g}x{1:<4.2g}'.format(m.s5xgap,m.s5ygap)
    detsl = '{0:4.2g}x{1:<4.2g}'.format(m.s6xgap,m.s6ygap)
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
    fig = plt.figure(figsize=[10,8])
    
    plt.plot(nprx,npry/flipping_ratio,'-ob',linewidth=2,label=nprlab)
    plt.plot(prx,pry,'-og',linewidth=2,label=prlab)
    plt.plot(xb,yb*bkg_scale,'-ok',linewidth=2,label=bkglab)
    plt.plot(prx,DIF,'-or',linewidth=2,label=diflab)
    
    plt.legend(loc=0, fontsize=16)
    plt.ylabel(vary2, fontsize=18)
    plt.xlabel('Energy (keV)', fontsize=18)
    #plttl = ttl2+'\n'+cmd+'\npsi={}, ss ={}, ds ={}, atten = {}'.format(sampsl,detsl,atten1)
    plt.title(ttl, fontsize=14)
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('PolEng_{}_{}_{}.png'.format(hkl,T,psival))

def scanabscor(num=0,u=1,eta_offset=0.0,chi_offset=0.0):
    """
    Calculate absorption correction
     A = abscor(num,u)
    """
    
    " Get data"
    try: 
        d = readscan(num)
    except:
        d = num
    
    " Get angles"
    eta = d.metadata.eta - eta_offset
    chi = d.metadata.chi - chi_offset
    delta = d.metadata.delta
    #mu = np.deg2rad(d.metadata.mu)
    #gam = np.deg2rad(d.metadata.gam)
    
    return abscor(eta,chi,delta,u)

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
            print( '{:>20} : {:<20}'.format(k,d1.metadata[k]) )
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
                print( '{:>20} : {:10} : {:<10} {}'.format(k,d1.metadata[k],d2.metadata[k],diff) )
            except:
                print( '{} does not exist in #{}'.format(k,d2.metadata.SRSRUN) )

def checkpeaks(num,test=1,vary=''):
    "Check multiple runs for peaks"
    
    # Turn num into a list
    try:
        N = len(num)
    except:
        num=[num]
    
    for n in num:
        x,y,dy = getdata(n,vary=vary)[:3]
        rat = ispeak(y,dy,return_rat=True)
        print( '{} {:8.2f} {}'.format(n,rat,rat>test) )

def savescan(num=None,varx='',vary='',abscor=None):
    "Save scan as .dat file"
    
    # load data
    x,y,dy,varx,vary,ttl,d = getdata(num,varx,vary,abscor)
    
    # save data to file
    savefile = os.path.join(savedir, '{}_{}.dat'.format(num,saveable(vary)))
    head = '{}\n{}\n{}, {}, {}'.format(ttl,d.metadata.cmd,varx,vary,'error_'+vary)
    np.savetxt(savefile,(x,y,dy),header=head)
    print( 'Scan #{} has been saved as {}'.format(num,savefile) )

def loadscan(num,vary=''):
    "Load a scan.dat file saved with savescan"
    
    filename = os.path.join(savedir, '{}_{}.dat'.format(num,saveable(vary)))
    
    # Get Header Data
    with open(filename,'r') as ff:
        # Line 1 = title
        ttl = ff.readline().strip('# ')
        # Line 2 = scan command
        cmd = ff.readline().strip('# ')
        # Line 3 = variable names
        varx,vary,dvary = ff.readline().strip('# ').split(',')
        
    # Get x, y, dy data
    x,y,dy = np.loadtxt(filename)
    
    return x,y,dy,varx,vary,ttl

def create_analysis_file(runs,depvar='Ta',vary='',varx='',fit_type = 'pVoight',bkg_type='flat',peaktest=1,
                  abscor=None,plot='all',show_fits=True,mask_cmd=None,estvals=None,xrange=None,sortdep=True,
                  Nloop=10, Binit=1e-5, Tinc=2, change_factor=0.5, converge_max=100, min_change=0.01,
                  save=False,saveFIT=False):
    "Creates a new python analysis script in the anaysis folder"
    
    ttl = saveable(exp_title)
    filename = savedir + '/I16_Analysis_{}_{}.py'.format(ttl,depvar)
    
    n = 1
    while os.path.isfile(filename):
        n += 1
        filename = savedir + '/I16_Analysis_{}_{}_{}.py'.format(ttl,depvar,n)
    
    with open(filename,'w') as f:
        # Write comments at top
        f.write('# Diamond I16 Analysis Script\n')
        f.write('# {}\n'.format(exp_title))
        f.write('# Analysis of {} dependence\n'.format(depvar))
        f.write('# \n')
        f.write('# By User\n')
        f.write('# {}\n\n'.format(datetime.datetime.strftime(datetime.datetime.now(),'%d/%m/%Y')))
        
        # Import stuff
        f.write('import sys,os\n')
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
        
        # Normalisation Options and Pilatus Centre
        f.write('# Experiment Parameters\n')
        f.write('dp.exp_ring_current = {}\n'.format(exp_ring_current))
        f.write('dp.exp_monitor = {}\n'.format(exp_monitor))
        f.write('dp.normby = \'{}\'\n'.format(normby))
        f.write('dp.pil_centre = {}\n'.format(pil_centre))
        f.write('dp.exp_title = \'{}\'\n\n'.format(exp_title))
        
        # Scan numbers
        f.write('# Scan numbers\n')
        try:
            if np.max(np.diff(np.diff(runs))) == 0:
                f.write('scans = range({},{},{})\n\n'.format(runs[0],runs[-1],runs[1]-runs[0]))
            else:
                f.write('scans = {} \n\n'.format(str(list(runs))))
        except:
            f.write('scans = {} \n\n'.format(str(list(runs))))
        
        if mask_cmd is not None:
            try:
                mask_cmd = str( mask_cmd.tolist() )
            except:
                mask_cmd = str( list(mask_cmd) )
        
        # Fit
        f.write('# Fitting\n')
        f.write('mask_cmd = {}\n'.format(mask_cmd))
        f.write('estvals = {}\n'.format(str(estvals)))
        f.write('fitopt = dict(depvar=\'{}\',vary=\'{}\',varx=\'{}\',fit_type = \'{}\',bkg_type=\'{}\',peaktest={},\n'.format(depvar,vary,varx,fit_type,bkg_type,peaktest))
        f.write('              abscor={},plot=\'{}\',show_fits={},mask_cmd=mask_cmd,estvals=estvals,xrange={},sortdep={},\n'.format(abscor,plot,show_fits,xrange,sortdep))
        f.write('              Nloop={}, Binit={}, Tinc={}, change_factor={}, converge_max={}, min_change={},\n'.format(Nloop,Binit,Tinc,change_factor,converge_max,min_change))
        f.write('              save={},saveFIT={})\n\n'.format(save,saveFIT))
        f.write('fit,err = dp.fit_scans(scans,**fitopt)\n\n')
        # Load
        f.write('# Load fitted data:\n')
        f.write('#fit,err = dp.load_fits([{},{}],depvar=\'{}\',fit_type=\'{}\')\n\n'.format(runs[0],runs[-1],depvar,fit_type))
    
    print( 'New Analysis file written to ',filename )

def get_all_scannos():
    "Returns the scan numbers available in filedir"
    
    # Get all data files in folder
    ls=glob.glob(filedir+os.path.sep+'*.dat')
    ls = np.sort(ls)
    
    # Convert to int
    scannos = [np.int(os.path.split(file)[-1][:-4]) for file in ls]
    return scannos


"----------------------------Analysis Functions---------------------------"


def fit_scans(runs,depvar='Ta',vary='',varx='',fit_type = 'pVoight',bkg_type='flat',peaktest=1,
                  abscor=None,plot='all',show_fits=True,mask_cmd=None,estvals=None,xrange=None,sortdep=True,
                  Nloop=10, Binit=1e-5, Tinc=2, change_factor=0.5, converge_max=100, min_change=0.01,
                  save=False,saveFIT=False):
    """ 
     Automated routine to fit peaks and plot results within runs of multiple scans dependent  
     on another variable, such as temperature, energy or azimuthal depdences.
     
         fit_scans(runs,depvar='Ta',vary='',fit_type = 'pVoight',plot='all')
     OR
         val,err = fit_scans(runs,depvar='Ta',vary='',fit_type = 'pVoight',plot=None)
     
     INPUTS:
             runs : list/array of scan numbers to use
           depvar : ('Ta')    : Independent variable, e.g. 'Ta','Energy','psi'
             vary : ('')      : Dependent variable e.g. 'APD', 'roi1_sum', 'peakroi' or '' for automatic
             varx : ('')      : Independent axis variable e.g. 'eta','chi','delta', or '' for automatic
         fit_type : ('pVoight): Type of fit to perform e.g. 'Simple','pVoight','Lorentz','Gauss'
         bkg_type : ('flat')  : Type of background to use e.g. 'flat','slope','step'
         peaktest : (10)      : Parameter to determine whether a peak exists, see help(ispeak) 
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
             save : (None)    : Save a jpg of the plot in the default directory (True/'Name.jpg'/None)
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
             val,err = load_fits(runs,depvar='Ta',fit_type='pVoight')
    
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
        >    mask_cmd = [['x>1'],['x>2'],['x>1'],...] # Remove different regions in different scans, MUST be same length as runs
    
    ESTIMATES:
        > Estimate the initial parameters of each peak
        > input "estvals" takes list of lists, where each list contains initial parameters for each scan
        > E.G. For gaussian + slope:
        >    estvals = [height,cen,wid,bkg,slope] # same for all scans
        >    estvals = [[height1,cen1,wid1,bkg1,slope1],[height2,cen2,wid2,bkg2,slope2],...] # Different for each scan
    """
    
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
            mask_cmd *= len(runs) # 2D list with single element -> 2D list same length as runs
    
    "estvals should be a 2D list of peak values estimates"
    "There should be the same number of rows as scans"
    "Each row should contain initial estimates for each of the fit parameters"
    "E.G. if gause+slope: [height,cen,wid,bkg,slope]"
    if estvals is not None:
        if len(estvals) < len(runs):
            estvals = [estvals]
        if len(estvals) == 1:
            estvals *= len(runs)
    else:
        estvals = [None]*len(runs)
    
    # Define output dictionary names
    dict_names = ['Scan Number'] + depvar + ['Peak Height','Peak Centre','FWHM','Lorz frac','Background','Area','CHI2 per dof']
    
    # Pre-allocate variables
    valstore = np.zeros([len(runs),8+Ndep])
    errstore = np.zeros([len(runs),8+Ndep])
    x_exp,y_exp = [],[]
    x_fit,y_fit = [],[]
    
    "-----Loading-----"
    for n,run in enumerate(runs):
        d = readscan(run)
        if d is None: print( 'File for run #{} does not exist!'.format(run) ); return
        x,y,dy,labvarx,labvary,ttl = getdata(d,vary=vary,abscor=abscor)[:6]
        
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
        print( '{0:3.0f}/{1} {2} {3}: Peak={4:7.3g}  Amp={5:<8.0f}  Cen={6:<6.2f}  Wid={7: <5.2g}  Frac={8:<5.2g}  Bkg={9:<8.2g}  Int={10:<8.2g}    CHI{11:8.2g}'.format(n,len(runs)-1,run,depstr,peak_rat,*valstore[n,Ndep+1:]) )
    
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
    
    fttl = '{} #{:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
    if show_fits == True:
        # Plot each scan in an axis of a 16*16 figure
        Nplots = 25
        for n in range(int(np.ceil(len(runs)/float(Nplots)))):
            fig, axs = plt.subplots(5,5,figsize=[24,14])
            fig.subplots_adjust(hspace = 0.35,wspace=0.32,left=0.07,right=0.97)
            plt.suptitle(fttl,fontsize=14)
            #fig.subplots_adjust(wspace = 0.25)
            axs = axs.flatten()
            pltruns = runs[n*Nplots:(n+1)*Nplots]
            for axn,rn in enumerate(pltruns):
                calno = axn+n*Nplots
                axs[axn].plot(x_exp[calno],y_exp[calno],'b-+',linewidth=2)
                axs[axn].plot(x_fit[calno],y_fit[calno],'r-',linewidth=2)
                axs[axn].set_title('#{}: {}={}'.format(rn,depvar[0],valstore[calno,1]))
            "---Save Multiplot figure---"
            if save not in [None, False, '']:
                if type(save) is str:
                    saveplot('{} FITS {}'.format(save,n))
                else:
                    saveplot('{0} ScansFIT {1:1.0f}-{2:1.0f} {3} FITS {4}'.format(' '.join(depvar),runs[0],runs[-1],fit_type,n))
                
    # Sort by depvar
    if sortdep:
        idx = np.argsort(valstore[:,1])
        valstore = valstore[idx,:]
        errstore = errstore[idx,:]
    
    # Save valstore & errstore values in text files
    if saveFIT not in [None, False, '']:
        header = ','.join(dict_names)
        if type(saveFIT) is str:
            savefile = os.path.join(savedir, '{}.dat'.format(saveFIT))
            esavefile = os.path.join(savedir, '{}_errors.dat'.format(saveFIT))
            np.savetxt(savefile,valstore,header=header)
            np.savetxt(esavefile,errstore,header=header)
            print( 'Saved as {}'.format(savefile) )
        else:
            savefile = os.path.join(savedir, '{0} ScansFIT {1:1.0f}-{2:1.0f} {3}.dat'.format(' '.join(depvar),runs[0],runs[-1],fit_type))
            esavefile = os.path.join(savedir, '{0} ScansFIT {1:1.0f}-{2:1.0f} {3}_errors.dat'.format(' '.join(depvar),runs[0],runs[-1],fit_type))
            np.savetxt(savefile,valstore,header=header)
            np.savetxt(esavefile,errstore,header=header)
            print( 'Saved as {}'.format(savefile) )
            print( 'Reload this scan with:\n val,err = load_fits([{},{}],depvar={},fit_type=\'{}\')'.format(runs[0],runs[-1],depvar,fit_type) )
        
    
    "------Plotting------"
    try:
        plot.__iter__; # test if array
    except AttributeError:
        plot = [plot]
    for nplot in plot:
        if nplot == 'all':
            # 2D Plots of area, width and position
            fig = plt.figure(figsize=[18,12])
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
             
            fttl = '{} #{:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.suptitle(fttl,fontsize=14)
            fig.subplots_adjust(wspace = 0.25,hspace=0.3)
            #plt.tight_layout()
            # Stop offset in x tick labels
            #x_formatter = ticker.ScalarFormatter(useOffset=False)
            #ax3.xaxis.set_major_formatter(x_formatter)
        elif nplot == 'int':
            # 2D Plots of area
            fig = plt.figure(figsize=[8,8])
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,6],'-o',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+6],errstore[:,Ndep+6],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Integrated Sum('+labvary+')',fontsize=18)
            fttl = '{0} #{1:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            plt.ylim([0,max(valstore[:,Ndep+6])*1.1])
            fig.subplots_adjust(left=0.15)
        elif nplot == 'cen':
            # 2D Plots of centre
            fig = plt.figure(figsize=[8,8])
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,2],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+2],errstore[:,Ndep+2],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel(labvarx+' Centre',fontsize=18)
            fttl = '{0} #{1:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        elif nplot == 'wid':
            # 2D Plots of width
            fig = plt.figure(figsize=[8,8])
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,3],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+3],errstore[:,Ndep+3],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel(labvarx+' Width',fontsize=18)
            fttl = '{0} #{1:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        
        "---Save Figure---"
        if save not in [None, False, '']:
            if type(save) is str:
                saveplot('{} {}'.format(save,nplot))
            else:
                saveplot('{0} ScansFIT {1} {2:1.0f}-{3:1.0f} {4}'.format(' '.join(depvar),nplot,runs[0],runs[-1],fit_type))
    
    " Prepare output dicts"
    out_values = {}
    out_errors = {}
    for n,name in enumerate(dict_names):
        out_values[name] = valstore[:,n]
        out_errors[name] = errstore[:,n]
    
    return out_values,out_errors

def load_fits(runs=[0],depvar='Ta',plot=None,fit_type = 'pVoight',file=None,disp=False,save=False):
    """ 
     Load previously fitted data from fit_scans(), assuming it completed succesfully and was saved.
      > The inputs are used to determine the automatic filename
      > Both the values and error files are read
      > Plots can be re-generated
      > A typical filename is: 'Ta ScansFIT 00001-00010 pVoight.dat'
     
         val,err = load_fits(runs,depvar='Ta',fit_type = 'pVoight',plot=None)
     
     INPUTS:
             runs : list/array of scan numbers used in fit_scans(), only the first and final values are required
           depvar : 'Ta'      : Independent variable, e.g. 'Ta','Energy','psi'
         fit_type : 'pVoight  : Type of fit to perform e.g. 'Simple','pVoight','Lorentz','Gauss'
             plot : None      : Plot the loaded fit data, available: [None,'all','int','cen','wid']
             file : None      : If given, this filename will be used to open the data files
             disp : False     : If True, prints the data on the console
             save : False     : If true, saves images of the plots
        
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
    
    """
    
    # Turn depvar into a list
    if type(depvar) is str:
        depvar = [depvar]
    Ndep = len(depvar)
    
    if type(runs) is str:
        file = os.path.join(savedir,runs+'.dat')
        efile = os.path.join(savedir,runs+'_error.dat')
    
    if file is None:
        file = os.path.join(savedir, '{0} ScansFIT {1:1.0f}-{2:1.0f} {3}.dat'.format(' '.join(depvar),runs[0],runs[-1],fit_type))
        efile = os.path.join(savedir, '{0} ScansFIT {1:1.0f}-{2:1.0f} {3}_errors.dat'.format(' '.join(depvar),runs[0],runs[-1],fit_type))
    
    valstore = np.loadtxt(file)
    errstore = np.loadtxt(efile)
    with open(file,'r') as ff:
        first_line = ff.readline().strip()
    
    first_line = first_line.strip('#')
    names = first_line.split(',')
    #depvar = names[0]
    print( names )
    
    "-----Printing-----"
    if disp:
        for n in range(len(valstore)):
            print( '{0:3.0f}/{1} {2} {3}: Amp={5:<8.0f}  Cen={6:<6.2f}  Wid={7: <5.2g}  Frac={8:<5.2g}  Bkg={9:<8.2g}  Int={10:<8.2g}    CHI{11:8.2g}'.format(n,len(runs)-1,*valstore[n,Ndep+1:]) )
    
    "------Plotting------"
    try:
        plot.__iter__; # test if array
    except AttributeError:
        plot = [plot]
    for nplot in plot:
        if nplot == 'all':
            # 2D Plots of area, width and position
            fig = plt.figure(figsize=[18,12])
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
             
            fttl = '{} #{:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.suptitle(fttl,fontsize=14)
            fig.subplots_adjust(wspace = 0.25,hspace=0.3)
            #plt.tight_layout()
            # Stop offset in x tick labels
            #x_formatter = ticker.ScalarFormatter(useOffset=False)
            #ax3.xaxis.set_major_formatter(x_formatter)
        elif nplot == 'int':
            # 2D Plots of area
            fig = plt.figure(figsize=[8,8])
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,6],'-o',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+6],errstore[:,Ndep+6],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel('Integrated Sum('+labvary+')',fontsize=18)
            fttl = '{0} #{1:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            plt.ylim([0,max(valstore[:,Ndep+6])*1.1])
            fig.subplots_adjust(left=0.15)
        elif nplot == 'cen':
            # 2D Plots of centre
            fig = plt.figure(figsize=[8,8])
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,2],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+2],errstore[:,Ndep+2],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel(labvarx+' Centre',fontsize=18)
            fttl = '{0} #{1:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        elif nplot == 'wid':
            # 2D Plots of width
            fig = plt.figure(figsize=[8,8])
            fig.add_subplot(111) # Area
            #plt.plot(valstore[:,0],valstore[:,3],'-+',linewidth=2)
            plt.errorbar(valstore[:,1],valstore[:,Ndep+3],errstore[:,Ndep+3],fmt='-o',linewidth=2)
            plt.xlabel(depvar[0],fontsize=18)
            plt.ylabel(labvarx+' Width',fontsize=18)
            fttl = '{0} #{1:1.0f}-'.format(fit_type,runs[0])+ttl+'\n'+d.metadata.cmd
            plt.title(fttl,fontsize=14)
            plt.xlim(xrange)
            fig.subplots_adjust(left=0.15)
        
        "---Save Figure---"
        if save not in [None, False, '']:
            if type(save) is str:
                saveplot('{} {}'.format(save,nplot))
            else:
                saveplot('{0} ScansFIT {1} {2:1.0f}-{3:1.0f} {4}'.format(' '.join(depvar),nplot,runs[0],runs[-1],fit_type))
    
    " Prepare output dicts"
    out_values = {}
    out_errors = {}
    for n,name in enumerate(names):
        out_values[name] = valstore[:,n]
        out_errors[name] = errstore[:,n]
    
    return out_values,out_errors

def pil_peaks(runs,depvar='Ta',ROIsize=[31,31],cax=None,save=False):
    "Show pilatus peak search results"
    
    # Load data
    x,y,z,varx,vary,varz,ttl = joindata(runs,vary=depvar)
    
    # Plot each scan in an axis of a 16*16 figure
    Nplots = 25
    ROIcenvals=np.zeros([len(runs),2])
    for n in range(int(np.ceil(len(runs)/float(Nplots)))):
        fig, axs = plt.subplots(5,5,figsize=[24,14])
        fig.subplots_adjust(hspace = 0.35,wspace=0.32,left=0.07,right=0.97)
        plt.suptitle(ttl,fontsize=14)
        #fig.subplots_adjust(wspace = 0.25)
        axs = axs.flatten()
        pltruns = runs[n*Nplots:(n+1)*Nplots]
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
                saveplot('{} PilatusPeaks {:1.0f}-{:1.0f} {}'.format(depvar,runs[0],runs[-1],n))
        
    plt.figure(figsize=[10,6])
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
            saveplot('{} PilatusPeaks {:1.0f}-{:1.0f} i-j'.format(depvar,runs[0],runs[-1]))


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
    
    Sum of several runs: plotscan([#1,#2],sum=True)
    
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
        if d is None: print( 'File for run #{} does not exist!'.format(num) ); return
    except:
        d = num
        num = d.metadata.SRSRUN
    
    " Get metadata"
    m = d.metadata
    cmd = m.cmd # Scan command
    sampsl = '{0:4.2g}x{1:<4.2g}'.format(m.s5xgap,m.s5ygap)
    detsl = '{0:4.2g}x{1:<4.2g}'.format(m.s6xgap,m.s6ygap)
    atten1 = '{0:1.0f}'.format(m.Atten)
    if labels is None:
        lbl = str(m.SRSRUN)
        if len(varys)>0: lbl=vary
    else:
        if labels in m.keys():
            lbl = '{}, {} = {}'.format(str(m.SRSRUN),labels,getattr(m,labels))
        elif labels in d.keys():
            lbl = '{}, {} = {}'.format(str(m.SRSRUN),labels,np.mean(getattr(d,lables)))
    
    " Get data" 
    x,y,dy,varx,varynew,ttl = getdata(d,vary=vary,varx=varx,norm=norm)[:6]
    
    "---Sum runs---"
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
    fig = plt.figure(figsize=[10,8])
    
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
                if labels in d.keys():
                    lab_val = np.mean(getattr(d,labels))
                elif labels in mn.keys():
                    lab_val = getattr(mn,labels)
                else:
                    labval = '?'
                lbl = '{}, {} = {}'.format(str(mn.SRSRUN),labels,lab_val)
            
            " Subtract ROIs"
            if len(varys) > 0 and subtract==True:
                dy2 = dy2**2 # add errors in quadrature
                for nvary in varys:
                    # Get data 
                    x3,y3,dy3 = getdata(dn,vary=nvary,varx=varx)[:3]
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
            
        plt.legend(loc='best')
        #plt.legend(loc='upper left')
        #plt.legend(loc='upper right')
        
    "---Plot multiple vary values---"
    if len(varys) > 0 and subtract==False:
        for n,nvary in enumerate(varys):
            " Get data" 
            x2,y2,dy2,varx2,vary2,ttl2,dn = getdata(num,vary=nvary,varx=varx)
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
        
        plt.legend(loc='best')
        #plt.legend(loc='upper left')
        #plt.legend(loc='upper right')
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot(ttl)
    return

def plotpil(num,cax=None,varx='',imnum=None,bkg_scan=None,ROIcen=None,ROIsize=[75,67],save=False):
    "Pilatus image viewer, plotpil(#)"
    
    " Load data file"
    d = readscan(num)
    if d is None: print( 'File for run #{} does not exist!'.format(num) ); return
    
    " Get data"
    x,y,dy,varx,varynew,ttl = getdata(d,varx=varx)[:6]
    
    " Load pilatus images as volume"
    vol = getvol(num)
    
    " Get scan information"
    Nframe = len(d.path)
    cmd = d.metadata.cmd
    
    " Subtract one frame from all the others for background subtraction"
    if bkg_scan is not None:
        bkgcut = vol[:,:,bkg_scan].copy()
        for n in range(len(xvals)):
            vol[:,:,n] = vol[:,:,n] - bkgcut
    
    " Default initial frame"
    if imnum is None:
        imnum = int(Nframe//2)
    
    " Default colour thresholds"
    if cax is None:
        md = np.median(vol[:,:,imnum])
        mx = np.max(vol[:,:,imnum])
        cmax = md + 10**(0.7*np.log10(mx-md))
        if cmax <= 0: cmax = 1
        cax = [0,cmax]
        print( 'caxis set at [{0:1.3g},{1:1.3g}]'.format(cax[0],cax[1]) )
    
    " Create figure & plot 1st image"
    fig = plt.figure(figsize=[10,6])
    ax = fig.add_subplot(111)
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
    ax.set_aspect('equal')
    ax.autoscale(tight=True)
    
    " Create slider on plot"
    axsldr = plt.axes([0.15, 0.15, 0.65, 0.03], axisbg='lightgoldenrodyellow')
    
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
        imgno = round(sldr.val)
        p.set_data(vol[:,:,imgno-1])
        txt.set_text('{0} = {1}'.format(varx,x[imgno-1]))
        p.set_clim(cax)
        plt.draw()
        #fig.canvas.draw()
    sldr.on_changed(update)
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot(ttl+'_'+str(imnum))

def plotdwn(num,save=None):
    "Default plot of I16 data, plotscan(#), or plotscan(#,save=1)"
    # Load data
    d = readscan(num)
    # Get principle variable data
    cmd = d.metadata.cmd # Scan command
    # Get variables
    varx = cmd.split()[1]
    vary = cmd.split()[-1]
    if vary[0:3] == 'roi':
        vary = vary + '_sum'
    x = getattr(d,varx)
    y = getattr(d,vary)
    
    #Create plot
    np.plot.plot(x,y,title='#{0}: {1}'.format(num,d.metadata.cmd),name = '#{0}'.format(num))
    return

def plotscans(scans=[],depvar=None,vary='',varx='',fit=None,norm=True,logplot=False,save=False):
    """
    Plot multiple scans of I16 data, plotscans([#1,#2,...])
    
    Plotting with fits: plotscans([#1,#2,...],fit='Gauss'), fit can also be 'Lorentz' or 'pVoight'
    
    """
    
    "---Create plot---"
    fig = plt.figure(figsize=[10,8])
    
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
        cmd = m.cmd # Scan command
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
        plotcol = plot_colors[(n)-(n)//len(plot_colors)*len(plot_colors)]
        plt.errorbar(x,y,dy,fmt='-o',c=plotcol,linewidth=2,label=lbl)
        
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
            ax.plot(xfit,yfit,':',linewidth=2,c=plotcol,label='#{} Fit'.format(str(m.SRSRUN)))
            
            " Print results on plot"
#             t0 = '   ---- {} Fit ----    '.format(fit)
#             t1 = ' Amp = {}'.format(stfm(amp,damp))
#             t2 = ' Cen = {}'.format(stfm(cen,dcen))
#             t3 = ' Wid = {}'.format(stfm(wid,dwid))
#             t4 = ' Bkg = {}'.format(stfm(bkg,dbkg))
#             t5 = ' Int = {}'.format(stfm(ara,dara))
#             
#             axs = ax.axis()
#             axw = axs[1]-axs[0]
#             axh = axs[3]-axs[2]
#             if cen > axs[0] + axw/2.:
#                 tx = axs[0] + axw*0.05
#             else:
#                 tx = axs[0] + axw*0.65
#             ty = min(yfit) + 0.9*(max(yfit)-min(yfit))
#             txt = t0+'\n'+t1+'\n'+t2+'\n'+t3+'\n'+t4+'\n'+t5
#             plt.text(tx,ty,txt,color=txtcol,
#              horizontalalignment = 'left',
#              verticalalignment   = 'top',
#              multialignment      = 'left')
    
    plt.xlabel(varxnew, fontsize=18)
    plt.ylabel(varynew, fontsize=18)
    plttl = ttl+'\n'+cmd+'\nss ={}, ds ={}, atten = {}'.format(sampsl,detsl,atten1)
    plttl = '#{} - #{}'.format(scans[0],scans[-1])
    plt.title(plttl, fontsize=16)
    
    if logplot: plt.gca().set_yscale(u'log')
    #if len(varys) > 0 and subtract==True: plt.legend([lgd],loc='best')
    fig.subplots_adjust(left=0.15)
    fig.subplots_adjust(top=0.85)
    
    # Change formats in x & y axes so they are nicer
    plt.gca().get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%8.3g'))
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        
    plt.legend(loc='best')
    #plt.legend(loc='upper left')
    #plt.legend(loc='upper right')
    
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot(ttl)
    return

def plotscans3D(runs,depvar='Ta',vary='',varx='',logplot=False,save=False):
    " Plot 3D eta scans of energy or temp dependence"
    
    "---Create Figure---"
    fig = plt.figure(figsize=[14,12])
    ax = fig.add_subplot(111,projection='3d')
    
    "-----Load & Plot-----"
    if type(runs) == int:
        "2D scan"
        "e.g. scan sx -1.0 1.0 0.05 sy -1.0 1.0 0.05 BeamOK pil100k 0.1 roi1 roi2"
        d = readscan(runs)
        cmds = d.metadata.cmd.split()
        if depvar not in d.keys():
            depvar = cmds[1]
        if varx not in d.keys():
            varx = cmds[5]
        x,y,dy,labvarx,labvary,ttl,d = getdata(runs,varx=varx,vary=vary)
        z,y2,dy2,labvarx2,labvary2,ttl2,d2 = getdata(runs,varx=depvar,vary=vary)
        ax.plot(z,x,y)
        runs = [runs]
    else:
        for n,run in enumerate(runs):
            x,y,dy,labvarx,labvary,ttl,d = getdata(runs[n],vary=vary,varx=varx)
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
    
    ttl = '#{}-'.format(runs[0])+ttl
    ax.set_title(ttl,fontsize=14)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} Scans {1:1.0f}-{2:1.0f} 3D'.format(depvar,runs[0],runs[-1]))

def plotscans2D(runs,depvar='Ta',vary='',varx='',logplot=False,save=False):
    " Plot pcolor of multiple scans"
    
    "-----Loading-----"
    x,y,z,varx,vary,varz,ttl = joindata(runs,varx,depvar,vary)
    
    # Create Plot
    fig = plt.figure(figsize=[14,12])
    ax = plt.subplot(1,1,1)
    plt.pcolor(x,y,z)
    plt.axis('tight')
    cb = plt.colorbar()
    # Axis labels
    ax.set_xlabel(varx, fontsize=18)
    ax.set_ylabel(vary, fontsize=18)
    cb.set_label(varz, fontsize=18)
    plt.suptitle(ttl,fontsize=14)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} Scans {1:1.0f}-{2:1.0f} 2D'.format(depvar,runs[0],runs[-1]))
    
def plotscansSURF(runs,depvar='Ta',vary='',varx='',logplot=False,save=False):
    " Plot surface of multiple scans"
    
    "-----Loading-----"
    x,y,z,varx,vary,varz,ttl = joindata(runs,varx,depvar,vary)
    
    # Create Plot
    fig = plt.figure(figsize=[14,12])
    ax = plt.subplot(1,1,1,projection='3d')
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
    plt.suptitle(ttl,fontsize=14)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} Scans {1:1.0f}-{2:1.0f} SURF'.format(depvar,runs[0],runs[-1]))

def plotpilSURF(num,varx='',ROIcen=None,wid=10,save=False):
    """ 
    Plot CHI scans from pilatus
    ***NOT COMPLETE***
    """
    
    x,y,dy,varx,vary,ttl,d = getdata(num,varx=varx)
    vol = getvol(num)
    if ROIcen is None:
        ROIcen,frame = pilpeak(vol,disp=True)
    
    Chivals = range(ROIcen[1]-wid//2,ROIcen[1]+wid//2)
    slice = np.sum(vol[:,Chivals,:],1)
    X,Y = np.meshgrid(x,range(vol.shape[0]))
    
    
    fig = plt.figure(figsize=[14,12])
    ax = plt.subplot(1,1,1,projection='3d')
    surf = ax.plot_surface(X,Y,slice,rstride=1,cstride=1, cmap=plt.cm.jet)
    #cmap=plt.cm.ScalarMappable(cmap=plt.cm.jet)
    #cmap.set_array(slice)
    plt.axis('tight')
    cb = plt.colorbar(surf)
    
    # Axis labels
    ax.set_xlabel(varx, fontsize=18)
    ax.set_ylabel('Pilatus vertical pixel', fontsize=18)
    ax.set_zlabel(vary, fontsize=18)
    cb.set_label(vary, fontsize=18)
    plt.suptitle(ttl,fontsize=14)
    
    "---Save---"
    if save not in [None, False, '']:
        if type(save) is str:
            saveplot(save)
        else:
            saveplot('{0} PILSURF'.format(num))


"-----------------------Peak Fitting Functions----------------------------"


def FWHM(x,y,interpolate=False):
    "Calculate a simple FWHM from a peak"
    
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

def straightline(x,grad=1.0,inter=0.0):
    "Staigh line"
    return grad*x + inter

def linefit(x,y,dy=None,disp=False):
    "Fit a line to data, y = mx + c"
    
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
    
    # Print Results
    if disp:
        print( ' ------Line Fit:----- ' )
        print( '  Gradient = {0:10.3G} +/- {1:10.3G}'.format(grad,dgrad) )
        print( ' Intercept = {0:10.3G} +/- {1:10.3G}'.format(inter,dinter) )
        print( '     CHI^2 = {0:10.3G}'.format(chi) )
        print( '  CHI^2 per free par = {0:10.3G}'.format(chinfp) )
    return grad,inter,dgrad,dinter,yfit

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
    srt = np.argsort(y)
    cen = np.average( x[ srt[ -len(x)//5: ] ] ,weights=y[ srt[ -len(x)//5: ] ])
    
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

def peakfit(x,y,dy=None,type='pVoight',bkg_type='flat',peaktest=1,estvals=None,
            Nloop=10,Binit=1e-5,Tinc=2,change_factor=0.5,converge_max = 100,
            min_change=0.01,interpolate=False,debug=False,disp=False):
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
    def create_peak_fun(text_fn,params):
        inputs = ','.join(params)
        funcstr = 'def func(x,{}):\n    return {}'.format(inputs,text_fn)
        exec funcstr in globals(),locals()
        return func
    
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
    elif type.lower() in ['simp','simple','basic','s']:
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
    elif type.lower() in ['simp','simple','basic','s']:
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
        output['Background'] = fitfunc(xold[len(xold)/2],*fitvals)
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
    chi = np.sum( (y-ycomp)**2 / dy)
    dof = len(y) - len(fitvals) # Number of degrees of freedom (Nobs-Npar)
    chinfp = chi/dof
    output['CHI**2'] = chi
    output['CHI2 per dof'] = chinfp
    
    
    # Print Results
    if disp:
        print( ' ------{} Fit:----- '.format(type) )
        for estn in range(len(fitvals)):
            print( '{0:10s} = {1:10.3G} +/- {2:10.3G}'.format(valnames[estn],fitvals[estn],errvals[estn]) )
        print( '      Area = {0:10.3G} +/- {1:10.3G}'.format(ara,dara) )
        print( '     CHI^2 = {0:10.8G}'.format(chi) )
        print( '  CHI^2 per free par = {0:10.3G}'.format(chinfp) )
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

def gauss2D((X,Y),height=1,cen_x=0,cen_y=0,FWHM_x=.5,FWHM_y=.5,bkg=0):
    "Define 2D Gaussian"
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

def abscor(eta=0,chi=90,delta=0,mu=0,gamma=0,u=1.0,disp=False,plot=False):
    """
    Calculate absorption correction
     A = abscor(num,u)
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
    A = np.sin(theta-phi) / (u*( np.sin(theta-phi) + np.sin(theta+phi) )) # IT-C Table 6.3.3.1
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
        fig = plt.figure(figsize=[12,14])
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

def labels(ttl=None,xvar=None,yvar=None,zvar=None):
    " Add good labels to current plot "
    " labels(title,xlabel,ylabel,zlabel)"
    
    if ttl != None:
        plt.gca().set_title(ttl,fontsize=20,fontweight='bold')
    
    if xvar != None:
        plt.gca().set_xlabel(xvar,fontsize=18)
    
    if yvar != None:
        plt.gca().set_ylabel(yvar,fontsize=18)
    
    if zvar != None:
        # Don't think this works, use ax.set_zaxis
        plt.gca().set_xlabel(zvar,fontsize=18)
    return

def saveplot(name,dpi=None):
    "Saves current figure as a png in the savedir directory"
    
    if type(name) is int:
        name = str(aa)
    
    gcf = plt.gcf()
    savefile = os.path.join(savedir, '{}.png'.format(saveable(name)))
    gcf.savefig(savefile,dpi=dpi)
    print( 'Saved Figure {} as {}.png'.format(gcf.number,savefile) )

def frange(start,stop=None,step=1):
    "Like np.arange but ends at stop, rather than stop-step"
    " A = frange(0,5,1) = [0.,1.,2.,3.,4.,5.]"
    
    if stop is None:
        stop = start
        start = 0
    
    return list(np.arange(start,stop+step,step,dtype=float))

def stfm(val,err):
    """
    Create standard form string from value and uncertainty"
     str = stfm(val,err)
     Examples:
          '35.25 (1)' = stfm(35.25,0.01)
          '110 (5)' = stfm(110.25,5)
          '0.0015300 (5)' = stfm(0.00153,0.0000005)
          '1.56(2)E+6' = stfm(1.5632e6,1.53e4)
    """
    if np.log10(np.abs(err)) >= 0.:
        sigfig = np.ceil(np.log10(np.abs(err)))-1
        dec = 0.
    elif err == 0.:
        return '{} (0)'.format(val)
    elif np.isnan(err):
        return '{} (-)'.format(val)
    else:
        sigfig = np.floor(np.log10(np.abs(err))+0.025)
        dec = -sigfig
    
    rval = round(val/(10.**sigfig))*(10.**sigfig)
    rerr = round(err/(10.**sigfig))*(10.**sigfig)
    if np.log10(rval) > 0:
        ln = np.ceil(np.log10(rval))
    else:
        ln = sigfig
    
    if np.log10(np.abs(err)) < 0:
        rerr = err/(10.**sigfig)
        ln = ln + dec + 1
    
    # Large numbers - exponential notation
    if np.log10(np.abs(rval)) > 4.:
        pw = np.floor(np.log10(np.abs(rval)))
        rval = rval/(10.**pw)
        rerr = rerr/(10.**sigfig)
        fmt = '{0:1.2f}({1:1.0f})E{2:+1.0f}'
        return fmt.format(rval,rerr,pw)
    
    fmt = '{'+'0:{0:1.0f}.{1:1.0f}f'.format(ln,dec+0)+'} ({1:1.0f})'
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

def findranges(scannos,sep=':'):
    "Convert a list of numbers to a simple string"
    
    scannos = np.sort(scannos).astype(int)
    
    dif = np.diff(scannos)
    
    stt,stp = [scannos[0]],[dif[0]]
    for n in range(1,len(dif)):
        if scannos[n+1] != scannos[n]+dif[n-1]:
            stt += [scannos[n]]
            stp += [dif[n]]
    stt += [scannos[-1]]
    
    out = []
    for x in range(0,len(stt),2):
        if stp[x] == 1:
            out += ['{}{}{}'.format(stt[x],sep,stt[x+1])]
        else:
            out += ['{}{}{}{}{}'.format(stt[x],sep,stp[x],sep,stt[x+1])]
    return ','.join(out)

def numbers2string(scannos,sep=':'):
    "Convert a list of numbers to a simple string"
    
    if type(scannos) is str or type(scannos) is int or len(scannos) == 1:
        return str(scannos)
    
    scannos = np.sort(scannos).astype(str)
    
    n = len(scannos[0])
    while np.all([scannos[0][:-n] == x[:-n] for x in scannos]): 
        n -= 1
    
    if n == len(scannos[0]):
        return '{}-{}'.format(scannos[0],scannos[-1])
    
    inistr = scannos[0][:-(n+1)]
    strc = [i[-(n+1):] for i in scannos]
    liststr = findranges(strc,sep=sep)
    return '{}[{}]'.format(inistr,liststr)

def maskvals(x,y,dy,mask_cmd):
    "Returns masked arrays"
    mask = np.ones(len(x),dtype=bool)
    for m in mask_cmd:
        mask[ eval(m) ] = False
    return x[mask],y[mask],dy[mask]

class dict2obj(dict):
    "Convert dictionary object to class instance"
    def __init__(self,dictvals):
        self.__dict__.update(**dictvals)
        self.update(**dictvals)