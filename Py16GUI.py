# -*- coding: utf-8 -*-
"""
Selection of graphical user interfaces (GUIs) that allow fast viewing and analysis of data 
on the beamline I16 at Diamond Light Source Ltd.

By Dan Porter, PhD
Diamond
2016

Usage: 
On an I16 workstation (as I16user)
    - Double click the shortcut on the desktop "I16_Data_Viewer"
    - Select "Run in Console"

On another Diamond workstation, the Diamond NX server or when remotely connecting:
    - Open a terminal
    - Type:
        >> cd /dls_sw/i16/software/python/Py16/
        >> module load python/ana
        >> ipython -i --matplotlib tk Py16GUI.py

On another system
    - Open a terminal/ console/ command prompt
    Type:
        >> cd /direcotry of Py16GUI.py
        >> ipython -i --matplotlib tk Py16GUI.py

Operation:
1. On running the program, the main GUI (I16_Data_Viewer) will appear
2. Browse for the experimental and analysis directories 
    - the experimental directory is where the scan data is stored
    - the analysis directory is where you wish to store any data you save or scripts you create
3. Click "Last" to load the latest data with automaticaly choosen axes

More Operations:
    - Click "Update Pilatus Plot" to load the images from any area detector
    - Click "Last Scans" to see a list of the last 100 scans that are selectable to plot
    - Select a function from the "Fit" dropdown menu to apply a fit the current plot
    - Click "Multiplot/ Peak Analysis" to open the I16_Peak_Analysis GUI
    - Click "Export Plot" to generate an figure of the current scan
    - Click "Print Current Scan" to send an image of the current scan to the default printer
    - Click "Print All Figures" to send all open figures to a single 2x3 page
    - Click "Meta" to see a full list of metadata for this scan
    - Tick the box "Live Mode" to automatically update the GUI every 10s
    - Click "Params" to change default parameters, such as error bars and temperature sensor

Console Commands:
If run in an interactive terminal, there is access to more functions via the command line, using the module "Py16progs.py"
It is already instantiated/ imported in the current session and is named "pp". For example:
    pp.plotscan(123456) # plots the scan number 123456
    d = pp.readscan(123456) # creates a data holder object 'd' with scan data
    x,y,dy,varx,vary,ttl,d = pp.getdata(123456) # returns x/y data with automatic choice of variables
    help(pp.getdata) # see the documentation on this function

Scripting:
Simple scipts can be automatically written to the analysis directory from the "MultiPlot/ Pleak Analysis" window
Pressing "Create File" will write a script including all necessary imports and parameters.
The script can the be altered for more complex analysis. Calls to the many functions in Py16Progs.py can be
made using the imported name "dp", for example:
    dp.plotscan(123456) # plots the scan number 123456
    d = dp.readscan(123456) # creates a data holder object 'd' with scan data
    x,y,dy,varx,vary,ttl,d = dp.getdata(123456) # returns x/y data with automatic choice of variables

Custom Regions of Interest:
It is possible to define new regions of interest on an area detector and plot these.
Note that manual regions of interest have hot and broken pixels set to zero. 
    - Select the scan
    - Click "Update Pilatus Plot"
    - Edit the "Centre" and "ROI" boxes to define the box centre and size, press enter to update the plot
    - The button "Find Peak" will find the largest pixel, the buttons "roi2" and "roi1" will generate standard ROIs
    - In the "Y" dropdown menu, choose "Custom ROI" to plot the sum of this ROI
    - Or: Choose "ROI - bkg" to plot the background subtracted sum (the background is defined by a region twice the ROIs size)

**********************
Main GUIs:
I16_Data_Viewer - Starts automatically, set the experiment directory to see scan details, plot scan and area detector data, access other GUIs
I16_Peak_Analysis - Plot and analyse multiple scans, including peak fitting and integration
I16_Advanced_Fitting - More fitting options, including masks
colour_cutoffs - A separate GUI that will interactively change the colormap max/min of the current figure.

Version 4.2
Last updated: 21/02/19

Version History:
07/02/16 0.9    Program created
29/02/16 1.0    I16_Data_Viewer and I16_Peak_Analysis Finished
03/03/16 1.1    Tested on Beamline, removed errors, cleaned up
28/04/16 1.2    Added Find Peak button
05/05/16 1.3    Added helper bar, made pilatus button clearer
09/05/16 1.4    Added buttons for find files and meta
10/05/16 1.5    Added print buffer + close all buttons
20/05/16 1.6    Added fit function to fit option menu and fixed reploting issue
12/07/16 1.7    Added Advanced Peak Fitting GUI
10/08/16 1.8    Added logplot and diffplot options
08/09/16 1.9    Added auto pilatus update checkbox
24/09/16 2.0    Removed requirement for SciSoftPi data loader
07/10/16 2.1    Some minor corrections, addition of parameters window
17/10/16 2.2    Addition of rem. Bkg for pilatus and pilatus peakregion/ background lines
14/12/16 2.3    New option menus for exp directories and X,Y variables, other minor improvements
20/12/16 2.4    Main app now resizes for screensize, Mac option added. Fixes for option menus. New Check buttons
08/02/17 2.5    Log of pilatus images added, other bugs fixed
25/02/17 2.6    Minor corrections and fixes, including multi-variable advanced fitting and persistence of custom ROIs
11/07/17 2.7    Added scan selector
24/07/17 2.8    Added colour_cutoffs and other bug fixes
01/08/17 2.9    Added check for large pilatus arrays, parameter "max array" in parameters window
02/10/17 2.9    Added ability to turn off normalisation in multi-plots
06/10/17 3.0    Added log plot to multiplots, plus other fixes
10/10/17 3.1    Added I16_Meta_Display, multiplotting from scan selector
23/10/17 3.2    Several minor improvements, including choice of plots for fitting
01/12/17 3.3    Update for matplotlib V2.1, reduce figure dpi, fix scaling of detector images
06/12/17 3.4    Added buttons for pixel2hkl and pixel2tth, plus other fixes
26/02/18 3.5    Added Custom Details plus improvements to image plotting, more fit functions
09/03/18 3.6    Updated to correct for python3.6 test errors
18/03/18 3.7    Default savedir directory, save Py16 parameter files to savedir directory
01/05/18 3.8    Updated for new PA + pilatus3
16/07/18 3.9    Upgraded Print_Buffer, added new plot options to multi-plot, added multi_plot input
20/11/18 4.0    Added plot toolbar, fast buttons for pilatus, external text viewer for log and checkscan
14/12/18 4.1    Some bug fixes, added windows printing
21/02/19 4.2    Some bug fixes, save scan corrected for custom rois, improved More Check Options

###FEEDBACK### Please submit your bug reports, feature requests or queries to: dan.porter@diamond.ac.uk

@author: Dan Porter
I16, Diamond Light Source
2016
"""

"""
Future Ideas:
- Convert to Qt
- Add hover boxes over buttons
- During live mode, show last line of log file in info bar
- Advanced plot allow different estimates
- Advanced plot multi peak fitting
"""

import sys,os,datetime,time,subprocess,tempfile,glob,re
import __main__ as main
import numpy as np
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox
else:
    import tkinter as tk
    from tkinter import filedialog
    from tkinter import messagebox
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.colors import Normalize, LogNorm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

"""
# Import scisoftpy - dnp.io.load is used to read #.dat files
try:
    import scisoftpy as dnp # Make sure this is in your python path
except ImportError:
    # Find scisoftpy
    print 'Hang on... just looking for scisoftpy...'
    for dirName, subdirList, fileList in os.walk(os.path.abspath(os.sep)):
        if 'uk.ac.diamond.scisoft.python' in dirName:
            print 'Found scisoftpy at: ',dirName
            sys.path.insert(0,dirName)
            break
    import scisoftpy as dnp
"""
# Import Py16progs - interprets data loaded and other handy functions
if os.path.dirname(__file__) not in sys.path:
    print('Adding to path: ''{}'''.format(os.path.dirname(__file__)) )
    sys.path.insert(0,os.path.dirname(__file__))
import Py16progs as pp

# Version
Py16GUI_Version = 4.2

# Print layout
default_print_layout = [3,2]
print_platform = 'windows' # different print commands for different os's

# App Fonts
BF= ["Times", 12]
SF= ["Times New Roman", 14]
LF= ["Times", 14]
HF= ['Courier',12]
# App Figure Sizes
WINDOWS = {'scan':[6,4], 'pilatus': [6,2.50]}
LINUX = {'scan':[6,4], 'pilatus': [6,2.90]}
MAC = {'scan':[5,3], 'pilatus': [5,2]}
# Window size needs to be different on Linux
if 'linux' in sys.platform:
    print('Linux Distribution detected - adjusting figure size accordingly')
    NORMAL = LINUX
    print_platform = 'unix'
elif 'darwin' in sys.platform:
    print('Mac OS detected - adjusting figure size accordingly')
    NORMAL = MAC
    print_platform = 'unix'
else:
    NORMAL = WINDOWS
    import win32api # requried for windows printing

"------------------------------------------------------------------------"
"---------------------------I16_Data_Viewer------------------------------"
"------------------------------------------------------------------------"
class I16_Data_Viewer():
    """
    Main GUI to view and analyse I16 data quickly during or after an experiment.
    
    OPERATION:
     - Once started, enter your experiment data directory, or press "Browse"
    
     - Enter the Scan number you wish to look at and press ENTER, or use "Last" to get the latest scan.
    
     - To plot the data, select the fitting and normalisation options you would like, then press "Plot"
    
     - To see a pilatus image, select the intensity-cutoffs and region of interest (ROI) and press "Pilatus"
    
     - If you would like to send the data to the console, press "Send to Console", this will send the variables
       x,y,dy,xvar,yvar,ttl,d to the console, where d is the dataholder with all the raw data.
        
     - Py16progs is imported to the console by default as pp, so you can use these functions easily. e.g. pp.plotscan(0)
       For more info on Py16progs and to see what functions are available, type: help(pp)
       
     - The buttons "Export Plot/Pilatus" will create figures of the current scan that can be saved or printed.
     
     - The button "Multiplot/ Peak Analysis" will take you to the Peak Analysis GUI.
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,figsize=NORMAL):
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('I16 Data Viewer [V{}] by D G Porter [dan.porter@diamond.ac.uk]'.format(Py16GUI_Version))
        self.root.minsize(width=640, height=480)
        self.root.maxsize(width=1920, height=1200)
        #self.root.maxsize(width = self.root.winfo_screenwidth(), height = self.root.winfo_screenheight())
        #print self.root.winfo_screenwidth(), self.root.winfo_screenmmwidth()
        #print self.root.winfo_screenheight(), self.root.winfo_screenmmheight()
        
        # Get initial parameters
        initial_dir = pp.filedir
        initial_num = '000000'
        initial_pilcen = pp.pil_centre
        initial_ROI = [75,67]
        initial_INT = [0,10]
        
        # Update default save location for exported plots
        plt.rcParams["savefig.directory"] = pp.savedir
        
        self.pilatus_scan = 0
        self.extra_scannos = [] # scan numbers for multiplots
        
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,anchor=tk.N)
        
        "------------------------------Help Window------------------------------"
        frm_help = tk.Frame(frame,background='white',borderwidth=1)
        frm_help.pack(fill=tk.X)
        
        self.helper = tk.StringVar(frm_help,'Welcome to I16 Data Viewer, start by Browsing for your data directory.')
        lbl_help = tk.Label(frm_help,textvariable=self.helper,font=HF,background='white')
        lbl_help.pack(side=tk.LEFT,fill=tk.X,pady=5)

        "----------------------------Data Directory-----------------------------"
        # Data Folder
        frm_fldr = tk.Frame(frame)
        frm_fldr.pack(fill=tk.X)
        self.filedir = tk.StringVar(frm_fldr,initial_dir)
        lbl_fldr = tk.Label(frm_fldr, text='Data Folder: ',width=15,font=SF)
        lbl_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        ety_fldr = tk.Entry(frm_fldr, textvariable=self.filedir, width=50)
        ety_fldr.pack(side=tk.LEFT,padx=0,pady=5)
        
        exp_list = [pp.filedir]+pp.exp_list_get()
        opt_fldr = tk.OptionMenu(frm_fldr, self.filedir, *exp_list)
        opt_fldr.config(width=1,height=1,fg=opt_fldr['menu']['bg'],activeforeground=opt_fldr['menu']['bg'])
        opt_fldr.pack(side=tk.LEFT,padx=0,pady=5)
        
        btn_fldr = tk.Button(frm_fldr, text='Browse',font=BF, command=self.f_fldr_browse)
        btn_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        btn2_fldr = tk.Button(frm_fldr, text='Search',font=BF, command=self.f_fldr_find)
        btn2_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Peak Analysis
        btn_anal = tk.Button(frm_fldr, text='Peak Analysis',font=BF, command=self.f_anal)
        btn_anal.pack(side=tk.RIGHT,padx=2)
        
        # Scripting
        btn_para = tk.Button(frm_fldr, text='Script',font=BF, command=self.f_script)
        btn_para.pack(side=tk.RIGHT,padx=2)

        # Parameters
        btn_para = tk.Button(frm_fldr, text='Params',font=BF, command=self.f_params)
        btn_para.pack(side=tk.RIGHT,padx=2)

        # Help
        btn_para = tk.Button(frm_fldr, text='Help',font=BF, command=self.f_help)
        btn_para.pack(side=tk.RIGHT,padx=2)
        
        # Analysis Folder
        frm_fldr = tk.Frame(frame)
        frm_fldr.pack(fill=tk.X)
        self.savedir = tk.StringVar(frm_fldr,pp.savedir)
        lbl_fldr = tk.Label(frm_fldr, text='Analysis Folder: ',width=15,font=SF)
        lbl_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        ety_fldr = tk.Entry(frm_fldr, textvariable=self.savedir, width=50)
        ety_fldr.pack(side=tk.LEFT,padx=0,pady=5)
        
        sav_list = [pp.savedir]+pp.sav_list_get()
        opt_fldr = tk.OptionMenu(frm_fldr, self.savedir, *sav_list)
        opt_fldr.config(width=1,height=1,fg=opt_fldr['menu']['bg'],activeforeground=opt_fldr['menu']['bg'])
        opt_fldr.pack(side=tk.LEFT,padx=0,pady=5)
        
        btn_fldr = tk.Button(frm_fldr, text='Browse',font=BF, command=self.f_fldr2_browse)
        btn_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Live Mode
        self.livemode = tk.IntVar(frm_fldr,0)
        chk_live = tk.Checkbutton(frm_fldr, text='Live Mode',font=BF, command=self.f_livemode, \
                                  variable=self.livemode,onvalue = 1, offvalue = 0)
        chk_live.pack(side=tk.RIGHT,padx=0)
        
        # Differentiate Plot
        self.diffplot = tk.IntVar(frm_fldr,0)
        chk_diff = tk.Checkbutton(frm_fldr, text='Differentiate',font=BF,variable=self.diffplot, command=self.update_plot)
        chk_diff.pack(side=tk.RIGHT,padx=0)
        
        # Log Plot
        self.logplot = tk.IntVar(frm_fldr,0)
        chk_log = tk.Checkbutton(frm_fldr, text='Log',font=BF,variable=self.logplot, command=self.f_fldr2_log)
        chk_log.pack(side=tk.RIGHT,padx=0)
        
        # Auto Pilatus Plot
        self.autopilplot = tk.IntVar(frm_fldr,0)
        chk_log = tk.Checkbutton(frm_fldr, text='Pilatus',font=BF,variable=self.autopilplot, command=self.update_pilatus)
        chk_log.pack(side=tk.RIGHT,padx=0)
        
        
        "----------------------------Scan Number-----------------------------"
        # Scan number
        frm_scan = tk.Frame(frame)
        frm_scan.pack(fill=tk.X)
        self.scanno = tk.IntVar(frm_scan,initial_num)
        lbl_scan = tk.Label(frm_scan, text='Scan No: ',font=SF)
        lbl_scan.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scan = tk.Entry(frm_scan, textvariable=self.scanno, width=10)
        ety_scan.bind('<Return>',self.update)
        ety_scan.bind('<KP_Enter>',self.update)
        ety_scan.pack(side=tk.LEFT,padx=1,pady=5)
        btn_scan_dn = tk.Button(frm_scan, text='<', font=BF, command=self.f_scan_dn)
        btn_scan_dn.pack(side=tk.LEFT)
        btn_scan_up = tk.Button(frm_scan, text='>', font=BF, command=self.f_scan_up)
        btn_scan_up.pack(side=tk.LEFT,padx=1)
        btn_scan_st = tk.Button(frm_scan, text='Last',font=BF, command=self.f_scan_st)
        btn_scan_st.pack(side=tk.LEFT,padx=1)
        btn_scan_ld = tk.Button(frm_scan, text='Load',font=BF, command=self.f_scan_ld)
        btn_scan_ld.pack(side=tk.LEFT)
        btn_scan_mt = tk.Button(frm_scan, text='Meta',font=BF, command=self.f_scan_mt)
        btn_scan_mt.pack(side=tk.LEFT)
        
        "----------------------------Plot Options-----------------------------"
        # Plot button packed next to scan number buttons
        frm_popt = tk.Frame(frm_scan)
        frm_popt.pack(side=tk.RIGHT)
        
        # Plot button
        btn_plot = tk.Button(frm_popt, text='Plot',font=BF, command=self.f_popt_plot)
        btn_plot.pack(side=tk.LEFT,padx=5,pady=5)
        
        # varx box
        self.varx = tk.StringVar(frm_popt,'Auto')
        lbl_varx = tk.Label(frm_popt, text='X',font=SF)
        lbl_varx.pack(side=tk.LEFT,padx=(5,2),pady=5)
        #ety_varx = tk.Entry(frm_popt, textvariable=self.varx, width=8)
        #ety_varx.bind('<Return>',self.update_plot)
        #ety_varx.bind('<KP_Enter>',self.update_plot)
        #ety_varx.pack(side=tk.LEFT,padx=(2,5),pady=5)
        xlist = ['Auto']
        self.opt_varx = tk.OptionMenu(frm_popt, self.varx, *xlist, command=self.f_popt_varx)
        self.opt_varx.config(width=8)
        self.opt_varx.pack(side=tk.LEFT,padx=(2,5),pady=5)
        
        # vary box
        self.vary = tk.StringVar(frm_popt,'Auto')
        lbl_vary = tk.Label(frm_popt, text='Y',font=SF)
        lbl_vary.pack(side=tk.LEFT,padx=(5,2),pady=5)
        #ety_vary = tk.Entry(frm_popt, textvariable=self.vary, width=20)
        #ety_vary.bind('<Return>',self.update_plot)
        #ety_vary.bind('<KP_Enter>',self.update_plot)
        #ety_vary.pack(side=tk.LEFT,padx=(2,3),pady=5)
        ylist = ['Auto']
        self.opt_vary = tk.OptionMenu(frm_popt, self.vary, *ylist,command=self.f_popt_vary)
        self.opt_vary.config(width=10)
        self.opt_vary.pack(side=tk.LEFT,padx=(2,3),pady=5)
        
        # normalise menu
        normopts = ['rc','ic1','none']
        self.normtype = tk.StringVar(frm_popt, normopts[0])
        lbl_nrm = tk.Label(frm_popt, text='Norm: ', font=SF)
        lbl_nrm.pack(side=tk.LEFT,padx=0,pady=5)
        opt_nrm = tk.OptionMenu(frm_popt, self.normtype, *normopts,command=self.f_popt_norm)
        opt_nrm.config(width=4)
        opt_nrm.pack(side=tk.LEFT,padx=0,pady=5)
        
        # fit menu
        fitopts = ['None','Gauss','Lorentz','pVoight','Max','Sum']
        self.fittype = tk.StringVar(frm_popt, fitopts[0])
        lbl_fit = tk.Label(frm_popt, text='Fit: ', font=SF)
        lbl_fit.pack(side=tk.LEFT,padx=0,pady=5)
        opt_fit = tk.OptionMenu(frm_popt, self.fittype, *fitopts,command=self.f_popt_fit)
        opt_fit.config(width=5)
        opt_fit.pack(side=tk.LEFT,padx=0,pady=5)
        
        "----------------------------Scan details-----------------------------"
        # Create frame just for long commands
        frm_cmd = tk.Frame(frame)
        frm_cmd.pack(fill=tk.X)
        
        # Create frame        
        frm_detl = tk.Frame(frame)
        frm_detl.pack(side=tk.LEFT,fill=tk.Y,anchor=tk.NW)
        
        # Initilise frame variables
        self.cmd = tk.StringVar(frm_detl,'')
        self.N = tk.StringVar(frm_detl,'')
        self.HKL = tk.StringVar(frm_detl,'')
        self.ENG = tk.StringVar(frm_detl,'')
        self.T = tk.StringVar(frm_detl,'')
        
        self.atten = tk.StringVar(frm_detl,'')
        self.trans = tk.StringVar(frm_detl,'')
        self.mm = tk.StringVar(frm_detl,'')
        self.do = tk.StringVar(frm_detl,'')
        
        self.pol = tk.StringVar(frm_detl,'tth =   thp =   pol = ')
        
        self.eta = tk.StringVar(frm_detl,'')
        self.chi = tk.StringVar(frm_detl,'')
        self.dlt = tk.StringVar(frm_detl,'')
        self.mu  = tk.StringVar(frm_detl,'')
        self.gam = tk.StringVar(frm_detl,'')
        self.azir= tk.StringVar(frm_detl,'')
        self.psi = tk.StringVar(frm_detl,'')
        self.phi = tk.StringVar(frm_detl,'')

        self.sx = tk.StringVar(frm_detl,'')
        self.sy = tk.StringVar(frm_detl,'')
        self.sz = tk.StringVar(frm_detl,'')
        self.spara = tk.StringVar(frm_detl,'')
        self.sperp = tk.StringVar(frm_detl,'')
        
        self.ss = tk.StringVar(frm_detl,'')
        self.ds = tk.StringVar(frm_detl,'')
        
        
        self.runtime = tk.StringVar(frm_detl,'')
        self.timetaken = tk.StringVar(frm_detl,'')
        
        # Write values to frame
        self.writeval(frm_cmd,'Command',self.cmd,wid=80)
        self.writeval(frm_detl,'Npoints',self.N)
        self.writeval(frm_detl,'HKL',self.HKL)
        self.writeval(frm_detl,'Energy',self.ENG)
        self.writeval(frm_detl,'Temp',self.T)
        self.writeval(frm_detl,'Atten',self.atten)
        self.writeval(frm_detl,'Minimirrors',self.mm)
        self.writeval(frm_detl,'Detector offset',self.do)
        
        # Analyser
        frm_pol = tk.Frame(frm_detl)
        frm_pol.pack(fill=tk.X, pady=(10,0))
        lbl_pol = tk.Label(frm_pol, textvariable=self.pol,
                           font=SF, width= 40)
        lbl_pol.pack(side=tk.LEFT)
        
        # Eta
        frm_eta = tk.Frame(frm_detl)
        frm_eta.pack(fill=tk.X, pady=(10,0))
        lbl_eta1 = tk.Label(frm_eta, text='eta:',
                            font=SF, width= 12, anchor=tk.E)
        lbl_eta1.pack(side=tk.LEFT)
        lbl_eta2 = tk.Label(frm_eta, textvariable=self.eta,
                            font=SF, width= 8, anchor=tk.W)
        lbl_eta2.pack(side=tk.LEFT)
        lbl_mu1 = tk.Label(frm_eta, text='mu:',
                           font=SF, width= 12, anchor=tk.E)
        lbl_mu1.pack(side=tk.LEFT)
        lbl_mu2 = tk.Label(frm_eta, textvariable=self.mu,
                           font=SF, width= 8, anchor=tk.W)
        lbl_mu2.pack(side=tk.LEFT,fill=tk.X)
        
        # Delta
        frm_del = tk.Frame(frm_detl)
        frm_del.pack(fill=tk.X)
        lbl_del1 = tk.Label(frm_del, text='delta:',
                            font=SF, width= 12, anchor=tk.E)
        lbl_del1.pack(side=tk.LEFT)
        lbl_del2 = tk.Label(frm_del, textvariable=self.dlt,
                            font=SF, width= 8, anchor=tk.W)
        lbl_del2.pack(side=tk.LEFT)
        lbl_gam1 = tk.Label(frm_del, text='gamma:',
                           font=SF, width= 12, anchor=tk.E)
        lbl_gam1.pack(side=tk.LEFT)
        lbl_gam2 = tk.Label(frm_del, textvariable=self.gam,
                           font=SF, width= 8, anchor=tk.W)
        lbl_gam2.pack(side=tk.LEFT,fill=tk.X)
        
        # Chi, Phi
        frm_chi = tk.Frame(frm_detl)
        frm_chi.pack(fill=tk.X)
        lbl_chi1 = tk.Label(frm_chi, text='chi:',
                            font=SF, width= 12, anchor=tk.E)
        lbl_chi1.pack(side=tk.LEFT)
        lbl_chi2 = tk.Label(frm_chi, textvariable=self.chi,
                            font=SF, width= 8, anchor=tk.W)
        lbl_chi2.pack(side=tk.LEFT)
        lbl_phi1 = tk.Label(frm_chi, text='phi:',
                            font=SF, width= 12, anchor=tk.E)
        lbl_phi1.pack(side=tk.LEFT)
        lbl_phi2 = tk.Label(frm_chi, textvariable=self.phi,
                            font=SF, width= 8, anchor=tk.W)
        lbl_phi2.pack(side=tk.LEFT)
        
        # Psi
        frm_psi = tk.Frame(frm_detl)
        frm_psi.pack(fill=tk.X, pady=(5,0))
        lbl_psi1 = tk.Label(frm_psi, text='psi:',
                           font=SF, width= 12, anchor=tk.E)
        lbl_psi1.pack(side=tk.LEFT)
        lbl_psi2 = tk.Label(frm_psi, textvariable=self.psi,
                           font=SF, width= 9, anchor=tk.W)
        lbl_psi2.pack(side=tk.LEFT)
        lbl_psi2 = tk.Label(frm_psi, textvariable=self.azir,
                           font=SF, width= 8, anchor=tk.W)
        lbl_psi2.pack(side=tk.LEFT,fill=tk.X)
        
        # sx,sy,sz
        frm_sx = tk.Frame(frm_detl)
        frm_sx.pack(fill=tk.X, pady=(10,0))
        lbl_sx1 = tk.Label(frm_sx, text='sx:',
                           font=SF, width= 6, anchor=tk.E)
        lbl_sx1.pack(side=tk.LEFT)
        lbl_sx2 = tk.Label(frm_sx, textvariable=self.sx,
                           font=SF, width= 7, anchor=tk.W)
        lbl_sx2.pack(side=tk.LEFT)
        lbl_sy1 = tk.Label(frm_sx, text='sy:',
                           font=SF, width= 6, anchor=tk.E)
        lbl_sy1.pack(side=tk.LEFT)
        lbl_sy2 = tk.Label(frm_sx, textvariable=self.sy,
                           font=SF, width= 7, anchor=tk.W)
        lbl_sy2.pack(side=tk.LEFT)
        lbl_sz1 = tk.Label(frm_sx, text='sz:',
                           font=SF, width= 6, anchor=tk.E)
        lbl_sz1.pack(side=tk.LEFT)
        lbl_sz2 = tk.Label(frm_sx, textvariable=self.sz,
                           font=SF, width= 7, anchor=tk.W)
        lbl_sz2.pack(side=tk.LEFT)
        # sperp, spara
        frm_sp = tk.Frame(frm_detl)
        frm_sp.pack(fill=tk.X)
        lbl_sp1 = tk.Label(frm_sp, text='sperp:',
                           font=SF, width= 6, anchor=tk.E)
        lbl_sp1.pack(side=tk.LEFT)
        lbl_sp2 = tk.Label(frm_sp, textvariable=self.sperp,
                           font=SF, width= 7, anchor=tk.W)
        lbl_sp2.pack(side=tk.LEFT)
        lbl_sr1 = tk.Label(frm_sp, text='spara:',
                           font=SF, width= 6, anchor=tk.E)
        lbl_sr1.pack(side=tk.LEFT)
        lbl_sr2 = tk.Label(frm_sp, textvariable=self.spara,
                           font=SF, width= 7, anchor=tk.W)
        lbl_sr2.pack(side=tk.LEFT)
        
        # Sample Slits, Detector Slits
        frm_ss = tk.Frame(frm_detl)
        frm_ss.pack(fill=tk.X, pady=(10,0))
        lbl_ss1 = tk.Label(frm_ss, text='Sample Slits: ',
                           font=SF, width= 16, anchor=tk.E)
        lbl_ss1.pack(side=tk.LEFT)
        lbl_ss2 = tk.Label(frm_ss, textvariable=self.ss,
                           font=SF, width= 12, anchor=tk.W)
        lbl_ss2.pack(side=tk.LEFT)
        frm_ds = tk.Frame(frm_detl)
        frm_ds.pack(fill=tk.X, pady=(0,10))
        lbl_ds1 = tk.Label(frm_ds, text='Detector Slits: ',
                           font=SF, width= 16, anchor=tk.E)
        lbl_ds1.pack(side=tk.LEFT)
        lbl_ds2 = tk.Label(frm_ds, textvariable=self.ds,
                           font=SF, width= 12, anchor=tk.W)
        lbl_ds2.pack(side=tk.LEFT)
        
        # Custom 1,2
        frm_cst1 = tk.Frame(frm_detl)
        frm_cst1.pack(fill=tk.X, pady=(0,0))
        self.custom1 = tk.StringVar(frm_cst1,'Custom')
        self.custom1_val = tk.StringVar(frm_cst1,'--')
        ylist = ['Custom']
        self.opt_cust1 = tk.Button(frm_cst1, textvariable=self.custom1, font=BF, command=self.f_detl_custom1)
        #self.opt_cust1 = tk.OptionMenu(frm_cst1, self.custom1, *ylist,command=self.f_detl_custom1)
        self.opt_cust1.config(width=15)
        self.opt_cust1.pack(side=tk.LEFT,padx=(10,3))
        lbl_cust = tk.Label(frm_cst1, textvariable=self.custom1_val,font=SF)
        lbl_cust.pack(side=tk.LEFT,padx=(5,2))
        
        frm_cst2 = tk.Frame(frm_detl)
        frm_cst2.pack(fill=tk.X, pady=(0,10))
        self.custom2 = tk.StringVar(frm_cst2,'Custom')
        self.custom2_val = tk.StringVar(frm_cst2,'--')
        ylist = ['Custom']
        self.opt_cust2 = tk.Button(frm_cst2, textvariable=self.custom2, font=BF, command=self.f_detl_custom2)
        #self.opt_cust2 = tk.OptionMenu(frm_cst2, self.custom2, *ylist,command=self.f_detl_custom2)
        self.opt_cust2.config(width=15)
        self.opt_cust2.pack(side=tk.LEFT,padx=(10,3))
        lbl_cust = tk.Label(frm_cst2, textvariable=self.custom2_val,font=SF)
        lbl_cust.pack(side=tk.LEFT,padx=(5,2))
        
        # Time Info
        self.writeval(frm_detl,'Ran on',self.runtime)
        self.writeval(frm_detl,'Time Taken',self.timetaken)
        
        "-----------------------------Check EXP buttons-----------------------"
        # Continue in frm_detl, from bottom
        # Check Log
        frm_log = tk.Frame(frm_detl)
        frm_log.pack(side=tk.BOTTOM,fill=tk.X, anchor=tk.SW)
        self.logmins = tk.IntVar(frm_detl,10)
        btn_log = tk.Button(frm_log, text='Check Log', width=10, 
                            font=BF, command=self.f_checklog) 
        btn_log.pack(side=tk.LEFT,padx=5)
        lbl_log = tk.Label(frm_log, text='Check last N mins',font=SF,width=18)
        lbl_log.pack(side=tk.LEFT,padx=5)
        ety_log = tk.Entry(frm_log, textvariable=self.logmins, width=4)
        ety_log.pack(side=tk.LEFT,padx=5)
        btn_log_dn = tk.Button(frm_log, text='<',font=BF, command=self.f_log_dn)
        btn_log_dn.pack(side=tk.LEFT)
        btn_log_up = tk.Button(frm_log, text='>',font=BF, command=self.f_log_up)
        btn_log_up.pack(side=tk.LEFT,padx=5)
        
        # Check Num
        frm_chk = tk.Frame(frm_detl)
        frm_chk.pack(side=tk.BOTTOM,fill=tk.X)
        self.chknum = tk.IntVar(frm_detl,10)
        btn_chk = tk.Button(frm_chk, text='Check Scans', width=10, 
                            font=BF, command=self.f_checknum) 
        btn_chk.pack(side=tk.LEFT,padx=5)
        lbl_chk = tk.Label(frm_chk, text='Check last N scans',font=SF,width=18)
        lbl_chk.pack(side=tk.LEFT,padx=5)
        ety_chk = tk.Entry(frm_chk, textvariable=self.chknum, width=4)
        ety_chk.pack(side=tk.LEFT,padx=5)
        btn_chk_dn = tk.Button(frm_chk, text='<',font=BF, command=self.f_chk_dn)
        btn_chk_dn.pack(side=tk.LEFT)
        btn_chk_up = tk.Button(frm_chk, text='>',font=BF, command=self.f_chk_up)
        btn_chk_up.pack(side=tk.LEFT,padx=5)
        
        # More Check options
        frm_mor = tk.Frame(frm_detl)
        frm_mor.pack(side=tk.BOTTOM,fill=tk.X)
        btn_mor = tk.Button(frm_mor, text='More Check Options',font=BF, command=self.f_chk_mor)
        btn_mor.pack(side=tk.LEFT,padx=5)
        btn_exp = tk.Button(frm_mor, text='Check Exp',font=BF, command=self.f_chk_exp)
        btn_exp.pack(side=tk.LEFT,padx=5)
        btn_lat = tk.Button(frm_mor, text='Check Latt',font=BF, command=self.f_chk_lat)
        btn_lat.pack(side=tk.LEFT,padx=5)
        
        # Scan Selector
        self.scan_selector_showval=''
        frm_ssl = tk.Frame(frm_detl,relief=tk.RAISED)
        frm_ssl.pack(side=tk.BOTTOM,fill=tk.X)
        btn_ssl = tk.Button(frm_ssl, text='Last scans window',font=BF, command=self.f_scn_sel)
        btn_ssl.pack(side=tk.LEFT,fill=tk.X,padx=5)
        self.select_all_scans = tk.IntVar(frm_ssl,0)
        chk_ssl = tk.Checkbutton(frm_ssl, text='All scans (slow)?',variable=self.select_all_scans)
        chk_ssl.pack(side=tk.LEFT)
        
        
        "----------------------------Plotting Window-----------------------------"
        # Create frame on right hand side
        frm_rgt = tk.Frame(frame)
        frm_rgt.pack(side=tk.RIGHT, fill=tk.X, anchor=tk.NE, expand = tk.YES)
        # Create frame for plot
        frm_plt = tk.Frame(frm_rgt)
        frm_plt.pack(fill=tk.X,expand=tk.YES)
        
        self.fig1 = plt.Figure(figsize=figsize['scan'],dpi=80)
        self.fig1.patch.set_facecolor('w')
        self.ax1 = self.fig1.add_subplot(111)
        self.ax1.set_autoscaley_on(True)
        self.ax1.set_autoscalex_on(True)
        self.plt1, = self.ax1.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5],'o-',c=pp.plot_colors[0],linewidth=2)
        self.plt2, = self.ax1.plot([],[],'g:') # marker point for pilatus
        self.pfit, = self.ax1.plot([],[],'-',c=pp.plot_colors[-1],linewidth=2) # fit line
        self.extra_plots = [] 
        self.ax1.set_xlabel('varx')
        self.ax1.set_ylabel('vary')
        self.ax1.set_title('Scan number',fontsize=16)
        
        self.fig1.subplots_adjust(left=0.25,bottom=0.2,right=0.95)
        
        # Change formats in x & y axes so they are nicer
        self.ax1.get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%8.3g'))
        self.ax1.get_xaxis().get_major_formatter().set_useOffset(False)
        
        canvas = FigureCanvasTkAgg(self.fig1, frm_plt)
        canvas.get_tk_widget().configure(bg='black')
        #canvas.bind("<Button-1>",figureclick) # this doesn't work - FigureCanvasTkAgg doesn't have bind
        #canvas.bind("<B1-Motion>",figurehold)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, anchor=tk.NE, expand=tk.YES)
        #canvas.get_tk_widget().pack()
        #self.update_plot()

        # Add matplotlib toolbar under plot
        self.toolbar = NavigationToolbar2TkAgg( canvas, frm_rgt ) 
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP)
        
        "----------------------------Pilatus Options-----------------------------"
        # Pilatus option buttons below plot figure
        frm_pilopt = tk.Frame(frm_rgt)
        frm_pilopt.pack()
        
        # Pilatus variables
        self.pilatus_active = False
        self.pilatus_scan = 0
        self.pilpos = 0
        self.pilstr = tk.StringVar(frm_pilopt,'(0) eta = 0')
        self.pilcen_i = tk.IntVar(frm_pilopt,initial_pilcen[0])
        self.pilcen_j = tk.IntVar(frm_pilopt,initial_pilcen[1])
        self.roisiz_i = tk.IntVar(frm_pilopt,initial_ROI[0])
        self.roisiz_j = tk.IntVar(frm_pilopt,initial_ROI[1])
        self.pilint_i = tk.DoubleVar(frm_pilopt,initial_INT[0])
        self.pilint_j = tk.DoubleVar(frm_pilopt,initial_INT[1])
        
        # Plot button
        frm_pilopt1 = tk.Frame(frm_pilopt)
        frm_pilopt1.pack(side=tk.LEFT,fill=tk.BOTH)
        
        btn_pilopt = tk.Button(frm_pilopt1, text='Update Pilatus Plot',wraplength=60,font=BF, command=self.f_pilopt_plot)
        btn_pilopt.pack(side=tk.LEFT,fill=tk.Y,padx=2)
        
        # Pilatus centre + ROI size
        frm_pilopt2 = tk.Frame(frm_pilopt)
        frm_pilopt2.pack(fill=tk.X)
        
        lbl_pilcen = tk.Label(frm_pilopt2, text='Centre:',font=SF,width=6)
        lbl_pilcen.pack(side=tk.LEFT,padx=1)
        ety_pilcen_i = tk.Entry(frm_pilopt2, textvariable=self.pilcen_i, width=4)
        ety_pilcen_i.bind('<Return>',self.update_pilatus)
        ety_pilcen_i.bind('<KP_Enter>',self.update_pilatus)
        ety_pilcen_i.pack(side=tk.LEFT,padx=1)
        ety_pilcen_j = tk.Entry(frm_pilopt2, textvariable=self.pilcen_j, width=4)
        ety_pilcen_j.bind('<Return>',self.update_pilatus)
        ety_pilcen_j.bind('<KP_Enter>',self.update_pilatus)
        ety_pilcen_j.pack(side=tk.LEFT,padx=1)
        
        lbl_roisiz = tk.Label(frm_pilopt2, text='ROI:',font=SF,width=5)
        lbl_roisiz.pack(side=tk.LEFT,padx=1)
        ety_roisiz_i = tk.Entry(frm_pilopt2, textvariable=self.roisiz_i, width=4)
        ety_roisiz_i.bind('<Return>',self.update_pilatus)
        ety_roisiz_i.bind('<KP_Enter>',self.update_pilatus)
        ety_roisiz_i.pack(side=tk.LEFT,padx=1)
        ety_roisiz_j = tk.Entry(frm_pilopt2, textvariable=self.roisiz_j, width=4)
        ety_roisiz_j.bind('<Return>',self.update_pilatus)
        ety_roisiz_j.bind('<KP_Enter>',self.update_pilatus)
        ety_roisiz_j.pack(side=tk.LEFT,padx=1)
        
        btn_peak = tk.Button(frm_pilopt2, text='Find Peak',font=BF,command=self.f_pilopt_peak)
        btn_peak.pack(side=tk.LEFT,padx=2)
        #btn_droi = tk.Button(frm_pilopt2, text='nroi',font=BF,command=self.f_pilopt_nroi)
        #btn_droi.pack(side=tk.LEFT,padx=2)
        btn_def2 = tk.Button(frm_pilopt2, text='roi2',font=BF,command=self.f_pilopt_def2)
        btn_def2.pack(side=tk.LEFT,padx=2)
        btn_def1 = tk.Button(frm_pilopt2, text='roi1',font=BF,command=self.f_pilopt_def1)
        btn_def1.pack(side=tk.LEFT,padx=2)
        
        # Pilatus Image intensity cutoffs + pos buttons
        frm_pilopt3 = tk.Frame(frm_pilopt)
        frm_pilopt3.pack(fill=tk.X)
        
        lbl_pilint = tk.Label(frm_pilopt3, text='CLim:',font=SF,width=4)
        lbl_pilint.pack(side=tk.LEFT,padx=2)
        ety_pilint_i = tk.Entry(frm_pilopt3, textvariable=self.pilint_i, width=4)
        ety_pilint_i.bind('<Return>',self.update_pilatus)
        ety_pilint_i.bind('<KP_Enter>',self.update_pilatus)
        ety_pilint_i.pack(side=tk.LEFT,padx=1)
        ety_pilint_j = tk.Entry(frm_pilopt3, textvariable=self.pilint_j, width=6)
        ety_pilint_j.bind('<Return>',self.update_pilatus)
        ety_pilint_j.bind('<KP_Enter>',self.update_pilatus)
        ety_pilint_j.pack(side=tk.LEFT,padx=1)
        
        lbl_pilpos = tk.Label(frm_pilopt3, textvariable=self.pilstr,font=SF,width=14)
        lbl_pilpos.pack(side=tk.LEFT,padx=2)
        btn_pilpos0 = tk.Button(frm_pilopt3, text='<<',font=BF, command=self.f_pilopt_posleftfast)
        btn_pilpos0.pack(side=tk.LEFT,padx=0)
        btn_pilpos1 = tk.Button(frm_pilopt3, text='<',font=BF, command=self.f_pilopt_posleft)
        btn_pilpos1.pack(side=tk.LEFT,padx=0)
        btn_pilpos2 = tk.Button(frm_pilopt3, text='>',font=BF, command=self.f_pilopt_posright)
        btn_pilpos2.pack(side=tk.LEFT,padx=0)
        btn_pilpos3 = tk.Button(frm_pilopt3, text='>>',font=BF, command=self.f_pilopt_posrightfast)
        btn_pilpos3.pack(side=tk.LEFT,padx=0)
        
        frm_pilopt3a = tk.Frame(frm_pilopt3)
        frm_pilopt3a.pack(side=tk.LEFT,fill=tk.X)
        
        # Remove background
        self.rembkg = tk.IntVar(frm_pilopt3,0)
        chk_rbkg = tk.Checkbutton(frm_pilopt3a, text='Rem. Bkg',font=('Times',8),borderwidth=0,variable=self.rembkg, command=self.f_pilopt_rembkg)
        chk_rbkg.pack(side=tk.TOP,padx=0,pady=0)
        
        self.remfrm = tk.IntVar(frm_pilopt3,0)
        chk_rfrm = tk.Checkbutton(frm_pilopt3a, text='Rem. Frm',font=('Times',8),borderwidth=0,variable=self.remfrm, command=self.f_pilopt_remfrm)
        chk_rfrm.pack(side=tk.TOP,padx=0,pady=0)
        
        "----------------------------Pilatus Window------------------------------"
        # Figure window 2, below pilatus options buttons
        frm_pil = tk.Frame(frm_rgt)
        frm_pil.pack(fill=tk.X,expand=tk.YES)
        
        self.fig2 = plt.Figure(figsize=figsize['pilatus'],dpi=80)
        self.fig2.patch.set_facecolor('w')
        self.ax2 = self.fig2.add_subplot(111)
        self.ax2.set_xticklabels([])
        self.ax2.set_yticklabels([])
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])
        self.ax2.set_autoscaley_on(True)
        self.ax2.set_autoscalex_on(True)
        self.ax2.set_frame_on(False)
        self.pilim = self.ax2.imshow(np.zeros([195,487]))
        #default_image = pp.misc.imread( os.path.dirname(__file__)+'\default_image.png')
        #self.pilim = self.ax2.imshow(default_image)
        self.ax2.set_position([0,0,1,1])
        
        #self.fig2.subplots_adjust(left=0.2,bottom=0.2)
        
        ROIcen = initial_pilcen
        ROIsize = [75,67]
        pil_centre = initial_pilcen
        pil_size = [195,487]
        idxi = np.array([ROIcen[0]-ROIsize[0]//2,ROIcen[0]+ROIsize[0]//2+1])
        idxj = np.array([ROIcen[1]-ROIsize[1]//2,ROIcen[1]+ROIsize[1]//2+1])
        self.pilp1, = self.ax2.plot(idxj[[0,1,1,0,0]],idxi[[0,0,1,1,0]],'k-',linewidth=2) # ROI
        self.pilp2, = self.ax2.plot([pil_centre[1],pil_centre[1]],[0,pil_size[0]],'k:',linewidth=2) # vertical line
        self.pilp3, = self.ax2.plot([0,pil_size[1]],[pil_centre[0],pil_centre[0]],'k:',linewidth=2) # Horizontal line
        self.pilp4, = self.ax2.plot([],[],'r-',linewidth=2) # ROI background
        self.pilp5, = self.ax2.plot([],[],'y-',linewidth=2) # Peak region
        self.ax2.set_aspect('equal')
        self.ax2.autoscale(tight=True)
        
        canvas2 = FigureCanvasTkAgg(self.fig2, frm_pil)
        canvas2.get_tk_widget().configure(bg='black')
        canvas2.draw()
        canvas2.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, anchor=tk.NE, expand=tk.YES)
        
        "-----------------------------Final Options------------------------------"
        # Buttons at bottom right
        frm_fnl = tk.Frame(frm_rgt)
        frm_fnl.pack()
        
        # Button to send data to Console
        btn_send = tk.Button(frm_fnl, text='Send to Console',font=BF, command=self.f_fnl_send)
        btn_send.pack(side=tk.LEFT,padx=2,pady=5)
        
        # Button to send data to Console
        btn_splot = tk.Button(frm_fnl, text='Export Plot',font=BF, command=self.f_fnl_splot)
        btn_splot.pack(side=tk.LEFT,padx=2,pady=5)
        
        # Button to send data to Console
        btn_spil = tk.Button(frm_fnl, text='Export Pilatus',font=BF, command=self.f_fnl_spil)
        btn_spil.pack(side=tk.LEFT,padx=2)
        
        # Button to send data to Console
        btn_spil = tk.Button(frm_fnl, text='pil2hkl',font=BF, command=self.f_fnl_phkl)
        btn_spil.pack(side=tk.LEFT,padx=2)
        
        # Button to send data to Console
        btn_spil = tk.Button(frm_fnl, text='pil2tth',font=BF, command=self.f_fnl_ptth)
        btn_spil.pack(side=tk.LEFT,padx=2)
        
        # Row 2
        frm_fnl2 = tk.Frame(frm_rgt)
        frm_fnl2.pack()
        
        # Button to save plot
        btn_send = tk.Button(frm_fnl2, text='Save Current Scan',font=BF, command=self.f_fnl_splotsave)
        btn_send.pack(side=tk.LEFT,padx=2,pady=1)
        
        # Button to print plot
        btn_send = tk.Button(frm_fnl2, text='Print Current Scan',font=BF, command=self.f_fnl_splotprint)
        btn_send.pack(side=tk.LEFT,padx=2,pady=1)
        
        # Button to print all figures
        btn_send = tk.Button(frm_fnl2, text='Print All Figures',font=BF, command=self.f_fnl_splotbuffer)
        btn_send.pack(side=tk.LEFT,padx=2,pady=1)
        
        # Button to close all figures
        btn_send = tk.Button(frm_fnl2, text='Close All',font=BF, command=self.f_fnl_splotclose)
        btn_send.pack(side=tk.LEFT,padx=2,pady=1)
        
        
        # Check widget size vs screen size
        if self.root.winfo_height() > self.root.winfo_screenheight()-80:
            print('Screen Height = ',self.root.winfo_screenheight())
            print('App Height = ',self.root.winfo_height())
            print('Oh... this is a small screen. I\'ll just make the viewer smaller...')
            # App Fonts
            BF[1] -= 2
            SF[1] -= 2
            LF[1] -= 2
            HF[1] -= 2
            # App Figure Sizes
            NORMAL['scan'] = [NORMAL['scan'][0]*0.9,NORMAL['scan'][1]*0.9]
            NORMAL['pilatus'] = [NORMAL['pilatus'][0]*0.9,NORMAL['pilatus'][1]*0.9]
            
            self.root.destroy()
            I16_Data_Viewer()
        
        # Load initial data
        if pp.latest() is not None:
            self.f_scan_st()
        if not hasattr(sys, 'ps1'):
            # If not in interactive mode, start mainloop
            self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
            self.root.mainloop()
        
    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_fldr_opt(self):
        "Load previous experiment"
        pass
    
    def f_fldr_browse(self):
        "Browse for data directory"
        inidir = self.filedir.get()
        dir = filedialog.askdirectory(initialdir=inidir)
        self.filedir.set(dir)
        self.savedir.set(dir +os.path.sep+'processing')
        #self.helper.set('Now set the analysis folder - saved images and scripts will be stored here')
    
    def f_fldr_find(self):
        "Search for data directory"
        inidir = self.filedir.get()
        num = self.scanno.get()
        
        if num > 1000:
            self.helper.set('Searching file directories for scan #{}'.format(num))
            newdir = pp.findfile(num,topdir=inidir)
            self.filedir.set(newdir)
        else:
            self.helper.set('Please enter a valid scan number so I can search for its directory')
    
    def f_fldr2_browse(self):
        "Browse for analysis directory"
        inidir = self.savedir.get()
        dir = filedialog.askdirectory(initialdir=inidir)
        self.savedir.set(dir)
        self.helper.set('Now click "Last" to load the latest scan, or enter a scan number and press Enter')
    
    def f_fldr2_log(self):
        "Activate log plot"
        
        if self.logplot.get():
            if self.pilint_i.get() < 0.1:
                self.pilint_i.set(0.1)
            self.pilim.norm = LogNorm()
        else:
            self.pilim.norm = Normalize()
        self.fig2.canvas.draw()
        self.update_plot()
    
    def f_help(self):
        "Launch help GUI"
        I16_Text_Display(__doc__)

    def f_script(self):
        "Launch python editor"
        self.set_files()
        scanno = self.scanno.get()
        outstr = pp.example_script(scanno)
        I16_Python_Editor(outstr)

    def f_params(self):
        "Launch parameters GUI"
        
        self.set_files()
        
        I16_Params()
    
    def f_anal(self):
        "Launch analysis GUI"
        
        self.set_files()
        
        setvarx = self.varx.get()
        setvary = self.vary.get()
        
        if self.remfrm.get():
            remfrm = '_sfm'
        else:
            remfrm = ''
        
        if setvarx == 'Auto':
            setvarx = ''
        if setvary == 'Auto':
            setvary = ''
        elif setvary == 'Custom ROI':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif setvary == 'ROI - bkg':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi_bkg{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        
        if len(self.extra_scannos) > 0:
            multirange = [self.scanno.get()]+list(self.extra_scannos)
        else:
            multirange=None
        
        I16_Peak_Analysis(self.scanno.get(),setvarx,setvary,multirange)
    
    def f_livemode(self):
        "Activate Live Mode"
        
        if self.livemode.get() == 1:
            self.helper.set('Live Mode activated - plot will update every 5 seconds.')
            self.scanno.set(pp.latest())
            self.update_details()
            self.update_plot()
            self.root.after(5000,self.f_livemode)
        else:
            self.helper.set('Live Mode deactivated')
    
    def f_scan_dn(self):
        "Decrease Scan number"
        num = self.scanno.get()
        self.scanno.set(num-1)
        self.update_details()
        self.update_plot()
    
    def f_scan_up(self):
        "Increase Scan number"
        num = self.scanno.get()
        last = pp.latest()
        if num < last:
            self.scanno.set(num+1)
        self.update_details()
        self.update_plot()
    
    def f_scan_st(self):
        "Latest Scan number"
        
        self.set_files()
        
        num = pp.latest()
        if num is None:
            num = '000000'
            self.helper.set('There are no scan files in this folder')
        self.scanno.set(num)
        self.update_details()
        self.update_plot()
    
    def f_scan_ld(self):
        "Update GUI for selected scan number"
        self.update_details()
    
    def f_scan_mt(self):
        "Send metadata to console"
        # Get the parameters
        self.set_files()
        scanno = self.scanno.get()
        
        self.helper.set('Displaying metadata in new window')
        I16_Meta_Display(scanno)
    
    def f_popt_plot(self):
        "Plot current scan"
        self.update_plot()
    
    def f_popt_varx(self,x):
        "Get varx and plot current scan"
        self.varx.set(x)
        self.pilatus_active = False
        self.update_plot()
    
    def f_popt_vary(self,x):
        "Get vary and plot current scan"
        self.vary.set(x)
        self.update_plot()
    
    def f_popt_norm(self,x):
        "Plot scan with fit"
        self.update_plot()
    
    def f_popt_fit(self,x):
        "Plot scan with fit"
        self.update_plot()
    
    def f_detl_custom1(self):
        "Asign value to custom detail"
        self.set_files()
        scanno = self.scanno.get()
        d = pp.readscan(scanno)

        metafields = d.metadata.keys()
        out = I16_Meta_Select(self, self.custom1, metafields)
    
    def f_detl_custom2(self):
        "Asign value to custom detail"
        self.set_files()
        scanno = self.scanno.get()
        d = pp.readscan(scanno)

        metafields = d.metadata.keys()
        out = I16_Meta_Select(self, self.custom2, metafields)
        #self.root.wait_window(out.root)
        self.update_details()
    
    def f_log_dn(self):
        "Decrease Scan number"
        num = self.logmins.get()
        self.logmins.set(num-1)
    
    def f_log_up(self):
        "Increase Scan number"
        num = self.logmins.get()
        self.logmins.set(num+1)
    
    def f_chk_dn(self):
        "Decrease Scan number"
        num = self.chknum.get()
        self.chknum.set(num-1)
    
    def f_chk_up(self):
        "Increase Scan number"
        num = self.chknum.get()
        self.chknum.set(num+1)
    
    def f_checknum(self):
        "checknum(-N,0)"
        num = self.chknum.get()
        outstr = pp.checkscans(-num,0)
        I16_Text_Display(outstr, 'Check Scans')
    
    def f_checklog(self):
        "checklog(mins=N)"
        num = self.logmins.get()
        outstr = pp.checklog(mins=num)
        I16_Text_Display(outstr, 'Check Log')
    
    def f_scn_sel(self):
        " Scan selector button launches scan selector"
        allscans = self.select_all_scans.get()
        showval = self.scan_selector_showval
        if allscans:
            scannos = pp.get_all_scannos()
        else:
            scannos = range(pp.latest()-100,pp.latest()+1)
        
        I16_Scan_Selector(self,scannos,showval)
    
    def f_chk_mor(self):
        "Check log file"
        I16_Check_Log(self.scanno.get())
    
    def f_chk_exp(self):
        "Check experiment button - opens a message box"
        
        " Code taken fro Py16progs.checkexp"
        
        self.set_files()
        filedir = pp.filedir
        
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
        fd = pp.readscan(fn)
        ld = pp.readscan(ln)
        
        # Get first and last file creation dates (not reliable)
        #ft = datetime.datetime.fromtimestamp(os.path.getctime(ls[0]))
        #lt = datetime.datetime.fromtimestamp(os.path.getctime(ls[-1]))
        
        # Use SRSDAT and SRSTIM - because date isn't in every file
        # Note that SRSDAT stores months and days without 0's, 
        # so 11 Feb (2015112) and 1 Dec (2015112) look the same! Not sure how to fix this!
        ft = datetime.datetime.strptime(str(fd.metadata.SRSDAT)+' '+str(fd.metadata.SRSTIM),'%Y%m%d %H%M%S')
        lt = datetime.datetime.strptime(str(ld.metadata.SRSDAT)+' '+str(ld.metadata.SRSTIM),'%Y%m%d %H%M%S')
        
        # Format times
        ft = ft.strftime('%Y-%m-%d   %H:%M')
        lt = lt.strftime('%Y-%m-%d   %H:%M')
        
        #if 'date' in fd.metadata.keys(): ft = fd.metadata.date
        #if 'date' in ld.metadata.keys(): lt = ld.metadata.date
        
        # Print data to messge box
        ttl = 'Experiment: '+filedir
        msg = ' Scans in directory: {}\n'.format(filedir)
        msg+= ' First scan: #{0}    {1}\n'.format(fn,ft)
        msg+= '  Last scan: #{0}    {1}\n'.format(ln,lt)
        msg+= '  No. scans: {}'.format(len(ls))
        messagebox.showinfo(ttl,msg)
    
    def f_chk_lat(self):
        "Check scan lattice parameters - opens a message box"
        
        scanno = self.scanno.get()
        
        d = pp.readscan(scanno)
        if d is None:
            return
        
        m = d.metadata
        ttl = 'Lattice Parameters #{}'.format(m.SRSRUN)
        msg = '---Lattice Parameters for scan #{}---\n'.format(m.SRSRUN)
        msg+= '    a = {:7.4f},     b = {:7.4f},      c = {:7.4f}\n'.format(m.a,m.b,m.c)
        msg+= 'alpha = {:7.2f},  beta = {:7.2f},  gamma = {:7.2f}'.format(m.alpha1,m.alpha2,m.alpha3)
        messagebox.showinfo(ttl,msg)
    
    def f_pilopt_plot(self):
        "Plot pilatus images"
        self.update_pilatus()
    
    def f_pilopt_peak(self):
        "Determine peak position in pilatus"
        
        if self.pilatus_active == False:
            self.helper.set('No pilatus image active, update the pilatus plot first')
            return
        
        self.helper.set('Find Peak: centres at the largest pixel, use "nroi" to define a new region of interest here')
        ROIsizei = self.roisiz_i.get()
        ROIsizej = self.roisiz_j.get()
        
        if self.rembkg.get() and self.remfrm.get():
            vary = 'nroi_peak_bkg_sfm[{},{}]'.format(ROIsizei,ROIsizej)
        elif self.rembkg.get():
            vary = 'nroi_peak_bkg[{},{}]'.format(ROIsizei,ROIsizej)
        elif self.remfrm.get():
            vary = 'nroi_peak_sfm[{},{}]'.format(ROIsizei,ROIsizej)
        else:
            vary = 'nroi_peak[{},{}]'.format(ROIsizei,ROIsizej)
        #self.vary.set(vary)
        
        [ROIceni,ROIcenj],frame = pp.pilpeak(self.vol,disp=True)
        
        self.pilpos = frame
        self.pilcen_i.set(str(ROIceni))
        self.pilcen_j.set(str(ROIcenj))
        
        # Show peakregion
        pry1,prx1,pry2,prx2 = pp.peakregion
        self.pilp5.set_xdata([prx1,prx2,prx2,prx1,prx1])
        self.pilp5.set_ydata([pry1,pry1,pry2,pry2,pry1])
        
        self.update_pilatus()
        self.update_plot()
    
    def f_pilopt_nroi(self):
        "Send pil values to vary"
        
        self.helper.set('nroi: Create a new Region of Interest, apply it by pressing "Plot"')
        ROIceni = self.pilcen_i.get()
        ROIcenj = self.pilcen_j.get()
        ROIsizei = self.roisiz_i.get()
        ROIsizej = self.roisiz_j.get()
        
        if self.rembkg.get() and self.remfrm.get():
            vary = 'nroi_bkg_sfm[{},{},{},{}]'.format(ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif self.rembkg.get():
            vary = 'nroi_bkg[{},{},{},{}]'.format(ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif self.remfrm.get():
            vary = 'nroi_sfm[{},{},{},{}]'.format(ROIceni,ROIcenj,ROIsizei,ROIsizej)
        else:
            vary = 'nroi[{},{},{},{}]'.format(ROIceni,ROIcenj,ROIsizei,ROIsizej)
       
        self.vary.set(vary)
    
    def f_pilopt_def2(self):
        "Reset the pilatus centre and ROI"
        
        self.pilcen_i.set(pp.pil_centre[0])
        self.pilcen_j.set(pp.pil_centre[1])
        self.roisiz_i.set(75)
        self.roisiz_j.set(67)
        self.update_pilatus()
    
    def f_pilopt_def1(self):
        "Reset the pilatus centre and ROI"
        
        self.pilcen_i.set(pp.pil_centre[0])
        self.pilcen_j.set(pp.pil_centre[1])
        self.roisiz_i.set(15)
        self.roisiz_j.set(13)
        self.update_pilatus()
    
    def f_pilopt_posleft(self):
        "Move to previous pilatus image"
        
        if self.pilatus_active == False: return
        if self.pilpos-1 < 0: return
        
        # Increment the image by -1
        self.pilpos -= 1
        imgno = self.pilpos
        xval = self.pilvals[imgno]
        self.pilstr.set('({}) {} = {:1.3f}'.format(imgno,self.pilvar,xval))
        
        # Update pilatus image
        self.pilim.set_data(self.vol[:,:,imgno])
        self.fig2.canvas.draw()
        
        # Update marker point on plot axis
        ylim = self.ax1.get_ylim()
        self.plt2.set_xdata([xval,xval])
        self.plt2.set_ydata(ylim)
        self.fig1.canvas.draw()
    
    def f_pilopt_posright(self):
        "Move to next pilatus image"
        
        if self.pilatus_active == False: return
        if self.pilpos+1 >= len(self.pilvals): return
        
        # Increment the image by +1
        self.pilpos += 1
        imgno = self.pilpos
        xval = self.pilvals[imgno]
        self.pilstr.set('({}) {} = {:1.3f}'.format(imgno,self.pilvar,xval))
        
        # Update pilatus image
        self.pilim.set_data(self.vol[:,:,imgno])
        self.fig2.canvas.draw()
        
        # Update marker point on plot axis
        ylim = self.ax1.get_ylim()
        self.plt2.set_xdata([xval,xval])
        self.plt2.set_ydata(ylim)
        self.fig1.canvas.draw()

    def f_pilopt_posleftfast(self):
        "Move back 1/10th of positions"
        
        if self.pilatus_active == False: return

        # Increment the image by -10%
        incval = len(self.pilvals)//10
        self.pilpos -= incval
        if self.pilpos < 0: self.pilpos = 0 # go to zero
        
        imgno = self.pilpos
        xval = self.pilvals[imgno]
        self.pilstr.set('({}) {} = {:1.3f}'.format(imgno,self.pilvar,xval))
        
        # Update pilatus image
        self.pilim.set_data(self.vol[:,:,imgno])
        self.fig2.canvas.draw()
        
        # Update marker point on plot axis
        ylim = self.ax1.get_ylim()
        self.plt2.set_xdata([xval,xval])
        self.plt2.set_ydata(ylim)
        self.fig1.canvas.draw()
    
    def f_pilopt_posrightfast(self):
        "Move forward 1/10th of position"
        
        if self.pilatus_active == False: return

        # Increment the image by +10%
        incval = len(self.pilvals)//10
        self.pilpos += incval
        if self.pilpos >= len(self.pilvals): self.pilpos = len(self.pilvals)-1 # go to end

        imgno = self.pilpos
        xval = self.pilvals[imgno]
        self.pilstr.set('({}) {} = {:1.3f}'.format(imgno,self.pilvar,xval))
        
        # Update pilatus image
        self.pilim.set_data(self.vol[:,:,imgno])
        self.fig2.canvas.draw()
        
        # Update marker point on plot axis
        ylim = self.ax1.get_ylim()
        self.plt2.set_xdata([xval,xval])
        self.plt2.set_ydata(ylim)
        self.fig1.canvas.draw()
    
    def f_pilopt_rembkg(self):
        "Remove background from ROI"
        
        if self.rembkg.get():
            self.helper.set('Set ROI with background removal. Backgound is average of area in RED')
            self.update_pilatus()
        else:
            self.update_pilatus()
    
    def f_pilopt_remfrm(self):
        "Remove background from ROI"
        
        if self.remfrm.get():
            self.helper.set('Remove first frame from all images in scan')
            self.pilatus_active = False
            self.update_pilatus()
        else:
            self.helper.set('Return to uncorrected images')
            self.pilatus_active = False
            self.update_pilatus()
    
    def f_fnl_send(self):
        "Send instructions to console"
        
        # Set globals
        global x,y,dy,varx,vary,ttl,d
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        setvarx = self.varx.get()
        setvary = self.vary.get()
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        # test the vary
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(scanno,vary=setvary,varx=setvarx,norm=norm)
        except AttributeError:
            setvarx = ''
            setvary = ''
        
        if self.remfrm.get():
            remfrm = '_sfm'
        else:
            remfrm = ''
        if setvary == 'Custom ROI':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif setvary == 'ROI - bkg':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi_bkg{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        
        # Send the command
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(scanno,vary=setvary,varx=setvarx,norm=norm)
            cmdstr = 'x,y,dy,varx,vary,ttl,d = pp.getdata({},varx=\'{}\',vary=\'{}\',norm={})'
            print( cmdstr.format(scanno,setvarx,setvary,norm) )
            self.helper.set('d = pp.readscan({})'.format(scanno) + ' has been sent to the console')
        except:
            print( 'd = pp.readscan({})'.format(scanno) )
            d = pp.readscan(scanno)
            self.helper.set('d = pp.readscan({})'.format(scanno) + ' has been sent to the console')
    
    def f_fnl_splot(self):
        "Send plotscan command to console"
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        logplot = self.logplot.get()
        diffplot = self.diffplot.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        # test the vary
        setvarx = self.varx.get()
        setvary = self.vary.get()
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(scanno,varx=setvarx,vary=setvary,norm=norm)
        except AttributeError:
            setvarx = ''
            setvary = ''
        
        if self.remfrm.get():
            remfrm = '_sfm'
        else:
            remfrm = ''
        if setvary == 'Custom ROI':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif setvary == 'ROI - bkg':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi_bkg{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        
        fittype = self.fittype.get()
        if fittype == 'None': fittype=None
        
        if len(self.extra_scannos) > 0:
            scanno = [scanno]+list(self.extra_scannos)
        
        # Send the command
        cmdstr = 'pp.plotscan({},varx=\'{}\',vary=\'{}\',fit=\'{}\',norm={})'
        print( cmdstr.format(scanno,setvarx,setvary,fittype,norm) )
        pp.plotscan(scanno,varx=setvarx,vary=setvary,fit=fittype,norm=norm,logplot=logplot,diffplot=diffplot)
        plt.show()
    
    def f_fnl_spil(self):
        "Send plotpil command to console"
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        LOG = self.logplot.get()
        setvarx = self.varx.get()
        
        ROIcen = [self.pilcen_i.get(),self.pilcen_j.get()]
        ROIsize = [self.roisiz_i.get(),self.roisiz_j.get()]
        cax = [self.pilint_i.get(),self.pilint_j.get()]
        imnum = self.pilpos
        ROIbkg = self.rembkg.get()
        if self.remfrm.get():
            bkg_img=0
        else:
            bkg_img=None
        
        # Send the command
        cmdstr = 'pp.plotpil({},cax={},varx={},imnum={},bkg_img={},ROIcen={},ROIsize={}, show_ROIbkg=ROIbkg)'
        print( cmdstr.format(scanno,cax,setvarx,imnum,bkg_img,ROIcen,ROIsize,ROIbkg) )
        pp.plotpil(scanno,cax=cax,varx=setvarx,imnum=imnum,bkg_img=bkg_img,ROIcen=ROIcen,ROIsize=ROIsize,show_ROIbkg=ROIbkg,log_colors=LOG)
        plt.show()
    
    def f_fnl_phkl(self):
        "Send plotpilhkl command to console"
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        ROIcen = [self.pilcen_i.get(),self.pilcen_j.get()]
        ijk = [ROIcen[0],ROIcen[1],int(self.N.get())//2]
        
        # Send the command
        cmdstr = 'pp.plotpilhkl_cuts({},hkl_centre=None,sum_tolarance=[0.05,0.05,0.05],max_points=301)'
        print( cmdstr.format(scanno) )
        pp.plotpilhkl_cuts(scanno,image_centre=ijk,sum_tolarance=[0.05,0.05,0.05],max_points=301)
        plt.show()
    
    def f_fnl_ptth(self):
        "Send plotpiltth command to console"
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        
        # Send the command
        cmdstr = 'pp.plotpiltth({},binsep=0.1)'
        print( cmdstr.format(scanno) )
        eval(cmdstr.format(scanno))
        #pp.plotpiltth(scanno)
        plt.show()
    
    def f_fnl_splotsave(self):
        "Send plotscan command to console, save result and close"
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        logplot = self.logplot.get()
        diffplot = self.diffplot.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        # test the vary
        setvarx = self.varx.get()
        setvary = self.vary.get()
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(scanno,varx=setvarx,vary=setvary,norm=norm)
        except AttributeError:
            setvarx = ''
            setvary = ''
        
        if self.remfrm.get():
            remfrm = '_sfm'
        else:
            remfrm = ''
        if setvary == 'Custom ROI':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif setvary == 'ROI - bkg':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi_bkg{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        
        fittype = self.fittype.get()
        if fittype == 'None': fittype=None
        
        # Send the command
        cmdstr = 'pp.plotscan({},varx=\'{}\',vary=\'{}\',fit=\'{}\',norm={})'
        print( cmdstr.format(scanno,setvarx,setvary,fittype,norm) )
        pp.plotscan(scanno,varx=setvarx,vary=setvary,fit=fittype,norm=norm,save=True,logplot=logplot,diffplot=diffplot)
        plt.close(plt.gcf())
        pp.savescan(scanno,varx=setvarx,vary=setvary,norm=norm)
        
        self.helper.set('Plot image saved in analysis folder as "S{} ... .png"'.format(scanno))
    
    def f_fnl_splotprint(self):
        "Send plotscan command to console, print result and close"
        
        # Get the parameters
        self.set_files()
        pp.normby = self.normtype.get()
        scanno = self.scanno.get()
        logplot = self.logplot.get()
        diffplot = self.diffplot.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        # test the vary
        setvarx = self.varx.get()
        setvary = self.vary.get()
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(scanno,varx=setvarx,vary=setvary,norm=norm)
        except AttributeError:
            setvarx = ''
            setvary = ''
        
        if self.remfrm.get():
            remfrm = '_sfm'
        else:
            remfrm = ''
        if setvary == 'Custom ROI':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        elif setvary == 'ROI - bkg':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            setvary = 'nroi_bkg{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
        
        fittype = self.fittype.get()
        if fittype == 'None': fittype=None
        
        # Send the command
        cmdstr = 'pp.plotscan({},varx=\'{}\',vary=\'{}\',fit=\'{}\',norm={},logplot={},diffplot={})'
        print( cmdstr.format(scanno,setvarx,setvary,fittype,norm,logplot,diffplot) )
        pp.plotscan(scanno,varx=setvarx,vary=setvary,fit=fittype,norm=norm,logplot=logplot,diffplot=diffplot)
        fig = plt.gcf()
        
        I16_Print_Buffer([fig.number],[3,2])
        plt.close(fig)
    
    def f_fnl_splotbuffer(self):
        " Send all open figures to the print buffer for 6-per-page page printing"
        
        "Get figure numbers automatically"
        fignos = plt.get_fignums()
        
        " Determine the number of print figures to generate"
        Nax = default_print_layout[0]*default_print_layout[1]
        Nfigs = np.ceil( len(fignos)/float(Nax) ).astype(np.int)
        
        self.helper.set('Sending {} figures to {} print buffers'.format(len(fignos),Nfigs))
        
        " Create the buffers"
        for nfig in range(Nfigs):
            buffer_fignos = fignos[nfig*Nax:nfig*Nax+Nax]
            I16_Print_Buffer(buffer_fignos,default_print_layout)
    
    def f_fnl_splotclose(self):
        " Close all figures"
        plt.close('all')

    
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    
    
    def set_files(self):
        "Set the filedir and savedir commands"
        
        if pp.filedir != self.filedir.get():
            pp.filedir = self.filedir.get()
            pp.exp_list_add()
        if pp.savedir != self.savedir.get():
            pp.savedir = self.savedir.get()
            pp.sav_list_add()
            pp.exp_parameters_load() # load Py16 parameters
    
    def writeval(self,mstr,label_str,label_var,wid=27):
        frame = tk.Frame(mstr)
        frame.pack(fill=tk.X)        
        str1 = tk.Label(frame, text='{:10s} : '.format(label_str),
                        font=SF, width= 14, anchor=tk.E)
        str1.pack(side=tk.LEFT)
        str2 = tk.Label(frame, textvariable=label_var,
                        font=SF, width= wid, anchor=tk.W)
        str2.pack(side=tk.LEFT)
    
    def on_closing(self):
        "End mainloop on close window"
        self.root.destroy()
    
    def update_details(self,event=None):
        "Load metadata for current scan"
        
        self.pilatus_active = False
        self.extra_scannos = []
        
        self.set_files()
        
        d = pp.readscan(self.scanno.get())
        
        if d is None:
            # Set frame variables
            self.cmd.set('')
            self.N.set('')
            self.HKL.set('')
            self.ENG.set('')
            self.T.set('')
            self.atten.set('')
            self.trans.set('')
            self.mm.set('')
            self.do.set('')
            self.pol.set('thp =   tthp =   pol = ')
    
            self.eta.set('')
            self.chi.set('')
            self.dlt.set('')
            self.mu .set('')
            self.gam.set('')
            self.azir.set('')
            self.psi.set('')
            self.phi.set('')
    
            self.sx.set('')
            self.sy.set('')
            self.sz.set('')
            self.spara.set('')
            self.sperp.set('')
            
            self.ss.set('')
            self.ds.set('')
            
            self.custom1_val.set('')
            self.custom2_val.set('')
            
            self.runtime.set('')
            self.timetaken.set('')
            return
        
        m = d.metadata
        mkeys = list(m.keys())
        keys = [x for x in d.keys() if type(d[x]) == np.ndarray ] # only keys linked to arrays
        keysx = ['Auto']+keys
        keys_all = ['Custom']+keys+mkeys
        
        # Add Cutsom ROI to keys for image files
        if hasattr(d,'path'):
            keysy = ['Auto','Custom ROI','ROI - bkg']+keys[::-1]
        else:
            keysy = ['Auto']+keys[::-1]
        
        if self.varx.get() not in keysx:
            self.varx.set('Auto')
        if self.vary.get() not in keysy:
            self.vary.set('Auto')
            self.fittype.set('None')
        
        # Set varx and vary optionmenu entrys
        self.opt_varx['menu'].delete(0, 'end')
        for key in keysx:
            #self.opt_varx['menu'].add_command(label=key, command=tk._setit(self.varx, key))
            self.opt_varx['menu'].add_command(label=key, command=lambda k=key: self.f_popt_varx(k))
        
        self.opt_vary['menu'].delete(0, 'end')
        for key in keysy:
            #self.opt_vary['menu'].add_command(label=key, command=tk._setit(self.vary, key))
            self.opt_vary['menu'].add_command(label=key, command=lambda k=key: self.f_popt_vary(k))
        
        # Set Custom1 and custom 2 option menu entrys
        """
        self.opt_cust1['menu'].delete(0, 'end')
        self.opt_cust2['menu'].delete(0, 'end')
        for key in keys_all:
            self.opt_cust1['menu'].add_command(label=key, command=lambda k=key: self.f_detl_custom1(k))
            self.opt_cust2['menu'].add_command(label=key, command=lambda k=key: self.f_detl_custom2(k))
        """
        
        # Initilise frame variables
        self.cmd.set(m.cmd_short)
        self.N.set(str(len(d[keys[0]])))
        self.HKL.set(m.hkl_str)
        self.ENG.set('{} keV'.format(m.Energy))
        self.T.set(m.temperature)
        
        self.do.set(m.delta_axis_offset)

        self.eta.set(m.eta)
        self.chi.set(m.chi)
        self.dlt.set(m.delta)
        self.mu .set(m.mu)
        self.gam.set(m.gam)
        self.azir.set('({0},{1},{2})'.format(m.azih,m.azik,m.azil))
        self.psi.set(np.round(m.psi,3))
        self.phi.set(np.round(m.phi,3))

        self.sx.set(m.sx)
        self.sy.set(m.sy)
        self.sz.set(m.sz)
        self.spara.set(m.spara)
        self.sperp.set(m.sperp)
        
        self.ss.set('{}'.format(pp.scanss(d)))
        self.ds.set('{}'.format(pp.scands(d)))
        
        # Minimirrors
        if m.m4pitch > 0.1: 
            mm = 'in' 
        else: 
            mm = 'out'
        self.mm.set('{} ({:4.2f} deg)'.format(mm,m.m4pitch))
        
        self.atten.set('{} ({}%)'.format(m.Atten,m.Transmission*100))
        self.trans.set(m.Transmission)
        
        # Analyser
        self.pol.set('thp: {}     tthp: {}     pol: {}'.format(m.thp,m.tthp,m.stoke))
        
        # Custom
        custom1 = self.custom1.get()
        custom2 = self.custom2.get()
        if custom1 in d.metadata.keys():
            self.custom1_val.set('{}'.format(getattr(d.metadata,custom1)))
        elif custom1 in keys:
            self.custom1_val.set('{}'.format(np.mean(getattr(d,custom1))))
        else:
            self.custom1_val.set('--')
        if custom2 in d.metadata.keys():
            self.custom2_val.set('{}'.format(getattr(d.metadata,custom2)))
        elif custom2 in keys:
            self.custom2_val.set('{}'.format(np.mean(getattr(d,custom2))))
        else:
            self.custom2_val.set('--')
        
        # Timing
        time = d.TimeSec[-1]-d.TimeSec[0]
        hours = np.int(np.floor(time/3600))
        mins = np.int(np.floor(np.remainder(time,3600)/60))
        secs = np.remainder(np.remainder(time,3600),60)
        self.runtime.set(m.date)
        self.timetaken.set('{} hours, {} mins, {} seconds'.format(hours,mins,secs))
    
    def update_plot(self,event=None):
        "Plot metadata for current scan"
        
        for eplt in self.extra_plots:
            eplt.remove()
        self.extra_plots = []
        
        self.set_files()
        pp.normby = self.normtype.get()
        d = pp.readscan(self.scanno.get())
        if d is None:
            return
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        yvar = self.vary.get()
        if self.remfrm.get():
            remfrm = '_sfm'
        else:
            remfrm = ''
        if yvar == 'Custom ROI':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            yvar = 'nroi{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
            self.helper.set('Plotting Region of Interest (ROI): '+yvar)
        elif yvar == 'ROI - bkg':
            ROIceni = self.pilcen_i.get()
            ROIcenj = self.pilcen_j.get()
            ROIsizei = self.roisiz_i.get()
            ROIsizej = self.roisiz_j.get()
            yvar = 'nroi_bkg{}[{},{},{},{}]'.format(remfrm,ROIceni,ROIcenj,ROIsizei,ROIsizej)
            self.helper.set('Plotting Region of Interest (ROI): '+yvar)
        
        
        xvar = self.varx.get()
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(d,varx=xvar,vary=yvar,norm=norm)
        except AttributeError:
            x,y,dy,varx,vary,ttl,d = pp.getdata(d,varx='',vary='',norm=norm)
        
        if self.diffplot.get():
            y = np.abs(np.gradient(y))
            vary = 'd {} / d {}'.format(vary,varx)
        
        if self.logplot.get():
            y = np.log(y)
            vary = 'log({})'.format(vary)
        
        # Update main Plot
        self.plt1.set_xdata(x)
        self.plt1.set_ydata(y)
        self.plt2.set_xdata([]) # remove pilatus marker
        self.plt2.set_ydata([])
        self.ax1.relim()
        self.ax1.autoscale_view()
        self.toolbar.update() # update toolbar home position
        
        self.ax1.set_xlabel(varx)
        self.ax1.set_ylabel(vary)
        self.ax1.set_title('#{}'.format(self.scanno.get()),fontsize=16)
        
        # Multiplots
        for n,enum in enumerate(self.extra_scannos):
            try:
                ex,ey,edy,evarx,evary,ettl,ed = pp.getdata(enum,varx=xvar,vary=yvar,norm=norm)
            except AttributeError:
                ex,ey,edy,evarx,evary,ettl,ed = pp.getdata(enum,varx='',vary='',norm=norm)
            
            if self.diffplot.get():
                ey = np.abs(np.gradient(ey))
            if self.logplot.get():
                ey = np.log(ey)
            col = pp.plot_colors[(n+1) %len(pp.plot_colors)]
            epl, = self.ax1.plot(ex,ey,'-o',c=col,lw=2)
            self.extra_plots += [epl]
        
        # Reset pilatus plot
        if self.autopilplot.get():
            self.update_pilatus()
        elif self.pilatus_scan == self.scanno.get():
            self.update_pilatus()
        else:
            self.pilatus_active = False
            self.pilatus_scan = 0
        
        # Perform Fit
        fit_type = self.fittype.get()
        if fit_type != 'None':
            out,err = pp.peakfit(x,y,type=fit_type,Nloop=100,peaktest=-100,interpolate=True,disp=True)
        
            xfit,yfit=out['x'],out['y']
            amp,damp = out['Peak Height'],err['Peak Height']
            cen,dcen = out['Peak Centre'],err['Peak Centre']
            wid,dwid = out['FWHM'],err['FWHM']
            bkg,dbkg = out['Background'],err['Background']
            ara,dara = out['Area'],err['Area']
            
            " add to plot"
            self.pfit.set_xdata(xfit)
            self.pfit.set_ydata(yfit)
        else:
            " remove from plot"
            self.pfit.set_xdata([])
            self.pfit.set_ydata([])
        
        self.ax1.relim()
        self.ax1.autoscale(True)
        self.ax1.autoscale_view()
        self.fig1.canvas.draw()
        self.toolbar.update() # update toolbar home position
    
    def update(self,event=None):
        "Update on <Enter>"
        self.update_details()
        self.update_plot()
    
    def update_pilatus(self,event=None):
        "Loads Pilatus data, plots data"
        
        if self.pilatus_active == False:
            self.helper.set('Loading pilatus images for scan #{}, use the arrows to scroll'.format(self.scanno.get()))
            self.pilatus_active = True
            self.pilatus_scan = self.scanno.get()
            self.set_files()
            
            
            # get data
            xvar = self.varx.get()
            x,y,dy,varx,vary,ttl,d = pp.getdata(self.scanno.get(),varx=xvar)
            # Check if scan has pilatus data
            if hasattr(d,'path'):
                self.vol = pp.getvol(self.scanno.get())
            else:
                self.pilatus_active = False
                self.helper.set('This isnt a pilatus scan!')
                return
            
            # Check volume size
            pil_size = self.vol.shape
            if pil_size[0]*pil_size[1]*len(x) > pp.pil_max_size:
                self.helper.set('The array size is too large! Returning only the first image')
                x = x[:1]
            
            # Subtract first frame from all images
            if self.remfrm.get():
                bkgcut = self.vol[:,:,0].copy()
                for n in range(len(x)):
                    self.vol[:,:,n] = self.vol[:,:,n] - bkgcut
            
            # Set pilatus position
            if self.livemode.get():
                self.pilpos = len(x)-1
            else:
                self.pilpos = int(len(x)//2)
            
            # Set pilatus variables
            imgno = self.pilpos
            xval = x[imgno]
            #self.pilpos = imgno
            self.pilvar = varx
            self.pilvals = x
            self.pilstr.set('({}) {}={:1.3f}'.format(imgno,varx,xval))
            
            # Automatically determine caxis
            md = np.median(self.vol[:,:,imgno])
            mx = np.max(self.vol[:,:,imgno])
            cmax = md + 10**(0.7*np.log10(mx-md))
            if cmax <= 0: cmax = 1
            self.pilint_i.set(0)
            self.pilint_j.set(int(cmax))
            
            if self.logplot.get():
                self.pilim.norm = LogNorm()
                self.pilint_i.set(0.1)
                self.pilint_j.set(mx)
            else:
                self.pilim.norm = Normalize()
            
            # Update pilatus image
            self.pilim.remove() # remove old image, stops memory buildup
            self.ax2.set_ylim(self.vol.shape[0],-0.5)
            self.ax2.set_xlim(-0.5,self.vol.shape[1])
            #self.pilim.set_data(self.vol[:,:,imgno])
            self.pilim = self.ax2.imshow(self.vol[:,:,imgno])
            
            #print(self.ax2.get_xlim(),self.ax2.get_ylim())
            #print(self.ax2.get_position())
            
            # Remove peakregion
            self.pilp5.set_xdata([])
            self.pilp5.set_ydata([])
        
        # Set intensity cutoffs
        self.pilim.set_clim([self.pilint_i.get(),self.pilint_j.get()])

        # update lines
        ROIcen = [self.pilcen_i.get(),self.pilcen_j.get()]
        ROIsize = [self.roisiz_i.get(),self.roisiz_j.get()]
        pil_centre = [self.pilcen_i.get(),self.pilcen_j.get()]
        idxi = np.array([ROIcen[0]-ROIsize[0]//2,ROIcen[0]+ROIsize[0]//2+1])
        idxj = np.array([ROIcen[1]-ROIsize[1]//2,ROIcen[1]+ROIsize[1]//2+1])
        
        pil_size = self.vol.shape
        self.pilp1.set_xdata(idxj[[0,1,1,0,0]])
        self.pilp1.set_ydata(idxi[[0,0,1,1,0]])
        self.pilp2.set_xdata([pil_centre[1],pil_centre[1]])
        self.pilp2.set_ydata([0,pil_size[0]])
        self.pilp3.set_xdata([0,pil_size[1]])
        self.pilp3.set_ydata([pil_centre[0],pil_centre[0]])
        
        if self.rembkg.get():
            "---Background ROI---" "FROM fun pilroi()"
            " The background region is twice the area of the required region"
            idxbi = np.array([ROIcen[0]-ROIsize[0],ROIcen[0]+ROIsize[0]])
            idxbj = np.array([ROIcen[1]-ROIsize[1],ROIcen[1]+ROIsize[1]])
            if idxbi[0] < 0: idxbi[1] = idxbi[1] - idxbi[0]; idxbi[0] = 0
            if idxbi[1] > pil_size[0]: idxbi[0] = idxbi[0] - (idxbi[1]-pil_size[0]); idxbi[1] = pil_size[0]
            if idxbj[0] < 0: idxbj[1] = idxbj[1] - idxbj[0]; idxbj[0] = 0
            if idxbj[1] > pil_size[1]: idxbj[0] = idxbj[0] - (idxbj[1]-pil_size[1]); idxbj[1] = pil_size[1]
            self.pilp4.set_xdata(idxbj[[0,1,1,0,0]])
            self.pilp4.set_ydata(idxbi[[0,0,1,1,0]])
        else:
            self.pilp4.set_xdata([])
            self.pilp4.set_ydata([])
        
        # Update pilatus image
        self.pilim.set_data(self.vol[:,:,self.pilpos])
        self.ax2.set_aspect('equal')
        #self.ax2.autoscale(tight=True)
        self.fig2.canvas.draw()
        
        # Update marker point on plot axis
        ylim = self.ax1.get_ylim()
        self.plt2.set_xdata([self.pilvals[self.pilpos],self.pilvals[self.pilpos]])
        self.plt2.set_ydata(ylim)
        self.fig1.canvas.draw()


"------------------------------------------------------------------------"
"--------------------------I16_Scan_Selector-----------------------------"
"------------------------------------------------------------------------"
class I16_Scan_Selector:
    """
    Displays list of previous scans, selecting them will change the currently
    displayed scan in the main window.
    Opens on click of the "Last Scans window" button in I16_Data_Viewer
        - If "All scans" is ticked in the main window, this window will display
          all scans from this experiment
        - Otherwise, the last 100 scans will be displayed
    
    The "Update" button adds the most recent scans to the top of the list
    The "Reload" button reloads the scans, but also displays data associated with
    any metadata defined in the "show" box.
    
    E.G. 
        Show: Ta >> Reload - scans are shown with temperatures
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,parent,scanrange,showval=''):
        # Create new window based on main window instance
        self.parent = parent
        self.root = tk.Toplevel(parent.root)
        self.root.wm_title('I16 Scan Selector by D G Porter [dan.porter@diamond.ac.uk]')
        viewerheight=parent.root.winfo_height()
        def_width = 200
        self.root.minsize(width=def_width, height=viewerheight)
        self.root.maxsize(width=800, height=viewerheight)
        
        # parent.root.winfo_geometry()
        #mainx = parent.root.winfo_x()
        #mainy = parent.root.winfo_y()
        geom = parent.root.geometry()
        mainx,mainy = [int(x) for x in geom.split('+')[1:]]
        self.root.geometry('%dx%d+%d+%d' % (def_width,viewerheight,mainx-def_width,mainy))
        
        frame = tk.Frame(self.root)
        frame.pack(fill=tk.BOTH,expand=tk.YES)
        
        # Create scan info
        self.scan_nos = np.asarray(scanrange[::-1],dtype=int)
        #self.scan_nos = pp.get_all_scannos()
        #fmt = '{cmd} | {date} | {mode:4s} {ss} {ds} {energy:5.3g}keV {temp:5.3g}K {hkl:17s} {show}{equal}{val} | {cmd}'
        if showval is None or len(showval) == 0:
            fmt = '{:6.0f}: {}'
        else:
            fmt = '{:6.0f}:  {}  {}'
        
        "-----------------------------Scan ListBox---------------------------"
        # Eval box with scroll bar
        frm_scan = tk.Frame(frame)
        frm_scan.pack(side=tk.TOP,fill=tk.BOTH,expand=tk.YES)
        
        scl_scanx = tk.Scrollbar(frm_scan,orient=tk.HORIZONTAL)
        scl_scanx.pack(side=tk.BOTTOM, fill=tk.BOTH)
        
        scl_scany = tk.Scrollbar(frm_scan)
        scl_scany.pack(side=tk.RIGHT, fill=tk.BOTH)
        
        self.lst_scan = tk.Listbox(frm_scan, font=HF, selectmode=tk.EXTENDED,#tk.SINGLE,#width=40,height=30,
                                xscrollcommand=scl_scanx.set,yscrollcommand=scl_scany.set)
        self.lst_scan.configure(exportselection=False)
        
        # Populate list box
        for scanno in self.scan_nos:
            d = pp.readscan(scanno)
            if d is None:
                self.lst_scan.insert(tk.END,'{:6.0f}: No Scan'.format(scanno)) 
                continue
            if showval is None or len(showval) == 0:
                strval = fmt.format(scanno,d.metadata.cmd_short)
            else:
                if showval == 'hkl':
                    val = pp.scanhkl(scanno)
                elif showval in d.keys():
                    val = np.mean(getattr(d,showval))
                elif showval in d.metadata.keys():
                    val = getattr(d.metadata,showval)
                else:
                    val = 'N/A'
                strval = fmt.format(scanno,val,d.metadata.cmd_short)
            self.lst_scan.insert(tk.END,strval) 
        
        self.lst_scan.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES)
        self.lst_scan.select_set(0)
        self.lst_scan.bind("<<ListboxSelect>>", self.f_scan_select)
        #print( self.lst_scan.curselection()[0] )
        
        scl_scanx.config(command=self.lst_scan.xview)
        scl_scany.config(command=self.lst_scan.yview)
        
        #self.txt_scan.config(xscrollcommand=scl_scanx.set,yscrollcommand=scl_scany.set)
        
        "----------------------------Showval entry---------------------------"
        frm_show = tk.Frame(frame)
        frm_show.pack()
        
        lbl_show = tk.Label(frm_show,text='Show:')
        lbl_show.pack(side=tk.LEFT)
        
        self.new_showval = tk.StringVar(frm_show,showval)
        ety_show = tk.Entry(frm_show,textvariable=self.new_showval, width=8)
        ety_show.pack(side=tk.LEFT,fill=tk.X,expand=tk.YES)
        
        "----------------------------Reload Button---------------------------"
        frm_btn = tk.Frame(frame)
        frm_btn.pack()
        
        btn_reld = tk.Button(frm_btn,text='Reload',font=BF, command=self.f_reload)
        btn_reld.pack(side=tk.LEFT)
        
        btn_updt = tk.Button(frm_btn,text='Update',font=BF, command=self.f_update)
        btn_updt.pack(side=tk.LEFT)
        
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_scan_select(self,event):
        "Select scan number from list box and send plot command to main window"
        
        current = list(self.lst_scan.curselection())
        scanno = self.scan_nos[current]
        
        self.parent.scanno.set(scanno[0])
        self.parent.update()
        if len(scanno) > 1:
            # Multi-plots
            self.parent.extra_scannos = scanno[1:]
            self.parent.update_plot()
            """
            # Try to color the highlighted scans the same colors as the plot *doesn't work*
            col_conv={'k':"black", 'r':"red", 'g': "green", 'b':"blue", 'c':"cyan", 'y':"yellow", 'm':"magenta"}
            for n in range(len(current)):
                col=pp.plot_colors[n%len(pp.plot_colors)]
                
                self.lst_scan.config(selectbackground=col_conv[col])
                #self.lst_scan.itemconfig(current[n],bg=col_conv[col])
            """
    
    def f_reload(self):
        "Closes the current selector window and reloads it"
        self.parent.scan_selector_showval = self.new_showval.get()
        self.root.destroy()
        self.parent.f_scn_sel()
    
    def f_update(self):
        "Closes the current selector window and reloads it"
        
        latest = pp.latest()
        last = self.scan_nos[0]
        if latest <= last: return
        
        new_scannos = range(last+1,latest+1)
        self.scan_nos = np.append(new_scannos[::-1],self.scan_nos)
        
        showval = self.new_showval.get()
        if showval is None or len(showval) == 0:
            fmt = '{:6.0f}: {}'
        else:
            fmt = '{:6.0f}:  {}  {}'
        
        for scanno in new_scannos:
            d = pp.readscan(scanno)
            if d is None:
                self.lst_scan.insert(tk.END,'{:6.0f}: No Scan'.format(scanno)) 
                continue
            if showval is None or len(showval) == 0:
                strval = fmt.format(scanno,d.metadata.cmd_short)
            else:
                if showval == 'hkl':
                    val = pp.scanhkl(scanno)
                elif showval in d.keys():
                    val = np.mean(getattr(d,showval))
                elif showval in d.metadata.keys():
                    val = getattr(d.metadata,showval)
                else:
                    val = 'N/A'
                strval = fmt.format(scanno,val,d.metadata.cmd_short)
            self.lst_scan.insert(0,strval) 


"------------------------------------------------------------------------"
"--------------------------I16_Meta_Display------------------------------"
"------------------------------------------------------------------------"
class I16_Meta_Display:
    """
    Displays the full meta data of the current scan
    Activated from the "Meta" button in I16_Data_Viewer
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,scanno=0):
        # Create scan info
        self.d = pp.readscan(scanno)
        scanno = self.d.metadata.SRSRUN
        
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('Metadata: #{}'.format(scanno))
        self.root.minsize(width=300, height=400)
        self.root.maxsize(width=800, height=800)
        
        #Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES,anchor=tk.N)
        
        "---------------------------Metadata ListBox---------------------------"
        # Eval box with scroll bar
        frm_meta = tk.Frame(frame)
        frm_meta.pack(side=tk.TOP,fill=tk.BOTH,expand=tk.YES)
        
        scl_metax = tk.Scrollbar(frm_meta,orient=tk.HORIZONTAL)
        scl_metax.pack(side=tk.BOTTOM, fill=tk.BOTH)
        
        scl_metay = tk.Scrollbar(frm_meta)
        scl_metay.pack(side=tk.RIGHT, fill=tk.BOTH)
        
        self.lst_meta = tk.Listbox(frm_meta, font=HF, selectmode=tk.SINGLE,width=40,height=40,
                                xscrollcommand=scl_metax.set,yscrollcommand=scl_metay.set)
        self.lst_meta.configure(exportselection=True)
        
        # Populate list box
        for k in self.d.metadata.keys():
            if k[0] == '_': continue # Omit _OrderedDict__root/map
            strval = '{:>20} : {:<20}'.format(k,self.d.metadata[k])
            self.lst_meta.insert(tk.END,strval)
        
        self.lst_meta.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES)
        #self.lst_meta.select_set(0)
        #self.lst_meta.bind("<<ListboxSelect>>", self.f_meta_select)
        #print( self.lst_meta.curselection()[0] )
        
        scl_metax.config(command=self.lst_meta.xview)
        scl_metay.config(command=self.lst_meta.yview)
        
        #self.txt_meta.config(xscrollcommand=scl_metax.set,yscrollcommand=scl_metay.set)
        
        "----------------------------Exit Button------------------------------"
        frm_btn = tk.Frame(frame)
        frm_btn.pack()
        
        btn_exit = tk.Button(frm_btn,text='Exit',font=BF, command=self.f_exit)
        btn_exit.pack()
        
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    
    def f_exit(self):
        "Closes the current metadata window"
        self.root.destroy()


"------------------------------------------------------------------------"
"--------------------------I16_Meta_Display------------------------------"
"------------------------------------------------------------------------"
class I16_Meta_Select:
    """
    Displays all metadata fields and returns a selection
    I16_Meta_Select(root, widget, ['field1','field2','field3'])
    Making a selection returns the field string
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,parent, field, meta_fields):
        self.parent = parent # main root 
        self.field = field # custom1/2 widget
        self.meta_fields = meta_fields
        
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('Select Metadata')
        self.root.minsize(width=30, height=300)
        self.root.maxsize(width=800, height=800)
        self.output = 'Custom'
        
        #Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES,anchor=tk.N)
        
        "---------------------------Metadata ListBox---------------------------"
        # Eval box with scroll bar
        frm_meta = tk.Frame(frame)
        frm_meta.pack(side=tk.TOP,fill=tk.BOTH,expand=tk.YES)
        
        scl_metax = tk.Scrollbar(frm_meta,orient=tk.HORIZONTAL)
        scl_metax.pack(side=tk.BOTTOM, fill=tk.BOTH)
        
        scl_metay = tk.Scrollbar(frm_meta)
        scl_metay.pack(side=tk.RIGHT, fill=tk.BOTH)
        
        self.lst_meta = tk.Listbox(frm_meta, font=HF, selectmode=tk.SINGLE,width=15,height=20,
                                xscrollcommand=scl_metax.set,yscrollcommand=scl_metay.set)
        self.lst_meta.configure(exportselection=True)
        
        # Populate list box
        for k in self.meta_fields:
            #if k[0] == '_': continue # Omit _OrderedDict__root/map
            strval = '{}'.format(k)
            self.lst_meta.insert(tk.END,strval)
        
        self.lst_meta.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES)
        self.lst_meta.select_set(0)
        self.lst_meta.bind("<<ListboxSelect>>", self.f_meta_select)
        #print( self.lst_meta.curselection()[0] )
        
        scl_metax.config(command=self.lst_meta.xview)
        scl_metay.config(command=self.lst_meta.yview)
        
        #self.txt_meta.config(xscrollcommand=scl_metax.set,yscrollcommand=scl_metay.set)
        
        "----------------------------Exit Button------------------------------"
        frm_btn = tk.Frame(frame)
        frm_btn.pack()
        
        btn_exit = tk.Button(frm_btn,text='Exit',font=BF, command=self.f_exit)
        btn_exit.pack()

        "-------------------------Start Mainloop------------------------------"
        self.root.protocol("WM_DELETE_WINDOW", self.f_exit)
        #self.root.mainloop()
        
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    
    def f_meta_select(self, event):
        "Choose a metadata field"
        self.output = self.meta_fields[self.lst_meta.curselection()[0]]
        if self.output[0] == '_': return # Omit _OrderedDict__root/map
        self.field.set(self.output)
        self.parent.update_details()
        self.root.destroy()

    def f_exit(self):
        "Closes the current metadata window"
        self.root.destroy()


"------------------------------------------------------------------------"
"--------------------------I16_Text_Display------------------------------"
"------------------------------------------------------------------------"
class I16_Text_Display:
    """
    General box to display any long string
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self, disp_str='', title=''):
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('Py16: {}'.format(title))
        self.root.minsize(width=200, height=100)
        self.root.maxsize(width=1200, height=800)
        
        # Dynamic box size
        lines = disp_str.split('\n')
        n_lines = len(lines)
        line_len = [len(ln) for ln in lines]
        max_len = np.max(line_len)

        box_height = n_lines
        box_width = max_len

        if n_lines > 20:
            box_height = 20
        if max_len > 200:
            box_width = 200
        

        #Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES,anchor=tk.N)
        
        "---------------------------Metadata ListBox---------------------------"
        # Eval box with scroll bar
        frm_text = tk.Frame(frame)
        frm_text.pack(side=tk.TOP,fill=tk.BOTH,expand=tk.YES)
        
        scl_textx = tk.Scrollbar(frm_text,orient=tk.HORIZONTAL)
        scl_textx.pack(side=tk.BOTTOM, fill=tk.BOTH)
        
        scl_texty = tk.Scrollbar(frm_text)
        scl_texty.pack(side=tk.RIGHT, fill=tk.BOTH)
        
        self.lst_text = tk.Text(frm_text, font=HF, width=box_width, height=box_height, wrap=tk.NONE,
                                xscrollcommand=scl_textx.set,yscrollcommand=scl_texty.set)
        self.lst_text.configure(exportselection=True)
        
        # Populate text box
        self.lst_text.insert(tk.END, disp_str)
        self.lst_text.configure(state=tk.DISABLED)
        
        self.lst_text.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES)
        #self.lst_text.select_set(0)
        #self.lst_text.bind("<<ListboxSelect>>", self.f_text_select)
        #print( self.lst_text.curselection()[0] )
        
        scl_textx.config(command=self.lst_text.xview)
        scl_texty.config(command=self.lst_text.yview)
        
        #self.txt_text.config(xscrollcommand=scl_textx.set,yscrollcommand=scl_texty.set)
        
        "----------------------------Exit Button------------------------------"
        frm_btn = tk.Frame(frame)
        frm_btn.pack()
        
        btn_exit = tk.Button(frm_btn,text='Exit',font=BF, command=self.f_exit)
        btn_exit.pack()
        
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    
    def f_exit(self):
        "Closes the current text window"
        self.root.destroy()


"------------------------------------------------------------------------"
"--------------------------I16_Python Editor-----------------------------"
"------------------------------------------------------------------------"
class I16_Python_Editor:
    """
    A very simple python editor, load and edit python files, execute them in 
    current python shell.
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self, disp_str='', filename=''):
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('I16 Python Editor by D G Porter [dan.porter@diamond.ac.uk]')
        self.root.minsize(width=200, height=100)
        self.root.maxsize(width=1200, height=800)
        
        box_height = 30
        box_width = 100

        self.savelocation = ''

        if os.path.isfile(filename):
            self.root.wm_title(filename)
            self.savelocation = filename

        #Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES,anchor=tk.N)
        
        "---------------------------Metadata ListBox---------------------------"
        # Eval box with scroll bar
        frm_text = tk.Frame(frame)
        frm_text.pack(side=tk.TOP,fill=tk.BOTH,expand=tk.YES)
        
        scl_textx = tk.Scrollbar(frm_text,orient=tk.HORIZONTAL)
        scl_textx.pack(side=tk.BOTTOM, fill=tk.BOTH)
        
        scl_texty = tk.Scrollbar(frm_text)
        scl_texty.pack(side=tk.RIGHT, fill=tk.BOTH)
        
        self.text = tk.Text(frm_text, 
            font=HF, 
            width=box_width, 
            height=box_height, 
            wrap=tk.NONE,
            background='white',
            xscrollcommand=scl_textx.set,
            yscrollcommand=scl_texty.set)
        self.text.configure(exportselection=True)
        self.text.bind('<Control-s>',self.f_save)
        self.text.bind('<Control-b>',self.f_run)
        self.text.bind('<Control-r>',self.f_run)
        
        # Populate text box
        self.text.insert(tk.END, disp_str)
        
        self.text.pack(side=tk.LEFT,fill=tk.BOTH,expand=tk.YES)
        
        scl_textx.config(command=self.text.xview)
        scl_texty.config(command=self.text.yview)
        
        #self.txt_text.config(xscrollcommand=scl_textx.set,yscrollcommand=scl_texty.set)
        
        "----------------------------Exit Button------------------------------"
        frm_btn = tk.Frame(frame)
        frm_btn.pack(fill=tk.X)
        
        btn = tk.Button(frm_btn,text='RUN',font=BF, command=self.f_run)
        btn.pack(side=tk.LEFT)
        btn = tk.Button(frm_btn,text='Exit',font=BF, command=self.f_exit)
        btn.pack(side=tk.RIGHT)
        btn = tk.Button(frm_btn,text='Save As', font=BF, command=self.f_saveas)
        btn.pack(side=tk.RIGHT)
        btn = tk.Button(frm_btn,text='Save', font=BF, command=self.f_save)
        btn.pack(side=tk.RIGHT)
        btn = tk.Button(frm_btn,text='Open', font=BF, command=self.f_open)
        btn.pack(side=tk.RIGHT)
        
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_run(self, event=None):
        "Run the code"
        code = self.text.get(1.0, tk.END)
        exec(code)
        print('Finished')

    def f_open(self):
        "Open a new file"
        newsavelocation=filedialog.askopenfilename(
            title='Open your python script',
            initialdir=pp.savedir, 
            initialfile='script.py',
            defaultextension='.py',
            filetypes = (("python file","*.py"),("all files","*.*")))

        if newsavelocation == '':
            return
        with open(newsavelocation) as file:
            disp_str = file.read()
        I16_Python_Editor(disp_str, newsavelocation)

    def f_save(self, event=None):
        "Save the file"
        if self.savelocation == '':
            self.f_saveas()
            return

        code = self.text.get(1.0, tk.END)
        with open(self.savelocation, 'wt') as outfile:
            outfile.write(code)
        self.root.wm_title(self.savelocation)
        print('Saved as {}'.format(self.savelocation))

    def f_saveas(self):
        "Save the file"
        code = self.text.get(1.0, tk.END)
        self.savelocation=filedialog.asksaveasfilename(
            title='Save your python script',
            initialdir=pp.savedir, 
            initialfile='script.py',
            defaultextension='.py',
            filetypes = (("python file","*.py"),("all files","*.*")))
        if self.savelocation != '':
            self.f_save()

    def f_exit(self):
        "Closes the current text window"
        self.root.destroy()


"------------------------------------------------------------------------"
"--------------------------I16_Peak_Analysis-----------------------------"
"------------------------------------------------------------------------"
class I16_Peak_Analysis:
    """
    I16_Peak_Analysis
    
    GUI allowing you to quickly plot multiscan runs such as temperature or energy dependences.
    
    OPERATION:
     - As above, enter your experiment data director, or press "Browse"
     - Also enter the directory you would like files saved in.
     
     - Give the title and the variable that changes between scans (Dependent)
     
     - Define the variable to analyse (Y var) and the fit options (leave blank for defaults).
     
     - Define the scan numbers required by either: 
         1. Specifying the range with first, last and step (then click "Generate")
         2. Specifying scan numbers directly in the box. Note that this box will be evaluated so you must
            specify a python array.
    
     - Use the buttons to select what you want to do next:
         "Plot Scans" - plot multiple scans on the same axis, with "Dependent" as the legend
         "Plot 3D" - Plot multiple scans on a 3D axis, with "Dependent" as the 3rd axis
         "Plot 2D" - As above, except creates a map where the height becomes a colour axis
         "Plot Surf" - As above but generates a surface
         "Fit Peaks" - On each scan, perform a fitting routine, storing the area, width, centre, etc, 
                       plot the results and save them to a .dat file.
    
    Custom regions of interest can be set for scans with area detectors:
        In the "Y" box, use the commands:
        nroi[110,242,75,67] - creates roi [ceni,cenj,widi,widj]
        nroi[31,31] - creates centred roi with size [widi,widj]
        nroi_peak[31,31] - create roi centred around largest pixel
        nroi_bkg[31,31] - create roi and subtract average value of surrounding pixels
        nroi_peak_bkg[11,15] - add functions together
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,initial_num=0,ini_varx='',ini_vary='', initial_range=None):
        # Create Tk inter instance
        root = tk.Tk()
        root.wm_title('I16 Peak Analysis by D G Porter [dan.porter@diamond.ac.uk]')
        root.minsize(width=610, height=420)
        root.maxsize(width=700, height=500)
        
        # Get initial parameters
        initial_dir = pp.filedir
        initial_sav = pp.savedir
        
        # Initial scan string
        if initial_range is None:
            scanstr = '[{},{}]'.format(initial_num,initial_num+1)
        else:
            scanstr = str(initial_range)
        
        # Update default save location for exported plots
        plt.rcParams["savefig.directory"] = pp.savedir
        
        frame = tk.Frame(root)
        frame.pack(side=tk.LEFT,anchor=tk.N)

        "----------------------------Data Directory-----------------------------"
        # Data Folder
        frm_fldr = tk.Frame(frame)
        frm_fldr.pack(fill=tk.X)
        self.filedir = tk.StringVar(frm_fldr,initial_dir)
        lbl_fldr = tk.Label(frm_fldr, text='Data Folder: ', width=15, font=SF)
        lbl_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        ety_fldr = tk.Entry(frm_fldr, textvariable=self.filedir, width=61)
        ety_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        btn_fldr = tk.Button(frm_fldr, text='Browse', font=BF, command=self.f_fldr_browse)
        btn_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        # Save Folder
        frm_save = tk.Frame(frame)
        frm_save.pack(fill=tk.X)
        self.savedir = tk.StringVar(frm_save,initial_sav)
        lbl_save = tk.Label(frm_save, text='Analysis Folder: ', width=15, font=SF)
        lbl_save.pack(side=tk.LEFT,padx=5,pady=5)
        ety_save = tk.Entry(frm_save, textvariable=self.savedir, width=61)
        ety_save.pack(side=tk.LEFT,padx=5,pady=5)
        btn_save = tk.Button(frm_save, text='Browse',font=BF, command=self.f_save_browse)
        btn_save.pack(side=tk.LEFT,padx=5,pady=5)
        
        "-----------------------------Fit Options----------------------------"
        frm_titl = tk.Frame(frame)
        frm_titl.pack(fill=tk.X)
        
        # Title
        self.title = tk.StringVar(frm_titl,pp.exp_title)
        lbl_titl = tk.Label(frm_titl, text='Title: ',font=SF)
        lbl_titl.pack(side=tk.LEFT,padx=5,pady=5)
        ety_titl = tk.Entry(frm_titl, textvariable=self.title, width=35)
        ety_titl.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Dependent variable
        self.depvar = tk.StringVar(frm_titl,pp.default_sensor)
        lbl_depv = tk.Label(frm_titl, text='Dependent: ',font=SF)
        lbl_depv.pack(side=tk.LEFT,padx=5,pady=5)
        ety_depv = tk.Entry(frm_titl, textvariable=self.depvar, width=15)
        ety_depv.pack(side=tk.LEFT,padx=5,pady=5)
        """
        # Differentiate Plot
        self.diffplot = tk.IntVar(frm_titl,0)
        chk_diff = tk.Checkbutton(frm_titl, text='Diff.',font=BF,variable=self.diffplot)
        chk_diff.pack(side=tk.RIGHT,padx=0)
        
        # Log Plot
        self.logplot = tk.IntVar(frm_titl,0)
        chk_log = tk.Checkbutton(frm_titl, text='Log',font=BF,variable=self.logplot)
        chk_log.pack(side=tk.RIGHT,padx=0)
        
        # Normalise Plot
        self.normplot = tk.IntVar(frm_titl,0)
        chk_log = tk.Checkbutton(frm_titl, text='Normalise',font=BF,variable=self.normplot)
        chk_log.pack(side=tk.RIGHT,padx=0)
        """
        frm_opt = tk.Frame(frame)
        frm_opt.pack(fill=tk.X)
        
        # vary box
        self.varx = tk.StringVar(frm_opt,ini_varx)
        lbl_varx = tk.Label(frm_opt, text='X: ',font=SF)
        lbl_varx.pack(side=tk.LEFT,padx=2,pady=5)
        ety_varx = tk.Entry(frm_opt, textvariable=self.varx, width=8)
        ety_varx.pack(side=tk.LEFT,padx=2,pady=5)
        
        # vary box
        self.vary = tk.StringVar(frm_opt,ini_vary)
        lbl_vary = tk.Label(frm_opt, text='Y: ',font=SF)
        lbl_vary.pack(side=tk.LEFT,padx=3,pady=5)
        ety_vary = tk.Entry(frm_opt, textvariable=self.vary, width=14)
        ety_vary.pack(side=tk.LEFT,padx=3,pady=5)
        
        # normalise menu
        normopts = ['rc','ic1','none']
        self.normtype = tk.StringVar(frm_opt, normopts[0])
        lbl_nrm = tk.Label(frm_opt, text='Norm: ', font=SF)
        lbl_nrm.pack(side=tk.LEFT,padx=0,pady=5)
        opt_nrm = tk.OptionMenu(frm_opt, self.normtype, *normopts)
        opt_nrm.config(width=4)
        opt_nrm.pack(side=tk.LEFT,padx=5,pady=5)
        
        # fit menu
        fitopts = ['Simple','Gauss','Lorentz','pVoight','Sum','Max']
        self.fittype = tk.StringVar(frm_opt, fitopts[3])
        lbl_fit = tk.Label(frm_opt, text='Fit: ', font=SF)
        lbl_fit.pack(side=tk.LEFT,padx=0,pady=5)
        opt_fit = tk.OptionMenu(frm_opt, self.fittype, *fitopts)
        opt_fit.config(width=6)
        opt_fit.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Peak Test
        self.peaktest = tk.DoubleVar(frm_opt,1)
        lbl_peak = tk.Label(frm_opt, text='I/sig test: ',font=SF)
        lbl_peak.pack(side=tk.LEFT,padx=5,pady=5)
        ety_peak = tk.Entry(frm_opt, textvariable=self.peaktest, width=4)
        ety_peak.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Save check box
        self.saveopt = tk.IntVar(frm_opt,0)
        chk_save = tk.Checkbutton(frm_opt, text='Save',font=BF,variable=self.saveopt,
                                  onvalue = 1, offvalue = 0)
        chk_save.pack(side=tk.LEFT,padx=5,pady=5)
        
        
        "----------------------------Scan Numbers----------------------------"
        # Scan numbers
        frm_scan = tk.Frame(frame)
        frm_scan.pack(fill=tk.X)
        self.first = tk.IntVar(frm_scan,initial_num)
        self.last = tk.StringVar(frm_scan,'{}+1'.format(initial_num))
        self.step = tk.IntVar(frm_scan,1)
        lbl_first = tk.Label(frm_scan, text='First: ',font=SF)
        lbl_first.pack(side=tk.LEFT,padx=5,pady=5)
        ety_first = tk.Entry(frm_scan, textvariable=self.first, width=10)
        ety_first.pack(side=tk.LEFT,padx=5,pady=5)
        lbl_last = tk.Label(frm_scan, text='Last: ',font=SF)
        lbl_last.pack(side=tk.LEFT,padx=5,pady=5)
        ety_last = tk.Entry(frm_scan, textvariable=self.last, width=10)
        ety_last.pack(side=tk.LEFT,padx=5,pady=5)
        lbl_step = tk.Label(frm_scan, text='Step: ',font=SF)
        lbl_step.pack(side=tk.LEFT,padx=5,pady=5)
        ety_step = tk.Entry(frm_scan, textvariable=self.step, width=5)
        ety_step.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Buttons
        btn_scan_st = tk.Button(frm_scan, text='Last',font=BF, command=self.f_scan_st)
        btn_scan_st.pack(side=tk.LEFT,padx=2)                  
        btn_scan_ld = tk.Button(frm_scan, text='Generate',font=BF, command=self.f_scan_ld)
        btn_scan_ld.pack(side=tk.LEFT)
        btn_scan_ld = tk.Button(frm_scan, text='Check',font=BF, command=self.f_scan_check)
        btn_scan_ld.pack(side=tk.LEFT)
        btn_scan_ld = tk.Button(frm_scan, text='Create File',font=BF, command=self.f_scan_makefile)
        btn_scan_ld.pack(side=tk.LEFT)
        
        # Eval box with scroll bar
        frm_scan2 = tk.Frame(frame)
        frm_scan2.pack(fill=tk.X,pady=[5,20])
        lbl_or = tk.Label(frm_scan2, text='Scans = ',font=SF)
        lbl_or.pack(side=tk.LEFT,padx=2,pady=[5,20])
        
        self.txt_scan = tk.Text(frm_scan2,width=65,wrap=tk.WORD,height=5)
        self.txt_scan.insert(tk.INSERT,scanstr)
        scl_scan = tk.Scrollbar(frm_scan2)
        scl_scan.config(command=self.txt_scan.yview)
        self.txt_scan.config(yscrollcommand=scl_scan.set)
        self.txt_scan.pack(side=tk.LEFT,padx=0)
        scl_scan.pack(side=tk.LEFT,fill=tk.Y)
        
        btn_adv_fit = tk.Button(frm_scan2, text='Adv. Fitting',wraplength=42,
                                font=BF, command=self.f_adv_fit)
        btn_adv_fit.pack(side=tk.LEFT,fill=tk.Y,padx=10)
        
        
        "------------------------------Buttons-------------------------------"
        frm_space = tk.Frame(frame)
        
        
        frm_btn0 = tk.Frame(frame)
        frm_btn0.pack(fill=tk.X)
        
        btn_plt1 = tk.Button(frm_btn0, text='Plot Scans',font=BF, width=23, height=2, command = self.f_btn_plot)
        btn_plt1.pack(side=tk.RIGHT, padx=2,pady=2)
        
        frm_fpt = tk.Frame(frm_btn0,borderwidth=2,relief=tk.RIDGE)
        frm_fpt.pack(side=tk.LEFT,fill=tk.X,expand=tk.YES)
        
        lbl_fpt = tk.Label(frm_fpt, text='Scan Plots: ',font=SF)
        lbl_fpt.pack(side=tk.LEFT,padx=5,pady=5)
        
        self.plot_indv = tk.IntVar(frm_btn0,0)
        chk_fpt = tk.Checkbutton(frm_fpt, text='multi-figure',font=SF,variable=self.plot_indv)
        chk_fpt.pack(side=tk.LEFT,padx=0)
        
        # Differentiate Plot
        self.diffplot = tk.IntVar(frm_btn0,0)
        chk_diff = tk.Checkbutton(frm_fpt, text='Diff.',font=BF,variable=self.diffplot)
        chk_diff.pack(side=tk.LEFT,padx=0)
        
        # Log Plot
        self.logplot = tk.IntVar(frm_btn0,0)
        chk_log = tk.Checkbutton(frm_fpt, text='Log',font=BF,variable=self.logplot)
        chk_log.pack(side=tk.LEFT,padx=0)
        
        # Normalise Plot
        self.normplot = tk.IntVar(frm_btn0,0)
        chk_log = tk.Checkbutton(frm_fpt, text='Eq. Area',font=BF,variable=self.normplot)
        chk_log.pack(side=tk.LEFT,padx=0)
        
        # Fit plots
        self.fitplots = tk.IntVar(frm_btn0,0)
        chk_ftp = tk.Checkbutton(frm_fpt, text='Fit',font=BF,variable=self.fitplots)
        chk_ftp.pack(side=tk.LEFT,padx=0)
        
        
        frm_btn1 = tk.Frame(frame)
        frm_btn1.pack()
        
        btn_plt2 = tk.Button(frm_btn1, text='Plot 3D',font=BF, width=19, height=2, command = self.f_btn_3D)
        btn_plt2.pack(side=tk.LEFT,padx=4,pady=2)
        
        btn_plt3 = tk.Button(frm_btn1, text='Plot 2D',font=BF, width=19, height=2, command = self.f_btn_2D)
        btn_plt3.pack(side=tk.LEFT,padx=4,pady=2)
        
        btn_plt4 = tk.Button(frm_btn1, text='Plot Surf',font=BF, width=19, height=2, command = self.f_btn_surf)
        btn_plt4.pack(side=tk.LEFT,padx=4,pady=2)
        
        
        frm_btn2 = tk.Frame(frame)
        frm_btn2.pack(fill=tk.X)
        
        btn_plt5 = tk.Button(frm_btn2, text='Fit Peaks',font=BF, width=23, height=2, command = self.f_btn_fit)
        btn_plt5.pack(side=tk.RIGHT, padx=2,pady=2)
        
        frm_fpt = tk.Frame(frm_btn2,borderwidth=2,relief=tk.RIDGE)
        frm_fpt.pack(side=tk.LEFT,fill=tk.X,expand=tk.YES)
        
        lbl_fpt = tk.Label(frm_fpt, text='Fit Plots: ',font=SF)
        lbl_fpt.pack(side=tk.LEFT,padx=5,pady=5)
        
        self.fitplot_all = tk.IntVar(frm_btn2,1)
        chk_fpt = tk.Checkbutton(frm_fpt, text='All',font=SF,variable=self.fitplot_all)
        chk_fpt.pack(side=tk.LEFT,padx=0)
        
        self.fitplot_ara = tk.IntVar(frm_btn2,0)
        chk_fpt = tk.Checkbutton(frm_fpt, text='Area',font=SF,variable=self.fitplot_ara)
        chk_fpt.pack(side=tk.LEFT,padx=0)
        
        self.fitplot_wid = tk.IntVar(frm_btn2,0)
        chk_fpt = tk.Checkbutton(frm_fpt, text='Width',font=SF,variable=self.fitplot_wid)
        chk_fpt.pack(side=tk.LEFT,padx=0)
        
        self.fitplot_cen = tk.IntVar(frm_btn2,0)
        chk_fpt = tk.Checkbutton(frm_fpt, text='Centre',font=SF,variable=self.fitplot_cen)
        chk_fpt.pack(side=tk.LEFT,padx=0)
        
        self.fitplot_2d = tk.IntVar(frm_btn2,0)
        chk_fpt = tk.Checkbutton(frm_fpt, text='2D',font=SF,variable=self.fitplot_2d)
        chk_fpt.pack(side=tk.LEFT,padx=0)
        
    
    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_fldr_browse(self):
        "Browse for data directory"
        inidir = self.filedir.get()
        dir = filedialog.askdirectory(initialdir=inidir)
        if len(dir) > 0:
            pp.filedir = dir
            self.filedir.set(dir)
    
    def f_save_browse(self):
        "Browse for data directory"
        inidir = self.savedir.get()
        dir = filedialog.askdirectory(initialdir=inidir)
        if len(dir) > 0:
            pp.savedir = dir
            self.savedir.set(dir)
    
    def f_scan_st(self):
        "Latest Scan number"
        try:
            num = pp.latest()
        except ValueError: 
            num = '000000'
            
        self.last.set(num)
    
    def f_scan_ld(self):
        "Update GUI for selected scan number"
        first = self.first.get()
        last = eval(self.last.get())
        step = self.step.get()
        
        if last > pp.latest(): 
            last = pp.latest()
            self.last.set(last)
        
        scans = range(first,last+1,step)
        
        # Update text in evalbox
        scanstr = str(scans)
        self.txt_scan.delete('1.0',tk.END)
        self.txt_scan.insert(tk.INSERT,scanstr)
        
    def f_scan_check(self):
        "Check the scan numbers"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        depvar = self.getdepvar()
        
        # Run checkscans
        outstr = pp.checkscans(scanno,showval=depvar)
        I16_Text_Display(outstr, 'Check Scans')
    
    def f_scan_makefile(self):
        "Create python analysis file"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get fit options
        depvar = self.depvar.get().split()
        yvar = self.vary.get()
        xvar = self.varx.get()
        save = self.saveopt.get()
        fit = self.fittype.get()
        test = self.peaktest.get()
        
        # Get plot options
        plot_opt = []
        if self.fitplot_all.get(): plot_opt += ['all']
        if self.fitplot_ara.get(): plot_opt += ['int']
        if self.fitplot_wid.get(): plot_opt += ['wid']
        if self.fitplot_cen.get(): plot_opt += ['cen']
        if self.fitplot_2d.get(): plot_opt += ['surface']
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        if save == 0: save = None
        
        pp.create_analysis_file(scanno, depvar=depvar, varx=xvar, vary=yvar, fit_type = fit, peaktest=test, norm=norm, plot=plot_opt, savePLOT=save)
    
    def f_adv_fit(self):
        "Start Advanced Fitting GUI"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get fit options
        depvar = self.getdepvar()
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        fit = self.fittype.get()
        test = self.peaktest.get()
        
        
        I16_Advanced_Fitting(scan_nos=scanno,ini_dependent=depvar,
                             ini_X=xvar,ini_Y=yvar,ini_fit=fit,
                             ini_Isig=test,ini_save=save)
    
    def f_btn_plot(self):
        "Button: Plot scans"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get plot options
        depvar = self.getdepvar()
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        log = self.logplot.get()
        normarea = self.normplot.get()
        fitplot = self.fitplots.get()
        dif = self.diffplot.get()
        indv_plot = self.plot_indv.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        if fitplot:
            fit = self.fittype.get()
        else:
            fit = None
        
        if indv_plot == True:
            for scn in scanno:
                pp.plotscan(scn, varx=xvar, vary=yvar, norm=norm, fit=fit, logplot=log, diffplot=dif, normalise=normarea, labels=depvar, save=save)
                plt.show()
        elif dif == False and normarea == False and len(scanno) > 7:
            pp.plotscans(scanno, depvar=depvar[0], varx=xvar, vary=yvar, norm=norm, fit=fit, logplot=log, save=save)
        else:
            pp.plotscan(scanno, varx=xvar, vary=yvar, norm=norm, fit=fit, fits=fitplot, logplot=log, diffplot=dif, normalise=normarea, labels=depvar, save=save)
        plt.show()
    
    def f_btn_3D(self):
        "Button: Plot 3D"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get plot options
        depvar = self.getdepvar()[0]
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        log = self.logplot.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        pp.plotscans3D(scanno,varx=xvar,vary=yvar,depvar=depvar,norm=norm,logplot=log,save=save)
        plt.show()
    
    def f_btn_2D(self):
        "Button: Plot 2D"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get plot options
        depvar = self.getdepvar()[0]
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        log = self.logplot.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        pp.plotscans2D(scanno,varx=xvar,vary=yvar,depvar=depvar,norm=norm,logplot=log,save=save)
        plt.show()
    
    def f_btn_surf(self):
        "Button: Plot Surf"
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get plot options
        depvar = self.getdepvar()[0]
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        log = self.logplot.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        pp.plotscansSURF(scanno,varx=xvar,vary=yvar,depvar=depvar,logplot=log,norm=norm,save=save)
        plt.show()
    
    def f_btn_fit(self):
        "Button: Fit Peaks"
        
        global fit,err
        
        # Set file directories
        pp.filedir = self.filedir.get()
        pp.savedir = self.savedir.get()
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scanstr=self.txt_scan.get('1.0','end-1c')
        scanno = eval(scanstr)
        
        # Get fit options
        depvar = self.getdepvar()
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        fit = self.fittype.get()
        test = self.peaktest.get()
        
        # Get plot options
        plot_opt = []
        if self.fitplot_all.get(): plot_opt += ['all']
        if self.fitplot_ara.get(): plot_opt += ['int']
        if self.fitplot_wid.get(): plot_opt += ['wid']
        if self.fitplot_cen.get(): plot_opt += ['cen']
        if self.fitplot_2d.get(): plot_opt += ['surface']
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        if save == 0: save = None
        
        print( 'fit,err = pp.fit_scans(scannos,vary={},depvar={},peaktest={},fit_type={},norm={},plot={},saveFIT={},savePLOT={})'.format(yvar,depvar,test,fit,norm,plot_opt,save,save) )
        fit,err = pp.fit_scans(scanno,varx=xvar,vary=yvar,depvar=depvar,
                               peaktest=test,fit_type=fit,norm=norm,plot=plot_opt,
                               saveFIT=save,savePLOT=save)
        plt.show()
        
    def getdepvar(self):
        "Reads the dependant variable and returns the values given"
        
        depvar = self.depvar.get()
        depvar = depvar.strip('[]')
        depvar = [x.strip() for x in depvar.split(',')]
        return depvar


"------------------------------------------------------------------------"
"-------------------------I16_Advanced_Fitting---------------------------"
"------------------------------------------------------------------------"
class I16_Advanced_Fitting:
    """
    More Fitting options, including masking
    Activated from "Adv. Fitting" button in I16_Peak_Analysis
    
    Takes scan numbers from I16_Peak_Analysis and allows navigation through
    all scans to check peak shapes etc.
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,scan_nos=[],ini_dependent=['Ta'],ini_X='',ini_Y='',ini_fit='pVoight',ini_Isig=1,ini_save=0):
        # Create Tk inter instance
        root = tk.Tk()
        root.wm_title('I16 Advanced Fitting by D G Porter [dan.porter@diamond.ac.uk]')
        root.minsize(width=700, height=600)
        root.maxsize(width=1200, height=700)
        
        if type(ini_dependent) is str:
            ini_dependent = [ini_dependent]
        
        # Initialise active array and masks array
        self.scan_nos = scan_nos
        self.Active = np.ones(len(scan_nos))
        self.Masks = [[] for n in range(len(scan_nos))]
        self.Current_Scan = 0
        self.depvals = []
        for n in range(len(scan_nos)):
            d = pp.readscan(scan_nos[n])
            if d is None:
                self.depvals += [np.nan]
            else:
                if ini_dependent[0] in d.metadata.keys():
                    self.depvals += [getattr(d.metadata,ini_dependent[0])]
                elif ini_dependent[0] in d.keys():
                    self.depvals += [np.mean(getattr(d,ini_dependent[0]))]
                else:
                    self.depvals += [n]
        
        # Update default save location for exported plots
        plt.rcParams["savefig.directory"] = pp.savedir
        
        frame = tk.Frame(root)
        frame.pack(side=tk.LEFT,anchor=tk.N)

        "----------------------------Data Directory-----------------------------"
        # Data Folder
        frm_fldr = tk.Frame(frame)
        frm_fldr.pack(fill=tk.X)
        lbl_fldr = tk.Label(frm_fldr, text='Data Folder: ',font=SF)
        lbl_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        ety_fldr = tk.Label(frm_fldr, text='{}'.format(pp.filedir), width=32)
        ety_fldr.pack(side=tk.LEFT,padx=5,pady=5)
        # Save Folder
        lbl_save = tk.Label(frm_fldr, text='Analysis Folder: ',font=SF)
        lbl_save.pack(side=tk.LEFT,padx=5,pady=5)
        ety_save = tk.Label(frm_fldr, text='{}'.format(pp.savedir), width=32)
        ety_save.pack(side=tk.LEFT,padx=5,pady=5)
        
        "-----------------------------Fit Options----------------------------"
        frm_titl = tk.Frame(frame)
        frm_titl.pack(fill=tk.X)
        
        # Title
        self.title = tk.StringVar(frm_titl,pp.exp_title)
        lbl_titl = tk.Label(frm_titl, text='Title: ',font=SF)
        lbl_titl.pack(side=tk.LEFT,padx=5,pady=5)
        ety_titl = tk.Entry(frm_titl, textvariable=self.title, width=47)
        ety_titl.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Dependent variable
        self.depvar = tk.StringVar(frm_titl,' '.join(ini_dependent))
        lbl_depv = tk.Label(frm_titl, text='Dependent: ',font=SF)
        lbl_depv.pack(side=tk.LEFT,padx=5,pady=5)
        ety_depv = tk.Entry(frm_titl, textvariable=self.depvar, width=10)
        ety_depv.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Plot range
        self.xrange = tk.StringVar(frm_titl,'[{:1.2f},{:1.2f}]'.format(np.min(self.depvals),np.max(self.depvals)))
        lbl_depr = tk.Label(frm_titl, text='Range: ',font=SF)
        lbl_depr.pack(side=tk.LEFT,padx=5,pady=5)
        ety_depr = tk.Entry(frm_titl, textvariable=self.xrange, width=20)
        ety_depr.pack(side=tk.LEFT,padx=5,pady=5)
        
        # 1st options row
        frm_opt = tk.Frame(frame)
        frm_opt.pack(fill=tk.X)
        
        # vary box
        self.varx = tk.StringVar(frm_opt,ini_X)
        lbl_varx = tk.Label(frm_opt, text='X: ',font=SF)
        lbl_varx.pack(side=tk.LEFT,padx=2,pady=5)
        ety_varx = tk.Entry(frm_opt, textvariable=self.varx, width=8)
        ety_varx.pack(side=tk.LEFT,padx=2,pady=5)
        
        # vary box
        self.vary = tk.StringVar(frm_opt,ini_Y)
        lbl_vary = tk.Label(frm_opt, text='Y: ',font=SF)
        lbl_vary.pack(side=tk.LEFT,padx=3,pady=5)
        ety_vary = tk.Entry(frm_opt, textvariable=self.vary, width=14)
        ety_vary.pack(side=tk.LEFT,padx=3,pady=5)
        
        # normalise menu
        normopts = ['rc','ic1','none']
        self.normtype = tk.StringVar(frm_opt, pp.normby)
        lbl_nrm = tk.Label(frm_opt, text='Norm: ', font=SF)
        lbl_nrm.pack(side=tk.LEFT,padx=0,pady=5)
        opt_nrm = tk.OptionMenu(frm_opt, self.normtype, *normopts)
        opt_nrm.config(width=4)
        opt_nrm.pack(side=tk.LEFT,padx=5,pady=5)
        
        # fit menu
        fitopts = ['Simple','Gauss','Lorentz','pVoight','Sum','Max']
        self.fittype = tk.StringVar(frm_opt, ini_fit)
        lbl_fit = tk.Label(frm_opt, text='Fit: ', font=SF)
        lbl_fit.pack(side=tk.LEFT,padx=0,pady=5)
        opt_fit = tk.OptionMenu(frm_opt, self.fittype, *fitopts)
        opt_fit.config(width=6)
        opt_fit.pack(side=tk.LEFT,padx=5,pady=5)
        
        bkgopts = ['Flat','Slope','Step']
        self.bkgtype = tk.StringVar(frm_opt, 'Flat')
        lbl_bkg = tk.Label(frm_opt, text='Bkg: ', font=SF)
        lbl_bkg.pack(side=tk.LEFT,padx=0,pady=5)
        opt_bkg = tk.OptionMenu(frm_opt, self.bkgtype, *bkgopts)
        opt_bkg.config(width=6)
        opt_bkg.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Peak Test
        self.peaktest = tk.DoubleVar(frm_opt,ini_Isig)
        lbl_peak = tk.Label(frm_opt, text='I/sig test: ',font=SF)
        lbl_peak.pack(side=tk.LEFT,padx=5,pady=5)
        ety_peak = tk.Entry(frm_opt, textvariable=self.peaktest, width=4)
        ety_peak.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Second Options row
        frm_opt2 = tk.Frame(frame)
        frm_opt2.pack(fill=tk.X)
        
        # Nloop
        self.Nloop = tk.IntVar(frm_opt2,10)
        lbl_loop = tk.Label(frm_opt2, text='Nloop: ', font=SF)
        lbl_loop.pack(side=tk.LEFT)
        ety_Nloop = tk.Entry(frm_opt2, textvariable=self.Nloop, width=4)
        ety_Nloop.pack(side=tk.LEFT,padx=5)
        
        # Binit
        self.Binit = tk.DoubleVar(frm_opt2,1e-5)
        lbl_Binit = tk.Label(frm_opt2, text='Binit: ', font=SF)
        lbl_Binit.pack(side=tk.LEFT)
        ety_Binit = tk.Entry(frm_opt2, textvariable=self.Binit, width=4)
        ety_Binit.pack(side=tk.LEFT,padx=5)
        
        # Tinc
        self.Tinc = tk.DoubleVar(frm_opt2,2)
        lbl_loop = tk.Label(frm_opt2, text='Tinc: ', font=SF)
        lbl_loop.pack(side=tk.LEFT)
        ety_Tinc = tk.Entry(frm_opt2, textvariable=self.Tinc, width=4)
        ety_Tinc.pack(side=tk.LEFT,padx=5)
        
        # Change
        self.Change = tk.DoubleVar(frm_opt2,0.5)
        lbl_loop = tk.Label(frm_opt2, text='Change: ', font=SF)
        lbl_loop.pack(side=tk.LEFT)
        ety_Change = tk.Entry(frm_opt2, textvariable=self.Change, width=4)
        ety_Change.pack(side=tk.LEFT,padx=5)
        
        # Converge
        self.Converge = tk.IntVar(frm_opt2,10)
        lbl_loop = tk.Label(frm_opt2, text='Converge: ', font=SF)
        lbl_loop.pack(side=tk.LEFT)
        ety_Converge = tk.Entry(frm_opt2, textvariable=self.Converge, width=4)
        ety_Converge.pack(side=tk.LEFT,padx=5)
        
        # Debug
        self.Debug = tk.IntVar(frm_opt2,0)
        chk_debug = tk.Checkbutton(frm_opt2, text='Debug',font=BF,variable=self.Debug,
                                  onvalue = 1, offvalue = 0)
        chk_debug.pack(side=tk.LEFT,padx=5)
        
        btn_fit = tk.Button(frm_opt2, text='Fit', font=BF, command = self.update_plot)
        btn_fit.pack(side=tk.LEFT, padx=5)
        
        # Make file
        btn_makefile = tk.Button(frm_opt2, text='Make File', font=BF, command = self.f_scan_makefile)
        btn_makefile.pack(side=tk.LEFT, padx=5)
        
        # Save check box
        self.saveopt = tk.IntVar(frm_opt2,ini_save)
        chk_save = tk.Checkbutton(frm_opt2, text='Save',font=BF,variable=self.saveopt,
                                  onvalue = 1, offvalue = 0)
        chk_save.pack(side=tk.LEFT,padx=0,pady=5)
        
        
        "-----------------------------Scan ListBox---------------------------"
        # Eval box with scroll bar
        frm_scan = tk.Frame(frame)
        frm_scan.pack(side=tk.LEFT,pady=[5,20])
        lbl_or = tk.Label(frm_scan, text='Scan No | {:10s} | Active | Masks '.format(ini_dependent[0]),font=SF)
        lbl_or.pack(side=tk.TOP,anchor=tk.NW,padx=2,pady=[5,1])
        
        scl_scanx = tk.Scrollbar(frm_scan,orient=tk.HORIZONTAL)
        scl_scanx.pack(side=tk.BOTTOM, fill=tk.X)
        
        scl_scany = tk.Scrollbar(frm_scan)
        scl_scany.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.lst_scan = tk.Listbox(frm_scan,width=40,height=30, font=SF, selectmode=tk.SINGLE,
                                xscrollcommand=scl_scanx.set,yscrollcommand=scl_scany.set)
        self.lst_scan.configure(exportselection=False)
        
        # Populate list box
        for n in range(len(scan_nos)):
            self.lst_scan.insert(tk.END,'{} | {:10.4g} | {} | {}'.format(scan_nos[n],self.depvals[n],int(self.Active[n]),str(self.Masks[n])))
        
        self.lst_scan.pack(side=tk.TOP,padx=5)
        self.lst_scan.select_set(0)
        self.lst_scan.bind("<<ListboxSelect>>", self.f_scan_select)
        #print( self.lst_scan.curselection()[0] )
        
        scl_scanx.config(command=self.lst_scan.xview)
        scl_scany.config(command=self.lst_scan.yview)
        
        #self.txt_scan.config(xscrollcommand=scl_scanx.set,yscrollcommand=scl_scany.set)
        
        
        "-------------------------------Scan Plot----------------------------"
        # Create frame on right hand side
        frm_rgt = tk.Frame(frame)
        frm_rgt.pack(side=tk.RIGHT, fill=tk.X, anchor=tk.NE, expand = tk.YES)
        
        #Create plot buttons
        frm_pbtn = tk.Frame(frm_rgt)
        frm_pbtn.pack(side=tk.LEFT, anchor=tk.NW)
        
        btn_plt1 = tk.Button(frm_pbtn, text='/\\',font=BF,height=6, command=self.f_plt_up)
        btn_plt1.pack()
        btn_plt2 = tk.Button(frm_pbtn, text='\\/',font=BF,height=6, command=self.f_plt_dn)
        btn_plt2.pack()
        
        # Create frame for plot
        frm_plt = tk.Frame(frm_rgt)
        frm_plt.pack(fill=tk.X,expand=tk.YES)
        
        self.fig1 = plt.Figure(figsize=[4,3],dpi=80)
        self.fig1.patch.set_facecolor('w')
        self.ax1 = self.fig1.add_subplot(111)
        self.ax1.set_autoscaley_on(True)
        self.ax1.set_autoscalex_on(True)
        self.plt1, = self.ax1.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5],'c+-',linewidth=1)
        self.plt2, = self.ax1.plot([],[],'bo-',linewidth=2) # marker point for pilatus
        self.pfit, = self.ax1.plot([],[],'r-',linewidth=2) # fit line
        self.ax1.set_xlabel('varx')
        self.ax1.set_ylabel('vary')
        self.ax1.set_title('Scan number',fontsize=16)
        
        self.fig1.subplots_adjust(left=0.25,bottom=0.2,right=0.95)
        
        # Change formats in x & y axes so they are nicer
        self.ax1.get_yaxis().set_major_formatter(mtick.FormatStrFormatter('%8.3g'))
        self.ax1.get_xaxis().get_major_formatter().set_useOffset(False)
        
        canvas = FigureCanvasTkAgg(self.fig1, frm_plt)
        canvas.get_tk_widget().configure(bg='black')
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, anchor=tk.NE, expand=tk.YES)
        #canvas.get_tk_widget().pack()
        #self.update_plot()
        
        "------------------------------Scan Values---------------------------"
        frm_val = tk.Frame(frm_rgt)
        frm_val.pack(side=tk.TOP, fill=tk.X)
        
        # Dependent value
        self.depval = tk.StringVar(frm_val,'{} = --'.format(ini_dependent[0]))
        lbl_depval = tk.Label(frm_val,textvariable=self.depval, font=SF)
        lbl_depval.pack(side=tk.LEFT)
        
        # Save check box
        self.scan_active = tk.IntVar(frm_val,1)
        chk_save = tk.Checkbutton(frm_val, text='Use this scan?',font=SF, command = self.f_chk_save,
                                  variable=self.scan_active,onvalue = 1, offvalue = 0)
        chk_save.pack(side=tk.RIGHT,padx=5,pady=5)
        
        "-----------------------------Fitted Values--------------------------"
        frm_fit = tk.Frame(frm_rgt)
        frm_fit.pack(side=tk.TOP, fill=tk.X)
        
        # Initilise frame variables
        self.Isig = tk.StringVar(frm_fit,'I/sig = --') # I/Sigma
        self.Height = tk.StringVar(frm_fit,'Height = --')
        self.Centre = tk.StringVar(frm_fit,'Centre = --')
        self.FWHM = tk.StringVar(frm_fit,'FWHM = --')
        self.Background = tk.StringVar(frm_fit,'Background = --')
        self.Area = tk.StringVar(frm_fit,'Area = --')
        
        # 3 stacked Frames, 2 variables per frame
        frm_fit1 = tk.Frame(frm_fit)
        frm_fit2 = tk.Frame(frm_fit)
        frm_fit3 = tk.Frame(frm_fit)
        
        frm_fit1.pack(fill=tk.X)
        frm_fit2.pack(fill=tk.X)
        frm_fit3.pack(fill=tk.X)
        
        lbl_Isig = tk.Label(frm_fit1, textvariable=self.Isig, font=SF)
        lbl_Height = tk.Label(frm_fit1, textvariable=self.Height, font=SF)
        lbl_Centre = tk.Label(frm_fit2, textvariable=self.Centre, font=SF)
        lbl_FWHM = tk.Label(frm_fit2, textvariable=self.FWHM, font=SF)
        lbl_Background = tk.Label(frm_fit3, textvariable=self.Background, font=SF)
        lbl_Area = tk.Label(frm_fit3, textvariable=self.Area, font=SF)
        
        lbl_Isig.pack(side=tk.LEFT,padx=5)
        lbl_Height.pack(side=tk.RIGHT,padx=5)
        lbl_Centre.pack(side=tk.LEFT,padx=5)
        lbl_FWHM.pack(side=tk.RIGHT,padx=5)
        lbl_Background.pack(side=tk.LEFT,padx=5)
        lbl_Area.pack(side=tk.RIGHT,padx=5)
        
        
        "----------------------------------Masks-----------------------------"
        frm_mask = tk.Frame(frm_rgt)
        frm_mask.pack(side=tk.TOP, fill=tk.X)
        
        frm_maskA = tk.Frame(frm_mask) # LEFT
        frm_maskB = tk.Frame(frm_mask) # RIGHT
        frm_mask1 = tk.Frame(frm_maskB) # RIGHT row 1
        frm_mask2 = tk.Frame(frm_maskB) # RIGHT row 2
        frm_mask3 = tk.Frame(frm_maskB) # RIGHT row 3
        
        frm_maskA.pack(side=tk.LEFT, fill=tk.Y)
        frm_maskB.pack(side=tk.LEFT, fill=tk.BOTH)
        frm_mask1.pack(fill=tk.X)
        frm_mask2.pack(fill=tk.X)
        frm_mask3.pack(fill=tk.X)
        
        # LEFT
        lbl_mask = tk.Label(frm_maskA, text = '---Masks---', font=HF)
        lbl_mask.pack()
        
        self.mask_all = tk.IntVar(frm_maskA,0)
        chk_mask_all = tk.Checkbutton(frm_maskA, text='All scans?', font=SF,
                                      variable=self.mask_all,onvalue = 1, offvalue = 0)
        chk_mask_all.pack()
        
        btn_remove = tk.Button(frm_maskA, text='Remove Masks', font=BF, command=self.f_mask_remove)
        btn_remove.pack()
        
        btn_print = tk.Button(frm_maskA, text='Print Masks', font=BF, command=self.f_mask_print)
        btn_print.pack()
        
        # RIGHT
        # x < val, x>val
        self.mask1 = tk.DoubleVar(frm_mask1,0)
        lbl_m1 = tk.Label(frm_mask1, text='x <  ',font=SF)
        lbl_m1.pack(side=tk.LEFT)
        ety_m1 = tk.Entry(frm_mask1, textvariable=self.mask1, width=4)
        ety_m1.pack(side=tk.LEFT)
        btn_m1 = tk.Button(frm_mask1, text='Mask', font=BF, command=self.f_mask_1)
        btn_m1.pack(side=tk.LEFT, padx = 5, pady = 5)
        
        self.mask2 = tk.DoubleVar(frm_mask1,0)
        lbl_m2 = tk.Label(frm_mask1, text='x > ',font=SF)
        lbl_m2.pack(side=tk.LEFT)
        ety_m2 = tk.Entry(frm_mask1, textvariable=self.mask2, width=4)
        ety_m2.pack(side=tk.LEFT)
        btn_m2 = tk.Button(frm_mask1, text='Mask', font=BF, command=self.f_mask_2)
        btn_m2.pack(side=tk.LEFT, padx = 5, pady = 5)
        
        # |x-val| < val, |x-val| > val
        self.mask3a = tk.DoubleVar(frm_mask2,0)
        self.mask3b = tk.DoubleVar(frm_mask2,0.1)
        lbl_m3a = tk.Label(frm_mask2, text='|x - ',font=SF)
        lbl_m3a.pack(side=tk.LEFT)
        ety_m3a = tk.Entry(frm_mask2, textvariable=self.mask3a, width=4)
        ety_m3a.pack(side=tk.LEFT)
        lbl_m3b = tk.Label(frm_mask2, text='| < ',font=SF)
        lbl_m3b.pack(side=tk.LEFT)
        ety_m3b = tk.Entry(frm_mask2, textvariable=self.mask3b, width=4)
        ety_m3b.pack(side=tk.LEFT)
        btn_m3 = tk.Button(frm_mask2, text='Mask', font=BF, command=self.f_mask_3)
        btn_m3.pack(side=tk.LEFT, padx = 5, pady = 5)
        
        self.mask4a = tk.DoubleVar(frm_mask3,0)
        self.mask4b = tk.DoubleVar(frm_mask3,0.1)
        lbl_m4a = tk.Label(frm_mask3, text='|x - ',font=SF)
        lbl_m4a.pack(side=tk.LEFT)
        ety_m4a = tk.Entry(frm_mask3, textvariable=self.mask4a, width=4)
        ety_m4a.pack(side=tk.LEFT)
        lbl_m4b = tk.Label(frm_mask3, text='| > ',font=SF)
        lbl_m4b.pack(side=tk.LEFT)
        ety_m4b = tk.Entry(frm_mask3, textvariable=self.mask4b, width=4)
        ety_m4b.pack(side=tk.LEFT)
        btn_m4 = tk.Button(frm_mask3, text='Mask', font=BF, command=self.f_mask_4)
        btn_m4.pack(side=tk.LEFT, padx = 5, pady = 5)
        
        "-------------------------------Start Fit----------------------------"
        frm_start = tk.Frame(frm_rgt)
        frm_start.pack(side=tk.TOP, fill=tk.X)
        
        btn_start = tk.Button(frm_start, text='Run Fit', font=BF, command=self.f_start)
        btn_start.pack(fill=tk.X)
        
        self.update_plot()
    
    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_scan_select(self,event):
        "Scan selection"
        
        self.update_plot()
    
    def f_plt_up(self):
        "Next plot up"
        
        current = self.lst_scan.curselection()[0]
        if current > 0:
            self.lst_scan.selection_clear(0,tk.END)
            self.lst_scan.activate(current-1)
            self.lst_scan.selection_set(current-1)
            self.update_plot()
    
    def f_plt_dn(self):
        "Next plot down"
        
        current = self.lst_scan.curselection()[0]
        if current < len(self.scan_nos)-1:
            self.lst_scan.selection_clear(0,tk.END)
            self.lst_scan.activate(current+1)
            self.lst_scan.selection_set(current+1)
            self.update_plot()
    
    def f_chk_save(self):
        "Activate or deactivate scan"
        
        current = self.lst_scan.curselection()[0]
        active = self.scan_active.get()
        self.Active[current] = active
        self.update_scans()
    
    def f_mask_remove(self):
        "Remove masks"
        
        if self.mask_all.get() > 0:
            # Remove masks from all scans
            self.Masks = [[] for n in range(len(self.scan_nos))]
        else:
            # Remove masks from selected scan
            current = self.lst_scan.curselection()[0]
            self.Masks[current] = []
        self.update_scans()
        self.update_plot()
    
    def f_mask_1(self):
        "Mask 1: x < val"
        
        val = self.mask1.get()
        mask = 'x < {}'.format(val)
        
        self.update_mask(mask)
    
    def f_mask_2(self):
        "Mask 2: x > val"
        
        val = self.mask2.get()
        mask = 'x > {}'.format(val)
        
        self.update_mask(mask)
    
    def f_mask_3(self):
        "Mask 3: |x-val1| < val2"
        
        val1 = self.mask3a.get()
        val2 = self.mask3b.get()
        mask = 'np.abs(x - {}) < {}'.format(val1,val2)
        
        self.update_mask(mask)
    
    def f_mask_4(self):
        "Mask 4: |x-val| > val2"
        
        val1 = self.mask4a.get()
        val2 = self.mask4b.get()
        mask = 'np.abs(x - {}) > {}'.format(val1,val2)
        
        self.update_mask(mask)
    
    def f_mask_print(self):
        "Print scans and masks"
        
        # Masks
        masks = np.array(self.Masks)
        masks = masks[np.where(self.Active)]
        mask_str = str(list(masks))
        
        # Run the Fit
        scannos = np.array(self.scan_nos)
        scannos = scannos[np.where(self.Active)]
        scan_str = str(list(scannos))
        
        print('scans = {}'.format(scan_str))
        print('masks = {}'.format(mask_str))
    
    def f_start(self):
        "Begin the fit"
        
        # Get the folders, titles etc.
        pp.normby = self.normtype.get()
        
        yvar = self.vary.get()
        xvar = self.varx.get()
        
        # Set Dependent value
        dependent = self.depvar.get().split()
        
        # Get Fitting values
        fit_type = self.fittype.get()
        bkg_type = self.bkgtype.get()
        peaktest = self.peaktest.get()
        Nloop = self.Nloop.get()
        Binit = self.Binit.get()
        Tinc = self.Tinc.get()
        Change = self.Change.get()
        Converge = self.Converge.get()
        Debug = self.Debug.get()
        save = self.saveopt.get()
        try:
            xrange = eval(self.xrange.get())
        except SyntaxError:
            xrange = None
        
        # Masks
        masks = np.array(self.Masks)
        masks = masks[np.where(self.Active)]
        
        # Run the Fit
        scannos = np.array(self.scan_nos)
        scannos = scannos[np.where(self.Active)]
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        for n in range(len(scannos)):
            print( n,scannos[n],masks[n] )
        
        pp.fit_scans(scannos,dependent,yvar,xvar,fit_type,bkg_type,peaktest,abscor=None,norm=norm,plot='all',show_fits=True,
                     mask_cmd=self.Masks,estvals=None,xrange=xrange,sortdep=True,Nloop=Nloop, Binit=Binit, Tinc=Tinc, 
                     change_factor=Change, converge_max=Converge, min_change=0.01,savePLOT=save,saveFIT=save)
        plt.show()
    
    def f_scan_makefile(self):
        "Create python analysis file"
        
        # Set file directories
        pp.exp_title = self.title.get()
        pp.normby = self.normtype.get()
        
        # Get scan numbers
        scannos = np.array(self.scan_nos)
        scannos = scannos[np.where(self.Active)]
        
        # Masks
        masks = np.array(self.Masks)
        masks = masks[np.where(self.Active)]
        
        # Get fit options
        depvar = self.depvar.get().split()
        xvar = self.varx.get()
        yvar = self.vary.get()
        save = self.saveopt.get()
        fit = self.fittype.get()
        test = self.peaktest.get()
        try:
            xrange = eval(self.xrange.get())
        except:
            xrange = None
        
        # Get Fitting values
        fit_type = self.fittype.get()
        bkg_type = self.bkgtype.get()
        peaktest = self.peaktest.get()
        Nloop = self.Nloop.get()
        Binit = self.Binit.get()
        Tinc = self.Tinc.get()
        Change = self.Change.get()
        Converge = self.Converge.get()
        Debug = self.Debug.get()
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        if save == 0: save = None
        
        pp.create_analysis_file(scannos,depvar,yvar,xvar,fit_type,bkg_type,peaktest,abscor=None,norm=norm,plot='all',show_fits=True,
                                mask_cmd=masks,estvals=None,xrange=xrange,sortdep=True,Nloop=Nloop, Binit=Binit, Tinc=Tinc, 
                                change_factor=Change, converge_max=Converge, min_change=0.01,savePLOT=save,saveFIT=save)
    
    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"
    def update_mask(self,mask):
        "Update the mask array, mask = string"
        
        if self.mask_all.get() > 0:
            # Apply mask to all scans
            for M in self.Masks:
                if mask not in M:
                    M += [mask]
        else:
            # Apply mask to this scan
            current = self.lst_scan.curselection()[0]
            self.Masks[current] += [mask]
        self.update_scans()
        self.update_plot()
    
    def update_scans(self):
        "Update the scan list and masks"
        
        current = self.lst_scan.curselection()[0]
        
        for n in range(len(self.scan_nos)):
            self.lst_scan.delete(n)
            self.lst_scan.insert(n,'{} | {:10.4g} | {} | {}'.format(self.scan_nos[n],self.depvals[n],int(self.Active[n]),str(self.Masks[n])))
        
        self.lst_scan.activate(current)
        self.lst_scan.selection_set(current)
    
    def update_plot(self):
        "Generate the plot, run a test fit and display the results"
        
        # Get the folders, titles etc.
        dependent = self.depvar.get().split()[0]
        pp.normby = self.normtype.get()
        
        current = self.lst_scan.curselection()[0]
        scanno = self.scan_nos[current]
        d = pp.readscan(scanno)
        if d is None:
            return
        
        if pp.normby == 'none':
            norm = False
        else:
            norm = True
        
        yvar = self.vary.get()
        xvar = self.varx.get()
        try:
            x,y,dy,varx,vary,ttl,d = pp.getdata(d,varx=xvar,vary=yvar,norm=norm)
        except AttributeError:
            x,y,dy,varx,vary,ttl,d = pp.getdata(d,varx='',vary='',norm=norm)
        
        # Mask the data
        x_mask,y_mask,dy_mask = pp.maskvals(x,y,dy,self.Masks[current])
        
        # Set Active value
        self.scan_active.set(int(self.Active[current]))
        
        # Set Dependent value
        if dependent in d.keys():
            depval = np.mean(getattr(d,dependent))
        elif dependent in d.metadata.keys():
            depval = getattr(d.metadata,dependent)
        else:
            depval = 0
        self.depval.set('{:s} = {:1.4g}'.format(dependent,depval))
        
        # Set I/Sig value
        Isig = pp.ispeak(y_mask,dy_mask,return_rat=True)
        self.Isig.set('I/sig = {:1.3g}'.format(Isig))
        
        # Update main Plot
        self.plt1.set_xdata(x)
        self.plt1.set_ydata(y)
        self.plt2.set_xdata(x_mask) # remove pilatus marker
        self.plt2.set_ydata(y_mask)
        self.ax1.relim()
        self.ax1.autoscale_view()
        
        self.ax1.set_xlabel(varx)
        self.ax1.set_ylabel(vary)
        self.ax1.set_title('#{}'.format(scanno,fontsize=16))
        
        # Get Fitting values
        fit_type = self.fittype.get()
        bkg_type = self.bkgtype.get()
        peaktest = self.peaktest.get()
        Nloop = self.Nloop.get()
        Binit = self.Binit.get()
        Tinc = self.Tinc.get()
        Change = self.Change.get()
        Converge = self.Converge.get()
        Debug = self.Debug.get()
        
        #print( fit_type,peaktest,Nloop,Binit,Tinc,Change,Converge,Debug )
        
        # Perform Fit
        out,err = pp.peakfit(x_mask,y_mask,dy_mask,type=fit_type,bkg_type=bkg_type,
                             peaktest=peaktest,Nloop=Nloop,Binit=Binit,Tinc=Tinc,
                             change_factor=Change,converge_max=Converge,
                             interpolate=True,debug=Debug)
        
        xfit,yfit=out['x'],out['y']
        amp,damp = out['Peak Height'],err['Peak Height']
        cen,dcen = out['Peak Centre'],err['Peak Centre']
        wid,dwid = out['FWHM'],err['FWHM']
        bkg,dbkg = out['Background'],err['Background']
        ara,dara = out['Area'],err['Area']
        
        # add to plot
        self.pfit.set_xdata(xfit)
        self.pfit.set_ydata(yfit)
        
        self.ax1.relim()
        self.ax1.autoscale_view()
        self.fig1.canvas.draw()
        
        # Update values
        self.Height.set('Height = {}'.format(pp.stfm(amp,damp)))
        self.Centre.set('Centre = {}'.format(pp.stfm(cen,dcen)))
        self.FWHM.set('FWHM = {}'.format(pp.stfm(wid,dwid)))
        self.Background.set('Background = {}'.format(pp.stfm(bkg,dbkg)))
        self.Area.set('Area = {}'.format(pp.stfm(ara,dara)))


"------------------------------------------------------------------------"
"---------------------------I16_Print_Buffer-----------------------------"
"------------------------------------------------------------------------"
class I16_Print_Buffer():
    """
    Adds images to a series of A4 pages allowing multiple plots to be printed
    on one page.
    Activated from the "Print All Figures" button in I16_Data_Viewer
    
    Usage:
        I16_Print_Buffer([fignos], [3,2])
          fignos = list of figure numbers, e.g [1,2,3 etc] or plt.get_fignums()
          ax_layout = [vertical figures on page, horizontal figures on page]
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,fignos=None,ax_layout=None):
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('I16 Print Buffer by D G Porter [dan.porter@diamond.ac.uk]')
        self.root.minsize(width=300, height=300)
        self.root.maxsize(width=1200, height=1100)
        
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,anchor=tk.N)
        
        if ax_layout is None:
            ax_layout = default_print_layout
        
        if fignos is None:
            fignos = plt.get_fignums()
        
        "----------------------------Plot Window---------------------------------"
        # Create frame for plot
        frm_plt = tk.Frame(frame)
        frm_plt.pack(fill=tk.X)
        
        " Find the temp directory"
        self.tmpdir = pp.tmpdir
        #ax_layout = [3,2] # [vertical, horizontal]
        Nax = ax_layout[0]*ax_layout[1]
        Nfigs = np.ceil( len(fignos)/float(Nax) ).astype(np.int)
        
        " Recursion of class to cover more than Nax figures"
        for nfig in range(1,Nfigs):
            buffer_fignos = fignos[nfig*Nax:nfig*Nax+Nax]
            I16_Print_Buffer(buffer_fignos,default_print_layout)
        
        " Save figures to temp directory"
        fignos = fignos[:Nax]
        for fn in fignos:
            fig = plt.figure(fn)
            fname = os.path.join(self.tmpdir,'Py16_figure%d.png' % fn)
            fig.savefig(fname)
        
        " Create the buffers"
        self.fig1 = plt.Figure(figsize=[8.27,11.69],dpi=80) #[horizontal, vertical]
        self.fig1.patch.set_facecolor('w')
        self.fig1.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
        
        for n,fn in enumerate(fignos[:Nax]):
            self.ax1 = self.fig1.add_subplot(ax_layout[0],ax_layout[1],n+1)
            fname = os.path.join(self.tmpdir,'Py16_figure%d.png' % fn)
            im = plt.imread(fname)
            self.ax1.imshow(im)
            self.ax1.axis('off')
            # Remove image after use
            os.remove(fname)
        
        canvas = FigureCanvasTkAgg(self.fig1, frm_plt)
        canvas.get_tk_widget().configure(bg='black')
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, anchor=tk.NE, expand=tk.YES)
        
        "----------------------------Print Button---------------------------------"
        # Create frame for button
        frm_btn = tk.Frame(frame)
        frm_btn.pack()
        
        # Button to print 
        btn_send = tk.Button(frm_btn, text='Print',font=BF, command=self.f_print)
        btn_send.pack(side=tk.LEFT,padx=2)
        
        # Button to save
        btn_send = tk.Button(frm_btn, text='Save',font=BF, command=self.f_save)
        btn_send.pack(side=tk.LEFT,padx=2)
    
    def f_print(self):
        " Print the buffer figure"
        fname = os.path.join(self.tmpdir,'Py16_buffer.pdf')
        self.fig1.savefig(fname,format='pdf',papertype='a4',dpi=300)
        print('Temp file created: %s'%fname)
        
        if print_platform is 'unix':
            print('Calling: lpr -o fit-to-page -r %s'%fname)
            subprocess.call(['lpr','-o','fit-to-page','-r',fname])
        else:
            # Windows
            print('Calling: print %s'%fname)
            win32api.ShellExecute (0, "print", fname, None, ".", 0)
        print( "Buffer Printed!" )
        # Remove file after printing
        try:
            os.remove(fname)
            print('Temp file removed')
        except:
            print('Couldn''t find %s'%fname)
            pass
        self.root.destroy()
    
    def f_save(self):
        " Save the buffer figure"
        opts = {'defaultextension':'.pdf',
                'filetypes':[('PDF','.pdf'),('PNG','.png')],
                'initialdir':pp.savedir,
                'initialfile':'I16_Print_Buffer',
                'parent':self.root,
                'title':'Save the print buffer:'}
        filename = filedialog.asksaveasfilename(**opts)
        self.fig1.savefig(filename,papertype='a4',dpi=300)
        print( "Buffer Saved!" )


"------------------------------------------------------------------------"
"-----------------------------I16_Check_Log------------------------------"
"------------------------------------------------------------------------"
class I16_Check_Log:
    """
    Functions to check multiple scans, check the log file and predict the end of a run
    Activated from the "More Check Options" button in I16_Data_Viewer
    
    Panel 1: Check scans
        - input the list of scan required in "Scans", this box is evaluated as python
        - Use the "show" box to add any metadata to the output
        - Click "Check Scans" - the output can be viewed in the console
    
    Panel 2: Check Log
        - Set either the Time or Scan number
        - Set the number of minutes (or tick "All time")
        - Click "Check Log" to show log file entries from "Mins" before Time/Scan until Time/Scan in the console
        - Use "Find Sting" to only show entries including this string
        - Tick "Show commands only" to only display commands
    
    Panel 3: Predict End
        - Set the number of the "First Scan" in a series
        - Set the expected "Last Scan" of a series (maths allowed)
        - Click "Predict End" to display the expexted end of the scan series in the console
        - Assumes each scan will be the same length on average.
        - Very useful for determining the end of a script
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self,ini_scan=0):
        # Create Tk inter instance
        root = tk.Tk()
        root.wm_title('Log Options')
        root.minsize(width=220, height=200)
        root.maxsize(width=300, height=1000)
        
        #Frame
        frame = tk.Frame(root)
        frame.pack(side=tk.LEFT,anchor=tk.N)
        
        "---------------------------Check Scan----------------------------"
        #Frame 1
        frame1 = tk.Frame(frame,borderwidth=2,relief=tk.RIDGE)
        frame1.pack(fill=tk.X,expand=tk.TRUE)
        
        # Line 1 title
        frame1t = tk.Frame(frame1)
        frame1t.pack()

        lbl_title = tk.Label(frame1t,text='Check Scans', font=BF, bd=1)
        lbl_title.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        # Line 1a
        frame1a = tk.Frame(frame1)
        frame1a.pack()
        
        # Scans
        self.scans_cmd = tk.StringVar(frame1,'range({},{},1)'.format(ini_scan-10,ini_scan))
        lbl_scncmd = tk.Label(frame1a,text='Scans:',font=SF)
        lbl_scncmd.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scncmd = tk.Entry(frame1a,textvariable=self.scans_cmd, width=20)
        ety_scncmd.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 1b
        frame1b = tk.Frame(frame1)
        frame1b.pack()
        
        # Scans
        self.scans_val = tk.StringVar(frame1,'')
        lbl_scnval = tk.Label(frame1b,text='Show:',font=SF)
        lbl_scnval.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scnval = tk.Entry(frame1b,textvariable=self.scans_val, width=20)
        ety_scnval.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 1c
        frame1c = tk.Frame(frame1)
        frame1c.pack()
        
        self.scans_typ = tk.StringVar(frame1,'')
        lbl_scnval = tk.Label(frame1c,text='Type:',font=SF)
        lbl_scnval.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scnval = tk.Entry(frame1c,textvariable=self.scans_typ, width=20)
        ety_scnval.pack(side=tk.LEFT,padx=5,pady=5)

        # Line 1d
        frame1d = tk.Frame(frame1)
        frame1d.pack()

        # Check Scans Button
        btn_scncmd = tk.Button(frame1d, text='Check Scans',font=BF,command=self.f_checkscans)
        btn_scncmd.pack(side=tk.LEFT)
        # Plot Meta Button
        btn_pltcmd = tk.Button(frame1d, text='Plot Data',font=BF,command=self.f_plotmeta)
        btn_pltcmd.pack(side=tk.LEFT)
        
        "---------------------------Check Log-----------------------------"
        #Frame 2
        frame2 = tk.Frame(frame,borderwidth=2,relief=tk.RIDGE)
        frame2.pack(fill=tk.X,expand=tk.TRUE)
        
        # Initial Values
        try:
            date = pp.readscan(ini_scan).metadata.date
        except AttributeError:
            date = 'Mon Sep 26 09:00:00 2016'
        now = datetime.datetime.strptime(date,'%a %b %d %H:%M:%S %Y')
        str_now = datetime.datetime.strftime(now,'%Y-%m-%d %H:%M')
        
        # Line 2 title
        frame2t = tk.Frame(frame2)
        frame2t.pack()

        lbl_title = tk.Label(frame2t,text='Check Log', font=BF, bd=1)
        lbl_title.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        # Line 2a
        frame2a = tk.Frame(frame2)
        frame2a.pack(fill=tk.X)
        
        # Time
        self.time = tk.StringVar(frame2,str_now)
        self.log_input1 = tk.IntVar(frame2,1)
        self.log_input2 = tk.IntVar(frame2,0)
        self.log_input = 'Time'
        chk_time = tk.Checkbutton(frame2a, text='Time:',font=SF,variable=self.log_input1,command=self.f_input_time)
        chk_time.pack(side=tk.LEFT,padx=5,pady=5)
        ety_time = tk.Entry(frame2a, textvariable=self.time, width=20)
        ety_time.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2b
        frame2b = tk.Frame(frame2)
        frame2b.pack(fill=tk.X)
        
        # Scan number
        self.scan = tk.IntVar(frame2,ini_scan)
        chk_scan = tk.Checkbutton(frame2b, text='Scan:',font=SF,variable=self.log_input2,command=self.f_input_scan)
        chk_scan.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scan = tk.Entry(frame2b, textvariable=self.scan, width=10)
        ety_scan.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2c
        frame2c = tk.Frame(frame2)
        frame2c.pack(fill=tk.X)
        
        # Mins
        self.mins = tk.IntVar(frame2,2)
        lbl_mins = tk.Label(frame2c, text='Mins: ',font=SF)
        lbl_mins.pack(side=tk.LEFT,padx=5,pady=5)
        ety_mins = tk.Entry(frame2c, textvariable=self.mins, width=5)
        ety_mins.pack(side=tk.LEFT,padx=5,pady=5)
        self.log_AllTime = tk.IntVar(frame2c,0)
        chk_all = tk.Checkbutton(frame2c, text='< All time',font=SF,variable=self.log_AllTime,onvalue = 1, offvalue = 0)
        chk_all.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2d
        frame2d = tk.Frame(frame2)
        frame2d.pack(fill=tk.X)
        
        # Find String
        self.find = tk.StringVar(frame2,'')
        lbl_find = tk.Label(frame2d, text='Find String: ',font=SF)
        lbl_find.pack(side=tk.LEFT,padx=5,pady=5)
        ety_find = tk.Entry(frame2d, textvariable=self.find, width=15)
        ety_find.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2e
        frame2e = tk.Frame(frame2)
        frame2e.pack(fill=tk.X)
        
        # Commands
        self.cmds = tk.IntVar(frame2,0)
        chk_cmds = tk.Checkbutton(frame2e, text='Show commands only',font=SF,variable=self.cmds,
                                  onvalue = 1, offvalue = 0)
        chk_cmds.pack(side=tk.LEFT,padx=0,pady=5)
        
        # Line 2f
        frame2f = tk.Frame(frame2)
        frame2f.pack(fill=tk.X)
        
        # Checklog Button
        btn_log = tk.Button(frame2f, text='Check Log',font=BF,command=self.f_checklog)
        btn_log.pack()
        
        
        "---------------------------Predict End----------------------------"
        #Frame 3
        frame3 = tk.Frame(frame,borderwidth=2,relief=tk.RIDGE)
        frame3.pack(fill=tk.X,expand=tk.TRUE)
        
        # Line 3 title
        frame3t = tk.Frame(frame3)
        frame3t.pack()

        lbl_title = tk.Label(frame3t,text='Predict End', font=BF, bd=1)
        lbl_title.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        # Line 3a
        frame3a = tk.Frame(frame3)
        frame3a.pack()
        
        # First scan
        self.first_scan = tk.IntVar(frame3,ini_scan)
        lbl_scan1 = tk.Label(frame3a,text='First Scan:',font=SF)
        lbl_scan1.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scan1 = tk.Entry(frame3a,textvariable=self.first_scan, width=20)
        ety_scan1.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 3b
        frame3b = tk.Frame(frame3)
        frame3b.pack()
        
        # Last scan
        self.last_scan = tk.StringVar(frame3,str(ini_scan) + '+10')
        lbl_scan2 = tk.Label(frame3b,text='Last Scan:',font=SF)
        lbl_scan2.pack(side=tk.LEFT,padx=5,pady=5)
        ety_scan2 = tk.Entry(frame3b,textvariable=self.last_scan, width=20)
        ety_scan2.pack(side=tk.LEFT,padx=5,pady=5)        
        # Line 1c
        frame3c = tk.Frame(frame3)
        frame3c.pack()
        
        # Checklog Button
        btn_scncmd = tk.Button(frame3c, text='Predict End',font=BF,command=self.f_prend)
        btn_scncmd.pack()
    
    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"
    def getdepvar(self):
        "Reads the dependant variable and returns the values given"
        
        depvar = self.scans_val.get()
        depvar = depvar.strip('[]')
        depvar = [x.strip() for x in depvar.split(',')]
        return depvar

    def f_checkscans(self):
        scans_command = self.scans_cmd.get()
        scans = eval(scans_command)
        
        showval = self.getdepvar() # allows multiple entries
        if len(showval[0]) == 0: showval = None

        find_scan = self.scans_typ.get()
        if len(find_scan) == 0: find_scan = None
        
        self.scans_cmd.set(str(scans))
        outstr = pp.checkscans(scans,showval=showval,find_scans=find_scan)
        I16_Text_Display(outstr, 'Check Scans {}-{}'.format(scans[0],scans[-1]))

    def f_plotmeta(self):
        scans_command = self.scans_cmd.get()
        scans = eval(scans_command)
        
        showval = self.getdepvar() # allows multiple entries
        if len(showval[0]) == 0: showval = 'TimeSec'
        
        self.scans_cmd.set(str(scans))
        pp.plotmeta(scans, fields=showval, use_time=False)
    
    def f_checklog(self):
        time = self.time.get()
        scan = self.scan.get()
        mins = self.mins.get()
        find = self.find.get()
        cmds = self.cmds.get()
        
        if self.log_input == 'Time':
            time = datetime.datetime.strptime(time,'%Y-%m-%d %H:%M')
            scan = datetime.datetime.strftime(time,'%Y-%m-%d %H:%M:%S,%f')
        
        if self.log_AllTime.get():
            mins = 'all'
        
        if len(find) == 0:
            find = None
        
        outstr = pp.checklog(scan, mins, cmds, find)
        I16_Text_Display(outstr, 'Check Log')
    
    def f_input_time(self):
        self.log_input1.set(1)
        self.log_input2.set(0)
        self.log_input = 'Time'
    
    def f_input_scan(self):
        self.log_input1.set(0)
        self.log_input2.set(1)
        self.log_input = 'Scan'
    
    def f_prend(self):
        fs = self.first_scan.get()
        if self.last_scan.get():
            ls = eval(self.last_scan.get())
            self.last_scan.set(str(ls))
        else:
            ls=None
        
        outstr = pp.prend(fs,ls)
        I16_Text_Display(outstr, 'Predict End')


"------------------------------------------------------------------------"
"--------------------------------I16_Params------------------------------"
"------------------------------------------------------------------------"
class I16_Params:
    """
    Set the default experimental parameters
    Activated from the "Params" button in I16_Data_Viewer
        Ring Current    Standard ring current for current experiment for normalisation
        ic1             Standard ic1monitor value for current experiment for normalisation
        Temp Sensor     Temperature sensor used as default
        Errors          Error funciton used (select none to remove error bars)
        Central Pixel   Centre pixel [i,j] in area detectors
        Max Pixel       Any pixels larger than this are hot pixels and set to zero
        Max array size  Stops PC form running out of memory on large pilatus scans
        Peak Region     Peak finding (nroi_peak) searches within this box
        Colours         Set the default colour order for plots
        Title           Set the title that will appear in plot titles
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self):
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('Parameters')
        self.root.minsize(width=300, height=200)
        self.root.maxsize(width=400, height=500)
        
        dwid = 15
        
        #Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,anchor=tk.N)
        
        "---------------------------Normalisation----------------------------"
        #Frame 1
        frame1 = tk.Frame(frame,borderwidth=2,relief=tk.RIDGE)
        frame1.pack(fill=tk.X,expand=tk.TRUE)
        
        # Line 1a
        frame1a = tk.Frame(frame1)
        frame1a.pack()
        
        # Ring Current
        self.exp_ring_current = tk.DoubleVar(frame1,pp.exp_ring_current)
        lbl_ring = tk.Label(frame1a,text='Ring Current:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_ring.pack(side=tk.LEFT,padx=5,pady=5)
        ety_ring = tk.Entry(frame1a,textvariable=self.exp_ring_current, width=20)
        ety_ring.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 1b
        frame1b = tk.Frame(frame1)
        frame1b.pack()
        
        # Monitor
        self.exp_monitor = tk.DoubleVar(frame1,pp.exp_monitor)
        lbl_mon = tk.Label(frame1b,text='ic1:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_mon.pack(side=tk.LEFT,padx=5,pady=5)
        ety_mon = tk.Entry(frame1b,textvariable=self.exp_monitor, width=20)
        ety_mon.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 1c
        frame1c = tk.Frame(frame1)
        frame1c.pack()
        
        # Errors
        tempopts = ['Ta','Tb','Tc','Td','tgas']
        self.sensortype = tk.StringVar(frame1, pp.default_sensor)
        lbl_sen = tk.Label(frame1c,text='Temp Sensor:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_sen.pack(side=tk.LEFT,padx=5,pady=5)
        opt_sen = tk.OptionMenu(frame1c, self.sensortype, *tempopts)
        opt_sen.config(width=dwid)
        opt_sen.pack(side=tk.LEFT,padx=0,pady=5)
        
        # Line 1d
        frame1d = tk.Frame(frame1)
        frame1d.pack()
        
        # Errors
        erroropts = ['sqrt(n)','std(N)','none']
        self.errortype = tk.StringVar(frame1, erroropts[0])
        lbl_err = tk.Label(frame1d,text='Errors:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_err.pack(side=tk.LEFT,padx=5,pady=5)
        opt_nrm = tk.OptionMenu(frame1d, self.errortype, *erroropts)
        opt_nrm.config(width=dwid)
        opt_nrm.pack(side=tk.LEFT,padx=0,pady=5)
        
        
        "---------------------------Pilatus Params-----------------------------"
        #Frame 2
        frame2 = tk.Frame(frame,borderwidth=2,relief=tk.RIDGE)
        frame2.pack(fill=tk.X,expand=tk.TRUE)
        
        # Line 2a
        frame2a = tk.Frame(frame2)
        frame2a.pack(fill=tk.X)
        
        # Pil_centre
        self.pil_centre = tk.StringVar(frame2,str(pp.pil_centre))
        lbl_cen = tk.Label(frame2a,text='Central Pixel:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_cen.pack(side=tk.LEFT,padx=5,pady=5)
        ety_cen = tk.Entry(frame2a,textvariable=self.pil_centre, width=20)
        ety_cen.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2b
        frame2b = tk.Frame(frame2)
        frame2b.pack(fill=tk.X)
        
        # Hot Pixel
        self.hot_pixel = tk.StringVar(frame2,str(pp.hot_pixel))
        lbl_hot = tk.Label(frame2b,text='Max Pixel:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_hot.pack(side=tk.LEFT,padx=5,pady=5)
        ety_hot = tk.Entry(frame2b,textvariable=self.hot_pixel, width=20)
        ety_hot.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2c
        frame2c = tk.Frame(frame2)
        frame2c.pack(fill=tk.X)
        
        # Max pilatus array size
        self.pil_max_size = tk.StringVar(frame2,'%1.3g'%pp.pil_max_size)
        lbl_max = tk.Label(frame2c,text='Max array size:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_max.pack(side=tk.LEFT,padx=5,pady=5)
        ety_max = tk.Entry(frame2c,textvariable=self.pil_max_size, width=20)
        ety_max.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 2d
        frame2d = tk.Frame(frame2)
        frame2d.pack(fill=tk.X)
        
        # Peak Region
        self.peakregion = tk.StringVar(frame2,str(pp.peakregion))
        lbl_pkrg = tk.Label(frame2d,text='Peak Region: [>y,>x,<Y,<X]',font=SF, justify=tk.RIGHT, width=dwid, wraplength=200)
        lbl_pkrg.pack(side=tk.LEFT,padx=5,pady=5)
        ety_pkrg = tk.Entry(frame2d,textvariable=self.peakregion, width=20)
        ety_pkrg.pack(side=tk.LEFT,padx=5,pady=5)
        
        
        "---------------------------Plotting----------------------------"
        #Frame 3
        frame3 = tk.Frame(frame,borderwidth=2,relief=tk.RIDGE)
        frame3.pack(fill=tk.X,expand=tk.TRUE)
        
        # Line 3a
        frame3a = tk.Frame(frame3)
        frame3a.pack()
        
        # Print layout
        self.print_layout = tk.StringVar(frame3,str(default_print_layout))
        lbl_pntlyt = tk.Label(frame3a,text='Print layout:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_pntlyt.pack(side=tk.LEFT,padx=5,pady=5)
        ety_pntlyt = tk.Entry(frame3a,textvariable=self.print_layout, width=20)
        ety_pntlyt.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 3b
        frame3b = tk.Frame(frame3)
        frame3b.pack()
        
        # Plot colours
        self.plot_colors = tk.StringVar(frame3,str(pp.plot_colors))
        lbl_pltcol = tk.Label(frame3b,text='Colours:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_pltcol.pack(side=tk.LEFT,padx=5,pady=5)
        ety_pltcol = tk.Entry(frame3b,textvariable=self.plot_colors, width=20)
        ety_pltcol.pack(side=tk.LEFT,padx=5,pady=5)
        
        # Line 3b
        frame3c = tk.Frame(frame3)
        frame3c.pack()
        
        # Last scan
        self.exp_title = tk.StringVar(frame3,pp.exp_title)
        lbl_titl = tk.Label(frame3c,text='Title:',font=SF, justify=tk.RIGHT, width=dwid)
        lbl_titl.pack(side=tk.LEFT,padx=5,pady=5)
        ety_titl = tk.Entry(frame3c,textvariable=self.exp_title, width=20)
        ety_titl.pack(side=tk.LEFT,padx=5,pady=5)
        
        "---------------------------Button----------------------------"
        #Frame 4
        frame4 = tk.Frame(frame,borderwidth=2)
        frame4.pack(fill=tk.X,expand=tk.TRUE)
        
        btn_update = tk.Button(frame4,text='Update',fg='red',font=BF,command=self.f_update)
        btn_update.pack(fill=tk.X,expand=tk.TRUE)
        
    
    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_update(self):
        "Update parameters and exit"
        
        pp.exp_ring_current = self.exp_ring_current.get()
        pp.exp_monitor = self.exp_monitor.get()
        pp.pil_centre = eval(self.pil_centre.get())
        pp.hot_pixel = eval(self.hot_pixel.get())
        pp.pil_max_size = eval(self.pil_max_size.get())
        pp.peakregion = eval(self.peakregion.get())
        pp.plot_colors = eval(self.plot_colors.get())
        pp.exp_title = self.exp_title.get()
        pp.default_sensor = self.sensortype.get()
        
        # Print layout
        new_print_layout = eval(self.print_layout.get())
        default_print_layout[0] = new_print_layout[0]
        default_print_layout[1] = new_print_layout[1]
        
        # Errors
        ch_error = self.errortype.get()
        if ch_error == 'sqrt(n)':
            pp.error_func = lambda x: np.sqrt(np.abs(x)+0.1)
        elif ch_error == 'std(N)':
            pp.error_func = pp.rolling_fun
        else:
            pp.error_func = lambda x: 0*x
        
        # Save parameter file
        pp.exp_parameters_save()
        
        # Close window
        self.root.destroy()

"------------------------------------------------------------------------"
"----------------------------------cutoffs-------------------------------"
"------------------------------------------------------------------------"
class colour_cutoffs:
    """
    Change the vmin/vmax colormap limits of the current figure
    Activate form the console by typing:
        colour_cutoffs()
    """
    "------------------------------------------------------------------------"
    "--------------------------GUI Initilisation-----------------------------"
    "------------------------------------------------------------------------"
    def __init__(self):
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('Change clim')
        #self.root.minsize(width=300, height=200)
        #self.root.maxsize(width=400, height=500)
        
        ini_vmin, ini_vmax = plt.gci().get_clim()
        
        #Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.LEFT,anchor=tk.N)
        
        # GCF button
        frm_btn = tk.Button(frame, text='Get Current Figure', font=BF, command=self.f_gcf)
        frm_btn.pack(fill=tk.X)
        
        # increment setting
        f_inc = tk.Frame(frame)
        f_inc.pack(fill=tk.X)
        
        inc_btn1 = tk.Button(f_inc, text='1', font=BF, command=self.f_but1)
        inc_btn1.pack(side=tk.LEFT)
        
        inc_btn2 = tk.Button(f_inc, text='100', font=BF, command=self.f_but2)
        inc_btn2.pack(side=tk.LEFT)
        
        inc_btn3 = tk.Button(f_inc, text='1000', font=BF, command=self.f_but3)
        inc_btn3.pack(side=tk.LEFT)
        
        self.increment = tk.DoubleVar(f_inc,1.0)
        inc_ety = tk.Entry(f_inc, textvariable=self.increment, width=6)
        inc_ety.pack(side=tk.LEFT)
        
        # Upper clim
        f_upper = tk.Frame(frame)
        f_upper.pack(fill=tk.X)
        
        up_left = tk.Button(f_upper, text='<', font=BF, command=self.f_upper_left)
        up_left.pack(side=tk.LEFT)
        
        self.vmin = tk.DoubleVar(f_upper,ini_vmin)
        up_edit = tk.Entry(f_upper, textvariable=self.vmin, width=12)
        up_edit.bind('<Return>',self.update)
        up_edit.bind('<KP_Enter>',self.update)
        up_edit.pack(side=tk.LEFT,expand=tk.YES)
        
        up_right = tk.Button(f_upper, text='>', font=BF, command=self.f_upper_right)
        up_right.pack(side=tk.LEFT)
        
        # Lower clim
        f_lower = tk.Frame(frame)
        f_lower.pack(fill=tk.X)
        
        lw_left = tk.Button(f_lower, text='<', font=BF, command=self.f_lower_left)
        lw_left.pack(side=tk.LEFT)
        
        self.vmax = tk.DoubleVar(f_lower,ini_vmax)
        lw_edit = tk.Entry(f_lower, textvariable=self.vmax, width=12)
        lw_edit.bind('<Return>',self.update)
        lw_edit.bind('<KP_Enter>',self.update)
        lw_edit.pack(side=tk.LEFT,expand=tk.YES)
        
        lw_right = tk.Button(f_lower, text='>', font=BF, command=self.f_lower_right)
        lw_right.pack(side=tk.LEFT)
        
        # Update button
        frm_btn = tk.Button(frame, text='Update', font=BF, command=self.update)
        frm_btn.pack(fill=tk.X)
        
        
    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"
    def f_gcf(self):
        #fig = plt.gcf()
        #fig.canvas.manager.window.raise_()
        
        new_vmin, new_vmax = plt.gci().get_clim()
        self.vmin.set(new_vmin)
        self.vmax.set(new_vmax)
    
    def f_but1(self):
        self.increment.set(1.0)
    
    def f_but2(self):
        self.increment.set(100)
    
    def f_but3(self):
        self.increment.set(1e3)
    
    def f_upper_left(self):
        inc = self.increment.get()
        cur_vmin = self.vmin.get()
        self.vmin.set( cur_vmin - inc)
        self.update()
    
    def f_upper_right(self):
        inc = self.increment.get()
        cur_vmin = self.vmin.get()
        self.vmin.set( cur_vmin + inc)
        self.update()
    
    def f_lower_left(self):
        inc = self.increment.get()
        cur_vmax = self.vmax.get()
        self.vmax.set( cur_vmax - inc)
        self.update()
    
    def f_lower_right(self):
        inc = self.increment.get()
        cur_vmax = self.vmax.get()
        self.vmax.set( cur_vmax + inc)
        self.update()
    
    def update(self,event=None):
        cur_vmin = self.vmin.get()
        cur_vmax = self.vmax.get()
        
        #fig = plt.gcf()
        #fig.canvas.manager.window.raise_()
        plt.clim(cur_vmin,cur_vmax)
        plt.show()


if __name__ == '__main__':
    print( 'I16_Data_Viewer By Dan Porter' )
    print( 'Launching GUI window...' )
    print( 'To restart the window, type: I16_Data_Viewer()' )
    print( 'Or to just load a scan, type: d = pp.readscan(scan_number)' )
    print( 'See "Py16Notes.txt" or type help(pp) for more info' )
    
    
    pp.recall_last_exp()
    dv=I16_Data_Viewer()
    #colour_cutoffs()
    
