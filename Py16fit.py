# -*- coding: utf-8 -*-
"""
Module: I16 Peak Fitting programs "Py16fits.py"

By Dan Porter, PhD
Diamond
2016

Usage: 
******In Script********
include the following lines at the top of your script

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/userscripts/python') # location of Py16Progs
import Py16Progs as p16

# Use functions:


******In Console*******
Run the file in the current console



**********************
Functions:
    
    

Version 0.9
Last updated: 30/06/16

Version History:
30/06/16 0.9    Program created from Py16progs.py V1.6

###FEEDBACK### Please submit your bug reports, feature requests or queries to: dan.porter@diamond.ac.uk

@author: Dan Porter
I16, Diamond Light Source
2016
"""

"""
Ideas for the future:
 - 

"""

import os
import numpy as np
from scipy.optimize import curve_fit # Peak fitting

# Parameter Names (will be used in output dicts)
Pheight = 'Peak Height'
Pcentre = 'Peak Centre'
Pwid = 'FWHM'
Pbkg = 'Background'
Pfrac = 'Lorentz Frac'
Pslope = 'Gradient'
Pstep = 'Step'
Psum = 'Area'

"-----------------------------------------------------------"
"-----------------------Generic Classes---------------------"
"-----------------------------------------------------------"
class fitfcn():
    ""
    names = []
    name = ''
    text_function = ''
    inputs = []
    params = []
    
    def func(self):
        pass
    def estimate(self,):

class simple():
    "Simple summation - no fit"
    names = ['simple','sum','s']
    name = 'Simple'
    text_function = ''
    inputs = []
    params = [Pheight,Pcentre,Pwid,Pbkg,Psum]
    
    def func(self):
        pass
    def area(self,values,errors):
        ara, dara = values[Psum],errors[Psum]
        return ara,dara

"-----------------------------------------------------------"
"-------------------------Peak Classes----------------------"
"-----------------------------------------------------------"
class gauss():
    "Define Gaussian"
    "From http://fityk.nieto.pl/model.html"
    names = ['gauss','gaussian','g']
    name = 'Gauss'
    text_function = 'height*np.exp(-np.log(2)*((x-cen)/(FWHM/2))**2)'
    inputs = ['height','cen','FWHM'] # inputs to text_funciton
    params = [Pheight,Pcentre,Pwid,Psum] # output parameters
    values = [np.nan,np.nan,np.nan,0] # Default values
    errors = [0,0,0,0] # Default errors
    
    def func(self,x,height=1,cen=0,FWHM=0.5,bkg=0):
        return height*np.exp(-np.log(2)*((x-cen)/(FWHM/2))**2) + bkg
    def area(self,values,errors):
        height, dheight = values[Pheight] , errors[Pheight]
        FWHM, dFWHM = values[Pwid] , errors[Pwid]
        
        sig = FWHM/(2*np.sqrt(2*np.log(2))) # Gaussian sigma
        dsig = dFWHM/((2*np.sqrt(2*np.log(2))))
        ara = np.abs(height*sig*np.sqrt(2*np.pi))
        dara = ara*np.sqrt( (dheight/height)**2 + (dsig/sig)**2 )
        return ara,dara

class lorentz():
    "Define Lorentzian"
    "From http://fityk.nieto.pl/model.html"
    names = ['lorentz','lorz','lorentzian','l']
    name = 'Lorentz'
    text_function = 'height/(1 + ((x-cen)/(FWHM/2))**2 )'
    inputs = ['height','cen','FWHM'] # inputs to text_funciton
    params = [Pheight,Pcentre,Pwid,Psum] # output parameters
    values = [np.nan,np.nan,np.nan,0] # Default values
    errors = [0,0,0,0] # Default errors
    
    def func(self,x,height=1,cen=0,FWHM=0.5,bkg=0):
        return height/(1 + ((x-cen)/(FWHM/2))**2 ) + bkg
    def area(self,values,errors):
        height, dheight = values[Pheight] , errors[Pheight]
        FWHM, dFWHM = values[Pwid] , errors[Pwid]
        
        ara = np.pi*height*FWHM/2
        dara = ara*np.sqrt( (dheight/height)**2 + (dFWHM/FWHM)**2 )
        return ara,dara

class pvoight():
    "Define pseudo-Voight"
    "From http://fityk.nieto.pl/model.html"
    names = ['pvoight','voight','pv','v']
    name = 'pVoight'
    text_function = 'height*( LorFrac/( 1.0 + (2.0*(x-cen)/FWHM)**2 ) + (1.0-LorFrac)*np.exp( -np.log(2)*(2.*(x-cen)/FWHM)**2 ) )'
    inputs = ['height','cen','FWHM' ,'LorFrac'] # inputs to text_funciton
    params = [Pheight,Pcentre,Pwid,Pfrac,Psum] # output parameters
    values = [np.nan,np.nan,np.nan,np.nan,0] # Default values
    errors = [0,0,0,0,0] # Default errors
    
    def func(self,x,height=1,cen=0,FWHM=0.5,LorFrac=0.5,bkg=0):
        HWHM = FWHM/2.0
        ln2 = 0.69314718055994529
        pos = x-cen
        L = LorFrac/( 1 + (pos/HWHM)**2 )
        G = (1-LorFrac)*np.exp( -ln2*(pos/HWHM)**2 )
        return height*(G + L) + bkg
    def area(self,values,errors):
        height, dheight = values[Pheight] , errors[Pheight]
        FWHM, dFWHM = values[Pwid] , errors[Pwid]
        LorFrac, dLorFrac = values[Pfrac] , errors[Pfrac]
        # Calculated Voight area = Gaussian + Voight
        sig = FWHM/(2*np.sqrt(2*np.log(2))) # Gaussian sigma
        dsig = dFWHM/((2*np.sqrt(2*np.log(2))))
        Gara = np.abs(height*sig*np.sqrt(2*np.pi))
        Lara = np.pi*height*FWHM/2
        ara = LorFrac*Lara + (1-LorFrac)*Gara
        # Error on area
        dGara = Gara*np.sqrt( (dheight/height)**2 + (dsig/sig)**2 )
        dLara = Lara*np.sqrt( (dheight/height)**2 + (dFWHM/FWHM)**2 )
        dVara1= (1-LorFrac)*Gara*np.sqrt( (dLorFrac/(1-LorFrac))**2 + (dGara/Gara)**2 )
        dVara2= LorFrac*Lara*np.sqrt( (dLorFrac/LorFrac)**2 + (dLara/Lara)**2 )
        dara = np.sqrt( dVara1**2 + dVara2**2 )
        return ara,dara

"-----------------------------------------------------------"
"---------------------Background Functions------------------"
"-----------------------------------------------------------"
class flat():
    names = ['flat','bkg','normal']
    name = 'flat'
    text_function = 'bkg'
    inputs = ['bkg'] # inputs to text_funciton
    params = [Pbkg] # output parameters
    values = [0] # Default values
    errors = [0] # Default errors
    
    def func(self,x,bkg):
        return bkg
    def area(self,values,errors):
        return 0,0

class slope():
    names = ['slope','sloping']
    name = 'slope'
    text_function = 'x*slope + bkg'
    inputs = ['bkg','slope'] # inputs to text_funciton
    params = [Pbkg,Pslope] # output parameters
    values = [0,0] # Default values
    errors = [0,0] # Default errors
    
    def func(self,x,bkg,slope=0):
        return x*slope + bkg
    def area(self,values,errors):
        return 0,0

class step():
    names = ['step']
    name = 'step'
    text_function = 'np.append(bkg-step*np.ones(len(x)/2),bkg+step*np.ones(len(x)/2))'
    inputs = ['bkg','step'] # inputs to text_funciton
    params = [Pbkg,Pstep] # output parameters
    values = [0,0] # Default values
    errors = [0,0] # Default errors
    
    def func(self,x,bkg,step=0):
        return np.append(bkg-step*np.ones(len(x)/2),bkg+step*np.ones(len(x)/2))
    def area(self,values,errors):
        return 0,0

all_functions = [simple(),gauss(),lorentz(),pvoight(),flat(),slope(),step()]


"-----------------------------------------------------------"
"---------------------Estimate Functions--------------------"
"-----------------------------------------------------------"
def estimate_FWHM(x,y,interpolate=False):
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
    
    wid = abs(hfpos2-hfpos1)
    dwid = abs(x[1]-x[0])
    return wid,dwid

def estimate_bkg(y):
    "Estimate background"
    
    #bkgrgn = np.concatenate( (y[:len(x)//5],y[-len(x)//5:]) ) # background method 1 - wrong if peak is off centre
    #bkgrgn = np.percentile(y,range(0,20)) # background method 2 - average lowest 5th of data
    #bkg = np.mean(bkgrgn) 
    h,bin = np.histogram(y,10)
    bincen = (bin[1:] + bin[:-1]) / 2.0
    bkg = bincen[np.argmax(h)]
    
    dbkg = np.sqrt(np.abs(bkg)+0.1)
    return bkg,dbkg

def estimate_height(y,bkg=0):
    "Estimate peak height"
    height = np.max(y) - bkg
    dheight = np.sqrt(np.abs(height)+0.1)
    return height,dheight

def estimate_cen(x,y):
    "Estimate peak centre"
    srt = np.argsort(y)
    cen = np.average( x[ srt[ -len(x)//5: ] ] ,weights=y[ srt[ -len(x)//5: ] ])
    dcen = abs(x[1]-x[0])
    return cen,dcen

def estimate_sum(x,y,bkg=0):
    "Estimate peak area"
    scanwid = abs(x[-1]-x[0])
    ara = np.sum(y-bkg)*scanwid/len(x)
    dara = np.sqrt(np.sum(y))*scanwid/len(x)
    return ara,dara
    
"-----------------------------------------------------------"
"-----------------------Parameter Table---------------------"
"-----------------------------------------------------------"
def estimate_params(x,y,inputs):
    """
    Return estimates of the given parameters
    """
    
    wid, dwid = estimate_FWHM(x,y)
    cen, dcen = estimate_cen(x,y)
    bkg, dbkg = estimate_bkg(y)
    amp, damp = estimate_height(y,bkg)
    ara, dara = estimate_area(x,y,bkg)
    
    tab=[]
    tab+=[{'names': ['height','peak height','amp','h',Pheight],
           'name': Pheight,
           'min_val': np.std(y),
           'max_val': 5*amp,
           'est_val': amp}]
    tab+=[{'names': ['cen','centre','pos','peak centre','c',Pcentre],
           'name': Pcentre,
           'min_val': np.min(x),
           'max_val': np.max(x),
           'est_val': cen}]
    tab+=[{'names': ['wid','fwhm','peak width','w',Pwid],
           'name': Pwid,
           'min_val': abs(x[1]-x[0]),
           'max_val': 2*(max(x)-min(x)),
           'est_val': wid}]
    tab+=[{'names': ['lorz_frac','frac','lorz frac','lf',Pfrac],
           'name': Pfrac,
           'min_val': -0.5,
           'max_val': 2,
           'est_val': 0.5}]
    tab+=[{'names': ['bkg','background','b',Pbkg],
           'name': Pbkg,
           'min_val': -np.inf,
           'max_val': np.inf,
           'est_val': bkg}]
    tab+=[{'names': ['slope','m','gradient','grad','g',Pslope],
           'name': Pslope,
           'min_val': -np.inf,
           'max_val': np.inf,
           'est_val': 0}]
    tab+=[{'names': ['step',Pstep],
           'name': Pstep,
           'min_val': -np.inf,
           'max_val': np.inf,
           'est_val': y[0]-y[-1] }]
    tab+=[{'names': ['area','sum','int',Parea],
           'name': Parea,
           'min_val': 0,
           'max_val': np.inf,
           'est_val': ara}]
    
    
    min_val = -np.inf*np.ones(len(inputs))
    max_val = np.inf*np.ones(len(inputs))
    est_val = np.zeros(len(inputs))
    val_name = list(inputs) # copy inputs
    
    for n,param in enumerate(inputs):
        for tb in tab:
            if param.lower() in tb['names']:
                min_val[n] = tb['min_val']
                max_val[n] = tb['max_val']
                est_val[n] = tb['est_val']
                val_name[n] = tb['name']
                continue
    return min_val,max_val,est_val




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
        print ' ------Simple Fit:----- '
        print ' Amplitude = {0:10.3G} +/- {1:10.3G}'.format(amp,damp)
        print '    Centre = {0:10.3G} +/- {1:10.3G}'.format(cen,dcen)
        print '      FWHM = {0:10.3G} +/- {1:10.3G}'.format(wid,dwid)
        print 'Background = {0:10.3G} +/- {1:10.3G}'.format(bkg,dbkg)
        print '      Area = {0:10.3G} +/- {1:10.3G}'.format(ara,dara)
    
    return amp,cen,wid,bkg,ara,damp,dcen,dwid,dbkg,dara

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
        print 'Zero detected - adding 0.001 to x values'
        offset = 0.001
        x = x + offset
    if any(np.isnan(dy)):
        print 'Ignoring errors due to NaNs'
        dy=np.ones(len(y))
    
    # Handle zero intensities
    y[ y<0.01 ] = 0.01
    dy[ dy<0.01 ] = 0.01
    
    # Starting parameters
    if Tc is None:
        Tc = x[len(x)//2]
    beta = 0.5
    amp = np.mean(y[:len(y)//10])
    print Tc,beta,amp
    
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
        print 'Fit didn''t work: oh dear'
        return
    
    # Print Results
    if disp:
        print ' ------Order Parameter Fit:----- '
        print '        Tc = {0:10.3G} +/- {1:10.3G}'.format(Tc,dTc)
        print '      Beta = {0:10.3G} +/- {1:10.3G}'.format(beta,dbeta)
        print '       Amp = {0:10.3G} +/- {1:10.3G}'.format(amp,damp)
        print '     CHI^2 = {0:10.3G}'.format(chi)
        print '  CHI^2 per free par = {0:10.3G}'.format(chinfp)
    return Tc,beta,amp,dTc,dbeta,damp,yfit

def ispeak(Y,dY=None,test = 1,disp=False,return_rat=False):
    "Determines whether a peak exists in the given dataset"
    
    if dY is None:
        dY = error_func(Y)
    
    "From Blessing, J. Appl. Cryst. (1997). 30, 421-426"
    "EQU: (1) + (6)"
    " Background estimation added by me"
    s = np.mean(Y)
    bkg = np.min(Y)
    wi = 1/dY
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
        print 'avg: ',s
        print 'bkg: ',bkg
        print 'signal: ',signal
        print 'error: ',err
        print 'rat: ',rat
    if return_rat:
        return rat
    return rat > test


def choose_fit(fit_name):
    """
    Returns peak function as text string and inputs
    """
    for dd in all_functions:
        if fit_name.lower() in dd.names:
            return dd
    return all_functions[0]

def create_fun(txt_fun,inputs):
    funcstr = 'def func(x,{}):\n    return {}'.format(','.join(inputs),txt_fun)
    exec funcstr
    return func

def create_peak_fun(peak,background):
    params = peak.params + background.params
    inp1 = ','.join(peak.inputs)
    inp2 = ','.join(background.inputs)
    funcstr = 'def func(x,{},{}):\n    return {} + {}'.format(inp1,inp2,peak.text_function,background.text_function)
    exec funcstr
    return func, params

def fcn_generator(*args):
    params = []
    inputs = []
    txt_fcn = []
    for peak in args:
        params += peak.params
        inputs += peak.inputs
        txt_fcn += [peak.text_function]
    
    inp = ','.join(inputs)
    code = '+'.join(txt_fcn)
    funcstr = 'def func(x,{}):\n    return {}'.format(inp,code)
    exec funcstr
    return func, params

def fit(x,y,dy,fitfunc,estvals):
    """
    Run SciPy.curve_fit and handle output
        x,y,dy should be arrays of data to fit
        fitfunc should be the function to fit
        estvals should be the initial estimates of fitfinc inputs
        
        y' = fitfunc(x,a,b,c)
        estvals = [ini_a, ini_b, ini_c]
        fitvals = [fnl_a, fnl_b, fnl_c]
        errvals - [err_a, err_b, err_c]
        
        If fit fails: 
            fitvals = estvals
            errvals = [nan, nan, nan]
    """
    
    try:
        fitvals, covmat = curve_fit(fitfunc,x,y,estvals,sigma=dy,absolute_sigma=True)
    except RuntimeError:
        fitvals = 1*estvals
        covmat = np.nan*np.eye(len(estvals))
    errvals = np.sqrt(np.diag(covmat))
    return fitvals, errvals

def chisq(yfit,y,dy=1):
    " Calculate CHI^2"
    return np.sum( (y-yfit)**2 / dy)

def peakfit2(x,y,dy=None,type='pVoight',bkg_type='flat',peaktest=1,
            Nloop=100,Binit=1e-5,Tinc=2,change_factor=0.2,conv_val = 10,
            min_change=0.01,interpolate=False,debug=False,disp=False):
    """ General Peak Fitting function to fit a profile to a peak in y = f(x)
    Allows several possible profiles to be used and can try to find the best estimates for 
    fitting parameters using an RMC-based least-squares routine.
    
    out,err = peakfit(x,y)
    out,err = peakfit(x,y,dy=None,type='pVoight',Nloop=0,Binit=1e-3,Tinc=2,change_factor=0.1,conv_val=5,interpolate=False,disp=False)
    
    Basic parameters:
        x = array of the dependent variable, e.g. eta, mu, phi
        y = array of the independent variable, e.g. maxval,roi1_sum, APD
        dy = errors on y (default = None)
        type = function type. Allowed: 'pVoight' (default), 'Gauss', 'Lorentz', 'Simple'*
    RMC options:
        Nloop = Number of iterations per temperature, default = 0 (RMC off)**
        Binit = Initial values of 1/kbT used for RMC, default = 1e-3 (lower = Higher temp)
        Tinc = After Nloop steps, the temperature is increased by factor Tinc, default = 2
        change_factor = Each parameter is multiplied by a normal distribution around 1 with width change_factor. (Default = 0.1)
    Output options:
        interpolate = True: The output fit will have interpolated (much finer) values in x and y. (Default = False)
        disp = True: The final fitted parameters will be displayed in the command line. (Dafault = False)
    
    Output:
        out = dict with fitted parameters
        err = dict with errors on fitted paramters
        
        out.keys() = ['Peak Height','Peak Centre','FWHM','Lorz frac','Background','Area','CHI**2','CHI2 per dof','x','y']
    
    * selecting type='simple' will not fit the data, just provide a very simple estimation.
    ** Nloop must be set > 0 for the RMC routine to be used, for Nloop=0, a simple gradient decend method from a simple estimation is used.
    
    Notes on the RMC routine:
     - see the code
    
    """
    
    # Set dy to 1 if not given
    if dy is None: dy=np.ones(len(y))
    
    # Remove zeros from x - causes errors in covariance matrix
    xold = 1.0*x
    offset = 0.
    if any(np.abs(x)<0.0001):
        print 'Zero detected - adding 0.0001 to x values'
        offset = 0.0001
        x = x + offset
    if any(np.isnan(dy)):
        print 'Ignoring errors due to NaNs'
        dy=np.ones(len(y))
    
    # Handle zero intensities
    y[ y<0.01 ] = y[ y<0.01 ]+0.01
    dy[ dy<0.01 ] = dy[ dy<0.01 ]+0.01
    
    # Select the peak and background functions
    peak = choose_fit(type)
    background = choose_fit(bkg_type)
    
    '-----------------------------------------------------------'
    '-------------------------FIT DATA--------------------------'
    '-----------------------------------------------------------'
    # Fitting not reuqired
    if peak.name == 'Simple':
        wid, dwid = estimate_FWHM(xold,y)
        cen, dcen = estimate_cen(xold,y)
        bkg, dbkg = estimate_bkg(y)
        amp, damp = estimate_height(y,bkg)
        ara, dara = estimate_area(xold,y,bkg)
        fitvals = [amp,cen,wid,bkg,ara]
        errvals = [damp,dcen,dwid,dbkg,dara]
        chi=0
        
    # Perform fitting
    else:
        # Check if a peak exists to fit
        peak_rat = ispeak(y,dy,test=peaktest,disp=False,return_rat=True)
        if debug: print 'Peak ratio: {:1.2g} ({:1.2g})'.format(peak_rat,peaktest)
        if peak_rat < peaktest:
            if debug: print 'No peak here (rat={:1.2g}). Fitting background instead!'.format(peak_rat)
            
            
            
            estvals = [0,bkg]
            minvals = [-np.inf,-np.inf]
            maxvals = [np.inf,np.inf]
        
        # Perform fitting
        # Create the fit function
        fitfunc, valnames = create_peak_fun(peak,background)
        
        # Estimate starting parameters
        minvals,maxvals,estvals = estimate_params(x,y,valnames)
        
        # Initial Fit (but don't update the estimators yet)
        try:
            fitvals, covmat = curve_fit(fitfunc,x,y,estvals,sigma=dy,absolute_sigma=True)
        except RuntimeError:
            if debug: print 'Initial fit failed!'
            fitvals = 1*estvals
            covmat = np.nan*np.eye(len(estvals))
        yfit = fitfunc(xold,*fitvals) # New curve
        chi = np.sum( (y-yfit)**2 / dy) # Calculate CHI^2
        if debug: print 'Initial Fit CHI**2 = ',chi
        
        # Check errors are reasonable
        errvals = np.sqrt(np.diag(covmat))
        if any(np.isnan(errvals)):
            chi = np.inf
        
        # Check new values are reasonable
        for n,val in enumerate(fitvals):
            if val < minvals[n] or val > maxvals[n]:
                if debug: print 'Initial value out of range: {} = {} ({}:{})'.format(valnames[n],val,minvals[n],maxvals[n])
                chi = np.inf # will not accept change if fitvalues fall out of range
        
        if debug: print 'Estimates: ',estvals
        if debug: print 'Initial Fit: ',fitvals,'Chi = ',chi
        
        changes = np.zeros(len(estvals))
        converge = 0
        Ntemp = 0
        while converge < conv_val:
            beta = Binit*Tinc**Ntemp
            if debug: print 'New Temperature: ',Ntemp,beta
            Ntemp += 1
            if Ntemp > Nloop:
                break
            for MCloop in range(Nloop):
                ini_estvals = 1*estvals # 1*estvals copies the array rather than links to it!
                # Loop over each estimator and randomly vary it
                for estn in range(len(estvals)):
                    inc_factor = np.random.normal(1,change_factor)
                    est_new = 1*estvals
                    est_new[estn] = est_new[estn]*inc_factor
                    try:
                        fitvals, covmat = curve_fit(fitfunc,x,y,est_new,sigma=dy,absolute_sigma=True)
                    except RuntimeError:
                        if debug: print beta,MCloop,estn,'Fit failed.'
                        continue
                    yfit = fitfunc(xold,*fitvals) # New curve
                    chi_new = np.sum( (y-yfit)**2 / dy) # Calculate CHI^2
                    
                    # Check errors are reasonable
                    errvals = np.sqrt(np.diag(covmat))
                    if any(np.isnan(errvals)):
                        chi_new = np.inf
                    
                    # Check new values are reasonable
                    for n,val in enumerate(fitvals):
                        #if debug: print beta,MCloop,estn,'CheckVal: ',n,val,minvals[n],maxvals[n]
                        if val < minvals[n] or val > maxvals[n]:
                            if debug: print 'Value out of range: {} = {} ({}:{})'.format(valnames[n],val,minvals[n],maxvals[n])
                            chi_new = np.inf # will not accept change if fitvalues fall out of range
                    if debug: print beta,MCloop,estn,'Vals:',fitvals,'Chi^2 = ',chi_new
                    
                    # Metropolis Algorithm
                    if chi_new < chi or np.exp(beta*(chi-chi_new)) > np.random.rand():
                        estvals = 1*fitvals # = 1*est_new
                        chi = 1*chi_new
                        changes[estn] += 1
                
                # Track changes
                chvals = np.divide(np.abs(np.subtract(estvals,ini_estvals)),ini_estvals)
                if np.any(chvals > min_change):
                    converge = 0
                else:
                    converge += 1
                
                if debug: print beta,MCloop,chi,'Changes: ',changes,chvals,converge
                
                # break the loop if the solution has converged
                if converge >= conv_val:
                    if debug: print 'Fit converged in {} temps!'.format(Ntemp-1)
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
        print 'Fit didnt work: use summation instead'
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
    ara,dara = peak.area(output,outerr)
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
        print ' ------{}+{} Fit:----- '.format(peak.name,background.name)
        for estn in range(len(fitvals)):
            print '{0:10s} = {1:10.3G} +/- {2:10.3G}'.format(valnames[estn],fitvals[estn],errvals[estn])
        print '      Area = {0:10.3G} +/- {1:10.3G}'.format(ara,dara)
        print '     CHI^2 = {0:10.8G}'.format(chi)
        print '  CHI^2 per free par = {0:10.3G}'.format(chinfp)
    return output,outerr

"""
x = np.arange(-5,5,0.1)
pk = choose_fit('pVoight')
bk = choose_fit('slope')
fitfunc, valnames = create_peak_fun(pk,bk)
y = fitfunc(x,5,0,1,0.8,20,0.2) + np.random.rand(len(x))
dy = np.sqrt(y)

out,err = peakfit2(x,y,dy,type='pVoight',bkg_type='slope',disp=True)

plt.figure()
plt.plot(x,y,label='Peak')
plt.plot(out['x'],out['y'],label='Fit')
plt.legend(loc=0)
"""