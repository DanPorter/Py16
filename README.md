# I16 Data Viewer

Py16 programs for data viewing and analysis of I16 data

By Dan Porter, 2020, Version 4.7

View scan and meta data from scan files generated at the beamline I16 at Diamond Light Source. Quickly look through scan data, plot the results with authomated normalisation and perform analysis such as peak fitting, generation of custom detector regions of interest (ROIs) and conversion of detector images to reciprocal space. Plus much more! Contains a feature rich graphical interface and a large module of automated functions for use in an interactive terminal or through scripting.

These programs workin work in almost any python environment, as long as you have numpy, scipy and matplotlib installed and matplotlib initialised in the tkinter backend. The software has been extensivley tested in Python 2.7+ and Python 3.7+.

Requirements: numpy, matplotlib, scipy, tkinter, imageio, h5py

See Py16Notes.txt for more info and Py16progs.html or Py16GUI.html for documentation

# Graphical interface
From a terminal/ command prompt containing Py16GUI.py and Py16progs.py:
```
ipython -i --matplotlib tk Py16GUI.py
```
Note that the terminal is still active meaning you can also interact manually with the data, using the Py16progs functions, accessed via 'pp.function'


# Interactive Console / Scripting
Example:
```
import Py16progs as pp
pp.filedir='directory/of/data'
pp.savedir='directory/of/analysis'

pp.plotscan(123456) \# plots the scan number 123456
d = pp.readscan(123456) \# creates a data holder object 'd' with scan data
x,y,dy,varx,vary,ttl,d = pp.getdata(123456) \# returns x/y data with automatic choice of variables
help(pp.getdata) \# see the documentation on this function
```
