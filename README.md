Py16 programs for data viewing and analysis of I16 data
By Dan Porter, 2016, Version 2.0

These programs are designed to be run in DAWN, though in principle should work in any python environment.

To run in a console:
	- Ensure you have Anaconda (or similar) Python 2.7 installed
	- Copy Py16progs.py and Py16GUI.py to a folder
	- Run a console/terminal/command prompt from this folder (linux: Right-click > console Windows: Shift+RightClick > command window  Mac: ?)
	- type: python py16GUI.py


To run these programs in DAWN on an I16 workstation (as I16user):
    - Open a conole
    - type: module load dawn
    - type: dawn
    - in DAWN, open the PyDev perspective (Window>Open Perspective>Pydev)
    - Open a PyDev interactive Console (Window>Show View>Console, then on the upper right of the console tab, press the button Open Console>PyDev Console 
    - Upon starting, the console should start the Py16 graphical user interface

To run these programs elsewhere:
    - You need Anaconda (or similar) Python 2.7
    - You must ensure you have tkinter set as the matplotlib backend (this is usually default)
    - Run the Py16GUI.py script.


Using Py16 programs:
You can either use the graphical user interface or use the functions from the Py16 module directly.

Py16 Data Viewer:
    - Start by adding your experiment data path and the path of the directory you would like to save files to
    - Type the scan number you would like to look at, or press the button "0" to get the latest scan
    - the x and y axis of plots are automatically choosen based on the scan command, to change the Y axis, type the required psudo device into the "Y axis" box and press "Plot"
    - Press "Pilatus" to load the set of pilatus images from the current scan, you can scroll through them using the arrows.
    - To save the current plot, you can export the plot to a new figure using "Export Plot" or you can save it directly to the analysis folder with "Save Plot"
    - Press "multiplot/ peak analysis" to go to the peak analysis window

Py16 functions:
    - Py16progs.py is a module of many useful functions, it is loaded automatically when you open a console in DAWN and named 'pp'
    - in the Console, type pp. and you should see a list of all the functions in the module
    - type help(pp.functionname) to see the allowed inputs, outputs and useage
Examples:
    d = readscan(534299) # Loads the data from this scan number into a python dictionary d
    x,y,dy,varx,vary,ttl,d = getdata(534299) # Automatically loads x,y data from scan, using the scan command to determine scanned values (x) and pseudo device (y), performs automatic normalisation
    plotscan(534299) # Automatically plots the data
    plotscan([534299,534300]) # plots multiple scans
    plotscan(534299,fit='pVoight') # plots with fitting
