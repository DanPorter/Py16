"""
GUI for viewing detector images

By Dan Porter
Feb 2021
"""

import sys, os
import numpy as np
from imageio import imread  # read Tiff images
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
try:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
except ImportError:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg

if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox
else:
    import tkinter as tk
    from tkinter import filedialog
    from tkinter import messagebox

# Fonts
TF = ["Times", 12]  # entry
BF = ["Times", 14]  # Buttons
SF = ["Times New Roman", 14]  # Title labels
MF = ["Courier", 8]  # fixed distance format
LF = ["Times", 14]  # Labels
HF = ["Courier", 12]  # Text widgets (big)
# Colours - background
bkg = 'snow'
ety = 'white'
btn = 'azure'  # 'light slate blue'
opt = 'azure'  # 'light slate blue'
btn2 = 'gold'
# Colours - active
btn_active = 'grey'
opt_active = 'grey'
# Colours - Fonts
txtcol = 'black'
btn_txt = 'black'
ety_txt = 'black'
opt_txt = 'black'
ttl_txt = 'black'
_figure_size = [8, 6]


class ImageGui:
    """
    A standalone GUI window that displays multiple images using a list of filenames
        ImageGui( file_list, name_list, title, initial_index)

    file_list: [] list of filenames
    name_list: [] list of strings for each filename (default=displays filename)
    title: str: title to display in GUI
    initial_index: int: first image to show
    """

    _increment = 1
    _increment_fast = 10

    def __init__(self, file_list, name_list=None, title='', initial_index=None, colormap=None, clim=None, logplot=False):
        """Initialise"""
        # Create Tk inter instance
        self.root = tk.Tk()
        self.root.wm_title('Image Display')
        # self.root.minsize(width=640, height=480)
        self.root.maxsize(width=self.root.winfo_screenwidth(), height=self.root.winfo_screenheight())
        self.root.tk_setPalette(
            background=bkg,
            foreground=txtcol,
            activeBackground=opt_active,
            activeForeground=txtcol)

        if initial_index is None:
            initial_index = len(file_list)//2
        if name_list is None:
            name_list = file_list
        if colormap is None:
            self._colormap = plt.get_cmap('viridis')
        else:
            self._colormap = colormap

        self.file_list = file_list
        self.name_list = name_list
        self._increment_fast = int(len(file_list)//self._increment_fast)

        image = imread(self.file_list[initial_index])

        if clim is None:
            if logplot:
                cmin = 0.1
                cmax = np.max(image)
            else:
                # Automatically determine caxis
                md = np.median(image)
                mx = np.max(image)
                cmax = md + 10 ** (0.7 * np.log10(mx - md))
                if cmax <= 0: cmax = 1
                cmax = int(cmax)
                cmin = 0
        else:
            cmin, cmax = clim

        # Create tkinter Frame
        frame = tk.Frame(self.root)
        frame.pack(side=tk.TOP, fill=tk.BOTH, expand=tk.YES)

        # variables
        self.name_text = tk.StringVar(frame, name_list[initial_index])
        self.index = tk.IntVar(frame, initial_index)
        self.logplot = tk.BooleanVar(frame, logplot)
        self.cmin = tk.DoubleVar(frame, cmin)
        self.cmax = tk.DoubleVar(frame, cmax)
        self.colormap = tk.StringVar(frame, self._colormap.name)
        all_colormaps = plt.colormaps()
        all_colormaps = ['viridis', 'Spectral', 'plasma', 'inferno', 'Greys', 'Blues', 'winter', 'autumn',
                         'hot', 'hot_r', 'hsv', 'rainbow', 'jet']

        # ---Image title---
        frm = tk.Frame(frame)
        frm.pack(fill=tk.X, expand=tk.YES, padx=3, pady=3)

        var = tk.Label(frm, text=title, font=SF, fg=ttl_txt)
        var.pack(pady=5)
        var = tk.Label(frm, textvariable=self.name_text, font=LF)
        var.pack(pady=3)

        # ---Figure window---
        frm = tk.Frame(frame)
        frm.pack(fill=tk.BOTH, expand=tk.YES)

        self.fig = plt.Figure(figsize=_figure_size, dpi=80)
        self.fig.patch.set_facecolor('w')
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_autoscaley_on(True)
        self.ax.set_autoscalex_on(True)
        self.ax.set_frame_on(False)
        self.current_image = self.ax.imshow(image, cmap=self._colormap)
        self.ax.set_position([0, 0, 1, 1])
        self.current_image.set_cmap(self._colormap)
        self.current_image.set_clim(self.cmin.get(), self.cmax.get())

        # Toolbar coordinates
        self.ax.format_coord = lambda x,y: "x:{:>5}, y:{:>5}".format(int(x), int(y))
        self.current_image.format_cursor_data = lambda z: "[I:{:>9.4g}]".format(z)

        """
        ROIcen = initial_pilcen
        ROIsize = [75, 67]
        pil_centre = initial_pilcen
        pil_size = [195, 487]
        idxi = np.array([ROIcen[0] - ROIsize[0] // 2, ROIcen[0] + ROIsize[0] // 2 + 1])
        idxj = np.array([ROIcen[1] - ROIsize[1] // 2, ROIcen[1] + ROIsize[1] // 2 + 1])
        self.pilp1, = self.ax2.plot(idxj[[0, 1, 1, 0, 0]], idxi[[0, 0, 1, 1, 0]], 'k-', linewidth=2)  # ROI
        self.pilp2, = self.ax2.plot([pil_centre[1], pil_centre[1]], [0, pil_size[0]], 'k:',
                                    linewidth=2)  # vertical line
        self.pilp3, = self.ax2.plot([0, pil_size[1]], [pil_centre[0], pil_centre[0]], 'k:',
                                    linewidth=2)  # Horizontal line
        self.pilp4, = self.ax2.plot([], [], 'r-', linewidth=2)  # ROI background
        self.pilp5, = self.ax2.plot([], [], 'y-', linewidth=2)  # Peak region
        self.ax.set_aspect('equal')
        self.ax.autoscale(tight=True)
        """

        canvas = FigureCanvasTkAgg(self.fig, frm)
        canvas.get_tk_widget().configure(bg='black')
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=tk.YES)

        #   ---Toolbar---
        frm = tk.Frame(frame)
        frm.pack(expand=tk.YES)

        # Add matplotlib toolbar under plot
        self.toolbar = NavigationToolbar2TkAgg(canvas, frm)
        self.toolbar.update()
        self.toolbar.pack(fill=tk.X, expand=tk.YES)

        #  ---Bottom part---
        #bottom = tk.Frame(frame)
        #bottom.pack(expand=tk.YES)

        #  ---Image Tools---
        side = tk.LabelFrame(frame, text='Image Tools', relief=tk.RIDGE)
        side.pack(side=tk.TOP, fill=tk.X, expand=tk.YES, padx=10, pady=3)

        #   ---Slider---
        frm = tk.Frame(side)
        frm.pack(expand=tk.YES, fill=tk.X, padx=3, pady=3)

        var = tk.Scale(frm, from_=0, to=len(self.file_list)-1, variable=self.index, font=BF,
                       sliderlength=60, orient=tk.HORIZONTAL, command=self.but_scale)
        #var.bind("<ButtonRelease-1>", self.but_scale)
        var.pack(side=tk.LEFT, fill=tk.X, expand=tk.YES, padx=5, pady=3)

        #   ---Image move buttons---
        frm = tk.Frame(side)
        frm.pack(expand=tk.YES, padx=3, pady=3)

        var = tk.Button(frm, text='<<', font=BF, command=self.but_left_fast,
                        bg=btn, activebackground=btn_active)
        var.pack(side=tk.LEFT)
        var = tk.Button(frm, text='<', font=BF, command=self.but_left,
                        bg=btn, activebackground=btn_active)
        var.pack(side=tk.LEFT)
        var = tk.Button(frm, text='>', font=BF, command=self.but_right,
                        bg=btn, activebackground=btn_active)
        var.pack(side=tk.LEFT)
        var = tk.Button(frm, text='>>', font=BF, command=self.but_right_fast,
                        bg=btn, activebackground=btn_active)
        var.pack(side=tk.LEFT)

        #  ---Color Tools---
        side = tk.LabelFrame(frame, text='Color Tools', relief=tk.RIDGE)
        side.pack(side=tk.TOP, fill=tk.X, expand=tk.YES, padx=10, pady=3)

        #   ---Colormaps---
        frm = tk.Frame(side)
        frm.pack(expand=tk.YES, padx=3, pady=3)

        var = tk.Label(frm, text='Colormap:', font=SF)
        var.pack(side=tk.LEFT)
        var = tk.OptionMenu(frm, self.colormap, *all_colormaps, command=self.optn_colormap)
        var.config(font=SF, width=10, bg=opt, activebackground=opt_active)
        var["menu"].config(bg=opt, bd=0, activebackground=opt_active)
        var.pack(side=tk.LEFT)

        #   ---Color limits---
        frm = tk.Frame(side)
        frm.pack(expand=tk.YES, padx=3, pady=3)

        var = tk.Label(frm, text='Clim  Min:', font=SF)
        var.pack(side=tk.LEFT)
        var = tk.Entry(frm, textvariable=self.cmin, font=TF, width=6, bg=ety, fg=ety_txt)
        var.pack(side=tk.LEFT)
        var.bind('<Return>', self.ent_clim)
        var.bind('<KP_Enter>', self.ent_clim)

        var = tk.Label(frm, text='  Max:', font=SF)
        var.pack(side=tk.LEFT)
        var = tk.Entry(frm, textvariable=self.cmax, font=TF, width=6, bg=ety, fg=ety_txt)
        var.pack(side=tk.LEFT)
        var.bind('<Return>', self.ent_clim)
        var.bind('<KP_Enter>', self.ent_clim)

        var = tk.Label(frm, text='  Log:', font=SF)
        var.pack(side=tk.LEFT)
        var = tk.Checkbutton(frm, variable=self.logplot, font=SF, command=self.chk_log)
        var.pack(side=tk.LEFT, padx=6)

        "-------------------------Start Mainloop------------------------------"
        self.load_image(initial_index)
        if not hasattr(sys, 'ps1'):
            # If not in interactive mode, start mainloop
            self.root.protocol("WM_DELETE_WINDOW", self.f_exit)
            self.root.mainloop()

    "------------------------------------------------------------------------"
    "--------------------------General Functions-----------------------------"
    "------------------------------------------------------------------------"

    def load_image(self, index):
        self.name_text.set(self.name_list[index])
        image = imread(self.file_list[index])
        self.current_image.set_data(image)
        self.current_image.set_cmap(self._colormap)
        if self.logplot.get():
            self.current_image.norm = LogNorm()
            if self.cmin.get() <= 0:
                self.cmin.set(0.1)
        else:
            self.current_image.norm = Normalize()
        self.current_image.set_clim(self.cmin.get(), self.cmax.get())
        self.toolbar.update()
        self.fig.canvas.draw()

    "------------------------------------------------------------------------"
    "---------------------------Button Functions-----------------------------"
    "------------------------------------------------------------------------"

    def but_scale(self, event=None):
        """Move scroll bar"""
        index = self.index.get()
        self.load_image(index)

    def but_left_fast(self):
        """Decrease image index 10%"""
        index = self.index.get() - self._increment_fast
        if index < 0:
            index = 0
        self.load_image(index)
        self.index.set(index)

    def but_left(self):
        """Decrease image index by 1"""
        index = self.index.get() - self._increment
        if index < 0:
            index = 0
        self.load_image(index)
        self.index.set(index)

    def but_right(self):
        """Increase image index by 1"""
        index = self.index.get() + self._increment
        if index >= len(self.file_list):
            index = len(self.file_list) - 1
        self.load_image(index)
        self.index.set(index)

    def but_right_fast(self):
        """Increase image index 10%"""
        index = self.index.get() + self._increment_fast
        if index >= len(self.file_list):
            index = len(self.file_list) - 1
        self.load_image(index)
        self.index.set(index)

    def optn_colormap(self, event=None):
        """Colormap option menu"""
        name = self.colormap.get()
        self._colormap = plt.get_cmap(name)

        index = self.index.get()
        self.load_image(index)

    def ent_clim(self, event=None):
        """Change clim"""
        cmin = self.cmin.get()
        cmax = self.cmax.get()
        self.current_image.set_clim(cmin, cmax)
        self.fig.canvas.draw()

    def chk_log(self, event=None):
        """Turn on log scale"""
        if self.logplot.get():
            # Log on
            self.cmin.set(0.1)
            self.cmax.set(int(np.max(self.current_image.get_array())))
        else:
            # log off
            self.cmin.set(0.1)
            image = self.current_image.get_array()
            md = np.median(image)
            mx = np.max(image)
            cmax = md + 10 ** (0.7 * np.log10(mx - md))
            if cmax <= 0: cmax = 1
            self.cmax.set(int(cmax))
        index = self.index.get()
        self.load_image(index)

    def f_exit(self):
        """Closes the current data window"""
        self.root.destroy()
