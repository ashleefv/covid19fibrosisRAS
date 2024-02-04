# SVG (Cell Plot) Tab
import os
from ipywidgets import Layout, Label, Text, Checkbox, Button, HBox, VBox, Box, \
    FloatText, BoundedIntText, BoundedFloatText, HTMLMath, Dropdown, interactive, Output
from collections import deque
import xml.etree.ElementTree as ET
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as mplc
import numpy as np
import zipfile
import glob
import platform

from pyMCDS import pyMCDS


hublib_flag = True
if platform.system() != 'Windows':
    try:
        from hublib.ui import Download
    except:
        hublib_flag = False
else:
    hublib_flag = False


class SVGTab(object):

    def __init__(self):

        self.output_dir = 'output/'

        constWidth = '180px'

        max_frames = 1
        self.svg_plot = interactive(self.plot_svg, frame=(0, max_frames), continuous_update=False)
        plot_size = '500px'  # small: controls the size of the tab height, not the plot (rf. figsize for that)
        plot_size = '700px'  # medium
        plot_size = '750px'  # medium
        self.svg_plot.layout.width = plot_size
        self.svg_plot.layout.height = plot_size
        self.use_defaults = True
        self.show_nucleus = 0  # 0->False, 1->True in Checkbox!
        self.show_edge = 1  # 0->False, 1->True in Checkbox!
        self.scale_radius = 1.0
        self.axes_min = 0.0
        self.axes_max = 2000   # hmm, this can change (TODO?)

        self.max_frames = BoundedIntText(
            min=0, max=99999, value=max_frames,
            description='Max',
            layout=Layout(width='160px'),
        )
        self.max_frames.observe(self.update_max_frames)

        self.show_nucleus_checkbox= Checkbox(
            description='nucleus', value=False, disabled=False,
            layout=Layout(width=constWidth),
        )
        self.show_nucleus_checkbox.observe(self.show_nucleus_cb)

        self.show_edge_checkbox= Checkbox(
            description='edge', value=True, disabled=False,
            layout=Layout(width=constWidth),
        )
        self.show_edge_checkbox.observe(self.show_edge_cb)

        items_auto = [Label('select slider: drag or left/right arrows'), 
            self.max_frames, 
            self.show_nucleus_checkbox,  
            self.show_edge_checkbox, 
         ]

        box_layout = Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    width='70%')
        row1 = Box(children=items_auto, layout=box_layout)

        if (hublib_flag):
            self.download_button = Download('svg.zip', style='warning', icon='cloud-download', 
                                            tooltip='You need to allow pop-ups in your browser', cb=self.download_cb)
            download_row = HBox([self.download_button.w, Label("Download all cell plots (browser must allow pop-ups).")])

            self.tab = VBox([row1, self.svg_plot, download_row])
        else:
            self.tab = VBox([row1, self.svg_plot])

    def update(self, rdir=''):

        if rdir:
            self.output_dir = rdir

        all_files = sorted(glob.glob(os.path.join(self.output_dir, 'snapshot*.svg')))
        if len(all_files) > 0:
            last_file = all_files[-1]
            self.max_frames.value = int(last_file[-12:-4])  # assumes naming scheme: "snapshot%08d.svg"


    def download_cb(self):
        file_str = os.path.join(self.output_dir, '*.svg')

        with zipfile.ZipFile('svg.zip', 'w') as myzip:
            for f in glob.glob(file_str):
                myzip.write(f, os.path.basename(f))   # 2nd arg avoids full filename path in the archive

    def show_nucleus_cb(self, b):
        global current_frame
        if (self.show_nucleus_checkbox.value):
            self.show_nucleus = 1
        else:
            self.show_nucleus = 0
        self.svg_plot.update()

    def show_edge_cb(self, b):
        if (self.show_edge_checkbox.value):
            self.show_edge = 1
        else:
            self.show_edge = 0
        self.svg_plot.update()


    def update_max_frames(self,_b):
        self.svg_plot.children[0].max = self.max_frames.value

    #-----------------------------------------------------
    def circles(self, x, y, s, c='b', vmin=None, vmax=None, **kwargs):
        """
        See https://gist.github.com/syrte/592a062c562cd2a98a83 

        Make a scatter plot of circles. 
        Similar to plt.scatter, but the size of circles are in data scale.
        Parameters
        ----------
        x, y : scalar or array_like, shape (n, )
            Input data
        s : scalar or array_like, shape (n, ) 
            Radius of circles.
        c : color or sequence of color, optional, default : 'b'
            `c` can be a single color format string, or a sequence of color
            specifications of length `N`, or a sequence of `N` numbers to be
            mapped to colors using the `cmap` and `norm` specified via kwargs.
            Note that `c` should not be a single numeric RGB or RGBA sequence 
            because that is indistinguishable from an array of values
            to be colormapped. (If you insist, use `color` instead.)  
            `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
        vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with `norm` to normalize
            luminance data.  If either are `None`, the min and max of the
            color array is used.
        kwargs : `~matplotlib.collections.Collection` properties
            Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
            norm, cmap, transform, etc.
        Returns
        -------
        paths : `~matplotlib.collections.PathCollection`
        Examples
        --------
        a = np.arange(11)
        circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
        plt.colorbar()
        License
        --------
        This code is under [The BSD 3-Clause License]
        (http://opensource.org/licenses/BSD-3-Clause)
        """

        if np.isscalar(c):
            kwargs.setdefault('color', c)
            c = None

        if 'fc' in kwargs:
            kwargs.setdefault('facecolor', kwargs.pop('fc'))
        if 'ec' in kwargs:
            kwargs.setdefault('edgecolor', kwargs.pop('ec'))
        if 'ls' in kwargs:
            kwargs.setdefault('linestyle', kwargs.pop('ls'))
        if 'lw' in kwargs:
            kwargs.setdefault('linewidth', kwargs.pop('lw'))


        zipped = np.broadcast(x, y, s)
        patches = [Circle((x_, y_), s_)
                for x_, y_, s_ in zipped]
        collection = PatchCollection(patches, **kwargs)
        if c is not None:
            c = np.broadcast_to(c, zipped.shape).ravel()
            collection.set_array(c)
            collection.set_clim(vmin, vmax)

        ax = plt.gca()
        ax.add_collection(collection)
        ax.autoscale_view()
        plt.draw_if_interactive()
        if c is not None:
            plt.sci(collection)
        return collection


    #-------------------------
    def plot_svg(self, frame):

        # global current_idx, axes_max
        global current_frame
        current_frame = frame
        fname = "snapshot%08d.svg" % frame
        full_fname = os.path.join(self.output_dir, fname)
        if not os.path.isfile(full_fname):
            print("Once output files are generated, click the slider.")   
            return

        xlist = deque()
        ylist = deque()
        rlist = deque()
        rgb_list = deque()
        tree = ET.parse(full_fname)
        root = tree.getroot()

        numChildren = 0
        for child in root:
            if self.use_defaults and ('width' in child.attrib.keys()):
                self.axes_max = float(child.attrib['width'])
            if self.use_defaults and ('width' in child.attrib.keys()):
                self.axes_max = float(child.attrib['width'])
            if child.text and "Current time" in child.text:
                svals = child.text.split()
                title_str = svals[2] + "d, " + svals[4] + "h"

            if ('id' in child.attrib.keys()):
                tissue_parent = child
                break

        cells_parent = None

        for child in tissue_parent:
            if (child.attrib['id'] == 'cells'):
                cells_parent = child
                break
            numChildren += 1

        num_cells = 0
        for child in cells_parent:
            for circle in child:  # two circles in each child: outer + nucleus
                xval = float(circle.attrib['cx'])

                s = circle.attrib['fill']
                if (s[0:3] == "rgb"):  # if an rgb string, e.g. "rgb(175,175,80)" 
                    rgb = list(map(int, s[4:-1].split(",")))  
                    rgb[:] = [x / 255. for x in rgb]
                else:     # otherwise, must be a color name
                    rgb_tuple = mplc.to_rgb(mplc.cnames[s])  # a tuple
                    rgb = [x for x in rgb_tuple]

                # test for bogus x,y locations (rwh TODO: use max of domain?)
                too_large_val = 10000.
                if (np.fabs(xval) > too_large_val):
                    print("bogus xval=", xval)
                    break
                yval = float(circle.attrib['cy'])
                if (np.fabs(yval) > too_large_val):
                    print("bogus xval=", xval)
                    break

                rval = float(circle.attrib['r'])
                xlist.append(xval)
                ylist.append(yval)
                rlist.append(rval)
                rgb_list.append(rgb)

                # For .svg files with cells that *have* a nucleus, there will be a 2nd
                if (self.show_nucleus == 0):
                    break

            num_cells += 1

        xvals = np.array(xlist)
        yvals = np.array(ylist)
        rvals = np.array(rlist)
        rgbs = np.array(rgb_list)


        self.fig = plt.figure(figsize=(7, 6))

        ############# plot SUBSTRATES  #############################

        str_name = 'output{:08d}.xml'.format(frame)
        substrate_name = "collagen"

        mcds = pyMCDS(str_name, 'output/')

        x = np.arange(self.axes_min, self.axes_max, 20)
        y = np.arange(self.axes_min, self.axes_max, 20)
        X, Y = np.meshgrid(x, y)
        substrate= mcds.get_concentrations( substrate_name )

        plt.contourf(X, Y, substrate[:,:,0])  # plt.cm.BuPu_r,  plt.cm.viridis
        cbar1 = plt.colorbar()
        plt.xlabel('X-distance ($μm$)')
        plt.ylabel('Y-distance ($μm$)')
        cbar1.set_label('Collagen $μg/μm^3$')
     

        ############# plot SVG  #############################

        if (self.show_edge):
            try:
                self.circles(xvals,yvals, s=rvals, color=rgbs, edgecolor='black', linewidth=0.5)
            except (ValueError):
                pass
        else:
            self.circles(xvals,yvals, s=rvals, color=rgbs)

        plt.xlim(self.axes_min, self.axes_max)
        plt.ylim(self.axes_min, self.axes_max)
        plt.title(title_str)




