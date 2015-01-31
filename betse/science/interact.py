#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME make the selection of points a blocking function...somehow...

import matplotlib.pyplot as plt

class PickObject(object):

    def __init__(self,cells,p):

        self.picked_indices = []
        self.plotPatch(cells,p)
      #  self.fig.canvas.mpl_connect('pick_event', self.get_inds)
        self.fig.canvas.mpl_connect('pick_event', self)

    def plotPatch(self,cells,p):

        self.fig = plt.figure()
        self.ax = plt.subplot(111)

        for i, poly in enumerate(cells.cell_verts):
            plygn = plt.Polygon(p.um*poly,facecolor='b')
            plygn.set_picker(True)
            plygn.set_alpha(0.25)
            plygn.ind = i
            self.ax.add_patch(plygn)
            self.ax.set_xlabel('Spatial distance x [um]')
            self.ax.set_ylabel('Spatial distance y [um]')
            self.ax.set_title('Select Cells with Mouse Click')

        self.ax.axis('equal')
        self.ax.autoscale_view()

        #self.fig.canvas.draw()
        plt.show(block=False)


    def get_inds(self,event):

        patch = event.artist
        patch.set_alpha(1.0)
        self.index = patch.ind
        self.fig.canvas.draw()
        self.picked_indices.append(self.index)
        print(self.picked_indices)

    def __call__(self, event):
        self.get_inds(event)

