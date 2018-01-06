#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import matplotlib.pyplot as plt

class PolyPicker(object):
    """
    Allows the user to interactively select cell polygons from a graph.
    Once selected, polygons turn bright blue.
    The process is complete when the figure window is closed.
    Returns a list of indices to cells, which will have properties
    changed (e.g. ion channels added) depending on how the polypicker
    is used in the main script.

    """

    def __init__(self,cells,p,resume_after_picking):

        self.resume_after_picking = resume_after_picking
        self.picked_indices = []
        self.plotPatch(cells,p)
      #  self.fig.canvas.mpl_connect('pick_event', self.get_inds)
        self.cidpick = self.fig.canvas.mpl_connect('pick_event', self)
        self.cidclose = self.fig.canvas.mpl_connect('close_event', self.on_close)

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

        plt.show(block=False)

    def on_pick(self,event):

        patch = event.artist
        patch.set_alpha(1.0)
        self.index = patch.ind
        self.fig.canvas.draw()
        self.picked_indices.append(self.index)

    def on_close(self, event):
        self.fig.canvas.mpl_disconnect(self.cidclose)
        self.fig.canvas.mpl_disconnect(self.cidpick)
        self.resume_after_picking(self)

    def __call__(self, event):
            self.on_pick(event)



