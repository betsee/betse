#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import numpy as np
from matplotlib import colors


# Method that generates a custom colormap (seems to work to produce mpl cmap object)

def make_cmap(color_set, cmap_name):

    new_cmap = colors.LinearSegmentedColormap.from_list(cmap_name, color_set, N=256)

    return new_cmap

# Custom colors in right format

magenta = np.array([239, 52, 236] ) /255 # nice magenta!

cyan = np.array([9, 232, 239] ) /255 # nice cyan!

orange = np.array([255, 164, 61] ) /255 # nice orange!

pale_blue = np.array([56, 132, 255] ) /255 # nice deep blue!

pale_red = np.array([244, 66, 66] ) /255 # nice gentle red!

gold = np.array([255, 231, 55] ) /255 # golden yellow!
yellow = np.array([255, 246, 130] ) /255

salmon = np.array([255, 111, 54] ) /255 # salmon orange/red!
salmon2 = np.array([255, 117, 71] ) /255 # salmon orange/red version 2!

aqua =  np.array([53, 255, 211] ) /255  # aqua!
aqua2 = np.array([71, 255, 218] ) /255 # aqua version 2!

light_green = np.array([184, 255, 104 ] ) /255 # light green!
light_purple = np.array([219, 104, 255] ) /255 # light purple!

# standard colours:
black = np.array([0, 0, 0])
green = np.array([0, 1, 0])
red = np.array([1, 0, 0])
blue = np.array([0, 0, 1])
grey = np.array([0.2, 0.2, 0.2])

# Definition of custom colormaps

cm_electric_green = make_cmap([black, green], 'betse_electric_green')
cm_electric_magenta = make_cmap([black, magenta], 'betse_electric_magenta')
cm_electric_orange = make_cmap([black, salmon2], 'betse_electric_orange')
cm_electric_gold = make_cmap([black, yellow], 'betse_electric_gold')
cm_electric_cyan = make_cmap([black, cyan], 'betse_electric_cyan')
cm_electric_blue = make_cmap([black, pale_blue], 'betse_electric_blue')

cm_electric_green2 = make_cmap([grey, green], 'betse_green_chalkboard')
cm_electric_magenta2 = make_cmap([grey, magenta], 'betse_magenta_chalkboard')
cm_electric_orange2 = make_cmap([grey, salmon], 'betse_orange_chalkboard')
cm_electric_gold2 = make_cmap([grey, gold], 'betse_gold_chalkboard')
cm_electric_cyan2 = make_cmap([grey, cyan], 'betse_cyan_chalkboard')
cm_electric_blue2 = make_cmap([grey, pale_blue], 'betse_blue_chalkboard')

cm_alien = make_cmap([salmon, black, aqua], 'betse_alien_solid')
cm_alien2 = make_cmap([salmon2, black, aqua2], 'betse_alien_pale')
cm_alien3 = make_cmap([salmon2, grey, aqua2], 'betse_alien_chalkboard')

cm_spring = make_cmap([light_purple, black, light_green], 'betse_purple_green_pale')
cm_spring2 = make_cmap([magenta, black, green], 'betse_purple_green_solid')
cm_spring3 = make_cmap([magenta, grey, green], 'betse_purple_green_chalkboard')

cm_rdbu = make_cmap([pale_blue, black, pale_red], 'betse_red_blue_pale')
cm_rdbu2 = make_cmap([blue, black, red], 'betse_red_blue_solid')
cm_rdbu3 = make_cmap([blue, grey, red], 'betse_red_blue_chalkboard')