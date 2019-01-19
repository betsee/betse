#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Post-simulation plot pipelining** (i.e., iteratively displaying and/or saving
plots produced after initialization and simulation) abstract base classes.
'''

#FIXME: I believe I've finally tracked down the issue relating to the following
#runtime "pyplot" warning:
#
#    pyplot.py:424: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
#
#The issue is that whenever we call a plotutil.plot* function below (e.g.,
#plotutil.plotSingleCellCData()), we localize the figure returned by that
#function. That figure will then be garbage collected on whichever of the
#following occurs last: (A) the corresponding local variable goes out of scope
#and (B) the corresponding plot window is closed by the user. Normally, neither
#would be a problem. Except that function is 1,400 lines long, which means that
#each figure's local variable effectively *NEVER* goes out of scope for the
#duration of plotting. Thus, figures will only be garbage collected *AFTER*
#this function terminates -- which is pretty much unacceptable.
#
#There are a couple solutions, thankfully. The simplest would simply be to stop
#localizing figures returned by plotutil.plot*() functions for all unused figure
#locals. The harder but probably more ideal solution would be to refactor all
#plotutil.plot*() functions to stop returning figures altogether. Since the current
#figure is *ALWAYS* accessible via the "matplotlib.pyplot.gcf()" getter (i.e.,
#[g]et[c]urrent[f]igure), there's no nead to explicitly return figures at all.
#
#In the unlikely event that we actually use a figure local, we need to either:
#
#* Extract all plotting logic related to that figure into a new helper function
#  of this module, ensuring that local will go out of scope immediately after
#  plotting that figure in that function.
#* Manually destroy that local with either "figure = None" or "del figure".
#
#Undomesticated unicorns running into the carefree sunset!
#FIXME: The above analysis for the following runtime "pyplot" warning is, sadly,
#completely wrong.
#
#    pyplot.py:424: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
#
#The underlying cause is simple. The pyplot.figure() function internally caches
#created figures into a global cache of such figures. There are many reasons why
#this is a bad idea, including the aforementioned warning as well as the
#probably non-thread-safety of this approach. The solution, of course, is to
#simply manually instantiate Figure() instances directly rather than call the
#pyplot.figure() function: e.g.,
#
#    # Instead of this...
#    pyplot.figure()
#    subplot = pyplot.subplot()
#
#    # ...do this.
#    fig = Figure()
#    subplot = fig.subplot()  # does this work? no idea.
#
#The obvious downside of the above approach, of course, is proper creation of
#non-blocking plots -- which will probably require creation and use of some
#sort of thread-safe cache. Assuming a one-to-one mapping is preserved between
#each non-blocking plot and each thread, the simplest mechanism would be to
#simply cache that plot's figure as an attribute of that thread. Sweet, no?

# ....................{ IMPORTS                           }....................
import matplotlib
from betse.lib.matplotlib import mplfigure, mplutil
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.science.phase.phasecls import SimPhase
from betse.science.pipe.export.pipeexpabc import SimPipeExportABC
from betse.util.io.log import logs
from betse.util.path import dirs, pathnames
from betse.util.type.types import type_check
from matplotlib import pyplot

# ....................{ SUBCLASSES                        }....................
class SimPipeExportPlotABC(SimPipeExportABC):
    '''
    Abstract base class of all **post-simulation plot pipelines** (i.e.,
    objects iteratively displaying and/or saving all plots produced after
    initialization and simulation enabled by the current simulation
    configuration).
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def init(self, phase: SimPhase) -> None:

        # Initialize our superclass for the current call to the _run() method.
        super().init(phase)

        # If saving post-simulation plots, create the top-level directory
        # containing these plots if needed.
        if phase.p.plot.is_after_sim_save:
            dirs.make_unless_dir(phase.export_dirname)

    # ..................{ SUPERCLASS                        }..................
    @type_check
    def _is_enabled(self, phase: SimPhase) -> bool:
        return phase.p.plot.is_after_sim

    # ..................{ PRIVATE ~ preparers               }..................
    @type_check
    def _export_prep(self, phase: SimPhase) -> None:
        '''
        Prepare to export the current plot for the passed simulation phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        #FIXME: DRY. This functionality perfectly duplicates the
        #AnimCellsWhileSolving.__enter__() method, which is bad. To resolve
        #this:
        #
        #* Shift that method into the "VisualCellsABC" superclass.
        #* Refactor all plots to subclass that superclass.

        # Id displaying this plot, do so in a non-blocking manner.
        if phase.p.plot.is_after_sim_show:
            # If the current matplotlib backend supports "true" non-blocking
            # behavior, prefer this non-deprecated approach.
            if mpl_config.is_backend_current_nonblockable():
                pyplot.show(block=False)
            # Else, fallback to the deprecated approach guaranteed to apply to
            # all matplotlib backends.
            else:
                # pass
                matplotlib.interactive(True)
                # pyplot.show()


    @type_check
    def _export(self, phase: SimPhase, basename: str) -> None:
        '''
        Export the current plot to the current screen if displaying plots
        and/or to a file with the passed basename if saving plots.

        Parameters
        -----------
        phase : SimPhase
            Current simulation phase.
        basename : str
            Basename excluding filetype of the plot to be exported. For
            convenience, this method internally prepends this basename by the
            identifying prefix ``fig_`` *before* writing this file.
        '''

        #FIXME: DRY. This functionality perfectly duplicates the
        #AnimCellsWhileSolving.__exit__() method, which is bad.

        # Id displaying this plot *AND* the current matplotlib backend fails to
        # support "true" non-blocking behavior...
        if (phase.p.plot.is_after_sim_show and
            not mpl_config.is_backend_current_nonblockable()):
            # Update all artists displayed by this plot.
            pyplot.draw()

            #FIXME: DRY. This functionality perfectly duplicates the
            #VisualCellsABC._show_frame(() method, which is also bad.

            # Temporarily yield the time slice for the smallest amount of time
            # required by the current matplotlib backend to handle queued
            # events in the GUI-specific event loop of the current process.
            with mplutil.deprecations_ignored():
                pyplot.pause(0.0001)

            # Disable the "fake" non-blocking behavior enabled by the prior
            # _export_prep() call.
            matplotlib.interactive(False)

        # If saving this plot...
        if phase.p.plot.is_after_sim_save:
            # Filetype and basename of the file to be saved.
            filetype = phase.p.plot.image_filetype
            basename = 'fig_{}.{}'.format(basename, filetype)

            # Absolute path of the file to be saved.
            filename = pathnames.join(phase.export_dirname, basename)

            # Log this saving attempt.
            logs.log_debug('Saving plot: %s', filename)

            # Save this plot to this file.
            pyplot.savefig(
                filename,
                dpi=phase.p.plot.image_dpi,
                format=filetype,
                transparent=True,
            )

        #FIXME: Non-ideal. Under a threading scenario, this will introduce race
        #conditions. Ideally, the figure associated with this plot should be
        #explicitly passed to the mplfigure.close_figure() function as is
        #currently done in the VisualCellsABC.close() method.

        # If *NOT* displaying this plot, close this plot to conserve resources.
        if not phase.p.plot.is_after_sim_show:
            mplfigure.close_figure_current()
