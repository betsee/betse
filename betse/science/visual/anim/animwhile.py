#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Mid-simulation animation** (i.e., animation produced *while* rather than
*after* solving a simulation) subclasses.
'''

# ....................{ IMPORTS                            }....................
import matplotlib
import numpy as np
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.science.export import expmath
from betse.science.phase.phasecls import SimPhase
from betse.science.visual.anim.animabc import AnimCellsABC
from betse.util.type.types import type_check, SequenceTypes
from matplotlib import pyplot

# ....................{ SUBCLASSES                         }....................
#FIXME: Rename "_cell_data_plot" to "_cell_body_plot".
#FIXME: Rename "_cell_edges_plot" to "_cell_edge_plot".

class AnimCellsWhileSolving(AnimCellsABC):
    '''
    Context manager animating a mid-simulation animation.

    This manager animates arbitrary cell data as a time series plotted over the
    cell cluster (e.g., cell membrane voltage as a function of time) *while*
    rather than *after* solving a simulation phase.

    Caveats
    ----------
    This animation *must* be the target clause of a ``with`` statement, ensuring
    safe setup and teardown of this animation's non-blocking behavior. Likewise,
    this animation's :meth:`plot_frame` method *must* be called within the body
    of this ``with`` statement after producing for each sampled time step the
    cell data visualized by this animation for the corresponding frame.

    Attributes
    ----------
    _cell_edges_plot : LineCollection
        Plot of the current or prior frame's cell edges.
    _cell_data_plot : Collection
        Plot of the current or prior frame's cell contents.
    _cell_verts_id : int
        Unique identifier for the array of cell vertices (i.e.,
        `cells.cell_verts`) when plotting the current or prior frame.
        Retaining this identifier permits the `_plot_frame_figure()` method to
        efficiently detect and respond to physical changes (e.g., deformation
        forces, cutting events) in the fundamental structure of the previously
        plotted cell cluster.
    _is_colorbar_autoscaling_telescoped : bool
        ``True`` if colorbar autoscaling is permitted to increase but not
        decrease the colorbar range *or* ``False`` otherwise (i.e., if
        colorbar autoscaling is permitted to both increase and decrease the
        colorbar range). Such telescoping assists in emphasizing stable
        long-term patterns in cell data at a cost of deemphasizing unstable
        short-term patterns. If colorbar autoscaling is disabled (i.e.,
        ``is_color_autoscaled`` is ``False``), this will be ignored.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        phase: SimPhase,

        # Optional parameters.

        #FIXME: Permit this option to be configured via a new configuration
        #option in the YAML file.
        #FIXME: Actually, just excise this everywhere. Ally dislikes this effect
        #and I certainly never require it; ergo, remove.

        is_colorbar_autoscaling_telescoped: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this in-simulation animation.

        Parameters
        ----------
        phase: SimPhase
            Current simulation phase.
        is_colorbar_autoscaling_telescoped : optional[bool]
            `True` if colorbar autoscaling is permitted to increase but _not_
            decrease the colorbar range _or_ `False` otherwise (i.e., if
            colorbar autoscaling is permitted to both increase and decrease the
            colorbar range). Such telescoping assists in emphasizing stable
            long-term patterns in cell data at a cost of deemphasizing unstable
            short-term patterns. If colorbar autoscaling is disabled (i.e.,
            `is_color_autoscaled` is `False`), this will be ignored. Defaults
            to `False`.

        See the superclass method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(
            *args,

            # Pass this simulation phase as is to our superclass.
            phase=phase,

            # Prevent the superclass from overlaying electric current or
            # concentration flux. Although this class does *NOT* animate a
            # streamplot, the time series required to plot this overlay is
            # unavailable until after the simulation ends.
            is_current_overlayable=False,

            # Save and show this mid-simulation animation only if this
            # configuration has enabled doing so.
            is_save=phase.p.anim.is_while_sim_save,
            is_show=phase.p.anim.is_while_sim_show,

            # Save this mid-simulation animation to a different parent
            # directory than that to which the corresponding post-simulation
            # animation is saved.
            save_dir_parent_basename='anim_while_solving',

            # Pass all remaining arguments as is to our superclass.
            **kwargs
        )
        # logs.log_debug('Showing mid-simulation animation: {}'.format(self._is_show))

        # Classify all remaining parameters.
        self._is_colorbar_autoscaling_telescoped = (
            is_colorbar_autoscaling_telescoped)

        # Unique identifier for the array of cell vertices. (See docstring.)
        self._cell_verts_id = id(self._phase.cells.cell_verts)

        #FIXME: This is a temp change until we get this right.
        #FIXME: Refactor to call the new
        #Cells.map_membranes_midpoint_to_cells_centre() method instead.

        # average the voltage to the cell centre
        vm_o = np.dot(
            self._phase.cells.M_sum_mems, self._phase.sim.vm) / (
                self._phase.cells.num_mems)

        # self._cell_time_series = self.sim.vm_time
        self._cell_time_series = self._phase.sim.vm_ave_time

        # cell_data_current = self.sim.vm
        cell_data_current = vm_o

        # Upscaled cell data for the first frame.
        cell_data = expmath.upscale_units_milli(cell_data_current)

        # Collection of cell polygons with animated voltage data.
        #
        # If *NOT* simulating extracellular spaces, only animate intracellular
        # spaces.
        self._cell_data_plot = self._plot_cells_sans_ecm(
            cell_data=cell_data)

        # Perform all superclass plotting preparation immediately *BEFORE*
        # plotting this animation's first frame.
        self._prep_figure(
            color_mappables=self._cell_data_plot,
            color_data=cell_data,
        )

    # ..................{ CONTEXTS                           }..................
    def __enter__(self) -> 'AnimCellsWhileSolving':
        '''
        Enter the runtime context for this context manager, returning the same
        context manager as the value bound to the identifier in the ``as``
        clause of the ``with`` block using this context manager.

        This special method temporarily enables non-blocking matplotlib behavior
        for the duration of this plot or animation.
        '''

        # If displaying this animation, do so in a non-blocking manner.
        if self._is_show:
            # If the current matplotlib backend supports "true" non-blocking
            # behavior, prefer this non-deprecated approach.
            if mpl_config.is_backend_current_nonblockable():
                pyplot.show(block=False)
            # Else, fallback to the deprecated approach guaranteed to apply to
            # all matplotlib backends.
            else:
                matplotlib.interactive(True)
                pyplot.show()

        # Bind this animation to the "as" clause of this "with" block.
        return self


    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        '''
        Exit the runtime context for this context manager given the passed
        exception metadata for an exception raised by the body of the ``with``
        block using this context manager if any, returning ``True`` only if that
        exception should be suppressed by this ``with`` block.

        This special method (in order):

        . Disables the non-blocking matplotlib behavior temporarily enabled by
          the prior :meth:`__enter__` call.
        . Safely closes this plot or animation.
        '''

        # Id displaying this animation *AND* the current matplotlib backend
        # fails to support "true" non-blocking behavior, disable the "fake"
        # non-blocking behavior enabled by the prior __enter__() call.
        if self._is_show and not mpl_config.is_backend_current_nonblockable():
            matplotlib.interactive(False)

        # Close this plot or animation.
        self.close()

        # Avoid suppressing exceptions raised by this "with" block.
        return False

    # ..................{ PLOTTERS                           }..................
    def _plot_frame_figure(self) -> None:

        # Upscaled cell data for the current time step.
        cell_data = expmath.upscale_units_milli(
            self._cell_time_series[self._time_step])

        #FIXME: Duplicated from above. What we probably want to do is define a
        #new _get_cell_data() method returning this array in a centralized
        #manner callable both here and above. Rays of deluded beaming sunspray!

        # If the unique identifier for the array of cell vertices has *NOT*
        # changed, the cell cluster has *NOT* fundamentally changed and need
        # only be updated with this time step's cell data.
        if self._cell_verts_id == id(self._phase.cells.cell_verts):
            # loggers.log_info(
            #     'Updating animation "{}" cell plots...'.format(self._type))
            self._update_cell_plots(cell_data)
        # Else, the cell cluster has fundamentally changed (e.g., due to
        # physical deformations or cutting events) and must be recreated.
        else:
            # loggers.log_info(
            #     'Reviving animation "{}" cell plots...'.format(self._type))

            # Prevent subsequent calls to this method from erroneously
            # recreating the cell cluster again.
            self._cell_verts_id = id(self._phase.cells.cell_verts)

            # Recreate the cell cluster.
            self._revive_cell_plots(cell_data)

        # Update the color bar with the content of the cell body plot *AFTER*
        # possibly recreating this plot above.
        if self._conf.is_color_autoscaled:
            cell_data_vm = cell_data

            # If autoscaling this colorbar in a telescoping manner and this is
            # *NOT* the first time step, do so.
            #
            # If this is the first time step, the previously calculated minimum
            # and maximum colors are garbage and thus *MUST* be ignored.
            if (self._is_colorbar_autoscaling_telescoped and
                not self._is_time_step_first):
                self._color_min = min(self._color_min, np.min(cell_data_vm))
                self._color_max = max(self._color_max, np.max(cell_data_vm))
            # Else, autoscale this colorbar in an unrestricted manner.
            else:
                self._color_min = np.min(cell_data_vm)
                self._color_max = np.max(cell_data_vm)

            # Autoscale the colorbar to these colors.
            self._rescale_color_mappables()

        # If displaying this frame...
        if self._is_show:
            # If the current event loop is idle, draw this frame; else, noop.
            # This is the OO-style equivalent to calling pyplot.draw().
            self._figure.canvas.draw_idle()
            # self._figure.canvas.draw()


    @type_check
    def _update_cell_plots(self, cell_data: SequenceTypes) -> None:
        '''
        Update *without* recreating all cell plots for this time step with the
        passed array of arbitrary cell data.

        This method is intended to be called *unless* physical changes
        (e.g., deformation forces, cutting events) in the underlying structure
        of the cell cluster have occurred for this simulation time step.

        Parameters
        -----------
        cell_data : SequenceTypes
            Arbitrary cell data defined on an environmental grid to be plotted.
        '''

        self._update_cell_plot_sans_ecm(
            cell_plot=self._cell_data_plot,
            cell_data=cell_data)


    @type_check
    def _revive_cell_plots(self, cell_data: SequenceTypes) -> None:
        '''
        Recreate all cell plots for this time step with the passed array of
        arbitrary cell data.

        This method is intended to be called in response to physical changes
        (e.g., deformation forces, cutting events) in the underlying structure
        of the cell cluster for this simulation time step. This method is both
        inefficient and destructive, and should be called only when needed.

        Parameters
        -----------
        cell_data : SequenceTypes
            Arbitrary cell data defined on an environmental grid to be plotted.
        '''

        self._cell_data_plot = self._revive_cell_plots_sans_ecm(
            cell_plot=self._cell_data_plot,
            cell_data=cell_data)
