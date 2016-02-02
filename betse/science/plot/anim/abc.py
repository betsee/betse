#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based animation classes.
'''

#FIXME: We should probably animate non-blockingly (e.g., by passing
#"block=False" to the plt.show() command. To do so, however, we'll probably have
#to implement an analogue to Matplotlib's "_pylab_helper.Gcf" global-like static
#class to ensure that object references are retained for all animations whose
#windows have not yet been closed. This implies, in turn, that we'll probably
#have to monkey-patch the appropriate Matplotlib event handlers on window close
#to release these object references. To do so, grep the Matplotlib codebase for
#"Gcf.destroy". Thunderous tales of woe!
#
#Actually, I believe that a simpler approach might be available. Rather
#than implementing yet another "_pylab_helper.Gcf"-like construct, we leverage
#the fact that animation objects should only live as long as their underlying
#figure objects by explicitly adding circular references between the two: e.g.,
#
#    # This is already done by the "PlotCells" superclass.
#    self._figure = pyplot.figure()
#
#    # Then just add this to the AnimCells.__init__() method *BEFORE* the
#    # self._figure.show(block=False) method is called.
#    self._figure.__BETSE_anim__ = self
#
#That said, we might not even need to do that much. Why? Because the
#"FuncAnimation" class is *ALWAYS* passed "self._plot_frame" -- which, being a
#bound method of the current animation object, should ensure that that object
#remains alive. Non-blocking animations may already work out of the box! Oh...
#no. The "FuncAnimation" class is instantiated into an attribute of this class,
#implying that when the last reference to instances of this class go away, they
#everything goes away. We probably will need to add circular references to the
#passed "_figure", as detailed above.

# ....................{ IMPORTS                            }....................
import numpy as np
from abc import abstractmethod
from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib.mpl import mpl_config
from betse.lib.matplotlib.anim import FileFrameWriter
from betse.science.plot.abc import PlotCells
from betse.util.io import loggers
from betse.util.path import dirs, paths
from betse.util.type import types
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation

#FIXME: Shift these functions into our superclass.
from betse.science.plot.plot import (env_stream, cell_stream)

# ....................{ BASE                               }....................
class AnimCells(PlotCells):
    '''
    Abstract base class of all animation classes.

    Instances of this class animate the spatial distribution of modelled
    variables (e.g., Vmem) over all time steps of the simulation.

    Attributes
    ----------
    _anim : FuncAnimation
        Low-level Matplotlib animation object instantiated by this high-level
        BETSE wrapper object.
    _current_density_magnitude_time_series : ndarray
        Time series of all current density magnitudes (i.e., `Jmag_M`) if the
        optional `_init_current_density()` method has been called for
        this animation _or_ `None` otherwise.
    _current_density_x_time_series : list
        Time series of all current density X components if the optional
        `_init_current_density()` method has been called for this
        animation _or_ `None` otherwise.
    _current_density_y_time_series : list
        Time series of all current density Y components if the optional
        `_init_current_density()` method has been called for this
        animation _or_ `None` otherwise.
    _current_density_stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of either electric current or concentration flux overlayed
        over this subclass' animation if `_is_plotting_current_overlay` is
        `True` _or_ `None` otherwise.
    _is_plotting_current_overlay : bool
        `True` if overlaying either electric current or concentration flux
        streamlines on this animation when requested by the current simulation
        configuration (as governed by the `p.I_overlay` parameter)_or_ `False`
        otherwise.
    _is_saving_plotted_frames : bool
        `True` if both saving and displaying animation frames _or_ `False`
        otherwise.
    _save_frame_template : str
        `str.format()`-formatted template which, when formatted with the 0-based
        index of the current frame, yields the absolute path of the image file
        to be saved for that frame.
    _writer_frames : MovieWriter
        Object writing frames from this animation to image files if enabled _or_
        `None` otherwise.
    _writer_savefig_kwargs : dict
        Dictionary of all keyword arguments to be passed to the
        `Figure.savefig()` method called to write each animation frame.
    _writer_video : MovieWriter
        Object encoding this animation to a video file if enabled _or_ `None`
        otherwise.
    '''

    # ..................{ PRIVATE ~ init                     }..................
    def __init__(
        self,
        is_current_overlayable: bool = False,
        is_ecm_required: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        is_current_overlayable : bool
            `True` if overlaying either electric current or concentration flux
            streamlines on this animation when requested by the current
            simulation configuration (as governed by the `p.I_overlay`
            parameter)_or_ `False` otherwise. All subclasses except those
            already plotting streamlines (e.g., by calling the superclass
            `_plot_stream()` method) should unconditionally enable this boolean.
            Defaults to `False`.
        is_ecm_required : bool
            `True` if this animation is specific to extracellular spaces or
            `False` otherwise. If `True` and extracellular spaces are currently
            disabled, an exception is raised. Defaults to `False`.
        '''
        assert types.is_bool(is_current_overlayable), (
            types.assert_not_bool(is_current_overlayable))
        assert types.is_bool(is_ecm_required), (
            types.assert_not_bool(is_ecm_required))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # If this subclass requires extracellular spaces but extracellular
        # spaces are currently disabled, raise an exception.
        if is_ecm_required and not self._p.sim_ECM:
            raise BetseExceptionParameters(
                'Animation "{}" requires extracellular spaces, which are '
                'disabled by the current simulation configuration.'.format(
                self._type))

        # Classify attributes to be possibly redefined below.
        self._current_density_magnitude_time_series = None
        self._current_density_x_time_series = None
        self._current_density_y_time_series = None
        self._current_density_stream_plot = None
        self._writer_frames = None
        self._writer_frames = None
        self._writer_video = None

        # If this subclass requests a current overlay, do so only if also
        # requested by the current simulation configuration.
        self._is_plotting_current_overlay = (
            is_current_overlayable and self._p.I_overlay)

        # True if both saving and displaying animation frames.
        self._is_saving_plotted_frames = (
            self._p.turn_all_plots_off is False and
            self._p.saveAnimations is True)

        # Type of animation attempt to be logged below.
        if self._p.turn_all_plots_off is False:
            animation_verb = 'Plotting'
        elif self._p.saveAnimations is True:
            animation_verb = 'Saving'
        # If neither displaying or saving this animation, this animation would
        # ideally reduce to a noop. Since this is a superclass method, however,
        # simply returning would have little effect; while raising an exception
        # would certainly have an effect, doing so would also require all
        # callers to explicitly catch and ignore that exception -- in which case
        # this object would hardly have reduced to a noop. In short, ignoring
        # this edge case is currently the only sane policy.

        # Log this animation as early as reasonably feasible.
        loggers.log_info(
            '{} animation "{}"...'.format(animation_verb, self._type))

        # If saving animations, prepare to do so.
        if self._p.saveAnimations is True:
            self._init_saving()


    def _init_saving(self) -> None:
        '''
        Initialize this animation for platform-compatible file saving if enabled
        by the current simulation configuration or noop otherwise.
        '''

        # Ensure that the passed directory and file basenames are actually
        # basenames and hence contain no directory separators.
        paths.die_unless_basename(self._type)

        # Path of the phase-specific parent directory of the subdirectory to
        # which these files will be saved.
        phase_dirname = None
        if self._p.plot_type == 'sim':
            phase_dirname = self._p.sim_results
        elif self._p.plot_type == 'init':
            phase_dirname = self._p.init_results
        else:
            raise BetseExceptionParameters(
                'Anim saving unsupported during the "{}" phase.'.format(
                    self._p.plot_type))

        #FIXME: Refactor all calls to os.makedirs() everywhere similarly.

        # Path of the subdirectory to which these files will be saved, creating
        # this subdirectory and all parents thereof if needed.
        save_dirname = paths.join(
            phase_dirname, 'animation', self._type)
        save_dirname = dirs.canonicalize_and_make_unless_dir(save_dirname)

        #FIXME: Pull the image filetype from the current YAML configuration
        #rather than coercing use of ".png".
        save_frame_filetype = 'png'

        #FIXME: This currently defaults to padding frames with six or seven
        #zeroes, on average. Let's make this a bit more aesthetic by padding
        #frames to only as many zeroes are absolutely required by the current
        #frame count. To do that, in turn, we'll probably need to shift
        #everything that follows in this method to the _animate() method, where
        #the actual frame count is finally passed.

        # Template yielding the basenames of frame image files to be saved.
        # The "{{"- and "}}"-delimited substring will reduce to a "{"- and "}"-
        # delimited substring after formatting, which subsequent formatting
        # elsewhere (e.g., in the "FileFrameWriter" class) will expand with the
        # 0-based index of the current frame number.
        save_frame_template_basename = '{}_{{:07d}}.{}'.format(
            self._type, save_frame_filetype)

        # Template yielding the absolute paths of frame image files to be saved.
        self._save_frame_template = paths.join(
            save_dirname, save_frame_template_basename)

        # Object writing frames from this animation to image files.
        self._writer_frames = FileFrameWriter()

        # Dictionary of keyword arguments to be passed to Figure.savefig().
        self._writer_savefig_kwargs = {
            # Plot the background of each animation frame as transparent rather
            # than pure-white.
            'transparent': True,
        }

        # If both saving and displaying animation frames, prepare for doing so.
        # See the _save_frame() method for horrid discussion.
        if self._is_saving_plotted_frames:
            self._writer_frames.setup(
                fig=self._figure,
                outfile=self._save_frame_template,

                #FIXME: Pass the actual desired "dpi" parameter.
                dpi=mpl_config.get_rc_param('savefig.dpi'),
            )

    # ..................{ PRIVATE ~ animate                  }..................
    def _animate(
        self,
        frame_count: int,
        *args, **kwargs
    ) -> None:
        '''
        Display and/or save this animation as requested by the current
        simulation configuration.

        This method is intended to be called as the last statement in the
        `__init__()` method of all subclasses of this superclass.

        Attributes
        ----------
        frame_count : int
            Number of frames to be animated.

        All remaining parameters will be passed to the superclass `__prep()`
        method as is.
        '''
        assert types.is_int(frame_count), types.assert_not_int(frame_count)

        # If plotting a current overlay, do so *AFTER* this subclass has already
        # performed all initial plotting but *BEFORE* our superclass overlays
        # its even more critical plotting data (e.g., cell labels).
        if self._is_plotting_current_overlay:
            self._plot_current_density()

        # Perform all superclass plotting preparation. This should typically be
        # performed immediately *BEFORE* creating this animation.
        self._prep(*args, **kwargs)

        #FIXME: For efficiency, we should probably be passing "blit=True," to
        #FuncAnimation(). Lemon grass and dill!

        # Create and assign an animation function to a local variable. If the
        # latter is *NOT* done, this function will be garbage collected prior
        # to subsequent plot handling -- in which case only the first plot will
        # be plotted without explicit warning or error. Die, matplotlib! Die!!
        self._anim = FuncAnimation(
            self._figure, self._plot_frame,

            # Number of frames to be animated.
            frames=frame_count,

            # Delay in milliseconds between consecutive frames.
            interval=100,

            # Indefinitely repeat this animation unless saving animations, as
            # doing so under the current implementation would repeatedly (and
            # hence unnecessarily) overwrite previously written files.
            repeat=not self._p.saveAnimations,
        )

        #FIXME: It appears that movies can be saved at this exact point via the
        #following lines:
        #
        #    video_filename = 'my_filename.mp4'
        #    video_encoder_name = 'ffmpeg'
        #    ani.save(filename=video_filename, writer=video_encoder_name)
        #
        #Both the "video_filename" and "video_encoder_name" variables used
        #above should be initialized from YAML configuration items. The latter
        #should probably be configured as a list of preferred encoders: e.g.,
        #
        #    # List of Matplotlib-specific names of preferred video encoders.
        #    # The first encoder available on the current $PATH will be used.
        #    video encoders: ['ffmpeg', 'avconv', 'mencoder']
        #
        #Given that, we would then need to iteratively search this list until
        #finding the first encoder available on the current $PATH. Doesn't look
        #too hard, actually. Grandest joy in the cornucopia of easy winter!
        #FIXME: Also note that movie encoders are configurable as follows,
        #which is probably a better idea than just passing a string above:
        #
        #    Writer = animation.writers[video_encoder_name]
        #    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        #    ani.save(filename=video_filename, writer=writer)
        #
        #Oh, wait. No, that's overkill. All of the above parameters are also
        #accepted by the ani.save() function itself, so string name it is!

        try:
            # If displaying animations, do so.
            if self._p.turn_all_plots_off is False:
                pyplot.show()
            # Else if saving animation frames, do so.
            elif self._p.saveAnimations is True:
                #FIXME: Pass the "dpi" parameter as well.
                self._anim.save(
                    filename=self._save_frame_template,
                    writer=self._writer_frames,
                    savefig_kwargs=self._writer_savefig_kwargs,
                )

                # For space efficiency, explicitly close this animation *AFTER*
                # saving this animation in a non-blocking manner.
                self._close()
        # plt.show() unreliably raises exceptions on window close resembling:
        #     AttributeError: 'NoneType' object has no attribute 'tk'
        # This error appears to ignorable and hence is caught and squelched.
        except AttributeError as exc:
            # If this is that exception, mercilessly squelch it.
            if str(exc) == "'NoneType' object has no attribute 'tk'":
                pass
            # Else, reraise this exception.
            else:
                raise

    # ..................{ PRIVATE ~ plot                     }..................
    def _plot_frame(self, frame_number: int) -> None:
        '''
        Iterate the current animation to the next frame.

        This method is iteratively called by Matplotlib's `FuncAnimation()`
        class instantiated by our `_animate()` method. The subclass
        implementation of this abstract method is expected to modify this
        animation's figure and axes so as to display the next frame. It is
        _not_, however, expected to save that figure to disk; frame saving is
        already implemented by this base class in a general-purpose manner.

        Specifically, this method (in order):

        . Calls the subclass `_plot_frame_figure()` method to update the current
          figure with this frame's data.
        . Updates the current figure's axes title with the current time.
        . Optionally writes this frame to disk if desired.

        Parameters
        ----------
        frame_number : int
            0-based index of the frame to be plotted.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Plot this frame onto this animation's figure.
        self._plot_frame_figure(frame_number)

        # If plotting a current overlay, do so *AFTER* this subclass has already
        # performed all plotting for this frame.
        if self._is_plotting_current_overlay:
            self._replot_current_density(frame_number)

        # Update this figure with the current time, rounded to three decimal
        # places for readability.
        self._axes.set_title('{} (time {:.3f}s)'.format(
            self._axes_title, self._sim.time[frame_number]))

        # If both saving and displaying animation frames, save this frame. Note
        # that if only saving but *NOT* displaying animations, this frame will
        # be handled by our _animate() method. Why? Because Matplotlib will fail
        # to iterate frames and hence call our _plot_next_frame() method calling
        # this method *UNLESS* our _animate() method explicitly calls the
        # FuncAnimation.save() method with the writer name "frame" signifying
        # our "FileFrameWriter" class to do so. (Look. It's complicated, O.K.?)
        # if self.p.saveAnimations is False:
        if self._is_saving_plotted_frames:
            self._writer_frames.grab_frame(**self._writer_savefig_kwargs)


    @abstractmethod
    def _plot_frame_figure(self, frame_number: int) -> None:
        '''
        Plot the frame with the passed 0-based index onto the current figure.

        Parameters
        ----------
        frame_number : int
            0-based index of the frame to be plotted.
        '''
        pass

    # ..................{ PRIVATE ~ plot : current           }..................
    def _init_current_density(self, is_gj_only: bool) -> None:
        '''
        Initialize all attributes pertaining to current density.

        Specifically, this method defines the `_current_density_x_time_series`,
        `_current_density_y_time_series`, and
        `_current_density_magnitude_time_series` attributes. These attributes
        are required both by this superclass for animating current overlays
        _and_ by current-specific subclasses.

        Parameters
        ----------
        is_gj_only : bool
            Either:
            * `True` if initializing these attributes to only the intracellular
              gap junction current density.
            * `False` if initializing these attributes to all current density
              (i.e., of both extracellular spaces _and_ intracellular gap
              junctions).
        '''
        assert types.is_bool(is_gj_only), types.assert_not_bool(is_gj_only)

        # Time series of all current density X and Y components.
        if is_gj_only is True:
            self._current_density_x_time_series = self._sim.I_gj_x_time
            self._current_density_y_time_series = self._sim.I_gj_y_time
        else:
            self._current_density_x_time_series = self._sim.I_tot_x_time
            self._current_density_y_time_series = self._sim.I_tot_y_time

        # Time series of all current density magnitudes (i.e., `Jmag_M`).
        self._current_density_magnitude_time_series = np.sqrt(
            np.array(self._current_density_x_time_series) ** 2 +
            np.array(self._current_density_y_time_series) ** 2) + 1e-30


    def _plot_current_density(self) -> None:
        '''
        Overlay the first frame of this subclass' animation with a streamplot of
        either electric current or concentration flux.
        '''

        # Boolean:
        #
        # * "True" if either extracellular spaces are disabled *OR* are enabled
        #   but only intracellular current is to be animated.
        # * "False" if extracellular spaces are enabled *AND* both intracellular
        #   and extracellular current is to be animated.
        is_gj_only = self._p.sim_ECM is False or self._p.IecmPlot is False

        # Initialize all attributes pertaining to current density.
        self._init_current_density(is_gj_only)

        # If animating only intracellular current, do so.
        if is_gj_only:
            self._axes_title = 'Gap Junction Current'
            self._current_density_stream_plot, self._axes = cell_stream(
                self._current_density_x_time_series[-1],
                self._current_density_y_time_series[-1],
                self._axes, self._cells, self._p)
        # If animating both extracellular and intracellular current, do so.
        else:
            self._axes_title = 'Total Current Overlay'
            self._current_density_stream_plot, self._axes = env_stream(
                self._current_density_x_time_series[-1],
                self._current_density_y_time_series[-1],
                self._axes, self._cells, self._p)


    def _replot_current_density(self, frame_number: int) -> None:
        '''
        Overlay the passed frame of this subclass' animation with a streamplot
        of either electric current or concentration flux.

        Parameters
        -----------
        frame_number : int
            0-based index of the frame to be plotted.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Current density magnitudes for this frame.
        Jmag_M = self._current_density_magnitude_time_series[frame_number]

        # Erase the prior frame's overlay and streamplot this frame's overlay.
        self._current_density_stream_plot = self._plot_stream(
            old_stream_plot=self._current_density_stream_plot,
            x=self._current_density_x_time_series[frame_number] / Jmag_M,
            y=self._current_density_y_time_series[frame_number] / Jmag_M,
            magnitude=Jmag_M,
        )

# ....................{ SUBCLASSES                         }....................
class AnimField(AnimCells):
    '''
    Abstract base class of all animations of electric field strength plotted on
    the current cell cluster.

    Attributes
    ----------
    _magnitude_time_series : list
        Electric field magnitudes as a function of time.
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's electric field magnitude.
    _stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of the current or prior frame's electric field.
    _x_time_series : list
        Electric field X components as a function of time.
    _y_time_series : list
        Electric field Y components as a function of time.
    _unit_x_time_series : list
        Electric field X unit components as a function of time. The resulting
        electric field vectors are **unit vectors** (i.e., have magnitude 1).
    _unit_y_time_series : list
        Electric field Y unit components as a function of time. The resulting
        electric field vectors are **unit vectors** (i.e., have magnitude 1).
    '''

    def __init__(
        self,
        x_time_series: list,
        y_time_series: list,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        x_time_series : list
            List of all electric field strength X components indexed by
            simulation time.
        y_time_series : list
            List of all electric field strength Y components indexed by
            simulation time.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(x_time_series), (
            types.assert_not_sequence_nonstr(x_time_series))
        assert types.is_sequence_nonstr(y_time_series), (
            types.assert_not_sequence_nonstr(y_time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            axes_x_label='Spatial distance [um]',
            axes_y_label='Spatial distance [um]',
            *args, **kwargs)

        # Classify all remaining parameters.
        self._x_time_series = x_time_series
        self._y_time_series = y_time_series

        # Electric field magnitudes and X and Y unit components.
        self._magnitude_time_series = []
        self._unit_x_time_series = []
        self._unit_y_time_series = []

        # Prefer an alternative colormap.
        self._colormap = self._p.background_cm


class AnimVelocity(AnimCells):
    '''
    Abstract base class of all animations of a velocity flow plotted on the
    current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            axes_x_label='Spatial distance [um]',
            axes_y_label='Spatial distance [um]',
            *args, **kwargs)

        # Prefer an alternative colormap.
        self._colormap = self._p.background_cm
