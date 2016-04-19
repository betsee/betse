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
#"FuncAnimation" class is *ALWAYS* passed "self.plot_frame" -- which, being a
#bound method of the current animation object, should ensure that that object
#remains alive. Non-blocking animations may already work out of the box! Oh...
#no. The "FuncAnimation" class is instantiated into an attribute of this class,
#implying that when the last reference to instances of this class go away, they
#everything goes away. We probably will need to add circular references to the
#passed "_figure", as detailed above.
#FIXME: Indeed, a (very minor) amount of research suggests that non-blocking
#animations should be trivially creatable by just ensuring that a reference to
#the instantiated FuncAnimation() object is retained, as intuited above: e.g.,
#
#    https://stackoverflow.com/questions/21099121/python-matplotlib-unable-to-call-funcanimation-from-inside-a-function

# ....................{ IMPORTS                            }....................
from abc import abstractmethod

import numpy as np
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation

from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib.anim import FileFrameWriter
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.science.plot.abc import PlotCells
from betse.util.io.log import logs
from betse.util.path import dirs, paths
from betse.util.type import types

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
        over this subclass' animation if `_is_overlaying_current` is
        `True` _or_ `None` otherwise.
    _time_step : int
        0-based index of the frame currently being plotted, corresponding to the
        0-based sampled time step currently being simulated.
    _is_overlaying_current : bool
        `True` if overlaying either electric current or concentration flux
        streamlines on this animation when requested by the current simulation
        configuration (as governed by the `p.I_overlay` parameter)_or_ `False`
        otherwise.
    _is_overlaying_current_gj_only : bool
        `True` if only overlaying intracellular current _or_ `False` if
        overlaying both intra- and extracellular current. Ignored if current is
        _not_ being overlayed at all (i.e., if `_is_overlaying_current` is
        `False`).
    _is_saving_shown_frames : bool
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

    # ..................{ LIFECYCLE                          }..................
    def __init__(
        self,

        #FIXME: I'm not terribly happy with defaulting these parameters here.
        #Ideally, we should define a new "AnimCellsAfterSolving" subclass simply
        #passing these defaults to this method as normal parameters. After doing
        #so, eliminate these defaults, thus forcing these parameters to
        #*ALWAYS* be passed.
        save_dir_parent_basename: str = 'anim',

        is_current_overlayable: bool = False,
        is_overlaying_current_gj_only: bool = None,
        is_ecm_required: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        save_dir_parent_basename : str
            Basename of the parent directory of the subdirectory to which this
            animation's frames will be saved when requested by the current
            simulation configuration. Defaults to a suitably generic basename.
        is_current_overlayable : bool
            `True` if overlaying either electric current or concentration flux
            streamlines on this animation when requested by the current
            simulation configuration (as governed by the `p.I_overlay`
            parameter) _or_ `False` otherwise. All subclasses except those
            already plotting streamlines (e.g., by calling the superclass
            `_plot_stream()` method) should unconditionally enable this boolean.
            Defaults to `False`.
        is_overlaying_current_gj_only : bool
            `True` if only overlaying intracellular current _or_ `False` if
            overlaying both intra- and extracellular current. Ignored if current
            is _not_ being overlayed at all (i.e., if `_is_overlaying_current`
            is `False`). If `None`, defaults to the following state:
            * `False` if extracellular spaces are enabled _and_ both
               intracellular and extracellular current is to be animated.
            * `True` if either extracellular spaces are disabled _or_ are
               enabled but only intracellular current is to be animated.
        is_ecm_required : bool
            `True` if this animation is specific to extracellular spaces or
            `False` otherwise. If `True` and extracellular spaces are currently
            disabled, an exception is raised. Defaults to `False`.
        '''
        assert types.is_str(save_dir_parent_basename), (
            types.assert_not_str(save_dir_parent_basename))
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

        # Default unpassed parameters.
        if is_overlaying_current_gj_only is None:
            is_overlaying_current_gj_only = not (
                self._p.sim_ECM and self._p.IecmPlot)

        # Classify defaulted parameters.
        self._is_overlaying_current_gj_only = is_overlaying_current_gj_only

        # Validate these attributes.
        assert types.is_bool(self._is_overlaying_current_gj_only), (
            types.assert_not_bool(self._is_overlaying_current_gj_only))

        # If this subclass requests a current overlay, do so only if also
        # requested by the current simulation configuration.
        self._is_overlaying_current = (
            is_current_overlayable and self._p.I_overlay)

        # True if both saving and displaying animation frames.
        self._is_saving_shown_frames = (
            self._is_showing and self._is_saving)

        # Type of animation attempt to be logged below.
        animation_verb = None
        if self._is_showing:
            animation_verb = 'Plotting'
        elif self._is_saving:
            animation_verb = 'Saving'
        # If neither displaying or saving this animation, this animation would
        # ideally reduce to a noop. Since this is a superclass method, however,
        # simply returning would have little effect; while raising an exception
        # would certainly have an effect, doing so would also require all
        # callers to explicitly catch and ignore that exception -- in which case
        # this object would hardly have reduced to a noop. In short, ignoring
        # this edge case is currently the only sane policy.

        # Log this animation as early as reasonably feasible.
        if animation_verb is not None:
            logs.log_info(
                '{} animation "{}"...'.format(animation_verb, self._type))

        # Classify attributes to be possibly redefined below.
        self._current_density_magnitude_time_series = None
        self._current_density_x_time_series = None
        self._current_density_y_time_series = None
        self._current_density_stream_plot = None
        self._time_step = 0
        self._writer_frames = None
        self._writer_frames = None
        self._writer_video = None

        # If saving animations, prepare to do so.
        if self._is_saving:
            self._init_saving(save_dir_parent_basename=save_dir_parent_basename)

        # If overlaying current, prepare to do so.
        if self._is_overlaying_current:
            self._init_current_density()


    def _init_saving(
        self,
        save_dir_parent_basename: str,
    ) -> None:
        '''
        Initialize this animation for platform-compatible file saving if enabled
        by the current simulation configuration or noop otherwise.

        Parameters
        ----------
        save_dir_parent_basename : str
            Basename of the parent directory of the subdirectory to which this
            animation's frames will be saved when requested by the current
            simulation configuration. Defaults to a suitably generic basename.
        '''
        assert types.is_str(save_dir_parent_basename), (
            types.assert_not_str(save_dir_parent_basename))

        # Ensure that the passed directory and file basenames are actually
        # basenames and hence contain no directory separators.
        paths.die_unless_basename(self._type)

        #FIXME: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #This is an utter travesty, but we have no choice but to hack detection
        #of the current loop type until we sort out just what is going on with
        #this boolean and/or string enumeration elsewhere.
        #FIXME: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if hasattr(self._p, 'plot_type'):
            plot_type = self._p.plot_type
        else:
            plot_type = 'sim' if self._p.run_sim else 'init'

        # Path of the phase-specific parent directory of the subdirectory to
        # which these files will be saved.
        loop_dirname = None
        if plot_type == 'sim':
            loop_dirname = self._p.sim_results
        elif plot_type == 'init':
            loop_dirname = self._p.init_results
        else:
            raise BetseExceptionParameters(
                'Anim saving unsupported during the "{}" loop.'.format(
                    plot_type))

        #FIXME: Refactor all calls to os.makedirs() everywhere similarly.

        # Path of the subdirectory to which these files will be saved, creating
        # this subdirectory and all parents thereof if needed.
        save_dirname = paths.join(
            loop_dirname, save_dir_parent_basename, self._type)
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
        if self._is_saving_shown_frames:
            self._writer_frames.setup(
                fig=self._figure,
                outfile=self._save_frame_template,

                #FIXME: Pass the actual desired "dpi" parameter.
                dpi=mpl_config.get_rc_param('savefig.dpi'),
            )


    # This method has been overriden to support subclasses that manually handle
    # animations rather than calling the _animate() method (e.g., the
    # "AnimCellsWhileSolving" subclass).
    def _prep_figure(self, *args, **kwargs) -> None:

        super()._prep_figure(*args, **kwargs)

        # If overlaying current, do so.
        if self._is_overlaying_current:
            self._plot_current_density()

    # ..................{ PROPERTIES                         }..................
    # The following testers are intended to be overridden by subclasses.
    #
    # The corresponding attributes (e.g., "_is_showing" for _is_showing())
    # *MUST* be defined via dynamic methods rather than static attributes passed
    # to this class' __init__() method (e.g., as an "is_saving" parameter.) Why?
    # Because chicken-and-the-egg constraints.  Specifically, the latter
    # approach prevents subclasses from passing a value dependent on the current
    # "Parameters" object to __init__(), as that object has yet to be classified
    # as the "_p" attribute yet. (Ugh.)

    @property
    def _is_showing(self) -> bool:
        '''
        `True` if interactively displaying this animation _or_ `False`
        otherwise.
        '''
        return not self._p.turn_all_plots_off


    @property
    def _is_saving(self) -> bool:
        '''
        `True` if non-interactively saving this animation as discrete frames
        and/or an encoded video _or_ `False` otherwise.
        '''
        return self._p.saveAnimations

    # ..................{ ANIMATORS                          }..................
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

        # Prepare for plotting immediately *BEFORE* plotting the first frame.
        self._prep_figure(*args, **kwargs)

        #FIXME: For efficiency, we should probably be passing "blit=True," to
        #FuncAnimation(). Lemon grass and dill!

        # Create and assign an animation function to a local variable. If the
        # latter is *NOT* done, this function will be garbage collected prior
        # to subsequent plot handling -- in which case only the first plot will
        # be plotted without explicit warning or error. Die, matplotlib! Die!!
        self._anim = FuncAnimation(
            self._figure, self.plot_frame,

            # Number of frames to be animated.
            frames=frame_count,

            # Delay in milliseconds between consecutive frames.
            interval=100,

            # Indefinitely repeat this animation unless saving animations, as
            # doing so under the current implementation would repeatedly (and
            # hence unnecessarily) overwrite previously written files.
            repeat=not self._is_saving,
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
            if self._is_showing:
                pyplot.show()
            # Else if saving animation frames, do so.
            elif self._is_saving:
                #FIXME: Pass the "dpi" parameter as well.
                self._anim.save(
                    filename=self._save_frame_template,
                    writer=self._writer_frames,
                    savefig_kwargs=self._writer_savefig_kwargs,
                )

                # For space efficiency, explicitly close this animation *AFTER*
                # saving this animation in a non-blocking manner.
                self.close()
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

    # ..................{ PLOTTERS                           }..................
    #FIXME: Rename "frame_number" to "time_step".
    def plot_frame(self, frame_number: int) -> None:
        '''
        Iterate this animation to the next frame.

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
        # loggers.log_info(
        #     'Updating animation "{}" frame number {}...'.format(
        #         self._type,
        #         len(self._sim.time) if frame_number == -1 else frame_number))

        # Classify the passed time step.
        self._time_step = frame_number

        # If plotting a current overlay, do so.
        if self._is_overlaying_current:
            #FIXME: Refactor _replot_current_density() to use
            #"self._time_step" rather than repass this parameter everywhere.
            self._replot_current_density(frame_number)

        # Update this figure with the current time, rounded to three decimal
        # places for readability.
        self._axes.set_title('{} (time {:.3f}s)'.format(
            self._axes_title, self._sim.time[frame_number]))

        #FIXME: Refactor _plot_frame_figure() to use "self._time_step" rather
        #than repass this parameter everywhere.

        # Plot this frame *AFTER* performing all superclass-specific plotting,
        # permitting the subclass to modify that plotting.
        self._plot_frame_figure(frame_number)

        # If both saving and displaying animation frames, save this frame. Note
        # that if only saving but *NOT* displaying animations, this frame will
        # be handled by our _animate() method. Why? Because Matplotlib will fail
        # to iterate frames and hence call our _plot_next_frame() method calling
        # this method *UNLESS* our _animate() method explicitly calls the
        # FuncAnimation.save() method with the writer name "frame" signifying
        # our "FileFrameWriter" class to do so. (Look. It's complicated, O.K.?)
        if self._is_saving_shown_frames:
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
    def _init_current_density(self) -> None:
        '''
        Initialize all attributes pertaining to current density.

        Specifically, this method defines the `_current_density_x_time_series`,
        `_current_density_y_time_series`, and
        `_current_density_magnitude_time_series` attributes. These attributes
        are required both by this superclass for animating current overlays
        _and_ by current-specific subclasses.
        '''

        # Time series of all current density X and Y components.
        if self._is_overlaying_current_gj_only is True:
            self._current_density_x_time_series = self._sim.I_gj_x_time
            self._current_density_y_time_series = self._sim.I_gj_y_time
        else:
            self._current_density_x_time_series = self._sim.I_tot_x_time
            self._current_density_y_time_series = self._sim.I_tot_y_time

        # Time series of all current density magnitudes (i.e., `Jmag_M`).
        self._current_density_magnitude_time_series = np.sqrt(
            np.asarray(self._current_density_x_time_series) ** 2 +
            np.asarray(self._current_density_y_time_series) ** 2) + 1e-30


    def _plot_current_density(self) -> None:
        '''
        Overlay the first frame of this subclass' animation with a streamplot of
        either electric current or concentration flux.
        '''

        # If animating only intracellular current, do so.
        if self._is_overlaying_current_gj_only:
            self._axes_title = 'Gap Junction Current'

            #FIXME: Is there any point to this? From what we can tell, the
            #"self._current_density_stream_plot" will simply be outright
            #replaced for the first and all subsequent frames. Galloping fish!
            self._current_density_stream_plot, self._axes = cell_stream(
                self._current_density_x_time_series[-1],
                self._current_density_y_time_series[-1],
                self._axes, self._cells, self._p)
        # If animating both extracellular and intracellular current, do so.
        else:
            self._axes_title = 'Total Current Overlay'

            #FIXME: Likewise.
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
        super().__init__(*args, **kwargs)

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
        super().__init__(*args, **kwargs)

        # Prefer an alternative colormap.
        self._colormap = self._p.background_cm
