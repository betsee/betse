#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based animation classes.
'''

#FIXME: Current overlays (as enabled by the "is_current_overlayable" boolean and
#animation-specific configuration options), appear to be broken. Panic station!
#FIXME: Actually, the current approach to implementing animation overlays is
#fundamentally flawed. We currently attempt to provide a crude form of plot
#composition (i.e., merging two or more types of plots together into a single
#plot) by adding new booleans to the "AnimCells" base class (e.g.,
#"is_current_overlayable") -- a fundamentally unwieldy and ultimately unworkable
#approach. By definition, you cannot provide true composability frow within a
#single class hierarchy. Instead, we need to split the specific process of
#plotting different types of artists (e.g., mesh plots, stream plots) from the
#general process of animating and saving frames and plots as follows:
#
#* Define a new "betse.science.plot.anim.plotter" submodule.
#* Define a new "CellsPlotterABC" abstract base class in this submodule. Plotter
#  classes encapsulate the plotting of a single type of plot (e.g., Vmem-style
#  cell cluster meshplot, electrical current streamplot). While they may
#  internally cache one or more time series required for this plotting, they do
#  *NOT* contain the axes, figure, plot, or animation currently being plotted.
#  Since these Matplotlib artists are shared between multiple plotters, the
#  existing "PlotCells" and "AnimCells" base classes retain ownership of these
#  artists.
#* Define the following abstract methods in "CellsPlotterABC":
#  * An init() method. Subclasses *MUST* redefine this method to initialize this
#    plotter instance (e.g., by internally caching one or more time series).
#  * A plot() method presumably accepting the index of the frame to be plotted.
#    Subclasses *MUST* redefine this method to plot the time series data for
#    this frame into the current plot or animation figure. To avoid circular
#    references (and hence unallocated figures), it wouldn't necessarily be a
#    bad idea to repass all objects required for plotting to each
#    CellsPlotterABC.plot() call rather than to the CellsPlotterABC.__init__()
#    constructor. Food for pleasant thought, anyway.
#  * *UHM.* Wait. There would clearly be *NO* circular references here, so
#    passing objects to CellsPlotterABC.__init__() once rather than repeatedly
#    to CellsPlotterABC.plot() probably makes the most sense. We absolutely
#    *MUST*, however, avoid passing the current "AnimCells" instance to
#    CellsPlotterABC.__init__(). Passing that instance to
#    CellsPlotterABC.plot(), however, should both be safe and desirable, as
#    doing so would then permit plotters to access relevant data on the current
#    animation (e.g., "AnimCells.time_step" providing the current frame number).
#* Add a new "plotters" parameter to the AnimCells.__init__() constructor,
#  classified as a new "AnimCells._plotters" instance variable. This parameter
#  and variable *MUST* be a list of "CellsPlotterABC" instances. The order of
#  plotters in this list defines the order in which these plotters are drawn and
#  hence overlaid one another (i.e., z-order).
#* Refactor AnimCells.__init__() or a method called by that method to iterate
#  over "self._plotters" and initialize each such plotter by calling
#  plotter.init().
#* Refactor AnimCells.plot_frame() or a related method to iterate over
#  "self._plotters" and draw each such plotter by calling plotter.draw().
#* Refactor all subclasses of "AnimCells" into one or more subclasses of
#  "CellsPlotterABC" instead, which may then be instantiated and composed
#  together into a new "plotters" list passed to CellsPlotterABC.__init__(). For
#  example:
#  * Split the existing "AnimGapJuncTimeSeries" subclass into:
#    * A new "CellsPlotterGapJunc" subclass plotting *ONLY* the gap junction
#      open state as a "LineCollection" overlay. This plotter subclass would
#      probably only be used for this specific purpose.
#    * A new "CellsPlotterTimeSeries" subclass plotting *ONLY* an arbitrary
#      time series for the cell cluster as a mesh plot underlay. Clearly, this
#      plotter subclass would be extensively reused elsewhere as well.
#* Replace all current overlay functionality in "AnimCells" with "plotters".
#* Refactor the configuration file from the current hard-coded non-composable
#  approach to a dynamic list-based approach permitting zero or more
#  user-defined animations, each consisting of one or more stock BETSE-defined
#  plotters, to be defined. Users would then be able to construct arbitrarily
#  simple or complex animations as required.
#
#So, yes. It's quite a bit of work. But it's absolutely essential as well,
#particularly for implementing a general-purpose BETSE GUI.

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
import numpy as np
from abc import abstractmethod
from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.lib.matplotlib.writer import video as video_writers
from betse.lib.matplotlib.writer.frames import FileFrameWriter
from betse.science.plot.plotabc import PlotCells
from betse.util.io.log import logs
from betse.util.path import dirs, paths
from betse.util.type.types import type_check, Sequence
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation
from scipy import interpolate

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
        streamlines on this animation _or_ `False` otherwise. By design, this
        boolean is `True` if and only if the following are all also `True`:
        * The `p.I_overlay` boolean, implying the current simulation
          configuration to request current overlays.
        * The `p.calc_J` boolean, implying the current simulation
          configuration to model such currents.
        * The `is_current_overlayable` boolean parameter passed to the
          `__init__()` method of this class, implying the current animation to
          support current overlays.
    _is_overlaying_current_gj_only : bool
        `True` only if overlaying intracellular current _or_ `False` otherwise
        (i.e., if overlaying both intra- and extracellular current). Ignored
        unless overlaying current (i.e., if `_is_overlaying_current` is `True`).
    _is_saving_shown_frames : bool
        `True` only if both saving _and_ displaying animation frames.
    _save_frame_template : str
        `str.format()`-formatted template which, when formatted with the 0-based
        index of the current frame, yields the absolute path of the image file
        to be saved for that frame.
    _writer_frames : MovieWriter
        Matplotlib object saving animation frames as images if doing so _or_
        `None` otherwise.
    _writer_savefig_kwargs : dict
        Dictionary of all keyword arguments to be passed to the
        `Figure.savefig()` method called to save each animation frame for both
        images and video.
    _writer_video : MovieWriter
        Matplotlib object saving animation frames as video if doing so _or_
        `None` otherwise.
    '''

    # ..................{ LIFECYCLE                          }..................
    @type_check
    def __init__(
        self,
        save_dir_parent_basename: str,
        is_current_overlayable: bool,

        #FIXME: For orthogonality, rename to "is_current_overlay_gj_only" and
        #the corresponding instance attributes similarly.
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
            simulation configuration.
        is_current_overlayable : bool
            `True` if overlaying either electric current or concentration flux
            streamlines on this animation when requested by the current
            simulation configuration (as governed by the `p.I_overlay` and
            `p.calc_J` parameters) _or_ `False` otherwise. All subclasses except
            those already plotting streamlines (e.g., by calling the superclass
            `_plot_stream()` method) should unconditionally enable this boolean.
        is_overlaying_current_gj_only : optional[bool]
            `True` if only overlaying intracellular current _or_ `False` if
            overlaying both intra- and extracellular current. Ignored if current
            is _not_ being overlayed at all (i.e., if `_is_overlaying_current`
            is `False`). If `None`, defaults to the following state:
            * `False` if extracellular spaces are enabled _and_ both
               intracellular and extracellular current is to be animated.
            * `True` if either extracellular spaces are disabled _or_ are
               enabled but only intracellular current is to be animated.
        is_ecm_required : optional[bool]
            `True` if this animation is specific to extracellular spaces or
            `False` otherwise. If `True` and extracellular spaces are currently
            disabled, an exception is raised. Defaults to `False`.
        '''

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # If this subclass requires extracellular spaces but extracellular
        # spaces are currently disabled, raise an exception.
        if is_ecm_required and not self._p.sim_ECM:
            raise BetseExceptionParameters(
                'Animation "{}" requires extracellular spaces, which are '
                'disabled by the current simulation configuration.'.format(
                self._label))

        # Default unpassed parameters.
        if is_overlaying_current_gj_only is None:
            is_overlaying_current_gj_only = not (
                self._p.sim_ECM and self._p.IecmPlot)

        # Classify defaulted parameters.
        self._is_overlaying_current_gj_only = is_overlaying_current_gj_only

        # If this subclass requests a current overlay, do so only if:
        #
        # * Requested by the current simulation configuration via "p.I_overlay".
        # * This configuration is modelling currents via "p.calc_J".
        self._is_overlaying_current = (
            is_current_overlayable and self._p.I_overlay)

        # True if both saving and displaying animation frames.
        self._is_saving_showing = self._is_showing and self._is_saving

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
                '%s animation "%s"...', animation_verb, self._label)

        # Classify attributes to be possibly redefined below.
        self._current_density_magnitude_time_series = None
        self._current_density_x_time_series = None
        self._current_density_y_time_series = None
        self._current_density_stream_plot = None
        self._time_step = 0
        self._writer_frames = None
        self._writer_video = None

        # If saving animations, prepare to do so.
        if self._is_saving:
            self._init_saving(save_dir_parent_basename=save_dir_parent_basename)

        # If overlaying current, prepare to do so.
        if self._is_overlaying_current:
            self._init_current_density()


    @type_check
    def _init_saving(
        self,
        save_dir_parent_basename: str,
    ) -> None:
        '''
        Initialize this animation for platform-compatible file saving if enabled
        by the current simulation configuration _or_ noop otherwise.

        Parameters
        ----------
        save_dir_parent_basename : str
            Basename of the parent directory of the subdirectory to which this
            animation's frames will be saved when requested by the current
            simulation configuration. Defaults to a suitably generic basename.
        '''

        # Ensure that the passed directory and file basenames are actually
        # basenames and hence contain no directory separators.
        paths.die_unless_basename(self._label)

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
                'Animation saving unsupported during the "{}" loop.'.format(
                    plot_type))

        #FIXME: Pull the image filetype from the current YAML configuration
        #rather than coercing use of ".png".
        #FIXME: To do so, let's:
        #
        #* Create a new "betse.science.plot.anim.config" submodule.
        #* Create a new "AnimConfig" class in that submodule.
        #* Define a AnimConfig.make() factory method ala the existing
        #  betse.science.tissue.picker.TissuePicker.make() factory method. This
        #  method should deserialize *ALL* settings in the
        #  "results options/save animations" YAML subsection.
        #* Refactor all of the following "save_*" locals into public
        #  "AnimConfig" attributes, which should then be accessed here directly.
        #* Call that method in the "Parameters" class, capturing the result to a
        #  new "Parameters.anim" instance variable.
        #FIXME: Ah, right-o. We already appear to have existing classes
        #encapsulating these configuration settings in the "saver" submodule.
        #Frankly, however, that submodule is somewhat too fine-grained. We
        #probably really only need a single "AnimConfig" class, as suggested
        #above, rather than a deep class hierarchy. Consider simplifying.

        # `True` only if saving all frames of this animation to disk as images.
        is_saving_frames = True

        # `True` only if saving all frames of this animation to disk as a video.
        # is_saving_video = True
        is_saving_video = False

        # Filetype of each frame image to be saved for this animation. Ignored
        # if `is_saving_frames` is `False`.
        save_frame_filetype = 'png'

        # Dots per inch (DPI) of each frame image to be saved for this
        # animation. Ignored if `is_saving_frames` is `False`.
        self._writer_frames_dpi = mpl_config.get_rc_param('savefig.dpi')

        # Filetype of the video to be saved for this animation. Ignored if
        # `is_saving_video` is `False`.
        save_video_filetype = 'mp4'

        # Dots per inch (DPI) of each frame of the video to be saved for this
        # animation. Ignored if `is_saving_frames` is `False`.
        self._writer_video_dpi = mpl_config.get_rc_param('savefig.dpi')

        #FIXME: Set a sane default. There appear to be two Matplotlib defaults
        #for FPS, depending on object instantiation:
        #
        #* For writers instantiated directly, FPS defaults to 5. (Bit low, no?)
        #* For writers instantiated indirectly, FPS defaults to a formula
        #  presumably intelligently depending on animation properties:
        #
        #      # Convert interval in ms to frames per second
        #      fps = 1000. / self._interval
        #
        #While we have no idea how the latter works, it's probably more useful.

        # Frames per second (FPS) of the video to be saved for this animation.
        # Ignored if `is_saving_video` is `False`.
        save_video_fps = None

        #FIXME: Set a sane default.
        # Dictionary mapping from each Matplotlib-specific metadata name to that
        # metadata's string value of the video to be saved for this animation.
        # Ignored if `is_saving_video` is `False`. Metadata names of common
        # interest include: "title", "artist", "genre", "subject", "copyright",
        # "srcform", and "comment".
        save_video_metadata = {}
        # save_video_metadata = {artist: 'Me'}

        # Bitrate of the video to be saved for this animation. Ignored if
        # `is_saving_video` is `False`.
        save_video_bitrate = mpl_config.get_rc_param('animation.bitrate')

        # Name of the codec with which to encode the video to be saved for this
        # animation. Ignored if `is_saving_video` is `False`.
        save_video_codec = mpl_config.get_rc_param('animation.codec')

        # List of the Matplotlib-specific names of all external encoders
        # supported by Matplotlib with which to encode this video (in descending
        # order of preference), automatically selecting the first encoder on the
        # current $PATH. Ignored if `is_saving_video` is `False`.
        save_video_writer_names = ['ffmpeg', 'avconv', 'mencoder']

        # If saving animation frames as either images or video, prepare to do so
        # in a manner common to both.
        if is_saving_frames or is_saving_video:
            # Dictionary of all keyword arguments to be passed to the
            # `Figure.savefig()` method called to save each animation frame for
            # both images and video.
            self._writer_savefig_kwargs = {
                # Plot the background of each animation frame as transparent
                # rather than pure white.
                'transparent': True,
            }

            #FIXME: Refactor all calls to os.makedirs() everywhere similarly.

            # Path of the subdirectory to which these files will be saved,
            # creating this subdirectory and all parents thereof if needed.
            save_dirname = paths.join(
                loop_dirname, save_dir_parent_basename, self._label)
            save_dirname = dirs.canonicalize_and_make_unless_dir(save_dirname)

        # If saving animation frames as images, prepare to do so.
        if is_saving_frames:
            #FIXME: This currently defaults to padding frames with six or seven
            #zeroes, on average. Let's make this a bit more aesthetic by padding
            #frames to only as many zeroes are absolutely required by the
            #current frame count. To do that, in turn, we'll probably need to
            #shift everything that follows in this method to the _animate()
            #method, where the actual frame count is finally passed.

            # Template expanding to the basename of each image to be saved.
            # The "FileFrameWriter" class subsequently expands the "{{"- and
            # "}}"-delimited substring to the 0-based index of the current
            # frame number.
            save_frame_template_basename = '{}_{{:07d}}.{}'.format(
                self._label, save_frame_filetype)

            # Template expanding to the absolute path of each image to be saved.
            self._writer_frames_template = paths.join(
                save_dirname, save_frame_template_basename)

            # Object writing animation frames as images.
            self._writer_frames = FileFrameWriter()

            #FIXME: Perform this for video as well.
            # If both saving and displaying animation frames, prepare to do so.
            # If only saving but *NOT* displaying animation frames, the setup()
            # method called below is already called by the MovieWriter.saving()
            # method called by the Anim.save() method called by the _animate()
            # method called below. (No comment on architectural missteps.)
            #
            # See the _save_frame() method for horrid discussion.
            if self._is_saving_showing:
                # Log this saving attempt.
                logs.log_debug(
                    'Saving animation frames "%s"...',
                    self._writer_frames_template)

                # Prepare to save these animation frames.
                self._writer_frames.setup(
                    fig=self._figure,
                    outfile=self._writer_frames_template,
                    dpi=self._writer_frames_dpi,
                )

        #FIXME: FFMpeg integration appears to be extremely fragile. To combat
        #this, consider enabling an FFMpeg-specific debug log with the
        #following (...although we have no idea where the log actually lives):
        #
        #    matplotlib.rcParams['animation.ffmpeg_args'] = '-report'

        # If saving animation frames as video, prepare to do so.
        if is_saving_video:
            # Class writing animation frames as video.
            VideoWriterClass = video_writers.get_first_class(
                save_video_writer_names)
            # print('found video writer: {}'.format(VideoWriterClass))

            # Basename of the video to be written.
            save_video_basename = '{}.{}'.format(
                self._label, save_video_filetype)

            # Absolute path of the video to be written.
            self._writer_video_filename = paths.join(
                save_dirname, save_video_basename)

            # Object writing animation frames as video.
            self._writer_video = VideoWriterClass(
                fps=save_video_fps,
                bitrate=save_video_bitrate,
                codec=save_video_codec,
                metadata=save_video_metadata,
            )

            # If both saving and displaying animation frames, prepare as above.
            if self._is_saving_showing:
                # Log this saving attempt.
                logs.log_debug(
                    'Saving animation video "%s"...',
                    self._writer_video_filename)

                # Prepare to save this animation video.
                self._writer_video.setup(
                    fig=self._figure,
                    outfile=self._writer_video_filename,
                    dpi=self._writer_video_dpi,
                )


    # This method has been overridden to support subclasses that manually handle
    # animations rather than calling the _animate() method (e.g., the
    # "AnimCellsWhileSolving" subclass).
    def _prep_figure(self, *args, **kwargs) -> None:

        super()._prep_figure(*args, **kwargs)

        # If overlaying current, do so.
        if self._is_overlaying_current:
            self._plot_current_density()

    # ..................{ PROPERTIES                         }..................
    @property
    def _is_saving(self) -> bool:
        return self._p.saveAnimations

    # ..................{ ANIMATORS                          }..................
    @type_check
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

        # Prepare for plotting immediately *BEFORE* plotting the first frame.
        self._prep_figure(*args, **kwargs)

        #FIXME: For efficiency, we should probably be passing "blit=True," to
        #FuncAnimation(). Lemon grass and dill!

        # Create and assign an animation function to a local variable. If the
        # latter is *NOT* done, this function will be garbage collected prior
        # to subsequent plot handling -- in which case only the first plot will
        # be plotted without explicit warning or error. Die, matplotlib! Die!!
        self._anim = FuncAnimation(
            # Figure to which the "func" callable plots each frame.
            fig=self._figure,

            # Callable plotting each frame.
            func=self.plot_frame,

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

                # If saving animation frames as images, finalize doing so.
                if self._writer_frames is not None:
                    self._writer_frames.finish()

                # If saving animation frames as video, finalize doing so.
                if self._writer_video is not None:
                    self._writer_video.finish()
            # Else if saving animation frames as...
            elif self._is_saving:
                #FIXME: If soving both frames and video, the current approach
                #does technically work but is highly inefficient. It appears to
                #recreate the animation in-memory twice! A sane alternative is:
                #
                #* If soving both frames and video:
                #  * Create only "self._writer_video". Hence,
                #    "self._writer_frames" should remain None.
                #  * Configure "self._writer_video" to write frames. This may
                #    require us to avoid the more efficient pipe-based method of
                #    writing video in favour of a temporary file-based method.
                #    If this is the case, we may also need to manually move all
                #    temporary files written by "self._writer_video" from the
                #    temporary directory to which they are written into the
                #    desired target directory.
                #
                #This should be feasible. Research is required, clearly.
                #FIXME: Actually, there's a far simpler way requiring no changes
                #to our existing methodology -- which is good. The core issue is
                #that the Anim.save() and Anim.saving() methods are overly
                #clumsy wrappers around core "MovieWriter" functions. Since we
                #already call the latter directly in our _plot_frame() method,
                #there's absolutely no incentive to calling the former here.
                #Instead, refactor this approach as follows:
                #
                #* Unconditionally call the MovieWriter.setup() method above.
                #  Specifically:
                #
                #      # Refactor conditionals like this...
                #      if self._is_saving_showing:
                #          self._writer_frames.setup(...)
                #
                #      # ...to merely this.
                #      self._writer_frames.setup(...)
                #
                #* Replace the calls to Anim.save() below with either:
                #  * Implicit iteration leveraging the existing pyplot.show()
                #    method call above by hiding this figure. I have no idea
                #    whether this actually works, but it seems worth a try.
                #    This is low-hanging fruit and hence preferable:
                #
                #      self._figure.canvas.Show(False)
                #      pyplot.show()
                #
                #    What happens? We have no idea. Simple to test, happily.
                #
                #  * If the implicit iteration approach fails, explicit
                #    iteration by essentially duplicating the iteration
                #    performed by the Anim.save() function below. While also
                #    trivial, the extra effort involved makes the above
                #    approach preferable.

                # ...images, do so. Due to non-orthogonality in the Matplotlib
                # API, all encoding parameters *EXCEPT* dots per inch (DPI) are
                # passable to the "self._writer_frames" constructor. DPI,
                # however, is *ONLY* passable to the save() method called below.
                if self._writer_frames is not None:
                    self._anim.save(
                        writer=self._writer_frames,
                        filename=self._writer_frames_template,
                        dpi=self._writer_frames_dpi,
                        savefig_kwargs=self._writer_savefig_kwargs,
                    )

                # ...video, do so.
                if self._writer_video is not None:
                    self._anim.save(
                        writer=self._writer_video,
                        filename=self._writer_video_filename,
                        dpi=self._writer_video_dpi,
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
    @type_check
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

        # Log this animation frame.
        logs.log_debug(
            'Plotting animation "%s" frame %d...',
            self._label,
            len(self._sim.time) if frame_number == -1 else frame_number)

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
        if self._is_saving_showing:
            # If saving animation frames as images, save this frame as such.
            if self._writer_frames is not None:
                self._writer_frames.grab_frame(**self._writer_savefig_kwargs)

            # If saving animation frames as video, save this frame as such.
            if self._writer_video is not None:
                self._writer_video.grab_frame(**self._writer_savefig_kwargs)


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
            I_grid_x_time = []
            I_grid_y_time = []

            # Interpolate data from cell centres to the xy-grid.
            cell_centres = (
                self._cells.cell_centres[:, 0], self._cells.cell_centres[:, 1])
            cell_grid = (self._cells.X, self._cells.Y)

            for i in range(0, len(self._sim.I_cell_x_time)):
                I_gj_x = self._cells.maskECM * interpolate.griddata(
                    cell_centres,
                    self._sim.I_cell_x_time[i],
                    cell_grid,
                    fill_value=0,
                    method=self._p.interp_type,
                )
                I_grid_x_time.append(I_gj_x)

                I_gj_y = self._cells.maskECM * interpolate.griddata(
                    cell_centres,
                    self._sim.I_cell_y_time[i],
                    cell_grid,
                    fill_value=0,
                    method=self._p.interp_type,
                )
                I_grid_y_time.append(I_gj_y)

            self._current_density_x_time_series = I_grid_x_time
            self._current_density_y_time_series = I_grid_y_time

        else:
            self._current_density_x_time_series = self._sim.I_tot_x_time
            self._current_density_y_time_series = self._sim.I_tot_y_time

        # Time series of all current density magnitudes (i.e., `Jmag_M`),
        # multiplying by 100 to obtain current density in units of uA/cm2.
        self._current_density_magnitude_time_series = 100*np.sqrt(
            np.asarray(self._current_density_x_time_series) ** 2 +
            np.asarray(self._current_density_y_time_series) ** 2) + 1e-15


    def _plot_current_density(self) -> None:
        '''
        Overlay the first frame of this subclass' animation with a streamplot of
        either electric current or concentration flux.
        '''

        # If animating only intracellular current, do so.
        if self._is_overlaying_current_gj_only:
            self._axes_title = 'Intracellular Current'

            # #FIXME: Is there any point to this? From what we can tell, the
            # #"self._current_density_stream_plot" will simply be outright
            # #replaced for the first and all subsequent frames. Galloping fish!
            # self._current_density_stream_plot, self._axes = cell_stream(
            #     self._current_density_x_time_series[-1],
            #     self._current_density_y_time_series[-1],
            #     self._axes, self._cells, self._p)
        # If animating both extracellular and intracellular current, do so.
        else:
            self._axes_title = 'Total Current Overlay'

            # #FIXME: Likewise.
            # self._current_density_stream_plot, self._axes = env_stream(
            #     self._current_density_x_time_series[-1],
            #     self._current_density_y_time_series[-1],
            #     self._axes, self._cells, self._p)

    @type_check
    def _replot_current_density(self, frame_number: int) -> None:
        '''
        Overlay the passed frame of this subclass' animation with a streamplot
        of either electric current or concentration flux.

        Parameters
        -----------
        frame_number : int
            0-based index of the frame to be plotted.
        '''

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
#FIXME: Refactor as follows:
#
#* Rename the existing "anim" submodule of this subpackage to "after".
#* Define a new "while" submodule of this subpackage.
#* Shift this and the following subclasses to the "after" submodule.
#* Shift the "AnimCellsWhileSolving" subclass to the "while" submodule.
class AnimCellsAfterSolving(AnimCells):
    '''
    Out-of-place animation of an arbitrary membrane-centric time series (e.g.,
    cell Vmem as a function of time), plotted over the cell cluster _after_
    rather than _during_ simulation modelling.
    '''

    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Save out-of-place animations to a different parent directory than
            # that to which in-place animations are saved.
            save_dir_parent_basename='anim',
            *args, **kwargs
        )


class AnimField(AnimCellsAfterSolving):
    '''
    Abstract base class of all animations of electric field strength plotted on
    the current cell cluster.

    Attributes
    ----------
    _magnitude_time_series : Sequence
        Electric field magnitudes as a function of time.
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's electric field magnitude.
    _stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of the current or prior frame's electric field.
    _x_time_series : Sequence
        Electric field X components as a function of time.
    _y_time_series : Sequence
        Electric field Y components as a function of time.
    _unit_x_time_series : Sequence
        Electric field X unit components as a function of time. The resulting
        electric field vectors are **unit vectors** (i.e., have magnitude 1).
    _unit_y_time_series : Sequence
        Electric field Y unit components as a function of time. The resulting
        electric field vectors are **unit vectors** (i.e., have magnitude 1).
    '''

    @type_check
    def __init__(
        self,
        x_time_series: Sequence,
        y_time_series: Sequence,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        x_time_series : Sequence
            Sequence (e.g., list, numpy array) of all electric field strength X
            components indexed by simulation time.
        y_time_series : Sequence
            Sequence (e.g., list, numpy array) of all electric field strength Y
            components indexed by simulation time.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class already plots a streamplot, prevent the
            # superclass from plotting another streamplot as an overlay.
            is_current_overlayable=False,
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


class AnimVelocity(AnimCellsAfterSolving):
    '''
    Abstract base class of all animations of a velocity flow plotted on the
    current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class already plots a streamplot, prevent the
            # superclass from plotting another streamplot as an overlay.
            is_current_overlayable=False,
            *args, **kwargs)

        # Prefer an alternative colormap.
        self._colormap = self._p.background_cm
