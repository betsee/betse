#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based animation classes.
'''

#FIXME: We should probably animate non-blockingly (e.g., by passing
#"block=False" to the plt.show() command. Since we don't, I suspect there are
#probably some subtle dragons here. Let's give it a go! Thunderous tales of woe!

# ....................{ IMPORTS                            }....................
import numpy as np
from abc import ABCMeta, abstractmethod  #, abstractstaticmethod
from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib.mpl import mpl_config
from betse.lib.matplotlib.anim import FileFrameWriter
from betse.util.io import loggers
from betse.util.path import dirs, paths
from betse.util.type import types
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation

# ....................{ BASE                               }....................
#FIXME: Privatize all public attributes declared below. Raving river madness!
#FIXME: Rename all following parameters and attributes:
#* "clrAutoscale" to "is_colorbar_autoscaled".
#* "clrMin" to "colorbar_min".
#* "clrMax" to "colorbar_max".

class Anim(object, metaclass=ABCMeta):
    '''
    Abstract base class of all animation classes.

    Instances of this class animate the spatial distribution of modelled
    variables (e.g., Vmem) over all time steps of the simulation.

    Attributes
    ----------
    sim : Simulation
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    _anim : FuncAnimation
        Low-level Matplotlib animation object instantiated by this high-level
        BETSE wrapper object.
    ax : FigureAxes
        Matplotlib figure axes providing the current animation frame data.
    _axes_bounds : list
        Spacial extent of the current 2D environment as a 4-element list
        conisting of (in order):
        1. Minimum value of the figure's X axis.
        2. Maximum value of the figure's X axis.
        3. Minimum value of the figure's Y axis.
        4. Maximum value of the figure's Y axis.
    _colorbar_mapping : object
        The Matplotlib mapping (e.g., `Image`, `ContourSet`) to which this
        animation's colorbar applies.
    _colorbar_title: str
        Text displayed above the figure colorbar.
    clrAutoscale : bool
        `True` if dynamically resetting the minimum and maximum colorbar values
        to be the corresponding minimum and maximum values for the current
        frame _or_ `False` if statically setting the minimum and maximum
        colorbar values to predetermined constants.
    clrMin : float
        Minimum colorbar value to be used. If `clrAutoscale` is `True`, the
        subclass is responsible for redefining this value as appropriate.
    clrMax : float
        Maximum colorbar value to be used. If `clrAutoscale` is `True`, the
        subclass is responsible for redefining this value as appropriate.
    colormap : Colormap
        Matplotlib colormap with which to create this animation's colorbar.
    fig : Figure
        Matplotlib figure providing the current animation frame.
    _figure_title : str
        Text displayed above the figure itself.
    _is_saving_plotted_frames : bool
        `True` if both saving and displaying animation frames _or_ `False`
        otherwise.
    _save_frame_template : str
        `str.format()`-formatted template which, when formatted with the 0-based
        index of the current frame, yields the absolute path of the image file
        to be saved for that frame.
    _axes_title : str
        Text displayed above the figure axes. If a non-`None` value for the
        `axes_title` parameter is passed to the `animate()` method, this is that
        value; else, this is the value of the `figure_title` parameter passed to
        the same method.
    _type : str
        Basename of the subdirectory in the phase-specific results directory
        to which all animation files will be saved _and_ the basename prefix of
        these files.
    _writer_frames : MovieWriter
        Object writing frames from this animation to image files if enabled _or_
        `None` otherwise.
    _writer_video : MovieWriter
        Object encoding this animation to a video file if enabled _or_ `None`
        otherwise.
    '''

    # ..................{ PRIVATE ~ init                     }..................
    def __init__(
        self,

        # Mandatory parameters.
        sim: 'Simulator',
        cells: 'Cells',
        p: 'Parameters',
        type: str,
        figure_title: str,
        colorbar_title: str,
        clrAutoscale: bool,
        clrMin: float,
        clrMax: float,

        # Optional parameters.
        colormap: 'Colormap' = None,
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        type : str
            Basename of the subdirectory in the phase-specific results directory
            to which all animation files will be saved _and_ the basename prefix
            of these files.
        figure_title : str
            Text displayed above the figure itself.
        colorbar_title: str
            Text displayed above the figure colorbar.
        clrAutoscale : bool
            `True` if dynamically resetting the minimum and maximum colorbar
            values to be the corresponding minimum and maximum values for the
            current frame _or_ `False` if statically setting the minimum and
            maximum colorbar values to predetermined constants.
        clrMin : float
            Minimum colorbar value to be used if `clrAutoscale` is `False`.
        clrMax : float
            Maximum colorbar value to be used if `clrAutoscale` is `False`.
        colormap : Colormap
            Matplotlib colormap to be used in this animation's colorbar or
            `None`, in which case the default colormap will be used.
        '''
        # Validate core parameters.
        assert types.is_simulator(sim), types.assert_not_simulator(sim)
        assert types.is_cells(cells), types.assert_not_parameters(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)

        # Default unpassed parameters.
        if colormap is None:
            colormap = p.default_cm

        # Validate all remaining parameters *AFTER* defaulting parameters.
        assert types.is_str_nonempty(type), (
            types.assert_not_str_nonempty(type, 'Animation type'))
        assert types.is_str_nonempty(figure_title), (
            types.assert_not_str_nonempty(figure_title, 'Figure title'))
        assert types.is_str_nonempty(colorbar_title), (
            types.assert_not_str_nonempty(colorbar_title, 'Colorbar title'))
        assert types.is_bool(clrAutoscale), types.assert_not_bool(clrAutoscale)
        assert types.is_numeric(clrMin), types.assert_not_numeric(clrMin)
        assert types.is_numeric(clrMax), types.assert_not_numeric(clrMax)
        assert types.is_matplotlib_colormap(colormap), (
            types.assert_not_matplotlib_colormap(colormap))

        # Classify *AFTER* validating parameters.
        self.sim = sim
        self.cells = cells
        self.p = p
        self._type = type
        self._figure_title = figure_title
        self._colorbar_title = colorbar_title
        self.clrAutoscale = clrAutoscale
        self.clrMin = clrMin
        self.clrMax = clrMax
        self.colormap = colormap

        # Classify attributes to be subsequently defined.
        self._axes_title = None
        self._colorbar_mapping = None
        self._writer_frames = None
        self._writer_video = None

        #FIXME: Abandon "pyplot", all who enter here!
        #FIXME: Call the
        #betse.lib.matplotlib.mpl.mpl_config.make_backend_figure() method rather
        #than pyplot.figure() here and everywhere.
        #FIXME: Actually, the
        #betse.lib.matplotlib.mpl.mpl_config.make_backend_figure() method
        #doesn't appear to be terribly helpful as one can't actually call show()
        #on that. The closest approximation is figure.canvas.print_figure(), but
        #that's intended principally for producing physical print copies (as the
        #name suggests) and mandates blocking in any case. Instead, what we
        #*REALLY* want to be doing is calling the current backend's
        #new_figure_manager() function. This creates a figure manager containing
        #the desired figure and figure canvas objects, all initialized in a sane
        #manner conducive to displaying and saving. We think, anyway!
        #
        #Given a "FigureManagerBase" instance, we can then call the current
        #backend's show() callable to do what we want. This callable accepts the
        #same exact parameters as the "pyplot.show()" function -- which is no
        #coincidence, as the latter defers to the former.
        #
        #Wait! We should just be able to call the show() method of the new
        #"FigureManagerBase" instance. The ShowBase.__call__() just defers to
        #that. By default, everything is non-blocking. If we want blocking
        #behavior (which we probably don't), then we'd probably want to call
        #the current backend's show() callable instead. Annnnnyway, let's just
        #go with the FigureManagerBase.show() approach for now.
        #FIXME: Since we're subverting the "pyplot" API, we'll want to make sure
        #that we're not opening too many figures: e.g.,
        #
        #    if len(list(six.iterkeys(_pylab_helpers.Gcf.figs))) >= 20:
        #         print('Uh oh!)
        #
        #The question is whether "_pylab_helpers.Gcf" applies if you only create
        #figure managers outside of the "pyplot" API. We suspect it does, but...
        #Oh, and note that the pyplot.get_fignums() function should probably be
        #called instead of the above code. It's likely "Gcf" will be moved away,
        #at some point.

        # Figure encapsulating this animation.
        self.fig = pyplot.figure()

        # Extent of the current 2D environment.
        self._axes_bounds = [
            self.cells.xmin * self.p.um,
            self.cells.xmax * self.p.um,
            self.cells.ymin * self.p.um,
            self.cells.ymax * self.p.um,
        ]

        # Figure axes scaled to the extent of the current 2D environment.
        self.ax = pyplot.subplot(111)
        self.ax.axis('equal')
        self.ax.axis(self._axes_bounds)

        # Initialize animation saving *AFTER* defining all attribute defaults.
        self._init_saving()


    def _init_saving(self) -> None:
        '''
        Initialize this animation for platform-compatible file saving if enabled
        by the current simulation configuration or noop otherwise.
        '''

        # True if both saving and displaying animation frames.
        self._is_saving_plotted_frames = (
            self.p.turn_all_plots_off is False and
            self.p.saveAnimations is True)

        # If animation saving is disabled, noop.
        if self.p.saveAnimations is False:
            return

        # Ensure that the passed directory and file basenames are actually
        # basenames and hence contain no directory separators.
        paths.die_unless_basename(self._type)

        # Path of the phase-specific parent directory of the subdirectory to
        # which these files will be saved.
        phase_dirname = None
        if self.p.plot_type == 'sim':
            phase_dirname = self.p.sim_results
        elif self.p.plot_type == 'init':
            phase_dirname = self.p.init_results
        else:
            raise BetseExceptionParameters(
                'Anim saving unsupported during the "{}" phase.'.format(
                    self.p.plot_type))

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

        # If both saving and displaying animation frames, prepare for doing so.
        # See the _save_frame() method for horrid discussion.
        if self._is_saving_plotted_frames:
            self._writer_frames.setup(
                fig=self.fig,
                outfile=self._save_frame_template,

                #FIXME: Pass the actual desired "dpi" parameter.
                dpi=mpl_config.get_rc_param('savefig.dpi'),
            )

    # ..................{ PRIVATE ~ animate                  }..................
    #FIXME: The "axes_x_label" and "axes_y_label" parameters should probably
    #just be passed to the __init__() method and then classified. Doing so
    #would simplify "AnimateField" subclasses by permitting that superclass to
    #define values for these parameters common to these subclasses. Bemusement!
    def _animate(
        self,
        frame_count: int,
        colorbar_mapping: object,
        axes_x_label: str,
        axes_y_label: str,
        axes_title: str = None,
        colorbar_values: np.ndarray = None,
    ) -> None:
        '''
        Display this animation if the current simulation configuration
        requests plots to be displayed or noop otherwise.

        This method is intended to be called as the last statement in the
        `__init__()` method of all subclasses of this superclass.

        Attributes
        ----------
        frame_count : int
            Number of frames to be animated.
        axes_title : str
            Optional text displayed above the figure axes but below the figure
            title (i.e., `_figure_title`) _or_ `None` if no such text is to be
            displayed. Defaults to `None`.
        axes_x_label : str
            Text displayed below the figure's X axis.
        axes_y_label : str
            Text displayed to the left of the figure's Y axis.
        colorbar_mapping : object
            The Matplotlib mapping (e.g., `Image`, `ContourSet`) to which this
            colorbar applies.
        colorbar_values : np.ndarray
            Optional multi-dimensional Numpy array containing all data values
            to be animated _or_ `None` if calculating this data during the
            animation initialization is infeasible or impractical (e.g., due to
            space and time constraints). If non-`None` _and_ colorbar
            autoscaling is requested (i.e., the initialization-time
            `clrAutoscale` parameter was `True`), the colorbar will be clipped
            to the maximum and minimum value in this matrix; else, the subclass
            is responsible for colorbar autoscaling. Defaults to `None`.
        '''
        assert types.is_int(frame_count), types.assert_not_int(frame_count)
        assert types.is_str_nonempty(axes_x_label), (
            types.assert_not_str_nonempty(axes_x_label, 'X axis label'))
        assert types.is_str_nonempty(axes_y_label), (
            types.assert_not_str_nonempty(axes_y_label, 'Y axis label'))

        #FIXME: No! Cease classifying this, please. Some animation subclasses
        #subsequently redefine the colorbar mapping. They shouldn't have to know
        #that they need to also reset this attribute to safely do so.

        # Classify passed parameters.
        self._colorbar_mapping = colorbar_mapping

        # If labelling each plotted cell with that cell's unique 0-based index,
        # do so.
        if self.p.enumerate_cells is True:
            for cell_index, cell_centre in enumerate(self.cells.cell_centres):
                self.ax.text(
                    self.p.um*cell_centre[0],
                    self.p.um*cell_centre[1],
                    cell_index,
                    va='center',
                    ha='center',
                )

        # If both a figure and axes title are defined, display the figure title
        # as such above the axes title.
        if axes_title is not None:
            self._axes_title = axes_title
            self.fig.suptitle(
                self._figure_title, fontsize=14, fontweight='bold')
        # Else, display the figure title as the axes title.
        else:
            self._axes_title = self._figure_title

        assert types.is_str_nonempty(self._axes_title), (
            types.assert_not_str_nonempty(self._axes_title, 'Axis title'))

        # Display the axes title and labels.
        self.ax.set_title(self._axes_title)
        self.ax.set_xlabel(axes_x_label)
        self.ax.set_ylabel(axes_y_label)

        # If a time series is passed *AND* colorbar autoscaling is requested,
        # clip the colorbar to the minimum and maximum values of this series.
        if colorbar_values is not None and self.clrAutoscale is True:
            assert types.is_sequence_nonstr(colorbar_values), (
                types.assert_not_sequence_nonstr(colorbar_values))

            # Flatten this two-dimensional matrix to a one-dimensional array,
            # providing efficient retrieval of minimum and maximum values.
            time_series_flat = np.ravel(colorbar_values)

            # Minimum and maximum values.
            self.clrMin = np.min(time_series_flat)
            self.clrMax = np.max(time_series_flat)

            # If these values are identical, coerce them to differ. This ensures
            # that the colorbar will be displayed in a sane manner.
            if self.clrMin == self.clrMax:
                self.clrMin = self.clrMin - 1
                self.clrMax = self.clrMax + 1

        # Set the colorbar range.
        colorbar_mapping.set_clim(self.clrMin, self.clrMax)

        # Display the colorbar.
        colorbar = self.fig.colorbar(colorbar_mapping)
        colorbar.set_label(self._colorbar_title)

        #FIXME: For efficiency, we should probably be passing "blit=True," to
        #FuncAnimation(). Lemon grass and dill!

        # Create and assign an animation function to a local variable. If the
        # latter is *NOT* done, this function will be garbage collected prior
        # to subsequent plot handling -- in which case only the first plot will
        # be plotted without explicit warning or error. Die, matplotlib! Die!!
        self._anim = FuncAnimation(
            self.fig, self._plot_frame,

            # Number of frames to be animated.
            frames=frame_count,

            # Delay in milliseconds between consecutive frames.
            interval=100,

            # Indefinitely repeat this animation unless saving animations, as
            # doing so under the current implementation would repeatedly (and
            # hence unnecessarily) overwrite previously written files.
            repeat=not self.p.saveAnimations,
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
            #FIXME: Refactor to *NOT* use the "pyplot" API. Juggling hugs!
            # If displaying animations, do so.
            if self.p.turn_all_plots_off is False:
                loggers.log_info(
                    'Plotting animation "{}"...'.format(self._type))

                # Display this animation.
                pyplot.show()
            # Else if saving animation frames, do so.
            elif self.p.saveAnimations is True:
                loggers.log_info(
                    'Saving animation "{}" frames...'.format(self._type))

                #FIXME: Pass the "dpi" parameter as well.
                self._anim.save(
                    filename=self._save_frame_template,
                    writer=self._writer_frames)
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

        . Calls the subclass `_plot_frame_figure()` method.
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

        #FIXME: Nope! This is terrible. Do this in subclasses if we do this
        #anywhere. Ideally, this will no longer be necessary *AFTER* we refactor
        #subclasses to precalculate the global "clrMin" and "clrMax" values.

        # For safety, rescale the colorbar regardless of whether or not the
        # minimum or maximum colorbar values have been modified by the above
        # call to the subclass _plot_frame_figure() method.
        self._colorbar_mapping.set_clim(self.clrMin, self.clrMax)

        # Update this figure with the current time, rounded to three decimal
        # places for readability.
        self.ax.set_title('{} (time {:.3f}s)'.format(
            self._axes_title, self.sim.time[frame_number]))

        # If both saving and displaying animation frames, save this frame. Note
        # that if only saving but *NOT* displaying animations, this frame will
        # be handled by our _animate() method. Why? Because Matplotlib will fail
        # to iterate frames and hence call our _plot_next_frame() method calling
        # this method *UNLESS* our _animate() method explicitly calls the
        # FuncAnimation.save() method with the writer name "frame" signifying
        # our "FileFrameWriter" class to do so. (Look. It's complicated, O.K.?)
        # if self.p.saveAnimations is False:
        if self._is_saving_plotted_frames:
            self._writer_frames.grab_frame()


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

# ....................{ BASE ~ field                       }....................
class AnimField(Anim):
    '''
    Abstract base class of all animations of an electric field plotted on the
    current cell cluster.

    Attributes
    ----------
    _Fx_time : list
        List of all electric field strength X components indexed by
        simulation time.
    _Fy_time : list
        List of all electric field strength Y components indexed by
        simulation time.
    '''

    def __init__(
        self,
        Fx_time: list,
        Fy_time: list,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        Fx_time : list
            List of all electric field strength X components indexed by
            simulation time.
        Fy_time : list
            List of all electric field strength Y components indexed by
            simulation time.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(Fx_time), (
            types.assert_not_sequence_nonstr(Fx_time))
        assert types.is_sequence_nonstr(Fy_time), (
            types.assert_not_sequence_nonstr(Fy_time))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Classify all remaining parameters.
        self._Fx_time = Fx_time
        self._Fy_time = Fy_time

        # Prefer an alternative colormap.
        self.colormap = self.p.background_cm
