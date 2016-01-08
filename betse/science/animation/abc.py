#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based animation classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod  #, abstractstaticmethod
from betse.exceptions import BetseExceptionParameters
# from betse.lib.matplotlib import mpl
from betse.util.path import dirs, paths
from betse.util.type import types
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation

# ....................{ BASE                               }....................
#FIXME: Privatize all public attributes declared below. Raving river madness!
#FIXME: Document the "clrAutoscale", "clrMin", and "clrMax" attributes. Sizzle!
#FIXME: Refactor the "savedAni" attribute into a template containing exactly one
#"{"- and "}"-delimited string, to satisfy our "FrameFileWriter" API. Currently,
#this attribute is merely the prefix for all frame image files to be saved.
#FIXME: Rename:
#* The "saveFolder" attribute to "_save_dir_basename".
#* The "saveFile" attribute to "_save_file_basename_prefix".
#Or, alternately (and probably preferably) combine the two attributes into a
#pathname template consumable by our new "FrameFileWriter" class.

class Animation(object, metaclass=ABCMeta):
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
    clrAutoscale : bool
        ???.
    clrMin : float
        ???.
    clrMax : float
        ???.
    colormap : Colormap
        Matplotlib colormap with which to create this animation's colorbar.
    savedAni : str
        Path prefix for all frame image files to be saved.
    saveFile : str
        Basename prefix of all frame image files to be saved.
    saveFolder : str
        Basename of the subdirectory in the phase-specific results directory
        to which all animation results will be saved.
    '''

    # ..................{ ABSTRACT ~ private                 }..................
    @abstractmethod
    def _plot_next_frame(self, frame_number: int) -> None:
        '''
        Iterate the current animation to the next frame.

        This method is iteratively called by Matplotlib's `FuncAnimation()`
        class instantiated by our `_animate()` method. The subclass
        implementation of this abstract method is expected to modify this
        animation's figure and axes so as to display the next frame. It is
        _not_, however, expected to save that figure to disk; frame saving is
        already implemented by this base class in a general-purpose manner.

        Parameters
        ----------
        frame_number : int
            0-based index of the frame to be plotted.
        '''
        pass

    # ..................{ CONCRETE ~ init                    }..................
    def __init__(
        self,

        # Mandatory parameters.
        sim: 'Simulator',
        cells: 'Cells',
        p: 'Parameters',
        clrAutoscale: bool,
        clrMin: float,
        clrMax: float,
        saveFolder: str,
        saveFile: str,

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
        clrAutoscale : bool
            ???.
        clrMin : float
            ???.
        clrMax : float
            ???.
        saveFolder : str
            Basename of the subdirectory in the phase-specific results directory
            to which all animation results will be saved.
        saveFile : str
            Basename prefix of all frame image files to be saved.
        colormap : Colormap
            Matplotlib colormap to be used in this animation's colorbar.
        '''
        # Validate core parameters.
        assert types.is_simulator(sim), types.assert_not_simulator(sim)
        assert types.is_cells(cells), types.assert_not_parameters(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)

        # Default unpassed parameters.
        if colormap is None:
            colormap = p.default_cm

        # Validate all remaining parameters *AFTER* defaulting parameters.
        assert types.is_bool(clrAutoscale), types.assert_not_bool(clrAutoscale)
        assert types.is_numeric(clrMin), types.assert_not_numeric(clrMin)
        assert types.is_numeric(clrMax), types.assert_not_numeric(clrMax)
        assert types.is_matplotlib_colormap(colormap), (
            types.assert_not_matplotlib_colormap(colormap))
        assert types.is_str_nonempty(saveFolder), (
            types.assert_not_str_nonempty(saveFolder, 'Save directory'))
        assert types.is_str_nonempty(saveFile), (
            types.assert_not_str_nonempty(saveFile, 'Save frame file prefix'))

        # Classify *AFTER* validating parameters.
        self.sim = sim
        self.cells = cells
        self.p = p
        self.clrAutoscale = clrAutoscale
        self.clrMin = clrMin
        self.clrMax = clrMax
        self.colormap = colormap
        self.saveFolder = saveFolder
        self.saveFile = saveFile

        # Initialize animation saving *AFTER* defining all attribute defaults.
        self._init_saving()


    def _init_saving(self) -> None:
        '''
        Initialize this animation for platform-compatible file saving if enabled
        by the current simulation configuration or noop otherwise.
        '''

        # If animation saving is disabled, noop.
        if self.p.saveAnimations is False:
            return

        # Ensure that the passed directory and file basenames are actually
        # basenames and hence contain no directory separators.
        paths.die_unless_basename(self.saveFolder)
        paths.die_unless_basename(self.saveFile)

        # Path of the phase-specific parent directory of the subdirectory to
        # which these files will be saved.
        phase_dirname = None
        if self.p.plot_type == 'sim':
            phase_dirname = self.p.sim_results
        elif self.p.plot_type == 'init':
            phase_dirname = self.p.init_results
        else:
            raise BetseExceptionParameters(
                'Animation saving unsupported during the "{}" phase.'.format(
                    self.p.plot_type))

        #FIXME: Refactor all calls to os.makedirs() everywhere similarly.

        # Path of the subdirectory to which these files will be saved, creating
        # this subdirectory and all parents thereof if needed.
        images_dirname = paths.join(phase_dirname, 'animation', self.saveFolder)
        images_dirname = dirs.canonicalize_and_make_unless_dir(images_dirname)

        # Path of the file to be saved.
        self.savedAni = paths.join(images_dirname, self.saveFile)

    # ..................{ CONCRETE ~ animate                 }..................
    def _animate(self, frame_count: int) -> None:
        '''
        Display this animation if the current simulation configuration
        requests plots to be displayed or noop otherwise.

        Attributes
        ----------
        frame_count : int
            Number of frames to be animated.
        '''
        assert types.is_int(frame_count), types.assert_not_int(frame_count)

        # Create and assign an animation function to a local variable. If the
        # latter is *NOT* done, this function will be garbage collected prior
        # to subsequent plot handling -- in which case only the first plot will
        # be plotted without explicit warning or error. Die, matplotlib! Die!!
        ani = FuncAnimation(
            self.fig, self._plot_next_frame,

            # Number of frames to be animated.
            frames=frame_count,

            # Delay in milliseconds between consecutive frames.
            interval=100,

            #FIXME: Ideally, we should be able to support both with a simple
            #boolean. Tinder sticks in a fragrant bonfire!

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
            if self.p.turn_all_plots_off is False:
                #FIXME: Refactor to *NOT* use the "pyplot" API. Juggling hugs!
                pyplot.show()
            elif self.p.saveAnimations:
                #FIXME: Uncomment when worky.
                # loggers.log_info(
                #     'Saving animation frames "{}"...'.format(video_filename))
                # ani.save(filename=video_filename, writer='frame')
                pass
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


    #FIXME: Refactor to just defer to FileFrameWriter.grab_frame().
    def _save_frame(self, frame_number: int) -> None:
        '''
        Write the frame with the 0-based passed index to disk.

        Parameters
        ----------
        frame_number : int
            0-based index of the frame to be written.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # If both saving and displaying animations, save this frame. Note that
        # if only saving but *NOT* displaying animations, this frame will be
        # handled by our _animate() method. Why? Because Matplotlib will fail to
        # iterate frames and hence call our _plot_next_frame() method calling
        # this method *UNLESS* our _animate() method explicitly calls the
        # FuncAnimation.save() method with the writer name "frame" signifying
        # our "FrameFileWriter" class to do so. (Look. It's complicated, O.K.?)
        if self.p.saveAnimations is False:
        #FIXME: Uncomment when _animate() actually calls FuncAnimation.save().
        # if not (
        #     self.p.saveAnimations is True and
        #     self.p.turn_all_plots_off is True):
            return

        self.fig.canvas.draw()

        # Filename of the file to be written.
        savename = self.savedAni + str(frame_number) + '.png'

        # #FIXME: Remove this debug statement later.
        # print('Saving animated frame: {}'.format(savename))
        pyplot.savefig(savename, format='png')

