#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based animation subclasses.
'''

#FIXME: A potential substantial speedup (albeit possibly non-portable to all
#possible Matplotlib backends) is as follows:
#
#    "One limitation of the methods presented above is that all figure elements
#     are redrawn with every call to draw, but we are only updating a single
#     element. Often what we want to do is draw a background, and animate just
#     one or two elements on top of it. As of matplotlib-0.87, GTKAgg, TkAgg,
#     WXAgg, and FLTKAgg support the methods discussed here.
#     The basic idea is to set the 'animated' property of the Artist you want
#     to animate (all figure elements from Figure to Axes to Line2D to Text
#     derive from the base class Artist). Then, when the standard canvas draw
#     operation is called, all the artists except the animated one will be
#     drawn. You can then use the method background =
#     canvas.copy_from_bbox(bbox) to copy a rectangular region (eg the axes
#     bounding box) into a a pixel buffer. In animation, you restore the
#     background with canvas.restore_region(background), and then call
#     ax.draw_artist(something) to draw your animated artist onto the clean
#     background, and canvas.blit(bbox) to blit the updated axes rectangle to
#     the figure. When I run the example below in the same environment that
#     produced 36 FPS for GTKAgg above, I measure 327 FPS with the techniques
#     below. See the caveats on performance numbers mentioned above. Suffice it
#     to say, quantitatively and qualitiatively it is much faster."
#
#For further details, see:
#
#    https://scipy.github.io/old-wiki/pages/Cookbook/Matplotlib/Animations.html#Animating_selected_plot_elements
#
#As the URL fragment "old-wiki" suggests, there probably exists a newer version
#of this technique ideally generalizable to all possible Matplotlib backends.

#FIXME: All animations should be displayed in a non-blocking rather than
#blocking manner, as required for parallelizing the animation pipeline. To
#minimize memory leaks while doing so, consider responding to animation window
#close events by explicitly closing the current animation on such events: e.g.,
#
#    def _hook_on_close(event) -> None:
#        print('Animation "{}" window closed.'.format(self._label))
#
#        #FIXME: Is this safe to call here? Presumably, but consider.
#        self._close()
#
#
#    def __init__(...):
#            .
#            -
#            .
#       # Register the close event handler defined above.
#       self._figure.canvas.mpl_connect('close_event', self._hook_on_close)
#            .
#            -
#            .

#FIXME: We should probably animate non-blockingly (e.g., by passing
#"block=False" to the plt.show() command). To do so, however, we'll probably
#have to implement an analogue to Matplotlib's "_pylab_helper.Gcf" global-like
#static class to ensure that object references are retained for all animations
#whose windows have not yet been closed. This implies, in turn, that we'll
#probably have to monkey-patch the appropriate Matplotlib event handlers on
#window close to release these object references. To do so, grep the Matplotlib
#codebase for "Gcf.destroy". Thunderous tales of woe!
#
#Actually, I believe that a simpler approach might be available. Rather
#than implementing yet another "_pylab_helper.Gcf"-like construct, we leverage
#the fact that animation objects should only live as long as their underlying
#figure objects by explicitly adding circular references between the two: e.g.,
#
#    # This is already done by the "VisualCellsABC" superclass.
#    self._figure = pyplot.figure()
#
#    # Then just add this to the AnimCellsABC.__init__() method *BEFORE* the
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
from betse.exceptions import BetseSimConfException
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.lib.matplotlib.writer import mplvideo
from betse.lib.matplotlib.writer.mplcls import ImageMovieWriter, NoopMovieWriter
from betse.science.phase.phaseenum import SimPhaseKind
from betse.science.visual.layer.vectorfield.lyrvecfldabc import (
    LayerCellsFieldColorlessABC)
from betse.science.visual.layer.vectorfield.lyrvecfldstream import (
    LayerCellsFieldStream)
from betse.science.visual.visabc import VisualCellsABC
from betse.util.io.log import logs
from betse.util.path import dirs, pathnames
from betse.util.type import iterables
from betse.util.type.types import type_check, BoolOrNoneTypes, IntOrNoneTypes
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation

# ....................{ BASE                               }....................
class AnimCellsABC(VisualCellsABC):
    '''
    Abstract base class of all animation subclasses.

    Instances of this class animate the spatial distribution of one or more
    modelled variables (e.g., membrane voltage) over all time steps of the
    current simulation.

    Attributes (Private)
    ----------
    _anim : FuncAnimation
        Low-level Matplotlib animation object instantiated by this high-level
        BETSE wrapper object.

    Attributes (Private: Time)
    ----------
    _time_step_count : int
        Number of frames to be plotted.
    _time_step_last : int
        0-based index of the last frame to be plotted, exactly equivalent to
        `self._time_step_count - 1`.
    _time_step : int
        0-based index of the current frame being plotted, corresponding to the
        0-based sampled time step currently being simulated.
    _time_step_absolute : int
        0-based index of the last frame to be plotted.

    Attributes (Private: Saving)
    ----------
    _save_frame_template : str
        :func:`str.format`-formatted template which, when formatted with the
        0-based index of the current frame, yields the absolute path of the
        image file to be saved for that frame.
    _writer_images : MovieWriter
        Matplotlib object saving animation frames as images if doing so *or*
        ``None`` otherwise.
    _writer_video : MovieWriter
        Matplotlib object saving animation frames as video if doing so *or*
        ``None`` otherwise.

    Attributes (Private: Current)
    ----------
    _is_current_overlayable : BoolOrNoneTypes
        ``True`` only if overlaying either electric current or concentration
        flux streamlines on this animation when requested by the current
        simulation configuration (as governed by the
        ``p.anim.is_overlay_current`` parameter).
    _is_current_overlay_only_gj : bool
        ``True`` only if overlaying intracellular current *or* ``False``
        otherwise (i.e., if overlaying both intra- and extracellular current).
        Ignored unless overlaying current (i.e., if :attr:`_is_current_overlay`
        is ``True``).
    '''

    # ..................{ LIFECYCLE                          }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        save_dir_parent_basename: str,

        #FIXME: Remove the "is_current_overlayable" and
        #"is_current_overlay_only_gj" parameters but *NOT* the corresponding
        #private attributes. Currently, the former parameter is only enabled
        #by subclasses that do *NOT* already plot streamlines. (Makes sense.)
        #Thanks to the new layer API, however, this base class can now
        #automatically detect when the subclass is plotting streamlines by:
        #
        #* Iteratively searching through the "_layers" list in the
        #  animate() function, permitting subclasses to append one or more
        #  layers that plot streamlines.
        #* If any layer in this list is an instance of the new
        #  "LayerCellsFieldStream" subclass, then the "_is_current_overlayable"
        #  attribute should be disabled.
        #* Else, that attribute should default to the corresponding setting(s)
        #  in the current simulation configuration.

        # Optional parameters.
        is_current_overlayable: BoolOrNoneTypes = None,
        is_current_overlay_only_gj: BoolOrNoneTypes = None,
        is_ecm_required: bool = False,
        time_step_count: IntOrNoneTypes = None,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Caveats
        ----------
        If neither displaying nor saving this animation, this animation would
        ideally reduce to a noop. Since this is a superclass method, however,
        simply returning would have little effect. While raising an exception
        would have an effect, doing so would also require all callers to
        explicitly catch and ignore that exception -- in which case this
        animation would hardly have reduced to a noop. In short, ignoring this
        edge case is currently the only sane policy.

        Parameters
        ----------
        save_dir_parent_basename : str
            Basename of the parent directory of the subdirectory to which this
            animation's frames will be saved when requested by the current
            simulation configuration.
        is_current_overlayable : optional[bool]
            ``True`` only if overlaying either electric current or concentration
            flux streamlines on this animation when requested by the current
            simulation configuration (as governed by the
            ``p.anim.is_overlay_current`` parameter). All subclasses except
            those already plotting streamlines (e.g., by calling the superclass
            :meth:`_plot_stream` method) should unconditionally enable this
            boolean.
        is_current_overlay_only_gj : optional[bool]
            ``True`` if only overlaying intracellular current *or* ``False`` if
            overlaying both intra- and extracellular current. Ignored if current
            is _not_ being overlayed at all (i.e., if
            :attr:`_is_current_overlay` is ``False``). If ``None``, defaults to
            the following state:
            * ``False`` if extracellular spaces are enabled _and_ both
               intracellular and extracellular current is being animated.
            * ``True`` if either extracellular spaces are disabled _or_ are
               enabled but only intracellular current is being animated.
        is_ecm_required : optional[bool]
            ``True`` only if this animation is specific to extracellular spaces.
            If ``True`` and extracellular spaces are currently disabled, an
            exception is raised. Defaults to ``False``.
        time_step_count : optional[int]
            Number of frames to be plotted. If `None`, defaults to the number
            of sampled time steps in the current tissue simulation.
        '''

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # If this subclass requires extracellular spaces but extracellular
        # spaces are currently disabled, raise an exception.
        if is_ecm_required and not self._phase.p.is_ecm:
            raise BetseSimConfException(
                'Animation "{}" requires extracellular spaces, which are '
                'disabled by the current simulation configuration.'.format(
                self._label))

        # Default unpassed parameters.
        #
        # If unpassed, overlay streamlines onto this animation only if all of
        # these conditions are satisfied:
        #
        # * The current simulation phase is any phase after the "seed" phase
        #   (i.e., any phase except the "seed" phase). Streamlines require
        #   simulation data defined during the "init" phase and are thus
        #   unavailable in earlier phases.
        # * The current simulation configuration enables this overlay.
        if is_current_overlayable is None:
            is_current_overlayable = (
                self._phase.kind >= SimPhaseKind.INIT and
                self._phase.p.anim.is_overlay_current)
        if is_current_overlay_only_gj is None:
            is_current_overlay_only_gj = not (
                self._phase.p.is_ecm and self._phase.p.IecmPlot)
        if time_step_count is None:
            time_step_count = len(self._phase.sim.time)

        # Classify all remaining parameters.
        self._is_current_overlayable = is_current_overlayable
        self._is_current_overlay_only_gj = is_current_overlay_only_gj
        self._time_step_count = time_step_count

        # 0-based index of the last frame to be plotted.
        self._time_step_last = self._time_step_count - 1

        # Classify attributes to be possibly redefined below.
        self._writer_images = None
        self._writer_video = None

        # 0-based index of the current frame.
        self._time_step = 0

        # If saving animations, prepare to do so.
        self._init_saving(save_dir_parent_basename)


    @type_check
    def _init_saving(self, save_dir_parent_basename: str) -> None:
        '''
        Initialize this animation for platform-compatible file saving if
        enabled by the current simulation configuration _or_ noop otherwise.

        Parameters
        ----------
        save_dir_parent_basename : str
            Basename of the parent directory of the subdirectory to which this
            animation's frames will be saved when requested by the current
            simulation configuration. Defaults to a suitably generic basename.
        '''

        # Animation configuration localized for convenience.
        anim_config = self._phase.p.anim

        # If it is *NOT* the case that...
        if not (
            # This animation is being saved...
            self._is_save and (
            # ...as either images or video.
            anim_config.is_images_save or anim_config.is_video_save)
        # Then this animation is unsaved. In this case, silently noop.
        ):
            return

        #FIXME: This is silly. Rather than prohibiting animation names
        #containing directory separators, simply sanitize this animation's name
        #by globally replacing all such separators with non-separator characters
        #guaranteed to be permitted in pathnames for all platforms (e.g.,
        #replacing "/" with "_" or "-" characters on POSIX systems).

        # If the human-readable name of this animation contains directory
        # separators and hence is *NOT* a valid basename, raise an exception.
        pathnames.die_unless_basename(self._label)

        # Path of the subdirectory to which these files will be saved,
        # creating this subdirectory and all parents thereof if needed.
        save_dirname = dirs.canonicalize_and_make_unless_dir(pathnames.join(
            self._phase.export_dirname, save_dir_parent_basename, self._label))

        # If saving animation frames as images, prepare to do so.
        if anim_config.is_images_save:
            #FIXME: This currently defaults to padding frames with six or seven
            #zeroes, on average. Let's make this a bit more aesthetic by
            #padding frames to only as many zeroes are absolutely required by
            #the current frame count. To do that, in turn, we'll probably need
            #to shift everything that follows in this method to the _animate()
            #method, where the actual frame count is finally passed.

            # Template expanding to the basename of each image to be saved.
            # The "ImageMovieWriter" class subsequently expands the "{{"- and
            # "}}"-delimited substring to the 0-based index of the current
            # frame number.
            save_frame_template_basename = '{}_{{:07d}}.{}'.format(
                self._label, anim_config.image_filetype)

            # Template expanding to the absolute path of each image to be
            # saved.
            writer_images_template = pathnames.join(
                save_dirname, save_frame_template_basename)

            # Object writing animation frames as images.
            self._writer_images = ImageMovieWriter()

            # Log this preparation.
            logs.log_debug(
                'Preparing to save animation frames "%s"...',
                writer_images_template)

            # Prepare to save these animation frames.
            self._writer_images.setup(
                fig=self._figure,
                outfile=writer_images_template,
                dpi=anim_config.image_dpi,
            )

        # If saving animation frames as video, prepare to do so.
        if anim_config.is_video_save:
            # Name of the first video encoder installed on the current system.
            video_writer_name = mplvideo.get_first_writer_name(
                anim_config.video_writer_names)
            # print('found video writer: {}'.format(VideoWriterClass))

            # Matplotlib animation writer class encapsulating this encoder.
            VideoWriterClass = mplvideo.get_writer_class(video_writer_name)

            # Name of the first video codec supported by both this video
            # encoder and the video container format corresponding to this
            # video's filetype.
            video_codec_name = mplvideo.get_first_codec_name(
                writer_name=video_writer_name,
                container_filetype=anim_config.video_filetype,
                codec_names=anim_config.video_codec_names,
            )

            # Basename of the video to be written.
            save_video_basename = '{}.{}'.format(
                self._label, anim_config.video_filetype)

            # Absolute path of the video to be written.
            writer_video_filename = pathnames.join(
                save_dirname, save_video_basename)

            # Object writing animation frames as video. Matplotlib animation
            # writer classes must *NOT* be manually instantiated, as doing so
            # fails to apply the monkey-patching applied by this function.
            self._writer_video = mplvideo.make_writer(
                cls=VideoWriterClass,
                bitrate=anim_config.video_bitrate,
                codec=video_codec_name,
                fps=anim_config.video_framerate,
                metadata=anim_config.video_metadata,
            )

            # Log this preparation.
            logs.log_debug(
                'Preparing to save animation video "%s"...',
                writer_video_filename)

            # Prepare to save this animation video. Matplotlib squelches
            # critical (technically non-fatal but effectively fatal)
            # warnings and errors emitted by the external command invoked
            # by this call to the MovieWriter.setup() method *UNLESS* the
            # current matplotlib-specific verbosity level is "debug".
            # Temporarily ensure this for the duration of this call.
            with mpl_config.reducing_log_level_to_debug_if_info():
                self._writer_video.setup(
                    fig=self._figure,
                    outfile=writer_video_filename,
                    dpi=anim_config.video_dpi,
                )

    # ..................{ PROPERTIES                         }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def time_step(self) -> int:
        '''
        0-based index of the current frame being plotted, corresponding to the
        0-based sampled time step currently being simulated.
        '''

        return self._time_step

    # ..................{ PREPARERS                          }..................
    # This method has been overridden to support subclasses that manually
    # handle animations rather than calling the _animate() method (e.g., the
    # "AnimCellsWhileSolving" subclass).
    def _prep_figure(self, *args, **kwargs) -> None:

        #FIXME: Refactor this method as follows:
        #
        #* Eliminate all current overlay-specific boolean attributes (e.g.,
        #  "_is_current_overlayable"), as detailed in an __init__() method
        #  comment above.
        #* Shift the _prep_layer_current() method into the superclass.
        #* Call the _prep_layer_current() method in the superclass
        #  _prep_figure() implementation.
        #* Remove this _prep_figure() implementation.

        # Append a layer overlaying current density *AFTER* all lower layers
        # (e.g., cell data) have been appended but *BEFORE* all higher layers
        # (e.g., cell labelling) have been appended.
        self._prep_layer_current()

        # Perform superclass figure preparation.
        super()._prep_figure(*args, **kwargs)


    def _prep_layer_current(self) -> None:
        '''
        Append a layer overlaying current density if this simulation
        configuration requests a current overlay and no layer in the layer
        sequence already plots streamlines _or_ noop otherwise.

        If an existing layer in the layer sequence already plots streamlines,
        this method avoids adding another layer also doing so. Since the
        current overlay also plots streamlines, attempting to plot streamlines
        over existing streamlines would produce an unintelligible plot or
        animation. (That would be bad.)
        '''

        # If either...
        if (
            # This simulation configuration requests no current overlay.
            not self._is_current_overlayable or

            # A layer in the layer sequence already plots a vector field.
            iterables.is_item_instance_of(
                iterable=self._layers, cls=LayerCellsFieldColorlessABC)
        # ...then silently noop.
        ):
            return
        # Else, overlay current.

        # Current density vector field to layer.
        field = None

        # If layering intracellular current, do so.
        if self._is_current_overlay_only_gj:
            logs.log_debug('Overlayering intracellular current...')
            field = self._phase.cache.vector_field.currents_intra
        # Else, layer extracellular current.
        else:
            logs.log_debug('Overlayering extracellular current...')
            field = self._phase.cache.vector_field.currents_extra

        # Append a layer overlaying this field.
        self._append_layer(LayerCellsFieldStream(field=field))

    # ..................{ ANIMATORS                          }..................
    @type_check
    def _animate(self, *args, **kwargs) -> None:
        '''
        Display and/or save this animation as requested by the current
        simulation configuration.

        This method is intended to be called as the last statement in the
        :meth:`__init__` method of all subclasses of this superclass.

        Parameters
        ----------
        All parameters are passed to the :meth:`_prep_figure` method.
        '''

        # Prepare for plotting immediately *BEFORE* plotting the first frame.
        self._prep_figure(*args, **kwargs)

        #FIXME: For efficiency, we should probably be passing "blit=True," to
        #FuncAnimation(). Unfortunately, doing so will necessitate
        #restructuring animations to conform to blitting-specific requirements,
        #including:
        #
        #* The definition of a new init_frame() method of this class, which
        #  plots all portions of this animation *NOT* changing between frames.
        #  Are there any? We have no idea. Axes ticks presumably never change
        #  between frames for any animation, so there should be *SOMETHING* to
        #  plot here.
        #* The passing of this method to the "FuncAnimation" instance defined
        #  below via the "init_func" parameter: e.g.,
        #      init_func=self.init_frame,
        #* Refactoring the plot_frame() method to return the tuple of all
        #  matplotlib artist objects modified for the current frame. This is
        #  probably infeasible in a generic manner given the current crusty
        #  design of this class, but should be trivial (...in theory) after
        #  redesigning this class to support composoble
        #
        #    http://devosoft.org/making-efficient-animations-in-matplotlib-with-blitting
        #
        #Lemon grass and dill!

        # Create and assign an animation function to a local variable. If the
        # latter is *NOT* done, this function will be garbage collected prior
        # to subsequent plot handling -- in which case only the first plot will
        # be plotted without explicit warning or error. Die, matplotlib! Die!!!
        self._anim = FuncAnimation(
            # Figure to which the "func" callable plots each frame.
            fig=self._figure,

            # Callable plotting each frame.
            func=self.plot_frame,

            # Number of frames to be animated.
            frames=self._time_step_count,

            #FIXME: The interval should, ideally, be synchronized with the FPS
            #used for video encoding. To guarantee this:
            #
            #* Generalize the FPS option in the configuration file to *ALL*
            #  animations. Currently, this option only applies to video
            #  encoding.
            #* Convert the currently configured FPS into this interval in
            #  milliseconds as follows:
            #
            #      interval = 1000.0 / fps

            # Delay in milliseconds between consecutive frames. To convert this
            # delay into the equivalent frames per second (FPS):
            #
            #      fps = 1000.0 / interval
            interval=200,

            #FIXME: This is a bit silly. Ideally, animations should *ALWAYS* be
            #repeatable. Once we've refactored away usage of the
            #Animation.save() method, refactor:
            #
            #* This parameter to unconditionally enable repeating: e.g.,
            #      repeat=True,
            #* The plot_frame() method to conditionally call MovieWriter
            #  methods (e.g., grab_frame(), finish()) *ONLY* if the current
            #  call to the plot_frame() method is the first such call for the
            #  current frame. While this state would be trivial for this class
            #  to record, perhaps matplotlib's "Animation" base class already
            #  records this state? Contemplate us up.

            # Indefinitely repeat this animation unless saving animations, as
            # doing so under the current implementation would repeatedly (and
            # hence unnecessarily) overwrite previously written files.
            repeat=not self._is_save,
        )

        try:
            # If displaying and optionally saving this animations, do so.
            if self._is_show:
                #FIXME: If the current backend is non-interactive (e.g.,
                #"Agg"), the following function call reduces to a noop. This is
                #insane, frankly. In this case, this animation's plot_frame()
                #is never called! No errors or warnings are logged, so it's
                #unclear who or what is the culprit here. If the pyplot.show()
                #function is indeed only supported by interactive backends, we
                #should do the following here:
                #
                #* Detect whether or not the current backend is
                #  non-interactive.
                #* If so, either:
                #  * Emit an explicit warning advising the user that this
                #    animation will almost certainly be silently ignored. This
                #    isn't terribly ideal, but it's better than zilch.
                #  * If this animation is currently being saved, simply perform
                #    the non-display save logic performed below -- which *DOES*
                #    behave as expected for non-interactive backends. Clearly,
                #    the culprit is the pyplot.show() function. Mournful sigh.

                # Display and optionally save this animations. Note that,
                # although this function is called in a blocking manner, the
                # GUI-driven event loops of some interactive backends appear to
                # ignore the request for blocking behavior and perform
                # non-blocking behaviour instead. This, in turn, prevents this
                # branch from reliably finalizing this animation by calling the
                # close() method. This differs from the non-interactive
                # saving-specific branch that follows, which is guaranteed to
                # behave in a blocking manner and hence *CAN* reliably call the
                # close() method. tl;dr: GUIs, so random.
                pyplot.show()
            # Else if only saving but not displaying this animation *AND* at
            # least one animation writer doing so is enabled, do so.
            elif self._is_save and (
                self._writer_images is not None or
                self._writer_video is not None
            ):
                # Save this animation by iteratively calling our plot_frame()
                # method to save each animation frame. Since this method
                # already manually saves each such frame for the case of both
                # displaying *AND* saving this animation via the above call to
                # the pyplot.show() function, that logic is reused here by
                # preventing this call to the Animation.save() method from
                # attempting to automatically save each such frame.
                #
                # By default, the Animation.save() method iteratively calls the
                # MovieWriter.grab_frame() method of the passed "writer" object
                # to save each such frame. If no such object is passed, this
                # object defaults to a new writer whose name is the current
                # value of the "animation.writer" rcparam. Hence, there exists
                # no means of preventing the Animation.save() method from
                # writing. However, there also exists no good alternative to
                # this method that iteratively calls our plot_frame() method
                # without also writing. For example:
                #
                # * The pyplot.show() function iterating frames silently
                #   reduces to a noop for non-interactive backends and is thus
                #   inapplicable as a general-purpose solution.
                # * The frame iteration automatically performed by the
                #   Animation.save() method is both non-trivial and requires
                #   calls to private methods of the matplotlib animation API.
                #   While this iteration could (and actually was, in the first
                #   implementation of this approach) be reduplicated here,
                #   doing so would be overly fragile and hence break under
                #   upstream changes to this private API.
                #
                # The robust solution is to instead pass the Animation.save()
                # method a writer reducing to a noop, circumventing conflicts
                # with the manual saving performed by our plot_frame() method.
                self._anim.save(
                    # Note that, since "NoopMovieWriter" maintains no state, a
                    # singleton "NoopMovieWriter" instance could technically be
                    # shared amongst all animation classes. However, since
                    # "NoopMovieWriter" construction is trivially fast, there are no
                    # demonstrable advantages and arguable disadvantages to
                    # doing so (e.g., code complexity, space consumption).
                    writer=NoopMovieWriter(),

                    # Pass an ignorable filename. To guarantee that an
                    # exception is raised on this method attempting to read or
                    # write this file, pass a filename guaranteed to be invalid
                    # on all supported platforms (e.g., containing null bytes).
                    # For understandable reasons, this parameter is mandatory.
                    filename=pathnames.INVALID_PATHNAME,
                )

                # Finalize saving this animation.
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

    # ..................{ CLOSERS                            }..................
    def close(self) -> None:
        '''
        Finalize displaying and/or saving this non-blocking animation.

        This method is intended to be called as the last statement in the
        `_animate()` method animating this non-blocking subclass. If this
        subclass is blocking, this method must _not_ be called.

        See Also
        ----------
        :meth:`VisualCellsABC.close`
            Superclass method destroying this animation's low-level plot and
            all memory associated with this plot.
        '''

        # Finalize this animation's low-level plot.
        super().close()

        # Finalize all writers saving this animation if any.
        self._close_writers()

        # Prevent this animation from being reused *AND* break hard cycles.
        self._anim = None


    def _close_writers(self) -> None:
        '''
        Finalize all writers saving this animation if any.
        '''

        # If saving animation frames as images...
        if self._writer_images is not None:
            # Finalize doing so.
            self._writer_images.finish()

            # Prevent this writer from being reused and break hard cycles.
            self._writer_images = None

        # If saving animation frames as video...
        if self._writer_video is not None:
            # Finalize doing so. For debuggability, temporarily escalate the
            # matplotlib-specific verbosity level.
            with mpl_config.reducing_log_level_to_debug_if_info():
                self._writer_video.finish()

            # Prevent this writer from being reused and break hard cycles.
            self._writer_video = None
