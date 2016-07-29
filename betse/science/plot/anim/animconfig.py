#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Animation serialization classes.
'''

#FIXME: Default the "copyright" entry of video metadata to
#"@ {}".format(current_year)".

#FIXME: Define saving-ordiented methods.

#FIXME: Redefine all animations-oriented "Parameters" booleans excluding
#"Parameters.turn_all_plots_off" (which applies to plots as well and hence is
#more general than merely animations) in terms of high-level instance variables
#of the "AnimConfig" class rather than in terms of low-level subdictionaries of
#the "Parameters.config" dictionary. For example, redefine:
#
#* "Parameters.createAnimations" in terms of "AnimConfig.is_postsim".

# ....................{ IMPORTS                            }....................
from betse.util.type import ints, types
from betse.util.type.types import type_check, MappingType, SequenceTypes

# ....................{ SUPERCLASS                         }....................
class AnimConfig(object):
    '''
    Object encapsulating both the configuration and writing of all animations
    (both mid- and post-simulation), as configured by the current tissue
    simulation's configuration file.

    This object saves (i.e., writes, serializes) in-memory animations to on-disk
    cache, image, and/or video files configured by this configuration.

    Attributes
    ----------
    is_withsim : bool
        `True` only if this configuration enables (but _not_ necessarily
        displays or saves) in-simulation animations.
    is_withsim_showing : bool
        `True` only if this configuration displays in-simulation animations.
        Ignored if `is_withsim` is `False`.
    is_withsim_saving : bool
        `True` only if this configuration saves in-simulation animations.
        Ignored if `is_midsim` is `False`.
    is_postsim : bool
        `True` only if this configuration enables (but _not_ necessarily
        displays or saves) post-simulation animations.
    is_postsim_showing : bool
        `True` only if this configuration displays post-simulation animations.
        Ignored if `is_postsim` is `False`.
    is_postsim_saving : bool
        `True` only if this configuration saves post-simulation animations.
        Ignored if `is_postsim` is `False`.
    is_images_saving : bool
        `True` only if this configuration saves animation frames as images.
    is_video_saving : bool
        `True` only if this configuration saves animation frames as video.
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        `is_images_saving` is `False`.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if `is_images_saving` is `False`.
    video_bitrate : int
        Bitrate in bits per second of all video files saved by this
        configuration. Ignored if `is_video_saving` is `False`.
    video_codec_names : Sequence
        List of the names of all encoder-specific codecs with which to encode
        animations (in descending order of preference), automatically:
        * Selecting the first codec supported by the selected writer.
        * Replacing any codec named `auto` by the name of a recommended codec
          specific to the selected writer and filetype. For details, see the
          `betse.lib.matplotlib.writer.mplvideo.get_first_codec_name()`
          function.
        Ignored if `is_video_saving` is `False`.
    video_dpi : int
        Dots per inch (DPI) of all frames of all video files saved by this
        configuration. Ignored if `is_images_saving` is `False`.
    video_filetype : str
        Filetype of all video files saved by this configuration. Ignored if
        `is_video_saving` is `False`. Supported filetypes include:
        * `mkv` (Matroska), an open-standard audio and video container
          supporting all relevant codecs and hence the default.
        * `avi`, Microsoft's obsolete proprietary audio and video container.
        * `gif`, a proprietary image format supporting video animations.
        * `ogv` (Theora), Xiph's open-standard audio and video container.
        * `mov` (QuickTime), Apple's proprietary audio and video container.
        * `mp4` (MPEG-4 Part 14), an open-standard audio and video container.
        * `webm` (WebM), Google's proprietary audio and video container.
    video_framerate : int
        Framerate in frames per second of all video files saved by this
        configuration. Ignored if `is_video_saving` is `False`.
    video_metadata : dict
        Dictionary mapping from the alphabetic lowercase name of video metadata
        supported by the active video encoder to that metadata's human-readable
        string to be embedded in all video files saved by this configuration.
        Ignored if `is_video_saving` is `False`. Supported names include:
        `title`, `artist`, `genre`, `subject`, `copyright`, `srcform`, and
        `comment`. If this dictionary does _not_ contain a `copyright` key, such
        a key will be automatically synthesized from the current year.
    video_writer_names : Sequence
        List of the names of all matplotlib animation writers with which to
        encode animations (in order of descending preference), automatically
        selecting the first writer installed on the current system. Ignored if
        `is_video_saving` is `False`. Supported names include:
        * `ffmpeg`, an open-source cross-platform audio and video encoder.
        * `avconv`, an open-source cross-platform audio and video encoder
          forked from (and largely interchangeable with) `ffmpeg`.
        * `mencoder`, an open-source cross-platform audio and video encoder
          associated with MPlayer, a popular media player.
        * `imagemagick`, an open-source cross-platform image manipulation
          suite supporting _only_ creation of animated GIFs.
    '''

    # ..................{ ABSTRACT ~ static                  }..................
    @staticmethod
    def make(params: 'Parameters') -> "AnimConfig":
        '''
        Factory method producing an instance of this class encapsulating the
        passed simulation configuration.

        Parameters
        ----------------------------
        params : Parameters
            Current simulation configuration.

        Returns
        ----------------------------
        AnimConfig
            Instance of this class encapsulating this configuration.
        '''
        assert types.is_parameters(params), types.assert_not_parameters(params)

        # For convenience, localize configuration subdictionaries.
        results = params.config['results options']
        while_solving = results['while solving']['animations']
        after_solving = results['after solving']['animations']
        export = results['export']['animations']
        images = export['images']
        video = export['video']

        # Create and return this instance.
        return AnimConfig(
            # In-simulation animations.
            is_withsim=while_solving['enabled'],
            is_withsim_showing=while_solving['show'],
            is_withsim_saving=while_solving['save'],

            # Post-simulation animations.
            is_postsim=after_solving['enabled'],
            is_postsim_showing=after_solving['show'],
            is_postsim_saving=after_solving['save'],

            # Image saving.
            is_images_saving=images['enabled'],
            image_filetype=images['filetype'],
            image_dpi=images['dpi'],

            # Video saving.
            is_video_saving=video['enabled'],
            video_bitrate=video['bitrate'],
            video_dpi=video['dpi'],
            video_filetype=video['filetype'],
            video_framerate=video['framerate'],
            video_metadata=video['metadata'],
            video_writer_names=video['writers'],
            video_codec_names=video['codecs'],
        )

    # ..................{ CONCRETE ~ public                  }..................
    @type_check
    def __init__(
        self,
        is_withsim: bool,
        is_withsim_showing: bool,
        is_withsim_saving: bool,

        # Post-simulation animations.
        is_postsim: bool,
        is_postsim_showing: bool,
        is_postsim_saving: bool,

        # Image saving.
        is_images_saving: bool,
        image_filetype: str,
        image_dpi: int,

        # Video saving.
        is_video_saving: bool,
        video_bitrate: int,
        video_dpi: int,
        video_filetype: str,
        video_framerate: int,
        video_metadata: MappingType,
        video_writer_names: SequenceTypes,
        video_codec_names: SequenceTypes,
    ) -> None:

        # Validate all passed integers as positive.
        ints.die_unless_positive(
            image_dpi, video_bitrate, video_dpi, video_framerate)

        #FIXME: Repetition is vile and demeaning. Design and leverage a new
        #@classify_params decorator here instead, please.

        # Classify the passed parameters.
        self.is_withsim = is_withsim
        self.is_withsim_showing = is_withsim_showing
        self.is_withsim_saving = is_withsim_saving
        self.is_postsim = is_postsim
        self.is_postsim_showing = is_postsim_showing
        self.is_postsim_saving = is_postsim_saving
        self.is_images_saving = is_images_saving
        self.image_filetype = image_filetype
        self.image_dpi = image_dpi
        self.is_video_saving = is_video_saving
        self.video_bitrate = video_bitrate
        self.video_dpi = video_dpi
        self.video_filetype = video_filetype
        self.video_framerate = video_framerate
        self.video_metadata = video_metadata
        self.video_writer_names = video_writer_names
        self.video_codec_names = video_codec_names
