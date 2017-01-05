#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Animation configuration and serialization classes.
'''

#FIXME: Default the "copyright" entry of video metadata to
#"@ {}".format(current_year)".

#FIXME: Define saving-ordiented methods.

# ....................{ IMPORTS                            }....................
from betse.util.type import ints, types
from betse.util.type.types import type_check, MappingType, SequenceTypes

# ....................{ SUPERCLASS                         }....................
class AnimConfig(object):
    '''
    Object encapsulating both the configuration and writing of all animations
    (both in- and post-simulation), parsed from the current configuration file.

    This object saves (i.e., writes, serializes) in-memory animations to on-disk
    cache, image, and/or video files configured by this configuration.

    Attributes
    ----------
    is_while_sim : bool
        `True` only if this configuration enables (but _not_ necessarily
        displays or saves) in-simulation animations.
    is_while_sim_show : bool
        `True` only if this configuration displays in-simulation animations.
        Ignored if `is_while_sim` is `False`.
    is_while_sim_save : bool
        `True` only if this configuration saves in-simulation animations.
        Ignored if `is_midsim` is `False`.
    is_after_sim : bool
        `True` only if this configuration enables (but _not_ necessarily
        displays or saves) post-simulation animations.
    is_after_sim_show : bool
        `True` only if this configuration displays post-simulation animations.
        Ignored if `is_after_sim` is `False`.
    is_after_sim_save : bool
        `True` only if this configuration saves post-simulation animations.
        Ignored if `is_after_sim` is `False`.
    is_images_save : bool
        `True` only if this configuration saves animation frames as images.
    is_video_save : bool
        `True` only if this configuration saves animation frames as video.
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        `is_images_save` is `False`.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if `is_images_save` is `False`.
    video_bitrate : int
        Bitrate in bits per second of all video files saved by this
        configuration. Ignored if `is_video_save` is `False`.
    video_codec_names : Sequence
        List of the names of all encoder-specific codecs with which to encode
        animations (in descending order of preference), automatically:
        * Selecting the first codec supported by the selected writer.
        * Replacing any codec named `auto` by the name of a recommended codec
          specific to the selected writer and filetype. For details, see the
          `betse.lib.matplotlib.writer.mplvideo.get_first_codec_name()`
          function.
        Ignored if `is_video_save` is `False`.
    video_dpi : int
        Dots per inch (DPI) of all frames of all video files saved by this
        configuration. Ignored if `is_images_save` is `False`.
    video_filetype : str
        Filetype of all video files saved by this configuration. Ignored if
        `is_video_save` is `False`. Supported filetypes include:
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
        configuration. Ignored if `is_video_save` is `False`.
    video_metadata : dict
        Dictionary mapping from the alphabetic lowercase name of video metadata
        supported by the active video encoder to that metadata's human-readable
        string to be embedded in all video files saved by this configuration.
        Ignored if `is_video_save` is `False`. Supported names include:
        `title`, `artist`, `genre`, `subject`, `copyright`, `srcform`, and
        `comment`. If this dictionary does _not_ contain a `copyright` key, such
        a key will be automatically synthesized from the current year.
    video_writer_names : Sequence
        List of the names of all matplotlib animation writers with which to
        encode animations (in order of descending preference), automatically
        selecting the first writer installed on the current system. Ignored if
        `is_video_save` is `False`. Supported names include:
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
        save = results['save']['animations']
        images = save['images']
        video = save['video']

        # Create and return this instance.
        return AnimConfig(
            # In-simulation animations.
            is_while_sim=while_solving['enabled'],
            is_while_sim_show=while_solving['show'],
            is_while_sim_save=while_solving['save'],

            # Post-simulation animations.
            is_after_sim=after_solving['enabled'],
            is_after_sim_show=after_solving['show'],
            is_after_sim_save=after_solving['save'],

            # Image saving.
            is_images_save=images['enabled'],
            image_filetype=images['filetype'],
            image_dpi=images['dpi'],

            # Video saving.
            is_video_save=video['enabled'],
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

        # In-simulation animations.
        is_while_sim: bool,
        is_while_sim_show: bool,
        is_while_sim_save: bool,

        # Post-simulation animations.
        is_after_sim: bool,
        is_after_sim_show: bool,
        is_after_sim_save: bool,

        # Image saving.
        is_images_save: bool,
        image_filetype: str,
        image_dpi: int,

        # Video saving.
        is_video_save: bool,
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
        self.is_while_sim = is_while_sim
        self.is_while_sim_show = is_while_sim_show
        self.is_while_sim_save = is_while_sim_save
        self.is_after_sim = is_after_sim
        self.is_after_sim_show = is_after_sim_show
        self.is_after_sim_save = is_after_sim_save
        self.is_images_save = is_images_save
        self.is_video_save = is_video_save
        self.image_filetype = image_filetype
        self.image_dpi = image_dpi
        self.video_bitrate = video_bitrate
        self.video_dpi = video_dpi
        self.video_filetype = video_filetype
        self.video_framerate = video_framerate
        self.video_metadata = video_metadata
        self.video_writer_names = video_writer_names
        self.video_codec_names = video_codec_names
