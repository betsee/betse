#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Animation serialization classes.
'''

#FIXME: Define saving-ordiented methods.

#FIXME: Redefine all animations-oriented "Parameters" booleans excluding
#"Parameters.turn_all_plots_off" (which applies to plots as well and hence is
#more general than merely animations) in terms of high-level instance variables
#of the "AnimConfig" class rather than in terms of low-level subdictionaries of
#the "Parameters.config" dictionary. For example, redefine:
#
#* "Parameters.createAnimations" in terms of "AnimConfig.is_postsim".

# ....................{ IMPORTS                            }....................
from betse.util.type import types
from betse.util.type.types import type_check, SequenceTypes

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
    is_midsim : bool
        `True` only if this configuration enables (but _not_ necessarily
        displays or saves) mid-simulation animations.
    is_midsim_showing : bool
        `True` only if this configuration displays mid-simulation animations.
        Ignored if `is_midsim` is `False`.
    is_midsim_saving : bool
        `True` only if this configuration saves mid-simulation animations.
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
    is_saving_images : bool
        `True` only if this configuration saves animation frames as images.
    is_saving_video : bool
        `True` only if this configuration saves animation frames as video.
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        `is_saving_images` is `False`.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if `is_saving_images` is `False`.
    video_bitrate : int
        Bitrate in bits per second of all video files saved by this
        configuration. Ignored if `is_saving_video` is `False`.
    video_dpi : int
        Dots per inch (DPI) of all frames of all video files saved by this
        configuration.  Ignored if `is_saving_images` is `False`.
    video_encoder_names : Sequence
        List of the matplotlib-specific names of all video encoders with which
        to encode animations (in order of descending preference), ignoring all
        video encoders _not_ installed on the current system. Ignored if
        `is_saving_video` is `False`.
    video_filetype : str
        Filetype of all video files saved by this configuration. Ignored if
        `is_saving_video` is `False`. Supported filetypes include:
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
        configuration. Ignored if `is_saving_video` is `False`.
    video_metadata : dict
        Dictionary mapping from the alphabetic lowercase name of video metadata
        supported by the active video encoder to that metadata's human-readable
        string to be embedded in all video files saved by this configuration.
        Ignored if `is_saving_video` is `False`. Supported names include:
        `title`, `artist`, `genre`, `subject`, `copyright`, `srcform`, and
        `comment`. If this dictionary does _not_ contain a `copyright` key, such
        a key will be automatically synthesized from the current year.
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

        anim = params.config['results options']
        save = anim['save animations']
        image = save['image']
        video = save['video']

        # Create and return this instance.
        return AnimConfig(
            # Mid-simulation animations.
            is_midsim=anim['plot while solving'],
            is_midsim_showing=not params.turn_all_plots_of,
            is_midsim_saving=anim['save solving plot'],

            # Post-simulation animations.
            is_postsim=anim['create all animations'],
            is_postsim_showing=not params.turn_all_plots_of,
            is_postsim_saving=anim['automatically save plots'],

            # Image saving.
            is_saving_images=image['enabled'],
            image_filetype=image['filetype'],
            image_dpi=image['dpi'],

            # Video saving.
            is_saving_video=video['enabled'],
            video_bitrate=video['bitrate'],
            video_dpi=video['dpi'],
            video_filetype=video['filetype'],
            video_framerate=video['framerate'],
            video_metadata=video['metadata'],
            video_encoder_names=video['encoders'],
        )

    # ..................{ CONCRETE ~ public                  }..................
    #FIXME: Implement us up.
    @type_check
    def __init__(
        self,
    ) -> None:

        # Classify the passed parameters.
        pass
