#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation animation subconfigurations.
'''

#FIXME: Default the "copyright" entry of video metadata to
#"@ {}".format(current_year)".

#FIXME: Define saving-ordiented methods.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMethodUnimplementedException
from betse.science.config.confabc import (
    SimConfABC, SimConfListableABC, SimConfList, conf_alias, conf_enum_alias)
from betse.util.type import ints
from betse.util.type.types import MappingType, NumericTypes, SequenceTypes
from enum import Enum

# ....................{ ENUMERATIONS                       }....................
SimConfAnimKind = Enum('SimConfAnimKind', (
    'CURRENTS_INTRA',
    'CURRENTS_TOTAL',
    'ELECTRIC_INTRA',
    'ELECTRIC_TOTAL',
    'VOLTAGES_INTRA',
    'VOLTAGES_TOTAL',
))
'''
Enumeration of all possible types of animations.

See the corresponding entry of the default YAML-based simulation configuration
file for further commentary.
'''

# ....................{ SUBCLASSES                         }....................
#FIXME: Rename to merely "SimConfAnimOne" *AFTER* eliminating the following
#filehandling._preserve_backward_importability() assignment:
#
#    sys.modules['betse.science.config.visual.confanim'].SimConfAnim = (
#        confanim.SimConfAnimAll)
class SimConfAnimOne(SimConfListableABC):
    '''
    YAML-backed simulation animation subconfiguration, encapsulating the
    configuration of a single animation (either in- or post-simulation) parsed
    from the list of all such animations in the current YAML-formatted
    simulation configuration file.

    Attributes (General)
    ----------
    kind : SimConfAnimKind
        Type of this animation.

    Attributes (Colorbar)
    ----------
    color_max : NumericTypes
        Maximum color value to be displayed by the colorbar. Ignored if
        :attr:`is_color_autoscaled` is ``True``.
    color_min : NumericTypes
        Minimum color value to be displayed by the colorbar. Ignored if
        :attr:`is_color_autoscaled` is ``True``.
    is_color_autoscaled : bool
        ``True`` if dynamically setting the minimum and maximum colorbar values
        for this animation to the minimum and maximum values flattened from the
        corresponding time series *or* ``False`` if statically setting these
        values to :attr:`color_min` and :attr:`color_max`.
    '''

    # ..................{ ALIASES                            }..................
    kind = conf_enum_alias("['type']", SimConfAnimKind)

    # ..................{ ALIASES ~ colorbar                 }..................
    is_color_autoscaled = conf_alias("['colorbar']['autoscale']", bool)
    color_min = conf_alias("['colorbar']['minimum']", NumericTypes)
    color_max = conf_alias("['colorbar']['maximum']", NumericTypes)

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Actually implement this method.
    def default(self) -> None:
        raise BetseMethodUnimplementedException()

# ....................{ SUBCLASSES ~ all                   }....................
class SimConfAnimAll(SimConfABC):
    '''
    YAML-backed simulation animations subconfiguration, encapsulating the
    configuration of all animations (both in- and post-simulation) parsed from
    the current YAML-formatted simulation configuration file.

    This subconfiguration saves (i.e., writes, serializes) in-memory animations
    to on-disk cache, image, and/or video files configured by this
    configuration.

    Attributes (While)
    ----------
    while_sim_pipeline : SimConfList
        List of all post-simulation animations to be animated.
    is_while_sim : bool
        ``True`` only if this configuration enables (but _not_ necessarily
        displays or saves) mid-simulation animations.
    is_while_sim_show : bool
        ``True`` only if this configuration displays mid-simulation animations.
        Ignored if `is_while_sim` is `False`.
    is_while_sim_save : bool
        ``True`` only if this configuration saves mid-simulation animations.
        Ignored if `is_midsim` is `False`.

    Attributes (After)
    ----------
    after_sim_pipeline : SimConfList
        List of all post-simulation animations to be animated.
    is_after_sim : bool
        ``True`` only if this configuration enables (but _not_ necessarily
        displays or saves) post-simulation animations.
    is_after_sim_show : bool
        ``True`` only if this configuration displays post-simulation animations.
        Ignored if `is_after_sim` is `False`.
    is_after_sim_save : bool
        ``True`` only if this configuration saves post-simulation animations.
        Ignored if `is_after_sim` is `False`.

    Attributes (Images)
    ----------
    is_images_save : bool
        ``True`` only if this configuration saves animation frames as images.
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        `is_images_save` is `False`.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if `is_images_save` is `False`.

    Attributes (Video)
    ----------
    is_video_save : bool
        ``True`` only if this configuration saves animation frames as video.
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
        ``is_video_save`` is ``False``. Supported names include:
        * ``ffmpeg``, an open-source cross-platform audio and video encoder.
        * ``avconv``, an open-source cross-platform audio and video encoder
          forked from (and largely interchangeable with) `ffmpeg`.
        * ``mencoder``, an open-source cross-platform audio and video encoder
          associated with MPlayer, a popular media player.
        * ``imagemagick``, an open-source cross-platform image manipulation
          suite supporting *only* creation of animated GIFs.
    '''

    # ..................{ ALIASES ~ while                    }..................
    is_while_sim = conf_alias(
        "['results options']['while solving']['animations']['enabled']", bool)
    is_while_sim_save = conf_alias(
        "['results options']['while solving']['animations']['save']", bool)
    is_while_sim_show = conf_alias(
        "['results options']['while solving']['animations']['show']", bool)

    # ..................{ ALIASES ~ after                    }..................
    is_after_sim = conf_alias(
        "['results options']['after solving']['animations']['enabled']", bool)
    is_after_sim_save = conf_alias(
        "['results options']['after solving']['animations']['save']", bool)
    is_after_sim_show = conf_alias(
        "['results options']['after solving']['animations']['show']", bool)

    # ..................{ ALIASES ~ save : images            }..................
    is_images_save = conf_alias(
        "['results options']['save']['animations']['images']['enabled']", bool)
    image_filetype = conf_alias(
        "['results options']['save']['animations']['images']['filetype']", str)
    image_dpi = conf_alias(
        "['results options']['save']['animations']['images']['dpi']", int)

    # ..................{ ALIASES ~ save : video             }..................
    is_video_save = conf_alias(
        "['results options']['save']['animations']['video']['enabled']", bool)
    video_bitrate = conf_alias(
        "['results options']['save']['animations']['video']['bitrate']", int)
    video_dpi = conf_alias(
        "['results options']['save']['animations']['video']['dpi']", int)
    video_filetype = conf_alias(
        "['results options']['save']['animations']['video']['filetype']", str)
    video_framerate = conf_alias(
        "['results options']['save']['animations']['video']['framerate']", int)
    video_metadata = conf_alias(
        "['results options']['save']['animations']['video']['metadata']",
        MappingType)
    video_writer_names = conf_alias(
        "['results options']['save']['animations']['video']['writers']",
        SequenceTypes)
    video_codec_names = conf_alias(
        "['results options']['save']['animations']['video']['codecs']",
        SequenceTypes)

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        self.after_sim_pipeline = SimConfList(
            confs=self._conf[
                'results options']['after solving']['animations']['pipeline'],
            conf_type=SimConfAnimOne,
        )

        #FIXME: Actually initialize to a valid "SimConfList", once the codebase
        #supports general-purpose in-simulation animations.
        self.while_sim_pipeline = []

        # Validate all configured integers to be positive.
        ints.die_unless_positive(
            self.image_dpi,
            self.video_bitrate,
            self.video_dpi,
            self.video_framerate,
        )
