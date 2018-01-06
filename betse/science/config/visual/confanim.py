#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation animation subconfigurations.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.yaml.yamlalias import yaml_alias, yaml_alias_int_positive
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.science.config.visual.confvisabc import (
    SimConfVisualCellsListItem, SimConfVisualCellsEmbedded)
from betse.util.type.types import type_check, MappingType, SequenceTypes

# ....................{ SUBCLASSES                         }....................
class SimConfAnimAll(YamlABC):
    '''
    YAML-backed simulation animations subconfiguration, encapsulating the
    configuration of all animations (both in- and post-simulation) parsed from
    the current YAML-formatted simulation configuration file.

    Attributes (General)
    ----------
    is_overlay_current : bool
        ``True`` only if this configuration overlays streamlines (e.g., electric
        current, concentration flux) onto appropriate animations.

    Attributes (While Solving)
    ----------
    is_while_sim : bool
        ``True`` only if this configuration displays and/or saves in-simulation
        animations.
    is_while_sim_save : bool
        ``True`` only if this configuration saves in-simulation animations.
    is_while_sim_show : bool
        ``True`` only if this configuration displays in-simulation animations.
    anim_while_sim : SimConfVisualCellsEmbedded
        Generic configuration applicable to all in-simulation animations. Ignored if
        :attr:``is_while_sim`` is ``False``.

    Attributes (After Solving)
    ----------
    is_after_sim : bool
        ``True`` only if this configuration displays and/or saves
        post-simulation animations.
    is_after_sim_save : bool
        ``True`` only if this configuration saves post-simulation animations.
    is_after_sim_show : bool
        ``True`` only if this configuration displays post-simulation animations.
    anims_after_sim : YamlList
        YAML-backed list of all post-simulation animations to be animated.
        Ignored if :attr:`is_after_sim` is ``False``.

    Attributes (Images)
    ----------
    is_images_save : bool
        ``True`` only if this configuration saves animation frames as images.
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        :attr:`is_images_save` is ``False``.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if :attr:`is_images_save` is ``False``.

    Attributes (Video)
    ----------
    is_video_save : bool
        ``True`` only if this configuration saves animation frames as video.
    video_bitrate : int
        Bitrate in bits per second of all video files saved by this
        configuration. Ignored if :attr:`is_video_save` is `False`.
    video_codec_names : Sequence
        List of the names of all encoder-specific codecs with which to encode
        animations (in descending order of preference), automatically:
        * Selecting the first codec supported by the selected writer.
        * Replacing any codec named `auto` by the name of a recommended codec
          specific to the selected writer and filetype. For details, see the
          `betse.lib.matplotlib.writer.mplvideo.get_first_codec_name()`
          function.
        Ignored if :attr:`is_video_save` is `False`.
    video_dpi : int
        Dots per inch (DPI) of all frames of all video files saved by this
        configuration. Ignored if :attr:`is_images_save` is `False`.
    video_filetype : str
        Filetype of all video files saved by this configuration. Supported
        filetypes include:
        * `mkv` (Matroska), an open-standard audio and video container
          supporting all relevant codecs and hence the default.
        * `avi`, Microsoft's obsolete proprietary audio and video container.
        * `gif`, a proprietary image format supporting video animations.
        * `ogv` (Theora), Xiph's open-standard audio and video container.
        * `mov` (QuickTime), Apple's proprietary audio and video container.
        * `mp4` (MPEG-4 Part 14), an open-standard audio and video container.
        * `webm` (WebM), Google's proprietary audio and video container.
        Ignored if :attr:`is_video_save` is `False`.
    video_framerate : int
        Framerate in frames per second of all video files saved by this
        configuration. Ignored if :attr:`is_video_save` is `False`.
    video_metadata : dict
        Dictionary mapping from the alphabetic lowercase name of video metadata
        supported by the active video encoder to that metadata's human-readable
        string to be embedded in all video files saved by this configuration.
        Supported names include: `title`, `artist`, `genre`, `subject`,
        `copyright`, `srcform`, and `comment`. If this dictionary does *not*
        contain a `copyright` key, such a key will be automatically synthesized
        from the current year. Ignored if :attr:`is_video_save` is `False`.
    video_writer_names : Sequence
        List of the names of all matplotlib animation writers with which to
        encode animations (in order of descending preference), automatically
        selecting the first writer installed on the current system. Supported
        names include:
        * ``ffmpeg``, an open-source cross-platform audio and video encoder.
        * ``avconv``, an open-source cross-platform audio and video encoder
          forked from (and largely interchangeable with) `ffmpeg`.
        * ``mencoder``, an open-source cross-platform audio and video encoder
          associated with MPlayer, a popular media player.
        * ``imagemagick``, an open-source cross-platform image manipulation
          suite supporting *only* creation of animated GIFs.
        Ignored if :attr:`is_video_save` is ``False``.
    '''

    # ..................{ ALIASES                            }..................
    is_overlay_current = yaml_alias(
        "['results options']['overlay currents']", bool)

    # ..................{ ALIASES ~ while                    }..................
    is_while_sim_save = yaml_alias(
        "['results options']['while solving']['animations']['save']", bool)
    is_while_sim_show = yaml_alias(
        "['results options']['while solving']['animations']['show']", bool)

    # ..................{ ALIASES ~ after                    }..................
    is_after_sim_save = yaml_alias(
        "['results options']['after solving']['animations']['save']", bool)
    is_after_sim_show = yaml_alias(
        "['results options']['after solving']['animations']['show']", bool)

    # ..................{ ALIASES ~ save : images            }..................
    is_images_save = yaml_alias(
        "['results options']['save']['animations']['images']['enabled']", bool)
    image_filetype = yaml_alias(
        "['results options']['save']['animations']['images']['filetype']", str)
    image_dpi = yaml_alias_int_positive(
        "['results options']['save']['animations']['images']['dpi']")

    # ..................{ ALIASES ~ save : video             }..................
    is_video_save = yaml_alias(
        "['results options']['save']['animations']['video']['enabled']", bool)
    video_bitrate = yaml_alias_int_positive(
        "['results options']['save']['animations']['video']['bitrate']")
    video_dpi = yaml_alias_int_positive(
        "['results options']['save']['animations']['video']['dpi']")
    video_filetype = yaml_alias(
        "['results options']['save']['animations']['video']['filetype']", str)
    video_framerate = yaml_alias_int_positive(
        "['results options']['save']['animations']['video']['framerate']")
    video_metadata = yaml_alias(
        "['results options']['save']['animations']['video']['metadata']",
        MappingType)
    video_writer_names = yaml_alias(
        "['results options']['save']['animations']['video']['writers']",
        SequenceTypes)
    video_codec_names = yaml_alias(
        "['results options']['save']['animations']['video']['codecs']",
        SequenceTypes)

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level dictionaries with high-level wrappers.
        self.anim_while_sim = SimConfVisualCellsEmbedded()

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        self.anims_after_sim = SimConfVisualCellsListItem.make_list()

    # ..................{ LOADERS                            }..................
    #FIXME: Default the "copyright" entry of video metadata to
    #"@ {}".format(current_year)".
    def load(self, *args, **kwargs) -> None:

        # Load our superclass with all passed arguments.
        super().load(*args, **kwargs)

        # Load all subconfigurations of this configuration.
        self.anim_while_sim.load(conf=self._conf[
            'results options']['while solving']['animations'])
        self.anims_after_sim.load(conf=self._conf[
            'results options']['after solving']['animations']['pipeline'])


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Unload all subconfigurations of this configuration.
        self.anim_while_sim.unload()
        self.anims_after_sim.unload()

    # ..................{ PROPERTIES ~ while                 }..................
    @property
    def is_while_sim(self) -> bool:
        return self.is_while_sim_save or self.is_while_sim_show


    @is_while_sim.setter
    @type_check
    def is_while_sim(self, is_while_sim: bool) -> None:
        self.is_while_sim_save = is_while_sim
        self.is_while_sim_show = is_while_sim

    # ..................{ PROPERTIES ~ after                 }..................
    @property
    def is_after_sim(self) -> bool:
        return self.is_after_sim_save or self.is_after_sim_show


    @is_after_sim.setter
    @type_check
    def is_after_sim(self, is_after_sim: bool) -> None:
        self.is_after_sim_save = is_after_sim
        self.is_after_sim_show = is_after_sim
