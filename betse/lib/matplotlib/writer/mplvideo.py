#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific classes writing animations as video.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMatplotlibException
from betse.util.path.command import pathables
from betse.util.type import strs
from betse.util.type.mappings import bidict
from betse.util.type.types import type_check, Sequence
from matplotlib.animation import writers

# ....................{ DICTS                              }....................
WRITER_NAME_TO_COMMAND_BASENAME = bidict(
    # AVConv-based video encoding with pipe-based writing.
    avconv='avconv',

    # AVConv-based video encoding with file-based writing.
    avconv_file='avconv',

    # FFMpeg-based video encoding with pipe-based writing.
    ffmpeg='ffmpeg',

    # FFMpeg-based video encoding with file-based writing.
    ffmpeg_file='ffmpeg',

    # Mencoder-based video encoding with pipe-based writing.
    mencoder='mencoder',

    # Mencoder-based video encoding with file-based writing.
    mencoder_file='mencoder',

    # ImageMagick-based animated GIF encoding with pipe-based writing.
    imagemagick='convert',

    # ImageMagick-based animated GIF encoding with file-based writing.
    imagemagick_file='convert',
)
'''
Bidirectional dictionary mapping the matplotlib-specific name of each video
encoder supported by matplotlib (e.g., `imagemagick`) to and from the basename
of that encoder's external command (e.g., `convert`).

These mappings are accessible as follows:

* `WRITER_NAME_TO_COMMAND_BASENAME`, the forward mapping from encoder names to
  command basenames.
* `WRITER_NAME_TO_COMMAND_BASENAME.reverse`, the reverse mapping from command
  basenames to encoder names.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_writer(writer_name: str) -> None:
    '''
    Raise an exception unless a video encoder with the passed
    matplotlib-specific name is recognized by both matplotlib and BETSE itself.

    See Also
    ----------
    is_writer
        Further details.
    '''

    if not is_writer(writer_name):
        raise BetseMatplotlibException(
            'Matplotlib animation video writer "{}" unrecognized.'.format(
                writer_name))

# ....................{ TESTERS                            }....................
@type_check
def is_writer(writer_name: str) -> bool:
    '''
    `True` only if a video encoder with the passed matplotlib-specific name
    (e.g., `ffmpeg`) is recognized by both matplotlib and BETSE itself.

    Specifically, this function returns `True` only if this name is a key of:

    * The global `matplotlib.animation.writers` dictionary, implying this
      encoder to have registered a matplotlib animation writer class (e.g.,
      `FFMpegWriter`).
    * The global `WRITER_NAME_TO_COMMAND_BASENAME` dictionary of this submodule,
      implying this encoder to be recognized by BETSE itself.
    '''

    return (
        # "True" only if this encoder is recognized by matplotlib.
        writers.is_available(writer_name) and
        # "True" only if this encoder is recognized by BETSE.
        writer_name in WRITER_NAME_TO_COMMAND_BASENAME
    )

# ....................{ GETTERS                            }....................

# ....................{ GETTERS ~ first                    }....................
@type_check
def get_first_writer_name(writer_names: Sequence) -> type:
    '''
    Get the matplotlib animation class (e.g., `ImageMagickWriter`) writing the
    first video encoder whose matplotlib-specific name (e.g., `imagemagick`) is
    in the passed list and whose corresponding command (e.g., `convert`) is in
    the current `${PATH}` if any _or_ raise an exception otherwise (i.e., if no
    such encoder's command is in the current `${PATH}`).

    This function iteratively searches for encoder commands in the same order as
    encoder names are listed in the passed list.

    Parameters
    ----------
    writer_names : Sequence
        List of the matplotlib-specific names of all encoders to search for.

    Returns
    ----------
    type
        Matplotlib animation writer class of the first such writer to be found.

    Raises
    ----------
    BetseMatplotlibException
        If either:
        * For any encoder in the passed list, this encoder's name is either:
          * _Not_ a key of the global `WRITER_NAME_TO_COMMAND_BASENAME`
            dictionary and hence unrecognized by BETSE.
          * _Not_ a key of the global `matplotlib.animation.writers` registry
            and hence unrecognized by matplotlib.
        * No corresponding command is found in the current `${PATH}`.
    '''

    # Basename of the first video encoder command in the ${PATH} if any or
    # "None" otherwise.
    writer_basename = None

    # For the matplotlib-specific name of each passed video encoder...
    for writer_name in writer_names:
        # If this encoder is unrecognized, raise an exception.
        die_unless_writer(writer_name)

        # Basename of this encoder's command.
        writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[
            writer_name]

        # If this command is in the ${PATH}, cease searching.
        if pathables.is_pathable(writer_basename):
            break
    # Else, no such command is in the ${PATH}. Raise an exception.
    else:
        # Human-readable string listing the names of all passed video encoders.
        writer_names_readable = (
            strs.join_as_conjunction_double_quoted(*writer_names))

        # Raise an exception containing this human-readable string.
        raise BetseMatplotlibException(
            'Matplotlib animation video writers {} '
            'not found in the current ${{PATH}}.'.format(
                writer_names_readable))

    # Return the matplotlib writer class for this encoder.
    return writers[writer_name]


@type_check
def get_first_codec_name(
    writer_name: str, container_filetype: str, codec_names: Sequence) -> str:
    '''
    Get the name of the first video codec (e.g., `libx264`) in the passed list
    supported by both the encoder with the passed matplotlib-specific name
    (e.g., `ffmpeg`) and the video container with the passed filetype (e.g.,
    `mkv`, `mp4`) if any _or_ raise an exception otherwise (i.e., if no such
    codecs are supported by both this encoder and container).

    This function iteratively searches for video codecs in the same order as
    listed in the passed list as follows:

    * If there are no remaining video codecs in this list to be examined, an
      exception is raised.
    * If the current video codec has the BETSE-specific name `auto`, the name of
      an intelligently selected codec supported by both this encoder and
      container if any is returned _or_ an exception is raised otherwise (i.e.,
      if no codecs are supported by both this encoder and container). Note that
      this codec's name rather than the BETSE-specific name `auto` is returned.
      See this function's body for further commentary.
    * Else if the current video codec is supported by both this encoder and
      container, this codec's name is returned.
    * Else the next video codec in this list is examined.

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific name of the video encoder to search codecs for.
    container_filetype: str
        Filetype of the video container to search codecs for.
    codec_names: Sequence
        List of the names of all codecs to search for.

    Returns
    ----------
    str
        Name of the first codec in the passed list supported by both this
        encoder and container.

    Raises
    ----------
    !!!!!!!!!!
    ! FIXME: Update us up, please!
    !!!!!!!!!!

    BetseMatplotlibException
        If either:
        * For any encoder in the passed list, this encoder's name is either:
          * _Not_ a key of the global `WRITER_NAME_TO_COMMAND_BASENAME`
            dictionary and hence unrecognized by BETSE.
          * _Not_ a key of the global `matplotlib.animation.writers` registry
            and hence unrecognized by Matplotlib.
        * No corresponding command is found in the current `${PATH}`.
    '''

    #FIXME: Terrible default. Instead, let's default to a video codec
    #conditionally dependent on the filetype of the desired video file.
    #Note that this codec is indeed a *VIDEO* rather than *AUDIO* codec, as
    #confirmed by matplotlib's use of the "-vcodec" CLI option to FFmpeg.
    #Makes sense. Matplotlib can hardly reasonably support audio, can it?
    #Note that the filetype uniquely specifies the container format. As a
    #starting place, consider:
    #
    #* If the desired codec is "auto" *AND*
    #* If the desired video encoder is "ffmpeg" and...
    #  * If this filetype is either ".mkv" (i.e., Matroska, a general-
    #    purpose container format that may contain relatively arbitrary
    #    video and audio streams, including WebM, MPEG-4, and H.26*) or
    #    ".webm" and...
    #    * If the "libvpx-vp9" encoder (i.e., WebM VP-9) is available,
    #      prefer that. All WebM encodings are unencumbered by onerous
    #      licensing or patent requirements, unlike both H.264 and H.265.
    #    * Else if the "libvpx" encoder (i.e., WebM VP-8) is available,
    #      fallback to that. VP-9 is newer and hence preferable, but VP-8
    #      suffices as well.
    #    * Else if this filetype is ".mkv", attempt both ".ogv" and ".mp4" codec
    #      dectection as well. Matroska supports a superset of all codecs
    #      supported by traditional MPEG-specific containers.
    #    * Else, raise an exception.
    #  * If this filetype is either ".mov", ".mp4", or ".avi" (three
    #    container formats which appear to be functionally synonymous for
    #    all extents and purposes) and...
    #    * If the "hevc" encoder (i.e., H.265 / HEVC / MPEG-H Part 2) is
    #      available, prefer that.
    #    * Else if the "libx264" encoder (i.e., H.264 / AVC / MPEG-4 Part
    #      10) is available, prefer that.
    #    * Else if the "mpeg4" encoder (i.e., MPEG-4 Part 2) is available,
    #      fallback to that.
    #    * Else if the "libxvid" encoder (i.e., MPEG-4 Part 2) is
    #      available, fallback to that. "libxvid" is an alternate MPEG-4
    #      encoder. I'm unclear whether that or "mpeg4" are preferable, but
    #      given the legacy nature of both... it probably doesn't matter.
    #    * Else if the "h263" encoder (i.e., H.263) is available, fallback
    #      to that. To quote Wikipedia: "MPEG-4 Part 2 is H.263 compatible
    #      in the sense that basic 'baseline' H.263 bitstreams are
    #      correctly decoded by an MPEG-4 Video decoder."
    #    * Else if the "mpeg2video" encoder (i.e., MPEG-2) is available,
    #      fallback to that. Sucks, but what can you do?
    #    * Else, raise an exception.
    #  * If this filetype is ".ogv" and...
    #    * If the "libtheora" encoder (i.e., Theora) is available, use
    #      that.
    #    * Else, raise an exception.
    #  * Else, raise an exception.
    #* If the desired video encoder is "imagemagick" and...
    #  * If this filetype is ".gif", coerce the codec to "None" to ensure
    #    that no edge-case functionality in the matplotlib codebase
    #    attempts to erroneously encode anything with a codec.
    #  * Else, raise an exception.
    #  * Note that a non-fatal warning should probably also be raised if
    #    the codec list is anything other than "['auto']", as all such
    #    codecs will be ignored in that case.
    #
    #To test whether a given encoding codec is available, we'll want to:
    #
    #* If the desired video encoder is "ffmpeg":
    #  * Create a new "betse.libs.ffmpeg" package containing a new
    #    "betse.libs.ffmpeg.ffmpegs" submodule providing an is_encoder()
    #    tester function returning "True" if this version of FFMpeg
    #    supports the encoding codec with the passed name: e.g.,
    #
    #        @type_check
    #        def is_encoder(encoder_name: str) -> bool:
    #
    #    There exist numerous possible implementations of this tester with
    #    the customary tradeoffs, including:
    #    * Capturing the stdout of the external "ffmpeg -encoders" command
    #      and grepping such stdout for a line whose second is the passed
    #      encoder name.
    #    * Capturing the stdout of the following external command, where
    #      ${encoder_name} is the passed encoder name:
    #
    #          $ ffmpeg -help encoder=${encoder_name}
    #
    #      Sadly, this command does *NOT* report non-zero exit status if
    #      the passed encoder is unavailable. It does, however, terminate
    #      its stdout with the following line:
    #
    #          Codec '${encoder_name}' is not recognized by FFmpeg.
    #
    #      Hence, merely grep such output via the str.endswith() builtin
    #      for the suffix "' is not recognized by FFmpeg."
    #  * Call betse.libs.ffmpeg.ffmpegs.is_encoder().
    #* If the desired video encoder is "avconv":
    #  * Equivalent to the "ffmpeg" approach except calling "avconv" rather
    #    than "ffmpeg" -- hopefully, anyway. To test this, we'll probably
    #    want to install and test "avconv" within a docker container.
    #* If the desired video encoder is "mencoder":
    #  * ...we have no idea, at the moment.
    #* If the desired video encoder is "imagemagick", print non-fatal
    #  warnings (as described above) for all codecs other than "auto".:

    # If this encoder is unrecognized, raise an exception.
    die_unless_writer(writer_name)

    #FIXME: Globalize us up. We'll probably need a module-level init() function
    #to do so sanely.
    WRITER_NAME_TO_CONTAINER_FILETYPE_TO_AUTO_CODEC_NAMES = {
        # FFmpeg.
        'ffmpeg': {
            # Audio Video Interleave.
            'avi': None,

            # Graphics Interchange Format.
            'gif': ['gif'],

            # Matroska.
            'mkv': None,

            # QuickTime.
            'mov': None,

            # MPEG-4 Part 14.
            'mp4': None,

            # Theora.
            'ogv': ['libtheora'],

            # WebM.
            'webm': [
                'libvpx-vp9',  # WebM VP-9
                'libvpx',   # WebM VP-8
            ],
        },

        # Libav.
        'avconv': None,

        #FIXME: Implement us up!
        # Mencoder.
        'mencoder': None,

        # ImageMagick, encoding videos only as old-school animated GIFs. Since
        # ImageMagick is *NOT* a general-purpose video encoder and hence fails
        # to support the notion of video codecs, a codec of "None" is used.
        'imagemagick': {
            'gif': [None],
        }
    }

    WRITER_NAME_TO_CONTAINER_FILETYPE_TO_AUTO_CODEC_NAMES
