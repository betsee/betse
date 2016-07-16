#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific classes writing animations as video.
'''

#FIXME: Consider contributing most or all of this submodule back to matplotlib.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMatplotlibException
from betse.util.path.command import pathables, runners
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

    # FFmpeg-based video encoding with pipe-based writing.
    ffmpeg='ffmpeg',

    # FFmpeg-based video encoding with file-based writing.
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


# Subsequently initialized by the init() function.
WRITER_NAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES = None
'''
Dictionary mapping from the matplotlib-specific name of each video encoder
supported by matplotlib (e.g., `ffmpeg`) to a nested dictionary mapping from
the filetype of each video container format recognized by BETSE (e.g.,
`mp4`) to a list of the names of all widely used video codecs supported by
that encoder (in descending order of widespread preference).

This dictionary is principally used by the `get_first_codec_name()` utility
function to obtain the preferred codec for a given combination of video
encoder and container format.
'''

# ....................{ INITIALIZERS                       }....................
# For simplicity, this function is called below on the first importation of this
# submodule rather than explicitly called by callers.
def init() -> None:
    '''
    Initialize all global variables declared by this submodule.
    '''

    # Globals initialized below.
    global WRITER_NAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES

    WRITER_NAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES = {
        # FFmpeg.
        'ffmpeg': {
            # Audio Video Interleave (AVI). AVI supports the smallest subset of
            # MPEG-centric codecs of all recognized container formats and hence
            # serves as the baseline for listing such codecs.
            'avi': [
                # H.264 / AVC / MPEG-4 Part 10.
                'libx264',
                # MPEG-4 Part 2 Advanced Simple Profile (ASP) via an external
                # shared library, typically providing better quality at lower
                # bitrates than the otherwise equivalent built-in codec for
                # MPEG-4 Part 2 (i.e., "mpeg4").
                'libxvid',
                # MPEG-4 Part 2 Advanced Simple Profile (ASP) support built-in
                # to all FFmpeg installations.
                'mpeg4',
                # H.263. Assuming videos leveraging this codec to be encoded as
                # basic "baseline" H.263 bitstreams, these videos are decodable
                # as is by conventional MPEG-4 decoders. To quote the eponymous
                # Wikipedia article on H.263: "MPEG-4 Part 2 is H.263 compatible
                # in the sense that basic 'baseline' H.263 bitstreams are
                # correctly decoded by an MPEG-4 Video decoder."
                'h263',
                # MPEG-2.
                'mpeg2video',
            ],

            # Graphics Interchange Format (GIF).
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
                'libvpx',      # WebM VP-8
            ],
        },

        # Libav.
        'avconv': None,

        # Mencoder.
        'mencoder': None,

        # ImageMagick, encoding videos only as old-school animated GIFs. Since
        # ImageMagick is *NOT* a general-purpose video encoder and hence fails
        # to support the notion of video codecs, a codec of "None" is used.
        'imagemagick': {
            'gif': [None],
        }
    }

    # Shorthand to preserve sanity below.
    codec_names = WRITER_NAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES

    # For FFmpeg, define the set of all codecs supported by the "mp4" (i.e.,
    # MPEG-4 Part 14) and "mov" (i.e., QuickTime) container formats to be the
    # same superset of those supported by the now-obsolete AVI container format.
    codec_names['ffmpeg']['mp4'] = codec_names['ffmpeg']['avi'] + [
        # H.265 / HEVC / MPEG-H Part 2.
        'hevc',
    ]
    codec_names['ffmpeg']['mov'] = codec_names['ffmpeg']['mp4']

    # For FFmpeg, define the set of all codecs supported by the "mkv" (i.e.,
    # Matroska) container format to the set of all codecs supported by all other
    # container formats, excluding animated GIFs. Matroska: it doth rocketh.
    codec_names['ffmpeg']['mkv'] = (
        codec_names['ffmpeg']['mp4'] +
        codec_names['ffmpeg']['ogv'] +
        codec_names['ffmpeg']['webm']
    )

    # Define Libav to support exactly all codecs supported by FFmpeg. Since the
    # two are well-synchronized forks of each other attempting (and mostly
    # succeeding) to preserve a common command-line API, assuming codec parity
    # is typically a safe assumption.
    codec_names['avconv'] = ['ffmpeg']

    # Define Mencoder to support exactly all codecs supported by FFmpeg. Since
    # Mencoder internally leverages the same "libavcodec" shared library
    # leveraged by FFmpeg, assuming codec parity is *ALWAYS* a safe assumption.
    # While Mencoder also provides a small handful of Mencoder-specific codecs
    # (see "mencoder -ovc help"), these codecs are commonly regarded as inferior
    # to their "libavcodec" counterparts. In either case, matplotlib internally
    # mandates use of "libavcodec"-provided codecs rather than Mencoder-specific
    # codecs in all Mencoder writer classes (e.g., "MencoderWriter").
    codec_names['mencoder'] = ['ffmpeg']

# ....................{ EXCEPTIONS                         }....................
def die_unless_writer(writer_name: str) -> None:
    '''
    Raise an exception unless a video encoder with the passed
    matplotlib-specific name is recognized by both matplotlib and BETSE itself.

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific alphanumeric lowercase name of the encoder to test.

    See Also
    ----------
    is_writer
        Further details.
    '''

    if not is_writer(writer_name):
        raise BetseMatplotlibException(
            'Matplotlib animation video writer "{}" unrecognized.'.format(
                writer_name))


def die_unless_writer_command(writer_name: str) -> None:
    '''
    Raise an exception unless a video encoder with the passed
    matplotlib-specific name is an external command in the current `${PATH}`.

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific alphanumeric lowercase name of the encoder to test.

    See Also
    ----------
    is_writer_command
        Further details.
    '''

    # If this writer is *NOT* in the ${PATH}...
    if not is_writer_command(writer_name):
        # Basename of this writer's command.
        writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[writer_name]

        # Raise an exception.
        raise BetseMatplotlibException(
            'Matplotlib animation video encoder "{}" '
            'not in the current ${{PATH}}.'.format(writer_basename))

# ....................{ TESTERS                            }....................
@type_check
def is_writer(writer_name: str) -> bool:
    '''
    `True` only if a video encoder with the passed matplotlib-specific name
    (e.g., `ffmpeg`) is recognized by both matplotlib and BETSE itself.

    Specifically, this function returns `True` only if this name is a key of:

    * The global `matplotlib.animation.writers` dictionary, implying this
      encoder to have registered a matplotlib animation writer class (e.g.,
      `FFmpegWriter`).
    * The global `WRITER_NAME_TO_COMMAND_BASENAME` dictionary of this submodule,
      implying this encoder to be recognized by BETSE itself.

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific alphanumeric lowercase name of the encoder to test.

    Returns
    ----------
    bool
        `True` only if this encoder is recognized by both matplotlib and BETSE.
    '''

    return (
        # "True" only if this encoder is recognized by matplotlib.
        writers.is_available(writer_name) and
        # "True" only if this encoder is recognized by BETSE.
        writer_name in WRITER_NAME_TO_COMMAND_BASENAME
    )


@type_check
def is_writer_command(writer_name: str) -> bool:
    '''
    `True` only if the video encoder with the passed matplotlib-specific name
    (e.g., `ffmpeg`) is installed as an external command in the current
    `${PATH}`, typically of the same name as the passed name.

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific alphanumeric lowercase name of the encoder to test.

    Returns
    ----------
    bool
        `True` only if this encoder is installed in the current `${PATH}`.

    Raises
    ----------
    BetseMatplotlibException
        If this writer is unrecognized by either matplotlib or BETSE itself.
    '''

    # If this writer is unrecognized, raise an exception.
    die_unless_writer(writer_name)

    # Basename of this encoder's command.
    writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[writer_name]

    # Return whether this command is in the current ${PATH}.
    return pathables.is_pathable(writer_basename)


@type_check
def is_writer_codec(writer_name: str, codec_name: str) -> bool:
    '''
    `True` only if the video encoder with the passed matplotlib-specific name
    (e.g., `ffmpeg`) supports the video codec with the passed encoder-specific
    name (e.g., `libx264`).

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific alphanumeric lowercase name of the encoder to test.
    codec_name : str
        Encoder-specific name of the codec to be tested for.

    Returns
    ----------
    bool
        `True` only if this encoder supports this codec.

    Raises
    ----------
    BetseMatplotlibException
        If this writer is either:
        * Unrecognized by matplotlib or BETSE itself.
        * Not found as an external command in the current `${PATH}`.

    See Also
    ----------
    is_writer_command
        Tester validating this writer.
    '''

    # If this encoder is unrecognized or not in the ${PATH}, raise an exception.
    die_unless_writer_command(writer_name)

    # Basename of this encoder's command.
    writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[writer_name]

    #FIXME: Implement support for all remaining encoders, including ImageMagick.

    # If this is FFmpeg, detect this codec by capturing help documentation
    # output by the "ffmpeg" command for this codec and parsing this output for
    # a string declaring this codec to be unrecognized. Sadly, this command
    # reports success if this codec is erroneously unrecognized! wut, FFmpeg?
    if writer_basename == 'ffmpeg':
        # Help documentation for this codec captured from an "ffmpeg" command.
        ffmpeg_codec_help = runners.run_with_stdout_captured(command_words=(
            'ffmpeg',
            '-help',
            'encoder=' + strs.shell_quote('${codec_name}'.format(codec_name)),
        ))

        # Return whether this documentation is suffixed by a string implying
        # this codec to be unrecognized or not. If this codec is unrecognized,
        # this documentation ends with the following line:
        #
        #     Codec '${codec_name}' is not recognized by FFmpeg.
        return ffmpeg_codec_help.endswith("' is not recognized by FFmpeg.")

# ....................{ GETTERS                            }....................
@type_check
def get_writer_class(writer_name: str) -> type:
    '''
    Get the matplotlib animation writer class (e.g., `ImageMagickWriter`)
    registered with the passed lowercase name (e.g., `imagemagick`).

    Parameters
    ----------
    writer_name : str
        Matplotlib-specific alphanumeric lowercase name of the class to obtain.

    Returns
    ----------
    type
        Matplotlib animation writer class registered with this name.

    Raises
    ----------
    BetseMatplotlibException
        If this writer is unrecognized by either matplotlib or BETSE itself.

    See Also
    ----------
    is_writer
        Tester validating this writer as recognized.
    '''

    # If this writer is unrecognized, raise an exception.
    die_unless_writer(writer_name)

    # Return this writer's class.
    return writers[writer_name]

# ....................{ GETTERS ~ first                    }....................
#FIXME: Revise docstring.
@type_check
def get_first_writer_name(writer_names: Sequence) -> str:
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
        List of the matplotlib-specific alphanumeric lowercase names of all
        encoders to search for.

    Returns
    ----------
    type
        Matplotlib animation writer class of the first such writer to be found.

    Raises
    ----------
    BetseMatplotlibException
        If either:
        * Any encoder in the passed list is unrecognized by matplotlib or BETSE.
        * No corresponding command is found in the current `${PATH}`.

    See Also
    ----------
    is_writer
        Tester validating this writer as recognized.
    '''

    # Basename of the first video encoder command in the ${PATH} if any or None.
    writer_basename = None

    # For the matplotlib-specific name of each passed video encoder...
    for writer_name in writer_names:
        # If this encoder is unrecognized, raise an exception.
        die_unless_writer(writer_name)

        # Basename of this encoder's command.
        writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[writer_name]

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
        Matplotlib-specific alphanumeric lowercase name of the video encoder to
        search for the passed codecs.
    container_filetype: str
        Filetype of the video container format to constrain this search to.
    codec_names: Sequence
        List of the encoder-specific names of all codecs to search for.

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
        If any of the following errors arise:
        * This encoder is either:
          * Unrecognized by matplotlib or BETSE itself.
          * Not found as an external command in the current `${PATH}`.

    See Also
    ----------
    is_writer
        Tester validating this writer.
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
    #    tester function returning "True" if this version of FFmpeg
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
    #  warnings (as described above) for all codecs other than "auto".

    # If this encoder is unrecognized or not in the ${PATH}, raise an exception.
    die_unless_writer_command(writer_name)

    #FIXME: Iteratively call is_writer_codec() here. Wowzer.

# ....................{ MAIN                               }....................
# Initialize all global variables declared by this submodule.
init()
