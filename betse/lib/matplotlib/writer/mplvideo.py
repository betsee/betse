#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific classes writing animations as video.
'''

#FIXME: Consider contributing most or all of this submodule back to matplotlib.

# ....................{ IMPORTS                            }....................
from matplotlib.animation import FFMpegBase, MovieWriter, writers
from betse.exceptions import BetseMatplotlibException
from betse.util.io.log import logs
from betse.util.os.shell import shellstr
from betse.util.path.command import cmdrun
from betse.util.type.cls import classes
from betse.util.type.mapping.mapcls import ReversibleDict
from betse.util.type.text import regexes, strs
from betse.util.type.types import (
    type_check, ClassType, NoneType, SequenceTypes, StrOrNoneTypes)

# ....................{ DICTS                              }....................
WRITER_NAME_TO_COMMAND_BASENAME = ReversibleDict(
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
WRITER_BASENAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES = None
'''
Dictionary mapping from the basename of the external command for each video
encoder recognized by matplotlib (e.g., ``ffmpeg``) to a nested dictionary
mapping from the filetype of each video container format recognized by BETSE
(e.g., ``mp4``) to a list of the names of all widely used video codecs supported
by that encoder (in descending order of general-purpose, subjective preference).

Since multiple matplotlib animation writers (e.g., ``ffmpeg``, ``ffmpeg_file``)
typically execute the same external command (e.g., ``ffmpeg``), this dictionary
maps from the latter rather than the former. Mapping from the former would
needlessly require redundant mappings for writers sharing the same command.

This dictionary is principally used by the
:func:`get_first_codec_name_supported` getter function to obtain the preferred
codec for a given combination of video writer and container format.
'''

# ....................{ INITIALIZERS                       }....................
# For simplicity, this function is called below on the first importation of this
# submodule rather than explicitly called by callers.
def _init() -> None:
    '''
    Initialize this submodule.
    '''

    _init_globals()


def _init_globals() -> None:
    '''
    Initialize all global variables declared by this submodule.
    '''

    # Globals initialized below.
    global WRITER_BASENAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES

    WRITER_BASENAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES = {
        # FFmpeg.
        'ffmpeg': {
            # Audio Video Interleave (AVI). AVI supports the smallest subset of
            # MPEG-centric codecs of all recognized container formats and hence
            # serves as the baseline for listing such codecs.
            'avi': (
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
            ),

            # Graphics Interchange Format (GIF).
            'gif': ('gif',),

            # Matroska.
            'mkv': None,

            # QuickTime.
            'mov': None,

            # MPEG-4 Part 14.
            'mp4': None,

            # Theora.
            'ogv': ('libtheora',),

            # WebM.
            'webm': (
                #FIXME: As of this writing (i.e., mid-2017), VP-9 is poorly
                #supported in media players (e.g., mpv, vlc) and hence disabled.
                #Reevaluate and consider enabling this otherwise sensible
                #default at some future date.
                # 'libvpx-vp9',  # WebM VP-9

                'libvpx',      # WebM VP-8
            ),
        },

        # Libav.
        'avconv': None,

        # Mencoder.
        'mencoder': None,

        # ImageMagick, encoding videos only as old-school animated GIFs. Since
        # ImageMagick is *NOT* a general-purpose video encoder and hence fails
        # to support the notion of video codecs, a codec of "None" is used.
        # This has beneficial side effects, including ensuring that no edge-case
        # functionality in the matplotlib codebase attempts to erroneously
        # encode anything with ImageMagick using a codec.
        'convert': {
            'gif': (None,),
        }
    }

    # Shorthand to preserve sanity below.
    codec_names = WRITER_BASENAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES

    # For FFmpeg, define the set of all codecs supported by the "mp4" (i.e.,
    # MPEG-4 Part 14) and "mov" (i.e., QuickTime) container formats to be the
    # those supported by the now-obsolete AVI container format plus all newer
    # codecs *NOT* supported by that format.
    codec_names['ffmpeg']['mp4'] = (
        # H.264 / AVC / MPEG-4 Part 10.
        'libx264',
        # H.265 / HEVC / MPEG-H Part 2. Since H.265 is significantly less
        # well-supported than H.264, the former takes precedence.
        'hevc',
    # All older codecs supported by the now-obsolete AVI container format. Since
    # the first codec is the aforementioned "libx264", this codec (and only this
    # codec) is ignored with a simple slice.
    ) + codec_names['ffmpeg']['avi'][1:]
    codec_names['ffmpeg']['mov'] = codec_names['ffmpeg']['mp4']

    # For FFmpeg, define the set of all codecs supported by the "mkv" (i.e.,
    # Matroska) container format to the set of all codecs supported by all other
    # container formats (excluding animated GIFs), giving preference to
    # open-standards codecs rather than proprietary codecs. Matroska: rock it!
    codec_names['ffmpeg']['mkv'] = (
        codec_names['ffmpeg']['webm'] +
        codec_names['ffmpeg']['ogv'] +
        codec_names['ffmpeg']['mp4']
    )

    # Define Libav to support exactly all codecs supported by FFmpeg. Since the
    # two are well-synchronized forks of each other attempting (and mostly
    # succeeding) to preserve a common command-line API, assuming codec parity
    # is typically a safe assumption.
    codec_names['avconv'] = codec_names['ffmpeg']

    # Define Mencoder to support exactly all codecs supported by FFmpeg. Since
    # Mencoder internally leverages the same "libavcodec" shared library
    # leveraged by FFmpeg, assuming codec parity is *ALWAYS* a safe assumption.
    # While Mencoder also provides a small handful of Mencoder-specific codecs
    # (see "mencoder -ovc help"), these codecs are commonly regarded as inferior
    # to their "libavcodec" counterparts. In either case, matplotlib internally
    # mandates use of "libavcodec"-provided codecs rather than Mencoder-specific
    # codecs in all Mencoder writer classes (e.g., "MencoderWriter").
    codec_names['mencoder'] = codec_names['ffmpeg']

# ....................{ EXCEPTIONS                         }....................
def die_unless_writer(writer_name: str) -> None:
    '''
    Raise an exception unless a matplotlib animation writer class registered
    with the passed name is recognized by both BETSE and matplotlib, implying
    the external command required by this class to be installed on the current
    system.

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to validate.

    See Also
    ----------
    :func:`is_writer`
        Further details.
    '''

    # For human-readable granularity in exception messages, call granular
    # testers rather than the catch-all is_writer() tester.
    #
    # If this writer is unrecognized by BETSE, raise an exception.
    die_unless_writer_betse(writer_name)

    # If this writer is unrecognized by matplotlib...
    if not is_writer_mpl(writer_name):
        # Basename of this writer's command.
        writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[writer_name]

        # Raise this exception.
        raise BetseMatplotlibException(
            'Matplotlib animation video writer "{}" '
            'not registered with matplotlib '
            '(i.e., command "{}" not found).'.format(
                writer_name, writer_basename))


def die_unless_writer_betse(writer_name: str) -> None:
    '''
    Raise an exception unless a matplotlib animation writer class registered
    with the passed name is recognized by BETSE.

    If an exception is *not* raised, this does *not* imply the external command
    required by this class to be runnable or even installed on this system.

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to validate.

    See Also
    ----------
    is_writer_betse
        Further details.
    '''

    if not is_writer_betse(writer_name):
        raise BetseMatplotlibException(
            'Matplotlib animation video writer "{}" '
            'unrecognized by BETSE.'.format(writer_name))

# ....................{ EXCEPTIONS ~ command               }....................
def die_unless_writer_command(writer_basename: str) -> None:
    '''
    Raise an exception unless at least one matplotlib animation writer class
    running the external command with the passed basename is recognized by both
    BETSE and matplotlib, implying this command to be installed on the current
    system.

    Parameters
    ----------
    writer_basename : str
        Basename of the external command of the writer to validate.

    See Also
    ----------
    is_writer_command
        Further details.
    '''

    if not is_writer_command(writer_basename):
        raise BetseMatplotlibException(
            'Matplotlib animation video writer command "{basename}" '
            'unrecognized by BETSE or not registered with matplotlib '
            '(i.e., command "{basename}" not found).'.format(
                basename=writer_basename))

# ....................{ TESTERS                            }....................
@type_check
def is_writer(writer_name: str) -> bool:
    '''
    ``True`` only if a matplotlib animation writer class (e.g., `FFMpegWriter`)
    registered with the passed name (e.g., `ffmpeg`) is recognized by both BETSE
    and matplotlib, implying the external command required by this class (e.g.,
    ``fmpeg`) to be installed on the current system.

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to test.

    Returns
    ----------
    bool
        ``True`` only if this writer is recognized by BETSE.

    See Also
    ----------
    :func:`is_writer_betse`, :func:`is_writer_mpl`
        Further details.
    '''

    return is_writer_betse(writer_name) and is_writer_mpl(writer_name)


@type_check
def is_writer_betse(writer_name: str) -> bool:
    '''
    ``True`` only if a matplotlib animation writer class (e.g., `FFMpegWriter`)
    with the passed name (e.g., `ffmpeg`) is recognized by BETSE.

    Specifically, this function returns ``True`` only if this name is a key of
    the global `WRITER_NAME_TO_COMMAND_BASENAME` dictionary of this submodule.
    Note this does *not* imply the external command required by this class
    (e.g., ``fmpeg`) to be installed on the current system.

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to test.

    Returns
    ----------
    bool
        ``True`` only if this writer is recognized by BETSE.
    '''

    return writer_name in WRITER_NAME_TO_COMMAND_BASENAME


@type_check
def is_writer_mpl(writer_name: str) -> bool:
    '''
    ``True`` only if a matplotlib animation writer class (e.g.,
    ``FFMpegWriter``) with the passed name (e.g., ``ffmpeg``) is registered with
    matplotlib, implying the external command required by this class (e.g.,
    ``fmpeg``) to be installed on the current system.

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to test.

    Returns
    ----------
    bool
        ``True`` only if this writer is recognized by matplotlib.
    '''

    return writers.is_available(writer_name)

# ....................{ TESTERS ~ command                  }....................
@type_check
def is_writer_command(writer_basename: str) -> bool:
    '''
    ``True`` only if at least one matplotlib animation writer class (e.g.,
    `MencoderWriter`) running the external command with the passed basename
    (e.g., `mencoder`) is recognized by both BETSE and matplotlib, implying this
    command to be installed on the current system.

    Parameters
    ----------
    writer_basename : str
        Basename of the external command of the writer to test.

    Returns
    ----------
    bool
        ``True`` only if this command is recognized by both BETSE and matplotlib.
    '''

    # If this command is run by at least one writer...
    if writer_basename in WRITER_NAME_TO_COMMAND_BASENAME.reverse:
        # Tuple of the names of all writers running this command.
        writer_names = WRITER_NAME_TO_COMMAND_BASENAME.reverse[writer_basename]

        # If at least one such writer is recognized, return True.
        for writer_name in writer_names:
            if is_writer(writer_name):
                return True

    # Else, no such writers exist. Return False.
    return False


@type_check
def is_writer_command_codec(
    writer_basename: str, codec_name: StrOrNoneTypes) -> bool:
    '''
    ``True`` only if the matplotlib animation writer class running the external
    command with the passed basename (e.g., ``ffmpeg``) supports the video codec
    with the passed encoder-specific name (e.g., ``libx264``).

    Specifically, this function returns ``True`` only if the passed basename is:

    * ``ffmpeg`` and the ``ffmpeg -help encoder={codec_name}`` command succeeds.
    * ``avconv`` and the ``avconv -help encoder={codec_name}`` command succeeds.
    * ``mencoder``, the ``mencoder -ovc help`` command lists the
      Mencoder-specific ``lavc`` video codec and either:
      * ``ffmpeg`` is in the current ``${PATH}`` and recursively calling this
        function as ``is_writer_codec('ffmpeg', codec_name)`` returns ``True``.
      * ``ffmpeg`` is *not* in the current ``${PATH}``, in which case this
        function assumes the passed codec to be supported and returns ``True``.
    * Any other passed basename (e.g., ``convert``, implying ImageMagick) *and*
      the passed codec is ``None``. These basenames are assumed to *not* actually
      be video encoders and thus support no video codecs.

    Parameters
    ----------
    writer_basename : str
        Basename of the external command of the video encoder to test.
    codec_name : str
        Encoder-specific name of the codec to be tested for.

    Returns
    ----------
    bool
        ``True`` only if this encoder supports this codec.

    Raises
    ----------
    BetseMatplotlibException
        If this basename is either:
        * Unrecognized by BETSE.
        * Unregistered with matplotlib, implying
        * Not found as an external command in the current `${PATH}`.
        * Mencoder and the `mencoder -ovc help` command fails to list the
          Mencoder-specific `lavc` video codec required by Matplotlib.
    '''

    # Log this detection attempt.
    logs.log_debug(
        'Detecting encoder "%s" codec "%s"...',
        writer_basename, codec_name)

    # Absolute path of this command.
    writer_filename = get_writer_command_filename(writer_basename)

    # For FFmpeg, detect this codec by capturing help documentation output by
    # the "ffmpeg" command for this codec and grepping this output for a string
    # stating this codec to be unrecognized. Sadly, this command reports success
    # rather than failure when this codec is unrecognized. (wut, FFmpeg?)
    if writer_basename == 'ffmpeg':
        # List of all shell words comprising the command to be tested.
        ffmpeg_command_words = (
            writer_filename,
            '-help',
            'encoder=' + shellstr.shell_quote(codec_name),
        )

        # Help documentation for this codec captured from "ffmpeg".
        ffmpeg_codec_help = cmdrun.get_stdout_or_die(ffmpeg_command_words)

        # Return whether this documentation is suffixed by a string implying
        # this codec to be unrecognized or not. If this codec is unrecognized,
        # this documentation ends with the following line:
        #
        #     Codec '${codec_name}' is not recognized by FFmpeg.
        return not ffmpeg_codec_help.endswith("' is not recognized by FFmpeg.")
    # For Libav, detect this codec in the same exact manner as for FFmpeg.
    elif writer_basename == 'avconv':
        # List of all shell words comprising the command to be tested.
        avconv_command_words = (
            writer_filename,
            '-help',
            'encoder=' + shellstr.shell_quote(codec_name),
        )

        # Help documentation for this codec captured from "avconv".
        avconv_codec_help = cmdrun.get_stdout_or_die(avconv_command_words)
        # print('avconv_command_words: {}'.format(avconv_command_words))
        # print('avconv_codec_help: {}'.format(avconv_codec_help))

        # Return whether this documentation is suffixed by an indicative string.
        return not avconv_codec_help.endswith("' is not recognized by Libav.")
    # For Mencoder, detect this codec by capturing help documentation output by
    # the "mencoder" command for *ALL* video codecs, grepping this output for
    # the "lavc" video codec required by matplotlib, and, if found, repeating
    # the above FFmpeg-specific logic to specifically detect this codec.
    elif writer_basename == 'mencoder':
        # Help documentation for all codecs captured from "mencoder".
        mencoder_codecs_help = cmdrun.get_stdout_or_die((
            writer_filename, '-ovc', 'help'))
        # print('mencoder codecs help: ' + mencoder_codecs_help)

        # If this output contains a line resembling the following, this
        # installation of Mencoder supports the requisite "lavc" codec:
        #     lavc     - libavcodec codecs - best quality!
        if regexes.is_match_line(
            text=mencoder_codecs_help, regex=r'^\s+lavc\s+'):
            # If the "ffmpeg" command is installed on the current system, query
            # that command for whether or not the passed codec is supported.
            # Note that the recursion bottoms out with this call, as the above
            # logic handling the FFmpeg writer does *NOT* recall this function.
            if is_writer_command('ffmpeg'):
                return is_writer_command_codec('ffmpeg', codec_name)
            # Else, "ffmpeg" is *NOT* in the ${PATH}. Since Mencoder implements
            # "lavc" codec support by linking against the "libavcodec" shared
            # library rather than calling the "ffmpeg" command, it's technically
            # permissible (albeit uncommon) for the "mencoder" but not "ffmpeg"
            # command to be in the ${PATH}. Hence, this does *NOT* indicate a
            # fatal error. This does, however, prevent us from querying whether
            # or not the passed codec is supported. In lieu of sensible
            # alternatives...
            else:
                # Log a non-fatal warning.
                logs.log_warning(
                    'Mencoder "libavcodec"-based video codec "{}" '
                    'possibly unavailable. Consider installing FFmpeg to '
                    'resolve this warning.'.format(codec_name))

                # Assume the passed codec to be supported.
                return True
        # Else, Mencoder fails to support the "lavc" codec. Raise an exception.
        else:
            raise BetseMatplotlibException(
                'Mencoder video codec "lavc" unavailable.')

    # For any other writer (e.g., ImageMagick), assume this writer to *NOT* be a
    # video encoder and hence support *NO* video codecs. In this case, return
    # True only if the passed codec is "None" -- signifying "no video codec."
    return codec_name is None

# ....................{ GETTERS                            }....................
@type_check
def get_writer_class(writer_name: str) -> ClassType:
    '''
    Matplotlib animation writer class (e.g., ``ImageMagickWriter``) registered
    with the passed name (e.g., ``imagemagick``).

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to obtain.

    Returns
    ----------
    ClassType
        Matplotlib animation writer class registered with this name.

    Raises
    ----------
    BetseMatplotlibException
        If this writer is unrecognized by either matplotlib or BETSE itself.
    '''

    # If this writer is unrecognized, raise an exception.
    die_unless_writer(writer_name)

    # Return this writer's class.
    return writers[writer_name]


@type_check
def get_writer_command_filename(writer_basename: str) -> str:
    '''
    Absolute path (e.g., ``/usr/bin/convert``) of the external command with the
    passed basename (e.g., ``convert``) run by a matplotlib animation writer
    class (e.g., ``ImageMagickWriter``).

    Parameters
    ----------
    writer_name : str
        Alphanumeric lowercase name of the writer to query.

    Returns
    ----------
    str
        Absolute path of this writer's command.

    Raises
    ----------
    BetseMatplotlibException
        If this writer is unrecognized by either matplotlib or BETSE itself.
    '''

    # If this command is unrecognized, raise an exception.
    die_unless_writer_command(writer_basename)

    # Tuple of the names of all writers running this command.
    writer_names = WRITER_NAME_TO_COMMAND_BASENAME.reverse[writer_basename]

    # Name of the first writer running this command. Since the absolute path of
    # this command is the same across all writers running this command, the
    # first writer is arbitrarily selected merely for simplicity.
    writer_name = writer_names[0]

    # This writer's class.
    writer_class = get_writer_class(writer_name)

    # Return the absolute path of this writer's command.
    return writer_class.bin_path()

# ....................{ GETTERS ~ first                    }....................
@type_check
def get_first_writer_name(writer_names: SequenceTypes) -> str:
    '''
    First name (e.g., ``imagemagick``) of the matplotlib animation writer class
    (e.g., ``ImageMagickWriter``) in the passed list recognized by both BETSE
    and matplotlib if any *or* raise an exception otherwise.

    This function iteratively searches for external commands in the same order
    as the passed list lists names.

    Parameters
    ----------
    writer_names : SequenceTypes
        List of the alphanumeric lowercase names of all writers to search for.

    Returns
    ----------
    str
        Alphanumeric lowercase name of the first such writer.

    Raises
    ----------
    BetseMatplotlibException
        If either:
        * Any writer in the passed list is unrecognized by BETSE.
        * No such writer is registered with matplotlib.
    '''

    # For the name of each passed writer...
    for writer_name in writer_names:
        # If this writer is unrecognized by BETSE, raise an exception. This
        # prevents BETSE from silently ignoring newly added writers not yet
        # recognized by BETSE.
        die_unless_writer_betse(writer_name)

        # If this writer is recognized by matplotlib, cease searching.
        if is_writer_mpl(writer_name):
            return writer_name

    # Else, no such command is in the ${PATH}. Raise an exception.
    raise BetseMatplotlibException(
        'Matplotlib animation video writers {} not found.'.format(
            strs.join_as_conjunction_double_quoted(*writer_names)))


@type_check
def get_first_codec_name(
    writer_name: str,
    container_filetype: str,
    codec_names: SequenceTypes,
) -> (str, NoneType):
    '''
    Name of the first video codec (e.g., ``libx264``) in the passed list
    supported by both the encoder with the passed matplotlib-specific name
    (e.g., ``ffmpeg``) and the container format with the passed filetype (e.g.,
    ``mkv``, ``mp4``) if that writer supports codecs (as most do), ``None`` if
    this writer supports no codecs (e.g., ``imagemagick``) and the passed list
    contains ``None``, *or* raise an exception otherwise (i.e., if no passed
    codecs are supported by both this writer and container format).

    Algorithm
    ----------
    This function iteratively searches for video codecs in the same order as
    listed in the passed list as follows:

    * If there are no remaining video codecs in this list to be examined, an
      exception is raised.
    * If the current video codec has the BETSE-specific name `auto`, the name of
      an intelligently selected codec supported by both this encoder and
      container if any is returned *or* an exception is raised otherwise (i.e.,
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
    codec_names: SequenceTypes
        Sequence of the encoder-specific names of all codecs to search for (in
        descending order of preference).

    Returns
    ----------
    str
        Name of the first codec in the passed list supported by both this
        encoder and container.

    Raises
    ----------
    BetseMatplotlibException
        If any of the following errors arise:
        * This writer is either:
          * Unrecognized by matplotlib or BETSE itself.
          * Not found as an external command in the current `${PATH}`.
        * This container format is unsupported by this writer.
        * No codec whose name is in the passed list is supported by both this
          writer and this container format.

    See Also
    ----------
    :func:`is_writer`
        Tester validating this writer.
    '''

    # If this writer is unrecognized, raise an exception.
    die_unless_writer(writer_name)

    # Basename of this writer's command.
    writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[writer_name]

    # Dictionary mapping from the filetype of each video container format to a
    # list of the names of all video codecs supported by this writer.
    container_filetype_to_codec_names = (
        WRITER_BASENAME_TO_CONTAINER_FILETYPE_TO_CODEC_NAMES[writer_basename])

    # If the passed container is unsupported by this writer, raise an exception.
    if container_filetype not in container_filetype_to_codec_names:
        raise BetseMatplotlibException(
            'Video container format "{}" unsupported by '
            'matplotlib animation video writer "{}".'.format(
                container_filetype, writer_name))

    # List of the names of all candidate codecs supported by both this encoder
    # and this container.
    codec_names_supported = container_filetype_to_codec_names[
        container_filetype]

    # List of the names of all candidate codecs to be detected below, with each
    # instance of "auto" in the original list of these names replaced by the
    # list of all codecs supported by both this encoder and container.
    codec_names_candidate = []

    # For the name of each candidate codec...
    for codec_name in codec_names:
        # If this is the BETSE-specific name "auto", append the names of all
        # codecs supported by both this encoder and container.
        if codec_name == 'auto':
            codec_names_candidate.extend(codec_names_supported)
        # Else, append only the name of this codec.
        else:
            codec_names_candidate.append(codec_name)

    # Log this detection attempt.
    logs.log_debug(
        'Detecting encoder "%s" codec from candidates: %r',
        writer_name, codec_names_candidate)

    # For the name of each preferred codec...
    for codec_name in codec_names_candidate:
        # If this encoder supports this codec *AND*...
        if is_writer_command_codec(writer_basename, codec_name):
            # Log this detection result.
            logs.log_debug(
                'Detected encoder "%s" codec "%s".', writer_name, codec_name)

            # If this container is not known to support this codec, log a
            # non-fatal warning. Since what this application thinks it knows
            # and what reality actually is need not coincide, this container
            # could actually still support this codec. Hence, this edge case
            # does *NOT* constitute a hard, fatal error.
            if codec_name not in codec_names_supported:
                logs.log_warning(
                    'Encoder "%s" container "%s" '
                    'not known to support codec "%s".',
                    writer_name, container_filetype, codec_name)

            # Return the name of this codec.
            return codec_name

    # Else, no passed codecs are supported by this combination of writer and
    # container format. Raise an exception.
    raise BetseMatplotlibException(
        'Codec(s) {} unsupported by '
        'encoder "{}" and/or container "{}".'.format(
            strs.join_as_conjunction_double_quoted(*codec_names),
            writer_name,
            container_filetype))

# ....................{ MAKERS                             }....................
#FIXME: Submit these fixes upstream, please. We'd rather not maintain them in
#perpuity, yet again.
@type_check
def make_writer(*args, cls: ClassType, **kwargs) -> MovieWriter:
    '''
    Create and return an instance of the passed matplotlib animation writer
    class, passing all remaining positional and keyword arguments to that
    class' ``__init__()`` method.

    This method wraps the ``__init__()`` method of the passed class with
    additional intelligence, effectively monkey-patching that method without
    unsafely doing so. Specifically, if:

    * A ``codec`` keyword argument is passed with value:
      * ``libx264``, the ``yuv420p`` pixel format is enabled. Note that the
        :meth:`FFMpegBase.output_args` method already does so onl* for the
        ``h264`` codec; since it should do so for both the ``libx264`` *and*
        ``h264`` codecs, this method coerces this behaviour.
      * ``libvpx``, the ``-auto-alt-ref 0`` command-line argument is passed.
        FFMpeg fails with the following fatal exception on failing to pass this
        parameter when selecting this popular codec:

        .. code:: bash

           MovieWriter.run: running command: ffmpeg -f rawvideo -vcodec rawvideo -s 1200x900 -pix_fmt rgba -r 15 -i pipe: -vcodec libvpx -b 1500k -metadata subject=Bioinformatics -metadata comment=Produced by BETSE. -metadata artist=BETSE -metadata genre=Bioinformatics -y /home/leycec/tmp/betse/sample_sim/RESULTS/init_1/anim_while_solving/Vmem/Vmem.mkv
           ffmpeg version 3.2.4 Copyright (c) 2000-2017 the FFmpeg developers
             built with gcc 4.9.4 (Gentoo 4.9.4 p1.0, pie-0.6.4)
           [betse] This run should take approximately 10.930000s to compute...
           [betse] Exporting "Vmem" frame 0 / 8...
           MovieWriter.grab_frame: Grabbing frame.
           Input #0, rawvideo, from 'pipe:':
             Duration: N/A, start: 0.000000, bitrate: 518400 kb/s
               Stream #0:0: Video: rawvideo (RGBA / 0x41424752), rgba, 1200x900, 518400 kb/s, 15 tbr, 15 tbn, 15 tbc
           Please use -b:a or -b:v, -b is ambiguous
           [libvpx @ 0x12bd960] v1.5.0
           [libvpx @ 0x12bd960] Transparency encoding with auto_alt_ref does not work
           Stream mapping:
             Stream #0:0 -> #0:0 (rawvideo (native) -> vp8 (libvpx))
           Error while opening encoder for output stream #0:0 - maybe incorrect parameters such as bit_rate, rate, width or height
           [betse] Exporting "Vmem" frame 1 / 8...
           MovieWriter.grab_frame: Grabbing frame.
           [betse] Exiting prematurely due to fatal error:

           BrokenPipeError: [Errno 32] Broken pipe

           Traceback (most recent call last):
             File "/home/leycec/py/betse/betse/cli/cliabc.py", line 163, in run
               profile_filename=self._profile_filename,
             File "<string>", line 65, in func_type_checked
             File "/home/leycec/py/betse/betse/util/py/profilers.py", line 130, in profile_callable
               profile_filename=profile_filename,
             File "/home/leycec/py/betse/betse/util/py/profilers.py", line 146, in _profile_callable_none
               return call(*args, **kwargs)
             File "/home/leycec/py/betse/betse/cli/climain.py", line 157, in _do
               return subcommand_method()
             File "/home/leycec/py/betse/betse/cli/climain.py", line 230, in _do_init
               return self._sim_runner.init()
             File "/home/leycec/py/betse/betse/science/simrunner.py", line 206, in init
               sim.run_sim_core(cells, p)
             File "/home/leycec/py/betse/betse/science/sim.py", line 929, in run_sim_core
               cells, p, time_steps, time_steps_sampled, anim_cells)
             File "<string>", line 109, in func_type_checked
             File "/home/leycec/py/betse/betse/science/sim.py", line 1238, in _run_sim_core_loop
               anim_cells.plot_frame(time_step=-1)
             File "<string>", line 15, in func_type_checked
             File "/home/leycec/py/betse/betse/science/visual/visabc.py", line 835, in plot_frame
               self._save_frame(time_step_absolute)
             File "/home/leycec/py/betse/betse/science/visual/visabc.py", line 1008, in _save_frame
               **self._writer_savefig_kwargs)
             File "/usr/lib64/python3.4/site-packages/matplotlib/animation.py", line 230, in grab_frame
               dpi=self.dpi, **savefig_kwargs)
             File "/usr/lib64/python3.4/site-packages/matplotlib/figure.py", line 1565, in savefig
               self.canvas.print_figure(*args, **kwargs)
             File "/usr/lib64/python3.4/site-packages/matplotlib/backend_bases.py", line 2232, in print_figure
               **kwargs)
             File "/usr/lib64/python3.4/site-packages/matplotlib/backends/backend_agg.py", line 519, in print_raw
               fileobj.write(renderer._renderer.buffer_rgba()

        For similar reports of the same underlying issue, see also:

        * https://stackoverflow.com/a/41311048/2809027
        * https://git.videolan.org/?p=ffmpeg.git;a=commitdiff;h=44bcb636c1a7052059330c3edbd900246e15ae60

    Parameters
    ----------
    cls : ClassType
        Matplotlib animation writer class to be instantiated.

    All remaining positional and keyword arguments are passed as is to this
    class' ``__init__()`` method.

    Returns
    ----------
    MovieWriter
        Matplotlib animation writer of this type instiated with these arguments.
    '''

    # If this class is *NOT* that of a Matplotlib animation writer, fail.
    classes.die_unless_subclass(subclass=cls, superclass=MovieWriter)

    # Name of the codec to be written if any or None otherwise.
    codec_name = kwargs.get('codec', None)

    # Additional command-line arguments to pass to the external command run by
    # this writer (e.g., "ffmpeg" for the "FFMpegWriter" class).
    writer_command_args = []

    # If this is an FFMpeg-specific writer passed a codec name of...
    if issubclass(cls, FFMpegBase) and codec_name is not None:
        # "libx264", force the "yuv420p" pixel format. See the docstring.
        if codec_name == 'libx264':
            writer_command_args.extend(['-pix_fmt', 'yuv420p'])
        # "libvpx", disable the "auto_alt_ref" feature. See the docstring.
        elif codec_name == 'libvpx':
            writer_command_args.extend(['-auto-alt-ref', '0'])

    # If passing this command additional command-line arguments, do so.
    if writer_command_args:
        kwargs['extra_args'] = writer_command_args

    # Create and return an instance of this class passed these arguments.
    return cls(*args, **kwargs)

# ....................{ MAIN                               }....................
# Initialize all global variables declared by this submodule.
_init()
