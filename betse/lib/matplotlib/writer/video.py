#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific classes writing animations as video.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMatplotlibException
from betse.util.command import pathables
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
Bidirectional dictionary mapping the Matplotlib-specific name of each video
encoder supported by Matplotlib (e.g., `imagemagick`) to and from the basename
of that encoder's external command (e.g., `convert`).

These mappings are accessible as follows:

* `WRITER_NAME_TO_COMMAND_BASENAME`, the forward mapping from encoder names to
  command basenames.
* `WRITER_NAME_TO_COMMAND_BASENAME.reverse`, the reverse mapping from command
  basenames to encoder names.
'''

# ....................{ FUNCTIONS                          }....................
@type_check
def get_first_class(video_writer_names: Sequence) -> type:
    '''
    Get the Matplotlib animation class (e.g., `ImageMagickWriter`) writing the
    first video encoder whose Matplotlib-specific name (e.g., `imagemagick`) is
    in the passed list and whose corresponding command (e.g., `convert`) is in
    the current `${PATH}` if any _or_ raise an exception otherwise (i.e., if no
    such encoder's command is in the current `${PATH}`).

    This function iteratively searches for encoder commands in the same order as
    encoder names are listed in the passed list.

    Parameters
    ----------
    video_writer_names : Sequence
        List of the Matplotlib-specific names of all encoders to search for.

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
            and hence unrecognized by Matplotlib.
        * No corresponding command is found in the current `${PATH}`.
    '''

    # Basename of the first video encoder command in the ${PATH} if any or
    # "None" otherwise.
    video_writer_basename = None

    # For the Matplotlib-specific name of each passed video encoder...
    for video_writer_name in video_writer_names:
        # If this encoder is unrecognized by Matplotlib itself, either upstream
        # changes have broken backward compatibility (unlikely) or the caller
        # specified an invalid encoder name (likely). Raise an exception.
        if not writers.is_available(video_writer_name):
            raise BetseMatplotlibException(
                'Matplotlib animation video writer "{}" '
                'unrecognized by Matplotlib.'.format(video_writer_name))

        # If this encoder is unrecognized by BETSE, upstream changes have
        # probably yet to be incorporated into BETSE. Raise an exception.
        if video_writer_name not in WRITER_NAME_TO_COMMAND_BASENAME:
            raise BetseMatplotlibException(
                'Matplotlib animation video writer "{}" '
                'unrecognized by BETSE.'.format(video_writer_name))

        # Basename of this encoder's command.
        video_writer_basename = WRITER_NAME_TO_COMMAND_BASENAME[
            video_writer_name]

        # If this command is in the ${PATH}, cease searching.
        if pathables.is_pathable(video_writer_basename):
            break
    # Else, no such command is in the ${PATH}. Raise an exception.
    else:
        # Human-readable string listing the names of all passed video encoders.
        video_writer_names_readable = (
            strs.join_as_conjunction_double_quoted(
                *video_writer_names))

        # Raise an exception containing this human-readable string.
        raise BetseMatplotlibException(
            'Matplotlib animation video writers {} '
            'not found in the current ${{PATH}}.'.format(
                video_writer_names_readable))

    # Return the Matplotlib writer class for this encoder.
    return writers[video_writer_name]
