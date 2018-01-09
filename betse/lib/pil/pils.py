#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Pillow, the Python Image Library
(PIL)-compatible fork underlying image I/O performed by this application.
'''

#FIXME: Revisit imageio when the following feature request is resolved in full:
#    https://github.com/imageio/imageio/issues/289

# ....................{ IMPORTS                            }....................
from PIL import Image
from betse.util.io.log import logs
from betse.util.type.types import IterableTypes  # type_check,

# ....................{ GETTERS                            }....................
def get_filetypes() -> IterableTypes:
    '''
    Set-like iterable (specifically, a dictionary view) of all ``.``-prefixed
    image filetypes supported by the current version of Pillow.

    Examples
    ----------
    >>> from betse.lib.pil import pils
    >>> pils.get_filetypes()
    dict_keys(['.flc', '.bmp', '.ppm', '.webp', '.j2k', '.jpf', '.jpe', '.pcd'])
    '''

    # Initialize Pillow if uninitialized.
    #
    # If Pillow is uninitialized, the "Image.EXTENSION" dictionary is empty.
    # Since the "betse.ignition" submodule already initializes Pillow,
    # explicitly doing so here should typically *NOT* be necessary. Since this
    # getter could technically be called from global scope prior to the
    # initialization performed by "betse.ignition" *AND* since this
    # initialization efficiently reduces to a noop if unnecessary, manually
    # initializing Pillow here is cost-free. (Cost-free is the way to be.)
    init()

    # Return this set-like iterable.
    return Image.EXTENSION.keys()

# ....................{ ENUMERATIONS                       }....................
def init() -> None:
    '''
    Initialize Pillow if uninitialized *or* reduce to a noop otherwise (i.e., if
    Pillow is already initialized).
    '''

    # Log this initialization.
    logs.log_debug('Initializing Pillow...')

    # If Pillow is already initialized, this function internally reduces to an
    # efficient noop in the expected manner (e.g., by accessing a private
    # boolean global).
    Image.init()
