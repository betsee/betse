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
from betse.util.path import pathnames
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import SetType  # type_check,

# ....................{ GETTERS                            }....................
@func_cached
def get_filetypes() -> SetType:
    '''
    Set of all image filetypes supported by the current version of Pillow.

    For generality, these filetypes are *not* prefixed by a ``.`` delimiter.

    Examples
    ----------
        >>> from betse.lib.pil import pils
        >>> pils.get_filetypes()
        {'flc', 'bmp', 'ppm', 'webp', 'j2k', 'jpf', 'jpe', 'pcd'}
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

    # Return a set of...
    return set(
        # This filetype strips of this prefixing ".".
        pathnames.undot_filetype(filetype_dotted)
        # For each "."-prefixed filetype supported by Pillow...
        for filetype_dotted in Image.EXTENSION.keys()
    )

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
