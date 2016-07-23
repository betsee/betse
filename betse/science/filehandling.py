#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level pickling facilities for saving and loading cell cluster and
simulation objects.
'''

#FIXME: For clarity, rename this module to "simsaver.py".

# ....................{ IMPORTS                            }....................
from betse.util.path.file import pickles
from betse.util.type.types import type_check
from collections.abc import Sequence
from scipy import interpolate as interp

# ....................{ SAVERS                             }....................
#FIXME: Consider replacing all calls to this function with calls to the
#pickles.save() function and then removing this function. It doesn't appear to
#serve any demonstrable point anymore. Moreover, the name of this function is
#no longer descriptive. It isn't simply used to save simulations anymore; it's
#been abused to save arbitrary collections of arbitrary objects!

@type_check
def saveSim(savePath: str, datadump: Sequence) -> None:
    '''
    Pickle the passed list of objects to the file with the passed path.

    For safety, any simulation object in this list should be pre-sanitized by
    calling the `safe_pickle()` function _before_ calling this function.

    Parameters
    ----------
    savePath : str
        Absolute or relative path to pickle to.
    datadump : Sequence
        List of all objects to be pickled.
    '''

    pickles.save(obj=datadump, filename=savePath)


#FIXME: This definitely works, but it's not quite ideal. Ideally, every object
#attribute should be pickled as is. There appear to be a number of solutions,
#each with corresponding tradeoffs:
#
#1. Implement the __getstate__() and __setstate__() methods in the "Cells",
#   "Parameters", and "Simulator" classes -- perhaps by inheriting a superclass
#   defining these methods. These methods would need to convert unpicklable to
#   picklable attributes (e.g., strings). See also:
#   https://docs.python.org/3/library/pickle.html#pickling-class-instances
#2. Use a third-party library instead of the built-in pickle functionality.
#   Alternatives supporting pickling of at least lambda functions include
#   "picloud" and "dill".
def safe_pickle(sim, p) -> None:
    '''
    Renders the the passed simulation object **pickle-able** (i.e., safely
    passable to the `saveSim()` function).

    This function destructively removes all interpolation functions, colormaps,
    and lambda functions from the passed simulation. For safety, this function
    should be passed a copy (either shallow or deep) of the current simulation
    object rather than this object itself.
    '''

    sim.gj_funk = None

    # Sanitize the simulation.
    for key, valu in vars(sim).items():
        if type(valu) == interp.interp1d or callable(valu):
            setattr(sim,key,None)

    for key, valu in vars(sim.dyna).items():
        if type(valu) == interp.interp1d or callable(valu):
            setattr(sim.dyna,key,None)

    # Sanitize the parameters.
    for key, valu in vars(p).items():
        if type(valu) == interp.interp1d or callable(valu):
            setattr(p,key,None)

# ....................{ LOADERS                            }....................
#FIXME: We should probably perform basic sanity checks on loaded objects --
#namely, that they were previously saved with the same version of BETSE. To do
#that, let's add a "_betse_version" attribute to the "Parameters" class: e.g.,
#
#    from betse import metadata
#    class Parameters:
#        def __init__(self):
#            self._betse_version = metadata.__version__
#
#We could then raise fatal exceptions if the "_betse_version" attribute of a
#loaded "Parameters" object is *NOT* equal to the current value of
#"metadata.__version__". Forests of sun flowers bask in ochre-burnished light!

#FIXME: For orthogonality, rename to load_init_or_sim().
def loadSim(loadPath) -> tuple:
    '''
    Unpickle the 3-tuple `(sim, cells, p)` describing a previously initialized
    or simulated cell cluster from the file with the passed path.

    For safety, the simulation object in this tuple has been sanitized by
    calling the `safe_pickle()` function and hence may _not_ be usable as is.

    Parameters
    ----------
    loadPath : str
        Absolute or relative path to unpickle from.

    Returns
    ----------
    (Simulator, Cells, Parameters)
        3-tuple `(sim, cells, p)` unpickled from this file.
    '''

    # Unpickle these objects.
    sim, cells, p = pickles.load(loadPath)

    #FIXME: Validate these objects.

    # Return these objects.
    return sim, cells, p


#FIXME: For orthogonality, rename to load_seed().
def loadWorld(loadPath) -> tuple:
    '''
    Unpickle the 2-tuple `(cells, p)` describing a previously seeded cell
    cluster from the file with the passed path.

    Parameters
    ----------
    loadPath : str
        Absolute or relative path to unpickle from.

    Returns
    ----------
    (Cells, Parameters)
        2-tuple `(cells, p)` unpickled from this file.
    '''

    # Unpickle these objects.
    cells, p = pickles.load(loadPath)

    #FIXME: Validate these objects.

    # Return these objects.
    return cells, p
