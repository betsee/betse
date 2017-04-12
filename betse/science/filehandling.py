#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level pickling facilities for saving and loading cell cluster and
simulation objects.
'''

#FIXME: For clarity, rename this module to "simpickler.py".

# ....................{ IMPORTS                            }....................
import sys
from betse.lib.pickle import pickles
from betse.util.type.types import type_check
from collections.abc import Sequence

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

    pickles.save(datadump, filename=savePath, is_overwritable=True)

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

    # Preserve backward importability with obsolete pickled objects.
    _preserve_backward_importability()

    # Unpickle these objects *AFTER* preserving backward importability.
    sim, cells, p = pickles.load(loadPath)

    #FIXME: Validate these objects.

    # Return these objects.
    return sim, cells, p


#FIXME: For orthogonality, rename to load_seed().
def loadWorld(loadPath) -> tuple:
    '''
    Unpickle the 2-tuple ``(cells, p)`` describing a previously seeded cell
    cluster from the file with the passed path.

    Parameters
    ----------
    loadPath : str
        Absolute or relative path to unpickle from.

    Returns
    ----------
    (Cells, Parameters)
        2-tuple ``(cells, p)`` unpickled from this file.
    '''

    # Preserve backward importability with obsolete pickled objects.
    _preserve_backward_importability()

    # Unpickle these objects *AFTER* preserving backward importability.
    cells, p = pickles.load(loadPath)

    #FIXME: Validate these objects.

    # Return these objects.
    return cells, p

# ....................{ LOADERS                            }....................
def _preserve_backward_importability() -> None:
    '''
    Preserve backward compatibility with pickled objects transitively referring
    to modified modules whose fully-qualified names at the time of pickling now
    differ from the current names for those these modules, typically due to
    these modules having since been moved, renamed, or outright removed.
    '''

    # Import all modules whose fully-qualified names have been modified.
    from betse.science import channels
    from betse.science.math import finitediff
    from betse.science.simulate import simphase
    from betse.science.tissue import tissuepick
    from betse.science.config import export
    from betse.science.config.export import confanim, confplot, confvis

    # Alias obsolete module names to current module objects.
    sys.modules['betse.science.config.export.confvisabc'] = confvis
    sys.modules['betse.science.config.visual'] = export
    sys.modules['betse.science.config.visual.confanim'] = confanim
    sys.modules['betse.science.config.visual.confplot'] = confplot
    sys.modules['betse.science.config.visual.confvisualabc'] = confvis
    sys.modules['betse.science.finitediff'] = finitediff
    sys.modules['betse.science.tissue.channels'] = channels
    sys.modules['betse.science.tissue.picker'] = tissuepick
    sys.modules['betse.science.plot.plotconfig'] = confplot
    sys.modules['betse.science.plot.anim.animconfig'] = confanim
    sys.modules['betse.science.visual.anim.animconfig'] = confanim
    sys.modules['betse.science.visual.plot.plotconfig'] = confplot

    # Alias obsolete to current class names.
    confanim.SimConfAnimOne = confvis.SimConfVisualCellsListItem
    confvis.SimConfVisualABC      = confvis.SimConfVisualCellsABC
    confvis.SimConfVisualMixin    = confvis.SimConfVisualCellsYAMLMixin
    confvis.SimConfVisualMolecule = confvis.SimConfVisualCellsNonYAML
    confvis.SimConfVisualGeneric  = confvis.SimConfVisualCellsEmbedded
    confvis.SimConfVisualListable = confvis.SimConfVisualCellsListItem
    confvis.SimConfVisual         = confvis.SimConfVisualCellsListItem
    confvis.SimConfListableVisual = confvis.SimConfVisualCellsListItem
    simphase.SimPhaseType = simphase.SimPhaseKind
    sys.modules['betse.science.config.visual.confanim'].SimConfAnim = (
        confanim.SimConfAnimAll)
    sys.modules['betse.science.config.visual.confplot'].SimConfPlot = (
        confplot.SimConfPlotAll)
    sys.modules['betse.science.visual.anim.animconfig'].AnimConfig = (
        confanim.SimConfAnimAll)
    sys.modules['betse.science.visual.plot.plotconfig'].PlotConfig = (
        confplot.SimConfPlotAll)
