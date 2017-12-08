#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Pickled simulation** (i.e., previously pickled seed, initialization, and
simulation files) backward compatibility facilities.
'''

# ....................{ IMPORTS                            }....................
import sys
from betse.util.io.log import logs
# from betse.util.type.types import type_check

# ....................{ UPGRADERS                          }....................
def upgrade_sim_imports() -> None:
    '''
    Upgrade the in-memory module and class structure of the active Python
    interpreter to reflect the newest structure of these modules and classes
    expected by the current version of this application.

    This function preserves backward importability (and hence compatibility)
    with all prior supported pickled simulation formats, converting the obsolete
    module and class names imported by these formats into their modern
    equivalents. Specifically, for each obsolete module or class name imported
    by a prior supported pickled simulation format, this function injects an
    in-memory alias mapping from the obselete to modern such name.

    Since the :mod:`pickle` API provides no explicit means of doing so, this
    function modifies the in-memory module and class structure of the active
    Python interpreter *before* the :mod:`pickle` API is invoked to deserialize
    pickled simulation files. Failing to call this function *before*
    deserializing simulation files pickled by older versions of this application
    reliably induces the :mod:`pickle` API to raise obscure and
    non-human-readable :class:`ImportError` exceptions.

    Ideally, leveraging the third-party :mod:`dill` dependency would
    automatically resolve such backward importability issues by transitively
    pickling all imports required by pickled files in those files. Sadly, this
    does *not* currently appear to be the case.

    See Also
    ----------
    :func:`betse.science.config.confcompat.upgrade_sim_conf`
        Further details on which prior formats exactly are supported.
    '''

    # Upgrade package imports to each successive format. For safety, each
    # upgrade is performed in strict chronological order.
    _upgrade_sim_imports_to_0_5_2()
    _upgrade_sim_imports_to_0_6_0()

# ....................{ UPGRADERS ~ 0.5.2                  }....................
def _upgrade_sim_imports_to_0_5_2() -> None:
    '''
    Upgrade the in-memory module and class structure of the active Python
    interpreter to reflect the newest structure of these modules and classes
    expected by version 0.5.2 (i.e., "Happiest Hodgkin") of this application.
    '''

    # Log this upgrade attempt.
    logs.log_debug('Upgrading simulation imports to 0.5.2 format...')

    # Import all modules whose fully-qualified names have been modified.
    from betse.lib.yaml.abc import yamlabc, yamllistabc
    from betse.science import channels
    from betse.science.math import finitediff
    from betse.science.simulate import simphase
    from betse.science.config.visual import confanim, confplot, confvisabc
    from betse.util.type.mapping import mapcls

    # Alias obsolete module names to current module objects.
    sys.modules['betse.science.config.confabc'] = yamlabc
    sys.modules['betse.science.config.visual.confvisualabc'] = confvisabc
    sys.modules['betse.science.finitediff'] = finitediff
    sys.modules['betse.science.tissue.channels'] = channels
    sys.modules['betse.science.plot.plotconfig'] = confplot
    sys.modules['betse.science.plot.anim.animconfig'] = confanim
    sys.modules['betse.science.visual.anim.animconfig'] = confanim
    sys.modules['betse.science.visual.plot.plotconfig'] = confplot
    sys.modules['betse.util.type.mappings'] = mapcls

    # Alias obsolete to current class names.
    yamlabc.SimConfList = yamllistabc.YamlList
    confanim.SimConfAnimOne       = confvisabc.SimConfVisualCellsListItem
    confvisabc.SimConfVisualABC      = confvisabc.SimConfVisualCellsABC
    confvisabc.SimConfVisualMixin    = confvisabc.SimConfVisualCellsYAMLMixin
    confvisabc.SimConfVisualMolecule = confvisabc.SimConfVisualCellsNonYAML
    confvisabc.SimConfVisualGeneric  = confvisabc.SimConfVisualCellsEmbedded
    confvisabc.SimConfVisualListable = confvisabc.SimConfVisualCellsListItem
    confvisabc.SimConfVisual         = confvisabc.SimConfVisualCellsListItem
    confvisabc.SimConfListableVisual = confvisabc.SimConfVisualCellsListItem
    simphase.SimPhaseType = simphase.SimPhaseKind
    sys.modules['betse.science.config.visual.confanim'].SimConfAnim = (
        confanim.SimConfAnimAll)
    sys.modules['betse.science.config.visual.confplot'].SimConfPlot = (
        confplot.SimConfPlotAll)
    sys.modules['betse.science.visual.anim.animconfig'].AnimConfig = (
        confanim.SimConfAnimAll)
    sys.modules['betse.science.visual.plot.plotconfig'].PlotConfig = (
        confplot.SimConfPlotAll)

# ....................{ UPGRADERS ~ 0.6.0                  }....................
def _upgrade_sim_imports_to_0_6_0() -> None:
    '''
    Upgrade the in-memory module and class structure of the active Python
    interpreter to reflect the newest structure of these modules and classes
    expected by version 0.6.0 of this application.
    '''

    # Log this upgrade attempt.
    logs.log_debug('Upgrading simulation imports to 0.6.0 format...')

    # Import all modules whose fully-qualified names have been modified.
    from betse.lib.yaml.abc import yamlabc  #, yamllistabc
    from betse.science.config.model import conftis
    from betse.science.config.visual import confanim, confplot, confvisabc
    from betse.science.tissue import tisprofile, tishandler
    from betse.science.tissue.event import tisevecut, tisevevolt
    from betse.science.tissue.picker import tispickcls, tispickimage
    from betse.science.tissue.picker.tispickcls import TissuePickerPercent
    from betse.science.tissue.picker.tispickimage import (
        TissuePickerImage, TissuePickerImageMask)
    from betse.science.tissue.tisprofile import CutProfile

    # Alias obsolete module names to current module objects.
    sys.modules['betse.lib.yaml.yamlabc'] = yamlabc
    sys.modules['betse.science.config.event.eventcut'] = tisevecut
    sys.modules['betse.science.config.event.eventvoltage'] = tisevevolt
    sys.modules['betse.science.config.export.confanim'] = confanim
    sys.modules['betse.science.config.export.confplot'] = confplot
    sys.modules['betse.science.config.export.confvis'] = confvisabc
    sys.modules['betse.science.config.tissue.conftis'] = conftis
    sys.modules['betse.science.tissue.bitmapper'] = tispickimage
    sys.modules['betse.science.tissue.handler'] = tishandler
    sys.modules['betse.science.tissue.tiscls'] = tisprofile
    sys.modules['betse.science.tissue.tissuecls'] = tisprofile
    sys.modules['betse.science.tissue.tissuepick'] = tispickcls

    # Alias obsolete to current class names.
    tispickcls.TissuePickerBitmap = TissuePickerImage
    tispickcls.TissuePickerRandom = TissuePickerPercent
    tispickimage.BitMapper = TissuePickerImageMask
    tisprofile.TissueCut = CutProfile
