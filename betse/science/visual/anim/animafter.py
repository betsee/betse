#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Post-simulation animation** (i.e., animation produced *after* rather than
*while* solving a simulation) subclasses.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config.export.visual.confvisabc import SimConfVisualCellsABC
from betse.science.phase.phasecls import SimPhase
from betse.science.visual.anim.animabc import AnimCellsABC
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUBCLASSES                         }....................
class AnimCellsAfterSolving(AnimCellsABC):
    '''
    Abstract base class of all post-simulation animation subclasses.

    Each subclass of this class animates arbitrary cell data as a time series
    plotted over the cell cluster (e.g., cell membrane voltage as a function of
    time) *after* rather than *while* solving simulations.
    '''

    @type_check
    def __init__(
        self,
        phase: SimPhase,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this post-simulation animation.

        Parameters
        ----------
        phase: SimPhase
            Current simulation phase.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize our superclass.
        super().__init__(
            *args,

            # Pass this simulation phase as is to our superclass.
            phase=phase,

            # Save and show this post-simulation animation only if this
            # configuration enables doing so.
            is_save=phase.p.anim.is_after_sim_save,
            is_show=phase.p.anim.is_after_sim_show,

            # Save all post-simulation animations to the same parent directory.
            save_dir_parent_basename='anim',

            # Pass all remaining arguments as is to our superclass.
            **kwargs
        )


#FIXME: Merge this class into the "AnimCellsAfterSolving" superclass *AFTER*
#refactoring all subclasses to leverage layers.
class AnimCellsAfterSolvingLayered(AnimCellsAfterSolving):
    '''
    Post-simulation animation of arbitrary cell data as a time series (e.g.,
    cell membrane voltage as a function of time), plotted over the cell
    cluster using a predefined sequence of layers.
    '''

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def __init__(
        self,
        conf: SimConfVisualCellsABC,
        layers: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        conf : SimConfVisualCellsABC
            Configuration for this visual, synchronized with the user-defined
            YAML-backed simulation configuration file for this phase.
        layers : SequenceTypes
            Sequence of all :class:`LayerCellsABC` instances collectively
            plotting each frame of this animation. See the
            :meth:`VisualCellsABC.__init__` method for further details.

        All remaining parameters are passed as is to the superclass
        :meth:`AnimCellsAfterSolving.__init__` method.
        '''

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            conf=conf,
            layers=layers,

            # Save this visual to a path containing this visual's name. Since
            # this name is guaranteed to be unique across all other visuals of
            # the same type, this path is guaranteed to be unique in at least
            # the directory containing all visuals of the same type.
            label=conf.name,
            **kwargs
        )

        # Display and/or save this animation.
        self._animate()

# ....................{ SUBCLASSES ~ obsolete              }....................
#FIXME: Replace use of all the following subclasses with the
#"AnimCellsAfterSolvingLayered" subclass *AFTER* refactoring all subclasses to
#leverage layers.
class AnimVelocity(AnimCellsAfterSolving):
    '''
    Abstract base class of all animations of a velocity flow plotted on the
    current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            *args,
            # Since this class already plots a streamplot, prevent the
            # superclass from plotting another streamplot as an overlay.
            is_current_overlayable=False,
            **kwargs
        )

        # Prefer an alternative colormap.
        self._colormap = self._phase.p.background_cm
