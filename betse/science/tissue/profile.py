#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to tissue and cut
profiles.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractstaticmethod
from betse.science.tissue.picker import TissuePicker
from betse.util.type import types

# ....................{ BASE                               }....................
class TissueProfileBase(object, metaclass = ABCMeta):
    '''
    Abstract base class of all tissue-centric profile classes.

    Instances of this class assign a subset of all cells matching
    subclass-specific criteria (e.g., explicit indexing, randomized selection,
    spatial location) to the corresponding tissue profile.

    Attributes
    ----------------------------
    name : str
        Unique name of the current profile.
    picker : TissuePicker
        Object assigning a subset of the cell population to this profile.
    '''

    # ..................{ ABSTRACT ~ static                  }..................
    @abstractstaticmethod
    def make(params: 'Parameters') -> 'TissueProfileBase':
        '''
        Factory method producing a concrete instance of this abstract base class
        from the passed tissue simulation configuration.

        Parameters
        ----------------------------
        params : Parameters
             Current tissue simulation configuration.

        Returns
        ----------------------------
        TissueProfileBase
            Concrete instance of this abstract base class.
        '''
        pass

    # ..................{ CONCRETE                           }..................
    def __init__(self, name: str, picker: 'TissuePicker') -> None:
        assert types.is_str(name), types.assert_not_str(name)
        assert types.is_tissue_picker(picker), \
            types.assert_not_tissue_picker(picker)

        self.name = name
        self.picker = picker

# ....................{ CUT                                }....................
class TissueProfile(TissueProfileBase):
    '''
    Profile identifying all cells to be initialized with the same parameters.
    '''

    @staticmethod
    def make(params: 'Parameters') -> 'TissueProfile':
        assert types.is_parameters(params), types.assert_not_parameters(params)
        tpd = params.config['tissue profile definition']

        #FIXME: Implement me!
        return None

# ....................{ CUT                                }....................
class CutProfile(TissueProfileBase):
    '''
    Profile identifying all cells to be permanently removed by a cutting
    event subsequently triggered during the tissue simulation.
    '''

    @staticmethod
    def make(params: 'Parameters') -> 'CutProfile':
        assert types.is_parameters(params), types.assert_not_parameters(params)
        tpd = params.config['tissue profile definition']

        # If cut profiles are enabled, return an instance of this class.
        if not tpd['profiles enabled']:
            cp = tpd['cut profile']
            picker = TissuePicker.make(cp['cell targets'], params)
            return CutProfile(cp['name'], picker)
        # Else, return the empty void of space.
        else:
            return None
