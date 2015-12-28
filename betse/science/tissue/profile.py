#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to tissue and cut
profiles.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractstaticmethod
from betse.exceptions import BetseExceptionParameters
from betse.science.tissue.picker import TissuePicker
from betse.util.type import types

# ....................{ BASE                               }....................
class Profile(object, metaclass = ABCMeta):
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
    def make(config: dict, params: 'Parameters') -> 'Profile':
        '''
        Factory method producing a concrete instance of this abstract base class
        from the passed tissue simulation configuration.

        Parameters
        ----------------------------
        config : dict
             Dictionary describing the contents of the profile to be created.
        params : Parameters
             Current tissue simulation configuration.

        Returns
        ----------------------------
        Profile
            Concrete instance of this abstract base class.
        '''
        assert types.is_mapping(config), types.assert_not_mapping(config)
        assert types.is_parameters(params), types.assert_not_parameters(params)

        tpd = params.config['tissue profile definition']

        # Object to be returned, defaulting to nothing.
        profile = None

        # If cut profiles are enabled, return an instance of this class.
        if tpd['profiles enabled']:
            profile_type = config['type']

            if profile_type == 'tissue':
                profile = TissueProfile(
                    name=config['name'],
                    picker=TissuePicker.make(config['cell targets'], params),
                    #FIXME: Finish us up the configuration bomb!
                )
            elif profile_type == 'cut':
                profile = CutProfile(
                    name=config['name'],
                    picker=TissuePicker.make(config['cell targets'], params),
                )
            else:
                raise BetseExceptionParameters(
                    'Profile type "{}"' 'unrecognized.'.format(profile_type))

        return profile

    # ..................{ CONCRETE                           }..................
    def __init__(self, name: str, picker: 'TissuePicker') -> None:
        assert types.is_str(name), types.assert_not_str(name)
        assert types.is_tissue_picker(picker), \
            types.assert_not_tissue_picker(picker)

        self.name = name
        self.picker = picker

# ....................{ CUT                                }....................
#FIXME: Implement us!
class TissueProfile(Profile):
    '''
    Parameters associated with a subset of the cell population.

    Since these parameters are typically derived from real-world observation of
    biological tissue, the spacialization of these parameters collectively
    defines a "tissue" and is hence referred to as a **tissue profile**.
    '''
    pass

# ....................{ CUT                                }....................
#FIXME: Implement us!
class CutProfile(Profile):
    '''
    Profile identifying all cells to be permanently removed by a cutting
    event subsequently triggered during the current tissue simulation.
    '''
    pass
