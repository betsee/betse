#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to tissue and cut
profiles.
'''

# ....................{ IMPORTS                            }....................
from betse.science.tissue.picker import TissuePicker
from betse.util.type import types

# ....................{ TISSUE                             }....................
class TissueProfile(object):
    '''
    Parameters associated with a subset of the cell population.

    Since these parameters are typically derived from real-world observation of
    biological tissue, the spacialization of these parameters collectively
    defines a "tissue" and is hence referred to as a **tissue profile**.

    Attributes
    ----------------------------
    name : str
        Unique profile name.
    picker : TissuePicker
        Object assigning a subset of the cell population to this profile.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @staticmethod
    def make(data: dict, params: 'Parameters') -> 'TissueProfile':
        '''
        Factory method producing an instance of this base class defined by the
        passed tissue simulation configuration.

        Parameters
        ----------------------------
        data : dict
             Dictionary describing the contents of the tissue profile to be
             created via the following key-value pairs:
             * `???`, a ???.
        params : Parameters
             Current tissue simulation configuration.

        Returns
        ----------------------------
        TissueProfile
            Instance of this class.
        '''
        assert types.is_parameters(params), types.assert_not_parameters(params)
        tpd = params.config['tissue profile definition']

        # Object to be returned, defaulting to nothing.
        profile = None

        #FIXME: Implement me!

        # If cut profiles are enabled, return an instance of this class.
        if tpd['profiles enabled']:
            profile = TissueProfile(
                name=data['name'],
                picker=TissuePicker.make(data['cell targets'], params),
            )

        return profile

    # ..................{ PUBLIC                             }..................
    def __init__(self, name: str, picker: 'TissuePicker') -> None:
        assert types.is_str(name), types.assert_not_str(name)
        assert types.is_tissue_picker(picker), \
            types.assert_not_tissue_picker(picker)

        self.name = name
        self.picker = picker

# ....................{ CUT                                }....................
#FIXME: Overkill identified. Specifically:
#
#* Fold this class into the new "EventCut" class elsewhere.
#* Fold the "Profile" class into the existing "TissueProfile" class above.
#* Eliminate both this and the "Profile" class.
class CutProfile(Profile):
    '''
    Profile identifying all cells to be permanently removed by a cutting
    event subsequently triggered during the current tissue simulation.
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
