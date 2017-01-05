#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to tissue and cut
profiles.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.exceptions import BetseSimConfigException
from betse.science.tissue.picker import (TissuePicker, TissuePickerBitmap)
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
    z_order : int
        1-based order in which this profile will be plotted with respect to all
        other plotted profiles. Specifically, profiles with larger z order will
        be plotted over profiles with smaller z order.
    picker : TissuePicker
        Object assigning a subset of the cell population to this profile.
    '''

    # ..................{ ABSTRACT ~ static                  }..................
    @staticmethod
    def make(
        profile_config: dict, params: 'Parameters', z_order: int) -> 'Profile':
        '''
        Factory method producing a concrete instance of this abstract base class
        from the passed tissue simulation configuration.

        Parameters
        ----------------------------
        profile_config : dict
            Dictionary describing the contents of the profile to be created.
        params : Parameters
            Current tissue simulation configuration.
        z_order : int
            1-based order in which this profile will be plotted with respect to
            all other plotted profiles.

        Returns
        ----------------------------
        Profile
            Concrete instance of this abstract base class.
        '''
        assert types.is_mapping(profile_config), \
            types.assert_not_mapping(profile_config)
        assert types.is_parameters(params), types.assert_not_parameters(params)
        assert types.is_int(z_order), types.assert_not_int(z_order)

        tpd = params.config['tissue profile definition']

        # Object to be returned, defaulting to nothing.
        profile = None

        # If cut profiles are enabled, return an instance of this class.
        if tpd['profiles enabled']:
            profile_type = profile_config['type']

            # If this is a tissue profile...
            if profile_type == 'tissue':
                #FIXME: Refactor to return the following class instead:
                # profile = TissueProfile(
                #     name=config['name'],
                #     picker=TissuePicker.make(config['cell targets'], params),
                # )

                profile = {
                    'type': profile_config['type'],
                    'name': profile_config['name'],
                    'z order': z_order,
                }

                profile['insular gj'] = profile_config['insular']
                profile['picker'] = TissuePicker.make(
                    profile_config['cell targets'], params)

                # For safety, coerce all diffusion constants to floats.
                profile['diffusion constants'] = {
                    key: float(value) for key, value in (
                        profile_config['diffusion constants'].items()) }

            # Else if this is a tissue profile...
            elif profile_type == 'cut':
                profile = ProfileCut(
                    name=profile_config['name'],
                    z_order=z_order,
                    picker=TissuePickerBitmap.make(
                        profile_config['bitmap'], params),
                )

            # Else, this profile is invalid. Raise an exception, matey!
            else:
                raise BetseSimConfigException(
                    'Profile type "{}"' 'unrecognized.'.format(profile_type))

        return profile

    # ..................{ CONCRETE                           }..................
    def __init__(self, name: str, z_order: int, picker: 'TissuePicker') -> None:
        '''
        Initialize this profile.

        Parameters
        ----------
        name : str
            Unique name of this profile.
        z_order : int
            1-based order in which this profile will be plotted with respect to
            all other plotted profiles. Specifically, profiles with larger z
            order will be plotted over profiles with smaller z order.
        picker : TissuePicker
            Object assigning a subset of the cell population to this profile.
        '''

        assert types.is_str(name), types.assert_not_str(name)
        assert types.is_int_ge(z_order, 1), types.assert_not_int_ge(z_order, 1)
        assert types.is_tissue_picker(picker), \
            types.assert_not_tissue_picker(picker)

        self.name = name
        self.z_order = z_order
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
class ProfileCut(Profile):
    '''
    Profile identifying all cells to be permanently removed by a cutting
    event subsequently triggered during the current tissue simulation.

    There exists a many-to-one relation between cut profiles and cutting events.
    That is to say, each cutting event references zero or more cut profiles.
    While unwieldy, this disconnection permits cutting event (but _not_ cut
    profile) parameters to be modified without requiring simulation reseeding or
    reinitialization. Welcome to the code jungle.
    '''
    pass
