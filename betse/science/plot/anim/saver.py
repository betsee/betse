#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Animation serialization classes.
'''

#FIXME: *UGH.* We overengineered this one, folks. Replace this entire class
#hierarchy with a single class encapsulating all saving functionality pertaining
#to frames and videos. Why? Consider the current approach. How do callers save
#frames and video under the current approach? Due to the current partitioning of
#frames and video logic into separate classes and hence objects, callers would
#presumably need to call something like:
#
#    self._anim_saver_frames.save()
#    self._anim_saver_video.save()
#
#Frankly, that's horrible. Such low-level details as whether or not frames
#and/or video saving is currently enabled should be encapsulated away from
#callers into a single high-level class.
#
#Indeed, let's go one step further and generalize this class to encapsulate
#*ALL* general-purpose animation functionality -- both saving and displaying.
#Hence, name this single class "AnimConfig".

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractstaticmethod
from betse.util.type import types
from betse.util.type.types import type_check, Sequence

# ....................{ SUPERCLASS                         }....................
#FIXME: Rename to "AnimSaver".
class AnimationSaver(object, metaclass=ABCMeta):
    '''
    Abstract base class of all animation serialization classes.

    Instances of this class serialize (i.e., save) in-memory simulation
    animations to on-disk cache, image, and/or video files.

    Attributes
    ----------
    is_enabled : bool
        `True` if this subclass-specific animation serialization is enabled or
        `False` otherwise.
    filetype : str
        Filetype of all files serialized by this subclass.
    '''

    # ..................{ ABSTRACT ~ static                  }..................
    @abstractstaticmethod
    def make(params: 'Parameters') -> "AnimationSaver":
        '''
        Factory method producing a concrete instance of this abstract base class
        from the passed simulation configuration.

        Parameters
        ----------------------------
        params : Parameters
            Current simulation configuration.

        Returns
        ----------------------------
        AnimationSaver
            Concrete instance of this abstract base class.
        '''

        pass

    # ..................{ CONCRETE ~ public                  }..................
    @type_check
    def __init__(self, is_enabled: bool, filetype: str) -> None:

        # Classify the passed parameters.
        self.is_enabled = is_enabled
        self.filetype = filetype

# ....................{ IMAGE                              }....................
#FIXME: Rename to "AnimSaverFrames".
class AnimationSaverFrames(AnimationSaver):
    '''
    Animation serialization class serializing each frame of an animation to an
    image file of some predefined filetype.

    Attributes
    ----------
    dpi : int
        Dots per inch (DPI) of all output image files.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @abstractstaticmethod
    def make(params: 'Parameters') -> "AnimationSaverFrames":
        assert types.is_parameters(params), types.assert_not_parameters(params)

        #FIXME: Non-ideal. Refactor in accordance with the corresponding
        #"sim_config.yaml" FIXME comment.

        saf = (
            params.config['results options']['save animations']['image frames'])
        return AnimationSaverFrames(
            is_enabled=saf,
            filetype='png',
            dpi=300,
        )

    # ..................{ PUBLIC                             }..................
    @type_check
    def __init__(self, is_enabled: bool, filetype: str, dpi: int) -> None:

        super().__init__(is_enabled, filetype)
        self.dpi = dpi

# ....................{ VIDEO                              }....................
#FIXME: Rename to "AnimSaverVideo".
class AnimationSaverVideo(AnimationSaver):
    '''
    Anim serialization class serializing an entire animation to an encoded
    video file of some predefined filetype.

    Attributes
    ----------
    encoder_names : list
        List of the Matplotlib-specific names of all supported video encoders
        (in order of descending preference). Unavailable video encoders will
        simply be ignored.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @abstractstaticmethod
    def make(params: 'Parameters') -> "AnimationSaverVideo":
        assert types.is_parameters(params), types.assert_not_parameters(params)

        sav = (
            params.config['results options']['save animations']['movie'])
        return AnimationSaverVideo(
            is_enabled=sav['enabled'],
            filetype=sav['filetype'],
            encoder_names=sav['encoders'],
        )

    # ..................{ PUBLIC                             }..................
    @type_check
    def __init__(
        self, is_enabled: bool, filetype: str, encoder_names: Sequence) -> None:

        super().__init__(is_enabled, filetype)
        self.encoder_names = encoder_names
