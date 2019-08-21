#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed wrapper subclasses.
'''

# ....................{ IMPORTS                           }....................
from abc import ABCMeta
from betse.exceptions import BetseYamlException
# from betse.util.io.log import logs
from betse.util.type.types import type_check, MappingOrSequenceTypes, NoneType

# ....................{ SUPERCLASSES                      }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: To avoid diamond inheritance issues in the "YamlListABC" subclass,
# neither this subclass nor subclasses of this subclass should define any
# methods already defined by the "MutableSequence" API. This includes but is
# *NOT* limited to the following "MutableSequence" methods: append(), clear(),
# count(), extend(), index(), insert(), mro(), pop(), register(), remove(), and
# reverse().
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class YamlABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **YAML-backed wrapper** (i.e., high-level object
    wrapping a low-level mapping or sequence of YAML-backed configuration
    settings) subclasses.

    Attributes
    ----------
    _conf : MappingOrSequenceOrNoneTypes
        Low-level container of related configuration settings loaded from and
        saved back to a YAML-formatted configuration file if the :meth:`load`
        method has been called *or* ``None`` otherwise.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self) -> None:
        '''
        Initialize this YAML wrapper.
        '''

        self._conf = None

    # ..................{ PROPERTIES ~ read-only            }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def is_loaded(self) -> bool:
        '''
        ``True`` only if this configuration is in the **loaded state** (i.e.,
        the :meth:`load` method has been called more recently than the
        :meth:`unload` method has been called).

        If ``True``, the :meth:`conf` property is safely accessible by callers.
        '''

        return self._conf is not None


    @property
    def conf(self) -> MappingOrSequenceTypes:
        '''
        Low-level container of related configuration settings loaded from and
        saved back to a YAML-formatted configuration file if the :meth:`load`
        method has been called more recently than the :meth:`unload` method has
        been called *or* raise an exception otherwise.

        Raises
        ----------
        BetseYamlException
             If the :meth:`load` method has yet be called.
        '''

        # If this property has yet to be set, raise an exception.
        self.die_unless_loaded()

        # Else, this property is a container. Return this container.
        return self._conf

    # ..................{ EXCEPTIONS                        }..................
    def die_unless_loaded(self) -> None:
        '''
        Raise an exception unless this configuration is in the loaded state.

        Equivalently, this method raises an exception if this configuration is
        in the unloaded state.

        Raises
        ----------
        BetseYamlException
            If this configuration is in the unloaded state.

        See Also
        ----------
        :meth:`is_loaded`
            Further details.
        '''

        if not self.is_loaded:
            raise BetseYamlException(
                'YAML configuration not loaded '
                '(i.e., load() method not called yet).')

    # ..................{ LOADERS                           }..................
    @type_check
    def load(self, conf: MappingOrSequenceTypes) -> None:
        '''
        Associate this configuration with the passed YAML-backed container.

        Parameters
        ----------
        conf : MappingOrSequenceTypes
            Low-level container of related configuration settings loaded from
            and saved back to a YAML-formatted configuration file.
        '''

        # Classify this container.
        self._conf = conf


    def unload(self) -> None:
        '''
        Deassociate this configuration from its YAML-backed container if any
        *or* reduce to a noop otherwise.

        This method serves as a low-level safety mechanism ensuring that the
        next access of the :meth:`conf` property will raise a human-readable
        exception.
        '''

        self._conf = None

# ....................{ TYPES                             }....................
# Intended for use in callable type validation.
YamlABCOrNoneTypes = (YamlABC, NoneType)
'''
Tuple of both the YAML-backed configuration type *and* the type of the
singleton ``None`` object.
'''
