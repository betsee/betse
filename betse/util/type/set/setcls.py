#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **set classes** (i.e., classes implementing set-like functionality,
typically by subclassing the builtin :class:`set` or :class:`frozenset`
container types or analogues thereof).
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, ClassType, IterableTypes, MappingType,)
from functools import wraps

# ....................{ GLOBALS                            }....................
_FROZENSET_METACLASS = type(frozenset)
'''
Metaclass of the builtin "frozenset" container type.
'''

# ....................{ METACLASSES                        }....................
class FrozenSetSubclassableMeta(ABCMeta, _FROZENSET_METACLASS):
    '''
    Metaclass of the abstract :class:`FrozenSetSubclassable` base class and all
    concrete subclasses thereof.

    This metaclass dynamically redefines *all* container-creating methods of the
    :class:`frozenset` superclass of the :class:`FrozenSetSubclassable` class
    (e.g., :meth:`frozenset.__or__`) within the currently declared concrete
    subclass of that class. Specifically, this metaclass redefines each such
    method to return an instance of the currently declared concrete subclass of
    :class:`FrozenSetSubclassable` rather than of the :class:`frozenset`
    superclass, preserving sane semantics and caller expectations.

    Design
    ----------
    The :class:`FrozenSetSubclassable` class is merely a placeholder subclass of
    the :class:`frozenset` type whose metaclass is this metaclass. Ideally, the
    work performed by this metaclass would directly reside in the
    :class:`FrozenSetSubclassable` class instead, in which case this metaclass
    would have *no* demonstrable reason to exist. However, this work requires
    access to the concrete subclass of the :class:`FrozenSetSubclassable` class
    being currently declared. Since each such subclass is accessible *only* from
    within the metaclass of the :class:`FrozenSetSubclassable` class rather than
    within that class itself, this work necessarily resides in this metaclass.

    For general-purpose usability, this metaclass subclasses both:

    * The metaclass of the :class:`frozenset` type (typically, the root
      metaclass :class:`type`), avoiding conflicts between this metaclass and
      the :class:`FrozenSetSubclassable` class subclassing the
      :class:`frozenset` type.
    * The :class:`ABCMeta` metaclass, avoiding conflicts between this metaclass
      and concrete subclasses of the :class:`FrozenSetSubclassable` class which
      additionally subclass one or more other abstract base classes which
      themselves leverage the :class:`ABCMeta` metaclass. (Subtle dragons.)

    For these and similar reasons, metaclass usage in Python should typically be
    kept to a minimum. The :class:`FrozenSetSubclassable` class violates this
    maxim because it absolutely must; in all other cases, alternate solutions
    *not* leveraging metaclasses should be implemented instead.

    See Also
    ----------
    https://stackoverflow.com/a/804973/2809027
        StackOverflow answer mildly inspiring this class.
    '''

    # ..................{ CONSTRUCTORS                       }..................
    def __new__(
        metacls,
        class_name: str,
        class_base_classes: IterableTypes,
        class_attrs: MappingType,
        **kwargs
    ) -> ClassType:
        '''
        Redefine all container-creating methods of the :class:`frozenset`
        superclass of the :class:`FrozenSetSubclassable` class (e.g.,
        :meth:`frozenset.__or__`) within the currently declared concrete
        subclass of that class identified by the passed parameters.
        '''

        # Tuple of the unqualified names of all container-creating methods defined
        # by the "frozenset" type, requiring redefinition in the
        # "FrozenSetSubclassable" subclass declared above.
        CREATION_METHOD_NAMES = (
            'copy',
            'difference',
            'intersection',
            'symmetric_difference',
            'union',
            '__and__',
            '__or__',
            '__rand__',
            '__ror__',
            '__rsub__',
            '__rxor__',
            '__sub__',
            '__xor__',
        )

        # Unsanitized "FrozenSetSubclassable" subclass.
        frozenset_subclass = super().__new__(
            metacls, class_name, class_base_classes, class_attrs)

        # For the name of each such method...
        for creation_method_name in CREATION_METHOD_NAMES:
            metacls._sanitize_creation_method(
                frozenset_subclass, creation_method_name)

        # Return this sanitized "FrozenSetSubclassable" subclass.
        return frozenset_subclass

    # ..................{ SANITIZERS                         }..................
    @staticmethod
    @type_check
    def _sanitize_creation_method(
        frozenset_subclass: ClassType, method_name: str) -> None:
        '''
        Redefine the container-creating method of the
        :class:`frozenset` superclass with the passed name in a cleverly
        automated manner circumventing all superclass issues.

        Design
        ----------
        This private static method is intended to be called only by the special
        static :meth:`__new__` method of this metaclass.

        This method is intentionally implemented as a discrete callable rather
        than inlined directly into the body of the :meth:`__new__` method, as
        the closure internally defined by this method expects the local
        ``frozenset_method`` variable captured by this closure to remain
        constant for the lifetime of that closure.

        Parameters
        ----------
        frozenset_subclass : ClassType
            :class:`FrozenSetSubclassable` subclass to redefine this method for.
        method_name : str
            Name of the container-creating method to be redefined.
        '''

        # Container-creating method with this name defined by "frozenset".
        frozenset_method = objects.get_method(
            obj=frozenset, method_name=method_name)

        # Closure sanitizing this method to return instances of the concrete
        # subclass inheriting from this class rather than of "frozenset",
        # wrapped in a manner propagating the name, docstring, and other
        # identifying metadata of the original method to this closure.
        @wraps(frozenset_method)
        def sanitized_creation_method(self, *args, **kwargs):

            # "frozenset" instance created by calling the superclass method with
            # all passed positional and keyword arguments.
            set_created = frozenset_method(self, *args, **kwargs)

            # Return either...
            return (
                # A new instance of this concrete subclass if the returned
                # object is in fact a "frozenset" instance.
                frozenset_subclass(set_created) if isinstance(
                    set_created, frozenset) else
                # This object as is otherwise. Technically, this should probably
                # never occur. Pragmatically, you know what they say about every
                # bad assumption we've ever made.
                set_created)

        # Override the superclass method with this closure.
        setattr(frozenset_subclass, method_name, sanitized_creation_method)

# ....................{ SUPERCLASSES                       }....................
class FrozenSetSubclassable(
    frozenset, metaclass=FrozenSetSubclassableMeta):
    '''
    Safely subclassable immutable set.

    Caveats
    ----------
    **Immutable set subclasses should always inherit from this base class.**
    Neither the builtin :class:`frozenset` container type nor the abstract
    :class:`collections.abc.Set` and :class:`collections.abc.Hashable` base
    classes should be inherited from.

    Subclasses should avoid declaring a metaclass or inheriting from another
    base class that declares a metaclass - excluding the standard
    :class:`ABCMeta` metaclass, which is compatible with this class by design.
    All other metaclasses should be considered incompatible. Violating this
    constraint typically raises the following runtime exception:

        TypeError: metaclass conflict: the metaclass of a derived class must be
        a (non-strict) subclass of the metaclasses of all its bases

    Subclasses must redefine the static ``__new__()`` method and *not* attempt
    to define the ``__init__()`` method. Immutable types are necessarily
    initialized at object creation time.

    By Python design, the ``__new__()`` method is static rather than a
    classmethod and hence *must* be manually redefined in *all* subclasses. This
    redefinition *must* internally call the superclass :func:`frozenset.__new__`
    method and return the result of doing so (i.e., an instance of the desired
    subclass type). This redefinition should ideally (but *not* necessarily)
    share a similar signature as the :func:`frozenset.__new__` method and pass
    all passed parameters to that method as is. In positional order, these are:

    #. **The type of the current subclass.** Since the ``__new__()`` method is
       static and hence *not* bound to the type of the current subclass, this
       type *must* be manually passed as the first argument to this method.
    #. **The optional iterable defining the contents of this immutable set.**
       If unpassed, this set reduces to the empty set.

    The trivial implementation of the ``__new__()`` method is as follows:

        def __new__(cls, *args):
            return super().__new__(cls, *args)

    Versus :class:`frozenset`
    ----------
    The :class:`frozenset` type is *not* safely subclassable, due to unfortunate
    design decisions baked into the C-based implementations of *all* builtin
    container types. The core issue pertains to the objects returned by
    **container-creating methods** (i.e., methods declared by these types that
    create and return instances of the same types). In the case of
    :class:`frozenset`, these methods include:

    * Binary set operations (e.g., the set union operator `|`, internally
      implemented by the :meth:`frozenset.__or__` special method).
    * Binary set methods (e.g., the set union method :meth:`frozenset.union`).
    * Copy set methods (e.g., the set copying method :meth:`frozenset.copy`).

    The specifics of the unsuitability of :class:`frozenset` depend on the major
    version of Python in use. Specifically, under:

    * Python 2.x, container-creating methods in both the :class:`frozenset` and
      :class:`set` types correctly created instances of subclasses inheriting
      from these types but incorrectly failed to call the `__init__` methods of
      these subclasses. This is referred to as `issue #1721812`_.

    .. _issue #1721812:
        https://mail.python.org/pipermail/python-bugs-list/2007-May/038471.html

    * Python 3.x, container-creating methods in both the :class:`frozenset` and
      :class:`set` types incorrectly resolved `issue #1721812`_ by creating
      instances of the corresponding base types (e.g., :class:`frozenset` or
      :class:`set`) in subclasses inheriting from these types rather than
      instances of these subclasses.

    Technically, the latter issue may be resolved by manually redefining *all*
    container-creating methods in :class:`frozenset` subclasses. Pragmatically,
    there exist at least 17 such methods, rendering such redefinition
    effectively infeasible: ``__ror__``, ``difference_update``, ``__isub__``,
    ``symmetric_difference``, ``__rsub__``, ``__and__``, ``__rand__``,
    ``intersection``, ``difference``, ``__iand__``, ``union``, ``__ixor__``,
    ``symmetric_difference_update``, ``__or__``, ``copy``, ``__rxor__``,
    ``intersection_update``, ``__xor__``, ``__ior__``, and ``__sub__``.

    Versus ``Set`` and ``Hashable``
    ----------
    Technically, the abstract :class:`collections.abc.Set` base class defining
    the official immutable set API *is* safely subclassable, but only under the
    following stipulations:

    * The abstract :class:`collections.abc.Hashable` base class should typically
      also be subclassed. Failing to do so raises exceptions on attempting to
      add subclass instances to builtin container types expecting hashable
      objects (e.g., :class:`dict`, :class:`set`).
    * All public methods defined by the :class:`frozenset` type (e.g.,
      :meth:`frozenset.union`, :meth:`frozenset.issubset`) should typically also
      be defined, both for usability *and* duck typing purposes.

    Satisfying these stipulations requires defining the following 13 methods:
    ``__contains__``, ``__eq__``, ``__hash__``, ``__init__``, ``__iter__``,
    ``__len__``, ``difference``, ``intersection``, ``symetric_difference``,
    ``union``, ``copy``, ``issubset``, and ``issuperset``.

    While feasible, doing so is surprisingly less trivial than redefining the
    17 container-creating :class:`frozenset` methods listed above. The latter
    all share the exact same semantics and hence are trivially redefined with
    automation, as this class demonstrates. The former, however, share *no*
    common semantics and hence each require manual redefinition.

    Moreover, the typical implementation of a :class:`collections.abc.Set`
    subclass encapsulates an internal :class:`frozenset` instance variable.
    Since both approaches require the :class:`frozenset` class *and* since
    directly subclassing that class is simpler than indirectly encapsulating
    that class in :class:`collections.abc.Set` subclasses, this class elects to
    directly subclass the :class:`frozenset` type instead.
    '''

    # ..................{ CONSTRUCTORS                       }..................
    # The entirety of the logic for this class resides in our metaclass..
    def __new__(cls, *args):
        return super().__new__(cls, *args)
