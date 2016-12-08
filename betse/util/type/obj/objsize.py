#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level wrappers around Pympler, an optional runtime dependency performing
memory profiling of Python applications.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.util.type import types
from betse.util.type.types import type_check, NumericTypes

# ....................{ CLASSES                            }....................
class SizeProfilableABC(object, metaclass=ABCMeta):
    '''
    Abstract base class signifying the subclass implementing this class to
    define a custom **size profile** (i.e., human-readable string synopsizing
    memory consumption by instances of this class).

    Classes _not_ implementing this class are provided a default size profile
    courtesy the general-purpose
    :func:`betse.lib.pympler.pymplers.print_vars_custom_size` function.
    '''

    # ..................{ SUBCLASS                           }..................
    @abstractmethod
    def get_size_profile(self, *args, **kwargs) -> str:
        '''
        Human-readable string synopsizing this object's memory consumption,
        accepting all optional arguments accepted by the :func:`get_sizes_vars`
        function.

        Parameters
        ----------
        All optional arguments accepted by the :func:`get_sizes_vars` function.

        Raises
        ----------
        BetseLibException
            If the optional runtime dependency Pympler is unsatisfied.
        '''

        pass

# ....................{ GETTERS                            }....................
@type_check
def get_size(obj: object, size_divisor: NumericTypes = 1) -> NumericTypes:
    '''
    Recursive in-memory size of the passed object, calculated with the external
    :mod:`pympler` package if importable _or_ raising an exception otherwise.

    This size is formally defined as the sum of:

    * The flat size of this object, defined as either:
      * If the type of this object provides the `__sizeof__()` special method,
        the return value of the :func:`sys.getsizeof` function passed this
        object.
      * Else, the sum of:
        * The basic size of this object, defined as the value of the
          undocumented `__basicsize__` attribute of the type of this object.
          Although undocumented, this appears to be the size of this object when
          **empty** (i.e., when this object contains no elements). For
          garbage-collected objects, this size "...includes the overhead for
          Python’s garbage collector (GC) as well as the space needed for
          refcounts (used only in certain Python builds)."
        * The item size of this object, defined as the value of the
          undocumented `__itemsize__` attribute of the type of this object.
          Although undocumented, this size appears to be:
          * If this object is a container, the sum of the sizes of the
            references referring to all objects contained by this object but
            _not_ the sizes of these objects themselves.
          * If this object is _not_ a container, 0.
    * The recursive size of each **referent** (i.e., object referred to by an
      instance variable) of this object, visited recursively up to some finite
      recursion depth. The size of a referent:
      * Visitable only by recursion exceeding this depth is silently ignored.
      * Visited multiple times is incorporated only once. Each subsequent
        visitation of this referent is silently ignored.

    Caveats
    ----------
    **This function requires the external :mod:`pympler` package,** an optional
    runtime dependency of this application. As the internal complexity of the
    :mod:`pympler.asizeof` submodule suggests, recursively calculating object
    sizes in a cross-platform and -interpreter portable manner generically
    applicable to _all_ possible objects is highly non-trivial. Since the only
    function defined by the stdlib to calculate object sizes (i.e.,
    :func:`sys.getsizeof`) behaves non-recursively, recursively calculating such
    sizes would effectively require wrapping that low-level function with a
    brute-force depth- or breadth-first search (DFS).

    Doing so would effectively reimplement the majority of the
    :mod:`pympler.asizeof` submodule, a lamentable time waste violating the
    Don't Repeat Yourself (DRY) principle in the best case and a futile failure
    in the worst. To preserve sanity, :mod:`pympler` is imported instead.

    Parameters
    ----------
    obj : object
        Object to iterate the variable sizes of.
    size_divisor: optional[NumericTypes]
        Number (i.e., integer, float) by which to divide the size in bytes of
        each variable of this object, hence converting from sizes in bytes to
        sizes in a larger denomination (e.g., :data:`betse.util.type.ints.KiB`,
        yielding sizes in kibibytes). Defaults to 1, yielding sizes in bytes.

    Returns
    ----------
    int
        Recursive in-memory size of this object in bytes, divided by the passed
        divisor if any.

    Raises
    ----------
    BetseLibException
        If the optional runtime dependency Pympler is unsatisfied.
    '''

    # Avoid circular import dependencies.
    from betse.lib import libs

    # If Pympler is unavailable, raise an exception.
    libs.die_unless_runtime_optional('Pympler')

    # Import Pympler *AFTER* validating its availability.
    from pympler import asizeof

    # Return the in-memory recursive size of this object, converted from bytes
    # to another denomination by the passed divisor.
    return asizeof.asizeof(obj) / size_divisor


#FIXME: Call this method as a new functional test detecting erroneously large
#attributes in these three principal objects. To do so, determine the largest
#attribute currently produced in these objects for the minified world
#environment used for all functional tests. Then raise an exception if any
#attribute of these objects is greater than or equal to 10 times this
#expected maximum.

@type_check
def get_sizes_vars(obj: object, size_divisor: NumericTypes = 1) -> tuple:
    '''
    2-tuple `(obj_size, vars_custom_size)` of this object's total recurusive
    size _and_ sequence of 2-tuples of the name and total recursive size of each
    **non-builtin variable** (i.e., variable whose name is _not_ both prefixed
    and suffixed by `__`) bound to the passed object (_in descending order of
    size_), calculated with the external :mod:`pympler` package if importable
    _or_ raising an exception otherwise.

    Parameters
    ----------
    obj : object
        Object to iterate the variable sizes of.
    size_divisor: optional[NumericTypes]
        Number (i.e., integer, float) by which to divide the size in bytes of
        each variable of this object, hence converting from sizes in bytes to
        sizes in a larger denomination (e.g., :data:`betse.util.type.ints.KiB`,
        yielding sizes in kibibytes). Defaults to 1, yielding sizes in bytes.

    Returns
    ----------
    (int, SequenceTypes)
        2-tuple `(obj_size, vars_custom_size)`, where
        * `obj_size` is this object's **total recursive size** (i.e., summation
          of the total non-recursive sizes of all non-builtin variables bound to
          this object).
        * `vars_custom_size` is a sequence of 2-tuples `(var_name, var_size)` of
          the name and recursive in-memory size of each variable bound to this
          object (_in descending order of size_).

    Raises
    ----------
    BetseLibException
        If the optional runtime dependency Pympler is unsatisfied.
    Exception
        If any variable of this object is a property whose
        :func:`property`-decorated method raises an exception.

    See Also
    ----------
    :func:`get_size`
        Further details on object size calculation.
    :func:`iter_vars_custom`
        Further details on non-builtin object variables.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import iterables
    from betse.util.type.obj import objs

    # List of all 2-tuples to be returned.
    vars_size = []

    # For the name and value of each non-builtin variable of this object (in
    # arbitrary order)...
    for var_name, var_value in objs.iter_vars_custom(obj):
        # Recursive in-memory size of this object.
        var_size = get_size(obj=var_value, size_divisor=size_divisor)

        # Append a 2-tuple of this variable's name and size.
        vars_size.append((var_name, var_size))

    # Total recursive size of this object, totalized over the second element of
    # each 2-tuple in this list.
    obj_size = iterables.sum_by_index(iterable=vars_size, subiterable_index=1)

    # This list sorted on the second element of each such 2-tuple (i.e., each
    # variable's size) in descending order.
    vars_size_sorted = iterables.sort_by_index_descending(
        iterable=vars_size, subiterable_index=1)

    # Return the expected 2-tuple.
    return obj_size, vars_size_sorted

# ....................{ GETTERS ~ profile                  }....................
#FIXME: Document the "vars_depth" parameter.
def get_size_profile(
    obj: object,
    *args,
    line_indent: str = '',
    vars_depth: int = 1,
    vars_max: int = None,
    **kwargs
) -> str:
    '''
    Human-readable synopsis of the size of the passed object with all remaining
    arguments passed as is to the :func:`get_sizes_vars` function, calculated
    with the external :mod:`pympler` package if importable _or_ raising an
    exception otherwise

    Specifically, this function returns a string tabulating:

    * The total recursive in-memory size in mebibytes of this object.
    * The names and total recursive in-memory sizes in mebibytes of all
      **non-builtin variables** (i.e., variables whose names are _not_ both
      prefixed and suffixed by `__`) bound to this object (_in descending order
      of size_).

    Motivation
    ----------
    Although Pympler _does_ define numerous classes and callables for profiling
    arbitrary objects, the resulting output is non-human-readable in the best
    case and ambiguous in the worst case. In either case, such output is mostly
    unusable. For example:

    >>> from betse.util.type.obj import objsize
    >>> from pympler.asizeof import asizesof
    >>> class ObjectLessons(object):
    ...     def __init__(self):
    ...         self.lessons = [
    ...             'Easter Island', 'St. Matthew Island', 'Greenland Norse',]
    ...         self.lesson_count = len(self.lessons)
    >>> object_lessons = ObjectLessons()
    >>> asizesof(object_lessons, stats=2.0)
     592 bytes:  <__main__.ObjectLessons object at 0x7fe937595208>
     592 bytes
       8 byte aligned
       8 byte sizeof(void*)
       1 object given
      10 objects sized
      10 objects seen
       3 recursion depth

       5 profiles:  total (% of grand total), average, and largest flat size:  largest object
       5 class str objects:  320 (54%), 64, 72:  'St. Matthew Island' leng 19!
       1 class dict object:  96 (16%), 96, 96:  {'lesson_count': 3, 'lessons': ['Easte..... Matthew Island', 'Greenland Norse']} leng 0
       1 class list object:  88 (15%), 88, 88:  ['Easter Island', 'St. Matthew Island', 'Greenland Norse'] leng 7!
       1 class __main__.ObjectLessons object:  56 (9%), 56, 56:  <__main__.ObjectLessons object at 0x7fe937595208>
       1 class int object:  32 (5%), 32, 32:  3 leng 1!

    Critically, note that:

    * The total size of the `ObjectLessons.lessons` list is actually 288 bytes
      rather than the significantly smaller 88 bytes reported above.
      Specifically:
      * The flat size of this list is 88 bytes, which Pympler erroneously
        reports to be the total size of this list.
      * The total flat size of all strings contained in this list is 200 bytes.
        Ergo, the total size of this list is 288 bytes -- as verified by
        explictly passing this list to the
        :func:`pympler.asizeof.asizeof` function, which then returns 288.
    * The total size of each variable of the passed `ObjectLessons` instance is
      _not_ reported. Only the total flat size of each type of object
      recursively visitable from this instance is reported, which is useless.

    Until Pympler addresses both concerns, this function profiles the memory
    consumption of arbitrary objects in a human-readable manner.

    Parameters
    ----------
    obj : object
        Object to synopsize the size of.
    vars_depth: optional[int]
        Defaults to 1, ...
    line_indent: optional[str]
        String to prefix each line of the returned string by, typically
        consisting of zero or more space characters. Defaults to the empty
        string.
    vars_max : optional[int]
        Maximum number of the largest non-builtin variables bound to this object
        to synopsize the sizes of. Defaults to `None`, in which case this
        function unconditionally synopsizes the sizes of _all_ non-builtin
        variables bound to this object (regardless of size).

    All remaining arguments are passed as is to the :func:`get_sizes_vars`
    function.

    Returns
    ----------
    str
        Human-readable synopsis of the size of this object.

    Raises
    ----------
    BetseLibException
        If the optional runtime dependency Pympler is unsatisfied.
    Exception
        If any variable of this object is a property whose
        :func:`property`-decorated method raises an exception.

    See Also
    ----------
    :func:`iter_vars_custom_size`
        Further details on object size calculation.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objs
    from betse.util.type.ints import MiB

    # If the caller failed to specify a size divisor, default to a divisor
    # denominating sizes in mebibytes for readability.
    kwargs.setdefault('size_divisor', MiB)

    # If either this is the bottom-most recursive call in the current subtree of
    # recursive calls *OR* this object is C-based (i.e., is either builtin or
    # defined by a C extension), halt the recursion by returning a single-line
    # profile rather than recursively introspecting this object.
    #
    # If this object is C-based, recursively introspecting this object would be
    # trivially feasible. Python code makes no distinction between pure-Python
    # and C-based objects. Indeed, differentiating between the two with Python
    # code is non-trivial! That said, recursively introspecting C-based objects
    # would be largely pointless; size profiles are principally intended to
    # profile the space consumed by variables bound to application-specific
    # pure-Python objects rather than those bound to external C-based objects.
    if vars_depth == 0 or objs.is_c_based(obj):
        # Total recursive size of only this object.
        obj_size = get_size(obj, *args, **kwargs)

        #FIXME: Generalize the format specifier "{:.04f} MiB" into a new local
        #constant declared above. DRY in all things.

        # Return a single-line profile.
        return '{:.04f} MiB'.format(obj_size)

    # Else, recursively introspect this object.
    #
    # Human-readable synopsis of the size of this object to be returned.
    size_profile = ''

    # Human-readable string to be appended to this synopsis on being returned.
    size_profile_suffix = ''

    # Total recursive size of this object *AND* a sequence of 2-tuples of the
    # name and recursive in-memory size of each non-builtin variable of this
    # object (in descending size order).
    obj_size, vars_name_size = get_sizes_vars(obj, *args, **kwargs)

    #FIXME: Convert the current "kwargs['size_divisor']" into a human-readable
    #denomenation (e.g., "GiB", "KiB", "B"). Since this is both non-trivial and
    #currently *NOT* required, we adopt a superior lazy approach.

    # Append this size to this profile.
    size_profile += '{}total: {:.04f} MiB'.format(line_indent, obj_size)

    # Increment the indentation level of object variables totalized below.
    line_indent += '  '

    # If a maximum number of the largest variables is requested *AND* the number
    # of variables bound to this object exceeds this maximum, ignore all
    # variables of smaller size. To ensure the synopsis suffix set below is of
    # the proper indentation level, this is done *AFTER* increasing this level.
    if vars_max is not None and len(vars_name_size) > vars_max:
        # Validate this number to be a positive integer.
        assert types.is_int_positive(vars_max), (
            types.assert_not_int_positive(vars_max))

        # Truncate this list to this maximum number of the first such 2-tuples.
        vars_name_size = vars_name_size[:vars_max]

        # Append an ellipses to this synopsis, notifying users that additional
        # variables exist but were truncated.
        size_profile_suffix = '\n{}...'.format(line_indent)

    #FIXME: Consider left-aligning the sizes printed by this iteration.
    #FIXME: Consider also printing the truncation of the value of each variable.

    # For the name and size of each of this object's non-builtin variables...
    for var_name, var_size in vars_name_size:
        # Append this size to this profile.
        size_profile += '\n{}{}: {:.04f} MiB'.format(
            line_indent, var_name, var_size)

    # Return this profile, suffixed by this suffix if any.
    return size_profile + size_profile_suffix
