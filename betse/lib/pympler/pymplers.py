#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level wrappers around Pympler, an optional runtime dependency performing
memory profiling of Python applications.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.objects import iter_vars_custom
from betse.util.type.types import type_check, NumericTypes, SequenceTypes

# ....................{ GETTERS                            }....................
@type_check
def get_object_size(obj: object, divisor: NumericTypes = 1) -> NumericTypes:
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
    divisor: optional[NumericTypes]
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
    return asizeof.asizeof(obj) / divisor

# ....................{ PRINTERS                           }....................
#FIXME: Call this method during profiling with our three principal objects:
#"Simulator", "Cells", and "Parameters".

def print_object_vars_custom_size(obj: object) -> None:
    '''
    Print a human-readable tabulation of the total recursive in-memory size in
    mebibytes of the passed object _and_ name and recursive in-memory size in
    mebibytes of each **non-builtin variable** (i.e., variable whose name is
    _not_ both prefixed and suffixed by `__`) bound to this object (_in
    descending order of size_), calculated with the external :mod:`pympler`
    package if importable _or_ raising an exception otherwise.

    Parameters
    ----------
    obj : object
        Object to print the variable sizes of.

    Raises
    ----------
    BetseLibException
        If the optional runtime dependency Pympler is unsatisfied.
    Exception
        If any variable of this object is a property whose
        :func:`property`-decorated method raises an exception.

    See Also
    ----------
    :func:`iter_object_vars_custom_size`
        Further details on object size calculation.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import iterables
    from betse.util.type.ints import MiB

    # Sequence of 2-tuples of the name and recursive in-memory size in mebibytes
    # of each non-builtin variable of to this object (in descending size order).
    vars_custom_size = iter_object_vars_custom_size(obj=obj, divisor=MiB)

    # Recursive in-memory size of this object, defined as the summation of the
    # sizes of these variables given by the second element of each 2-tuple.
    obj_size = iterables.sum_by_index(
        iterable=vars_custom_size, subiterable_index=1)

    # Print this object's size.
    print('total: {:.02f} MiB'.format(obj_size))

    #FIXME: Consider left-aligning the sizes printed by this iteration.
    #FIXME: Consider also printing the truncation of the value of each variable.

    # Print the name and size of this object's variables.
    for var_name, var_size in vars_custom_size:
        print('  {}: {:.02f} MiB'.format(var_name, var_size))

# ....................{ ITERATORS                          }....................
#FIXME: Call this method as a new functional test detecting erroneously large
#attributes in these three principal objects. To do so, determine the largest
#attribute currently produced in these objects for the minified world
#environment used for all functional tests. Then raise an exception if any
#attribute of these objects is greater than or equal to 10 times this
#expected maximum.

@type_check
def iter_object_vars_custom_size(
    obj: object, divisor: NumericTypes = 1) -> SequenceTypes:
    '''
    Sequence of 2-tuples of the name and recursive in-memory size of each
    **non-builtin variable** (i.e., variable whose name is _not_ both prefixed
    and suffixed by `__`) bound to the passed object (_in descending order of
    size_), calculated with the external :mod:`pympler` package if importable
    _or_ raising an exception otherwise.

    Parameters
    ----------
    obj : object
        Object to iterate the variable sizes of.
    divisor: optional[NumericTypes]
        Number (i.e., integer, float) by which to divide the size in bytes of
        each variable of this object, hence converting from sizes in bytes to
        sizes in a larger denomination (e.g., :data:`betse.util.type.ints.KiB`,
        yielding sizes in kibibytes). Defaults to 1, yielding sizes in bytes.

    Returns
    ----------
    SequenceTypes
        Each element of this sequence is a 2-tuple `(var_name, var_size)` of
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
    :func:`get_object_size`
        Further details on object size calculation.
    :func:`iter_vars_custom`
        Further details on non-builtin object variables.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import iterables

    # List of all 2-tuples to be returned.
    vars_size = []

    # For the name and value of each non-builtin variable of this object (in
    # arbitrary order)...
    for var_name, var_value in iter_vars_custom(obj):
        # Recursive in-memory size of this object.
        var_size = get_object_size(obj=var_value, divisor=divisor)

        # Append a 2-tuple of this variable's name and size.
        vars_size.append((var_name, var_size))

    # Return this list sorted on the second element of each such 2-tuple (i.e.,
    # each variable's size in bytes) in descending order.
    return iterables.sort_by_index_descending(
        iterable=vars_size, subiterable_index=1)
