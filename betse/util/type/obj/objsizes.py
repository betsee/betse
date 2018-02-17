#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level wrappers around Pympler, an optional runtime dependency performing
memory profiling of Python applications.
'''

#FIXME: Consider donating the get_size_profile() function back to Pympler, which
#could *REALLY* benefit from a human-readable CLI-based object graph profile.

# ....................{ IMPORTS                            }....................
from betse.util.type import types
from betse.util.type.types import type_check, NumericSimpleTypes

# ....................{ GETTERS                            }....................
@type_check
def get_size(obj: object, size_divisor: NumericSimpleTypes = 1) -> NumericSimpleTypes:
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
    size_divisor: optional[NumericSimpleTypes]
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
    libs.die_unless_runtime_optional('pympler')

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
def get_sizes_vars(obj: object, *args, **kwargs) -> tuple:
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

    All remaining positional and keyword arguments are passed as is to the
    :func:`get_size` function.

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
    from betse.util.type.numeric import numerics
    from betse.util.type.obj import objects

    # List of all 2-tuples to be returned.
    vars_size = []

    # For the name and value of each non-builtin variable of this object (in
    # arbitrary order)...
    for var_name, var_value in objects.iter_vars_custom(obj):
        # Recursive in-memory size of this object.
        var_size = get_size(obj=var_value, *args, **kwargs)

        # Append a 2-tuple of this variable's name and size.
        vars_size.append((var_name, var_size))

    # Total recursive size of this object, totalized over the second element of
    # each 2-tuple in this list.
    obj_size = numerics.sum_by_index(iterable=vars_size, subiterable_index=1)

    # This list sorted on the second element of each such 2-tuple (i.e., each
    # variable's size) in descending order.
    vars_size_sorted = iterables.sort_by_index_descending(
        iterable=vars_size, subiterable_index=1)

    # Return the expected 2-tuple.
    return obj_size, vars_size_sorted

# ....................{ GETTERS ~ profile                  }....................
#FIXME: Shift above.

_SIZE_FORMATTER = r'{:.04f} MiB'
'''
Format specifier converting the size in bytes interpolated into this string with
the :func:`format` builtin into a human-readable string.
'''


_TYPE_FORMATTER = r'<{.__class__.__name__}>: '
'''
Format specifier converting an object type interpolated into this string with
the :func:`format` builtin into a human-readable string.
'''


#FIXME: Revise docstring, removing reference to get_sizes_vars() in particular.

@type_check
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

    >>> from betse.util.type.obj import objsizes
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
    line_indent: optional[str]
        String to prefix each line of the returned string by, typically
        consisting of zero or more space characters. Defaults to the empty
        string.
    vars_depth: optional[int]
        Maximum depth of the recursion tree induced by calling this function.
        If this depth is:
        * Zero, a single-line synopsis of the total size of this object is
          returned _without_ totalling the sizes of any variables bound to this
          object.  This depth is appropriate for C-based objects (i.e., objects
          whose classes are either primitive builtins or defined by external C
          extensions), whose details are usually inconsequential when profiling
          the space consumed by pure-Python objects. Since depth is necessarily
          positive, recursion "bottoms out" at this depth.
        * Non-zero, a single-line synopsis of the total size of this object
          appended by the total sizes of all variables bound to this object is
          recursively returned. Specifically, for each such variable, this
          function recursively calls itself with:
          * `obj` parameter equal to the value of that variable.
          * `line_indent` parameter appended by two spaces.
          * `vars_depth` parameter decremented by one.
          * `vars_max` parameter divided by the total number of variables bound
            to this object.
        Defaults to 1, in which case a single-line synopsis of the total size of
        this object appended by the total sizes of all variables directly bound
        to this object is non-recursively returned.
    vars_max : optional[int]
        Maximum number of the largest non-builtin variables bound to this object
        to synopsize the sizes of. Defaults to `None`, in which case this
        function unconditionally synopsizes the sizes of _all_ non-builtin
        variables bound to this object (regardless of size).

    All remaining positional and keyword arguments are passed as is to the
    :func:`get_size` function.

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
    from betse.util.type import iterables
    from betse.util.type.text import strs
    from betse.util.type.numeric.ints import MiB
    from betse.util.type.obj import objects

    #FIXME: Convert the current "kwargs['size_divisor']" into a human-readable
    #denomenation (e.g., "GiB", "KiB", "B"). Since this is both non-trivial and
    #currently *NOT* required, we adopt a superior lazy approach.

    # If the caller failed to specify a size divisor, default to a divisor
    # denominating sizes in mebibytes for readability.
    kwargs.setdefault('size_divisor', MiB)

    # Set of the unique identifiers of all objects previously passed to a prior
    # recursive call by the closure defined below, preventing erroneous
    # reinspection of these objects.
    obj_ids = set()

    #FIXME: Revise docstring.

    # For efficiency, this function is intentionally *NOT* type-checked.
    def _get_size_profile_and_size(
        obj: object,
        line_indent: str,
        vars_depth: int,
        vars_max: int,
    ) -> tuple:
        '''
        Human-readable synopsis of the size of the passed object with all
        remaining arguments passed as is to the :func:`get_sizes_vars` function,
        calculated with the external :mod:`pympler` package if importable _or_
        raising an exception otherwise.

        Specifically, this function returns a string tabulating:

        * The total recursive in-memory size in mebibytes of this object.
        * The names and total recursive in-memory sizes in mebibytes of all
          **non-builtin variables** (i.e., variables whose names are _not_ both
          prefixed and suffixed by `__`) bound to this object (_in descending
          order of size_).

        Parameters
        ----------
        This function has the same signature with respect to passed parameters
        as the :func:`get_size_profile` function.

        Returns
        ----------
        (str, NumericOrNoneTypes)
            2-tuple `(obj_size_profile, obj_size)`, where:
            * `obj_size_profile` is a recursive size profile of this object.
            * `obj_size` is the total recursive size of this object in bytes if
              this object has _not_ already been passed to a prior recursive
              call of this closure or `None` otherwise. Note that, in this case,
              `obj_size_profile` is a relevant human-readable string indicating
              this circularity in the object graph.
        '''

        # If this object has already been passed to a prior recursive call of
        # this closure, skip this object and notify the caller of our having
        # done so by returning the empty 2-tuple.
        if id(obj) in obj_ids:
            # Total recursive size of this object, zeroed to prevent parent
            # recursive calls from erroneously re-summing with this size.
            obj_size = 0

            # Notify the end user of this object graph circularity.
            size_profile = (
                _TYPE_FORMATTER + '(duplicate object already sized above)'
            ).format(obj)

            # Return the expected 2-tuple up the recursive call stack.
            return size_profile, obj_size

        # Else, this object has *NOT* yet been passed to a prior recursive call
        # of this closure. To prevent recursive calls below from erroneously
        # reinspecting this object, record this object's unique identifier.
        obj_ids.add(id(obj))

        # True only if this function is immediately returning *WITHOUT*
        # recursivily calling itself below.
        #
        # Specifically, if any of the following conditions apply, halt this
        # recursive subtree by non-recursively returning a simple size profile:
        #
        # * This is the bottom-most recursive call in the current
        #   subtree of recursive calls.
        # * The parent recursive call requested that this child recursive call
        #   profile only this object itself (rather than both this object and
        #   all or some variables bound to this object).
        # * No variables are bound to this object.
        # * This object is C-based (i.e., is either builtin or defined by a C
        #   extension).
        #
        # Note that the second such test (i.e., "vars_max < 1") is technically
        # non-essential, as the list truncation performed below already
        # implicitly handles this edge cases without error. Doing so is highly
        # inefficient, however, and hence best avoided if feasible. It is
        # feasible, so it is.
        #
        # If this object is C-based, recursively introspecting this object would
        # be trivially feasible. Python code makes no distinction between
        # pure-Python and C-based objects. Indeed, differentiating between the
        # two with Python code is non-trivial! That said, recursively
        # introspecting C-based objects would be largely pointless; size
        # profiles are principally intended to profile the space consumed by
        # variables bound to application-specific pure-Python objects rather
        # than those bound to external C-based objects.
        is_nonrecursive = (
            vars_depth == 0 or vars_max < 1 or objects.is_c_based(obj))

        # If this function is possibly recursively calling itself below...
        if not is_nonrecursive:
            # Tuple of 2-tuples "(var_name, var_value)" for each non-builtin
            # variable bound to this object, permitting the number of 2-tuples
            # yielded by this generator to be summed below.
            #
            # For efficiency, this inefficient operation is performed only if
            # this function is possibly recursively calling itself below. If
            # this function is *NOT* doing so, this operation is avoided.
            vars_name_value = tuple(objects.iter_vars_custom(obj))

            # Number of such variables.
            vars_count = len(vars_name_value)

            # If no variables are bound to this object, return immediately.
            is_nonrecursive = not vars_count

        # If this function is returning *WITHOUT* recursivily calling itself...
        if is_nonrecursive:
            # Total recursive size of only this object.
            obj_size = get_size(obj, *args, **kwargs)

            # Total recursive size profile of only this object.
            size_profile = (
                _TYPE_FORMATTER + _SIZE_FORMATTER).format(obj, obj_size)

            # Return the expected 2-tuple up the recursive call stack.
            return size_profile, obj_size
        # Else, recursively return a size profile of both this object *AND* all
        # or some variables bound to this object. By the above logic, this
        # object is guaranteed to contain at least one such variable.

        # Human-readable string suffixing this synopsis if any or empty.
        size_profile_suffix = ''

        # Substring prefixing each line of variable-specific synopsis, indented
        # to the right by two additional spaces.
        var_line_indent = line_indent + '  '

        # Maximum depth of the recursion tree induced by calling this function
        # recursively below for each variable bound to this object. To ensure
        # that recursion bottoms out, the current maximum is decremented by one.
        var_vars_depth = vars_depth - 1

        # Maximum number of variables bound to each variable bound to this
        # object to synopsize the sizes of if any or None otherwise.
        var_vars_max = vars_max

        # If such a maximum is requested, distribute this maximum amongst all
        # available variables, rounding down to the nearest integer.
        #
        # Note that, if the number of such variables is larger than this
        # maximum, this maximum will be reduced to zero. This valid edge case is
        # explicitly handled by the conditional above.
        if var_vars_max is not None:
            var_vars_max = round(vars_max / vars_count)

        # List of 2-tuples "(var_size_profile, var_size)" for each variable
        # bound to this object, where:
        #
        # * "var_size_profile" is a recursive size profile of this variable.
        # * "var_size" is the total recursive size of this variable in bytes.
        #
        # Associating each profile with its size is required to subsequently
        # sort and possibly discard profiles by size.
        vars_size_profile = []

        # Total recursive size of this object, incremented by iteration below.
        obj_size = 0

        #FIXME: Consider left-aligning the sizes interpolated by this iteration.
        #FIXME: Consider interpolating the truncation of the value of each
        #variable in this iteration.

        # For the name and size of each of this object's variables...
        for var_name, var_value in vars_name_value:
            # Size profile and size of this variable, computed recursively.
            var_size_profile, var_size = _get_size_profile_and_size(
                obj=var_value,
                line_indent=var_line_indent,
                vars_depth=var_vars_depth,
                vars_max=var_vars_max,
            )

            # Prefix this variable's size profile by indentation and its name.
            var_size_profile = '{}{} {}'.format(
                var_line_indent, var_name, var_size_profile)

            # Append this size and size profile to this list as a 2-tuple.
            vars_size_profile.append((var_size_profile, var_size))

            # Increment the size of this parent object by this child variable.
            obj_size += var_size

        # This list sorted on the second element of each such 2-tuple (i.e.,
        # each variable's size) in descending order.
        vars_size_profile_sorted = iterables.sort_by_index_descending(
            iterable=vars_size_profile, subiterable_index=1)

        # If a maximum number of the largest variables is requested *AND* the
        # number of variables bound to this object exceeds this maximum, ignore
        # all variables of smaller size. To ensure the synopsis suffix set below
        # has the proper indentation, do so *AFTER* increasing this level.
        if vars_max is not None and vars_count > vars_max:
            # print('!!Truncating!!')
            # Validate this number to be a positive integer.
            assert types.is_int_positive(vars_max), (
                types.assert_not_int_positive(vars_max))

            # Truncate this list to this maximum number of the first 2-tuples.
            vars_size_profile_sorted = vars_size_profile_sorted[:vars_max]

            # Append an ellipses to this synopsis, notifying users that
            # additional variables exist but were truncated.
            size_profile_suffix = '{}...'.format(var_line_indent)

        # print('vars_size_profile_sorted: {}'.format(vars_size_profile_sorted))
        # print('vars_size_profile_sorted:')
        # for var_size_profile_sorted, var_size_sorted in (
        #     vars_size_profile_sorted):
        #     print('  str: "{}"; int: {}'.format(var_size_profile_sorted, var_size_sorted))

        #FIXME: For efficiency, define a new "_TYPE_SIZE_FORMATTER" global
        #concatenating "_TYPE_FORMATTER" and "_SIZE_FORMATTER". Then remove both
        #"_TYPE_FORMATTER" and "_SIZE_FORMATTER", which will then be extraneous.

        # Human-readable synopsis of...
        size_profile = (
            # The size of this object.
            (_TYPE_FORMATTER + _SIZE_FORMATTER).format(obj, obj_size) +

            # The size of each variable bound to this object.
            strs.join_by_index(
                iterable=vars_size_profile_sorted, subiterable_index=0) +

            # The substring suffixing this entire synopsis if any.
            size_profile_suffix +

            # A blank line, improving readability and aesthetics.
            '\n'
        )

        # Return the expected 2-tuple up the recursive call stack.
        return size_profile, obj_size

    # Size profile recursively computed from this object. The second element of
    # the tuple returned by this function is intended for use *ONLY* by
    # recursive calls by this function of itself and is hence ignorable here.
    size_profile, _ = _get_size_profile_and_size(
        obj,

        # Prefix the passed substring prefixing every output line by a newline,
        # conveniently simplifying the implementation of _get_size_profile().
        line_indent='\n' + line_indent,

        # Pass all other passed parameters as is.
        vars_depth=vars_depth,
        vars_max=vars_max,
    )

    # Munge this profile as follows:
    #
    # * Prefix this profile by this indentation. For inscrutable recursion
    #   reasons, the first element of the tuple returned by this function is a
    #   string guaranteed *NOT* to be prefixed by this indentation
    # * Remove all suffixing whitespace -- namely, the trailing blank line
    #   resulting again from inscrutable recursion reasons.
    size_profile = line_indent + strs.remove_whitespace_suffix(size_profile)

    # Return this munged profile.
    return size_profile
