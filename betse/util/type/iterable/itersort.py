#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-string iterable sorting** (i.e., sorting of non-string objects
implementing the abstract base class :class:`collections.abc.Iterable`)
facilities.
'''

#FIXME: Generalize to support mappings. Doing so may prove mildly non-trivial,
#as unordered mapping types (e.g., plain "dict" instances) will need to be
#converted into ordered mapping types (e.g., "OrderedDict" instances).

# ....................{ IMPORTS                           }....................
from betse.util.type.types import type_check, IterableTypes
from operator import itemgetter

# ....................{ SORTERS ~ ascending               }....................
@type_check
def sort_ascending(iterable: IterableTypes) -> IterableTypes:
    '''
    Iterable sorted from the passed iterable in ascending order.

    Each item of this iterable is compared to each other item of this iterable
    under the ``<`` operator, thus implicitly calling the ``__le__()`` special
    method of each such item. Each item should ideally but *not* necessarily be
    of the same type. If each item is:

    * A string, these strings are sorted in **ascending lexicographic order**
      (i.e., traditional order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **ascending numeric order.**

    Caveats
    ----------
    If the passed iterable is a generator, this function:

    * Internally coerces this generator into a list *before* attempting to
      sort the items yielded by this generator.
    * Externally returns a new list rather than a new generator.

    This constraint is mandated by the :func:`sorted` builtin, which raises the
    following exception on receiving a generator:

       TypeError: cannot create 'generator' instances

    Parameters
    ----------
    iterable : IterableTypes
        Unsorted iterable to be returned sorted. For generality, this iterable
        is *not* modified by this function.

    Returns
    ----------
    IterableTypes
        Iterable sorted from and of the same type as the passed iterable. For
        efficiency, this iterable is only a shallow rather than deep copy of
        the passed iterable. Note lastly that the class of the passed iterable
        *must* define an ``__init__()`` method accepting a list.
    '''

    # Defer to this lowor-level function. Since the default behaviour of the
    # sorted() builtin and thus this function is to sort in ascending
    # lexicographic order, no additional keyword arguments are required.
    return _sort_iterable(iterable)


@type_check
def sort_by_index_ascending(
    iterable: IterableTypes, subiterable_index: object) -> IterableTypes:
    '''
    Iterable of subiterables sorted from the passed iterable of subiterables in
    ascending order of the value of each element at the passed key or index of
    each subiterable of this iterable.

    Each item at the passed key or index of each subiterable of this iterable
    is compared to each other item at each other key or index of each other
    subiterable of this iterable via the ``<`` operator, implicitly calling the
    ``__le__()`` special method of these items. Each item is ideally but *not*
    necessarily of the same type. If each item is:

    * A string, these strings are sorted in **ascending lexicographic order**
      (i.e., traditional order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **ascending numeric order.**

    Parameters
    ----------
    iterable : IterableTypes
        Unsorted iterable of subiterables to be returned sorted. For
        generality, neither this iterable nor these subiterables are modified
        by this function.
    subiterable_index : object
        Object with which to index each subiterable of this iterable. The type
        of this object *must* be a type accepted by the ``__getitem__()``
        special method of each subiterable. Specifically, if each subiterable
        is a:

        * **Mapping** (e.g., :class:`dict`), this object *must* be hashable.
        * **Sequence** (e.g., :class:`list`, :class:`tuple`), this object
          *must* be either:

          * An integer.
          * A :func:`slice` object.

    Returns
    ----------
    IterableTypes
        Iterable of subiterables sorted from and of the same type as the passed
        iterable of subiterables. For efficiency, this iterable is only a
        shallow rather than deep copy of the passed iterable. Note lastly that
        the class of the passed iterable *must* define an ``__init__()`` method
        accepting a list.

    See Also
    ----------
    :func:`sort_ascending`
        Further details.
    '''

    # Defer to this lowor-level function.
    #
    # For efficiency, items of each subiterable of this iterable are retrieved
    # via the function internally created and called by calling an instance of
    # the standard "itemgetter" class. While technically an instance of that
    # class, that class is intended to be treated as a simple function passed
    # the desired index returning a function passed a subiterable returning the
    # element at that index: e.g.,
    #
    #     >>> penguins = [('adelie', 0xFEEDFACE), ('gentoo', 0xDEADBEEF)]
    #     >>> penguin_id = itemgetter(2)
    #     >>> penguin_id(penguins[0]) == 0xFEEDFACE
    #     True
    #
    # Despite the internal complexity and object overhead imposed by the
    # "itemgetter" class, this approach has been definitively profiled to be
    # 126% faster on average than the traditional
    # "lambda subiterable: subiterable[subiterable_index]" approach --
    # presumably due to hidden overhead imposed by capturing all local
    # variables into the closure context encapsulated by the lambda. See also
    # the following Stackoverflow answer exhibiting this profiling:
    #
    #     https://stackoverflow.com/a/17243726/2809027
    return _sort_iterable(iterable, key=itemgetter(subiterable_index))

# ....................{ SORTERS ~ descending              }....................
@type_check
def sort_descending(iterable: IterableTypes) -> IterableTypes:
    '''
    Iterable sorted from the passed iterable in descending order.

    Each item of this iterable is compared to each other item of this iterable
    via the ``>`` operator, implicitly calling the ``__ge__()`` special method
    of these items. Each item is ideally but *not* necessarily of the same
    type. If each item is:

    * A string, these strings are sorted in **descending lexicographic order**
      (i.e., reverse order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **descending numeric order.**

    See Also
    ----------
    :func:`sort_ascending`
        Further details.
    '''

    # Defer to this lowor-level function.
    return _sort_iterable(iterable, reverse=True)


@type_check
def sort_by_index_descending(
    iterable: IterableTypes, subiterable_index: object) -> IterableTypes:
    '''
    Iterable of subiterables sorted from the passed iterable of subiterables in
    descending order of the value of each element at the passed key or index of
    each subiterable of this iterable.

    Each element at the passed key or index of each subiterable of this
    iterable is compared to each other element at each other key or index of
    each other subiterable of this iterable via the ``>`` operator, implicitly
    calling the ``__ge__()`` special method of these elements. Each element is
    ideally but *not* necessarily of the same type. If each element is:

    * A string, these strings are sorted in **descending lexicographic order**
      (i.e., reverse order of dead-tree dictionaries and encyclopedias).
    * A number (i.e., integer or a float), these numbers are sorted in
      **descending numeric order.**

    See Also
    ----------
    :func:`sort_by_index_ascending`
        Further details.
    '''

    # Defer to this lowor-level function. See sort_by_index_descending() for
    # additional commentary, particularly for the "key" parameter.
    return _sort_iterable(
        iterable,
        key=itemgetter(subiterable_index),
        reverse=True,
    )

# ....................{ SORTERS ~ private                 }....................
@type_check
def _sort_iterable(iterable: IterableTypes, **kwargs) -> IterableTypes:
    '''
    Iterable sorted from the passed iterable according to the passed keyword
    arguments accepted by the :func:`sorted` builtin.

    Parameters
    ----------
    iterable : IterableTypes
        Unsorted iterable to be returned sorted. For generality, this iterable
        is *not* modified by this function.

    All remaining keyword arguments are passed as is to the :func:`sorted`
    builtin.

    Returns
    ----------
    IterableTypes
        Iterable sorted from and of the same type as the passed iterable. For
        efficiency, this iterable is only a shallow rather than deep copy of
        the passed iterable. Note lastly that the class of the passed iterable
        *must* define an ``__init__()`` method accepting a list.

    See Also
    ----------
    :func:`sort_ascending`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import generators

    # Input iterable to be sorted.
    #
    # If this iterable is a generator, coerce this generator into a list.
    # Generators are *NOT* safely sortable as is. See the function docsctring
    # for details. Note that coercing this generator into a list ensures that
    # the "list_sorted" generated below is efficiently returned as is without
    # requiring further inefficient type coercion (e.g., from a list to a
    # tuple).
    iterable_sortable = (
        list(iterable) if generators.is_generator(iterable) else iterable)

    # Type of both this iterable *AND* the output iterable to be returned.
    iterable_type = type(iterable_sortable)

    # Output list sorted from this iterable.
    list_sorted = sorted(iterable_sortable, **kwargs)

    # Return either:
    #
    # * If the input iterable is:
    #   * A list, efficiently return this output list as is, which already
    #     satisfies the exact type expected by the caller.
    # * Else, an iterable of the same type as that of the passed iterable.
    return list_sorted if iterable_type is list else iterable_type(list_sorted)
