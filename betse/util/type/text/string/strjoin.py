#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **string joining** (i.e., concatenation of two or more strings,
typically separated by a delimiting substring) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type import types
from betse.util.type.types import type_check, IterableTypes

# ....................{ JOINERS                           }....................
def join(*texts) -> str:
    '''
    Concatenation of the passed strings with no separating delimiter.

    This is a convenience function wrapping the standard ``"".join((...))``
    method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='')


@type_check
def join_by_index(iterable: IterableTypes, subiterable_index: object) -> str:
    '''
    Concatenation of each string at the passed key or index of each subiterable
    of the passed iterable, with no separating delimiter.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable of subiterables to be summed.
    subiterable_index : object
        Object with which to index each subiterable of this iterable. The type
        of this object *must* be a type accepted by the ``__getitem__()``
        special method of each subiterable. If each subiterable is a:

        * **Mapping** (e.g., :class:`dict`), this object *must* be hashable.
        * **Sequence** (e.g., :class:`list`, :class:`tuple`), this object
          *must* be either:

          * An integer.
          * A :func:`slice` object.

    Returns
    ----------
    str
        Concatenation of each string at this key or index of each subiterable
        of this iterable, with no separating delimiter.
    '''

    # Efficiency and simplicity combine here to form SuperHappyFunFunction.
    return ''.join(subiterable[subiterable_index] for subiterable in iterable)

# ....................{ JOINERS ~ on                      }....................
def join_on_newline(*texts) -> str:
    '''
    Join the passed strings with newline as the separating delimiter.

    This is a convnience function wrapping the standard
    ``"\n".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='\n')


def join_on_dot(*texts) -> str:
    '''
    Join the passed strings with a period as the separating delimiter.

    This is a convnience function wrapping the standard
    ``".".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter='.')


def join_on_space(*texts) -> str:
    '''
    Join the passed strings with a space as the separating delimiter.

    This is a convnience function wrapping the standard
    ``" ".join((...))`` method, whose syntax is arguably overly verbose.
    '''

    return join_on(*texts, delimiter=' ')


@type_check
def join_on(*texts: IterableTypes, delimiter: str) -> str:
    '''
    Join the passed strings with the passed separating delimiter.

    This is a convenience function wrapping the standard
    ``"...".join((...))`` method, whose syntax is arguably overly obfuscated.

    Parameters
    ----------
    texts : SequenceTypes
        Tuple of all strings to be joined, consisting of either:

        * One or more strings.
        * One iterable of strings.
    delimiter : str
        Substring to join each such string on.

    Returns
    ----------
    str
        String joined from the passed strings with the passed delimiter.
    '''

    # If only one object was passed and this object is a non-string iterable
    # (e.g, list, tuple), default the list of passed strings to this object.
    if len(texts) == 1 and types.is_iterable_nonstr(texts[0]):
        texts = texts[0]

    # Join these strings.
    return delimiter.join(texts)

# ....................{ JOINERS ~ as                      }....................
@type_check
def join_as(
    *texts,
    delimiter_if_two: str,
    delimiter_if_three_or_more_nonlast: str,
    delimiter_if_three_or_more_last: str
) -> str:
    '''
    Join the passed strings in a human-readable manner.

    If:

    * No strings are passed, the empty string is returned.
    * One string is passed, this string is returned as is without modification.
    * Two strings are passed, these strings are joined with the passed
      ``delimiter_if_two`` separator.
    * Three or more strings are passed:

      * All such strings except the last two are joined with the passed
        ``delimiter_if_three_or_more_nonlast`` separator.
      * The last two such strings are joined with the passed
        ``delimiter_if_three_or_more_last`` separator.

    Parameters
    ----------
    texts : tuple[str]
        Tuple of all strings to be joined.
    delimiter_if_two : str
        String separating each item of ``texts`` if ``len(texts) == 2``.
    delimiter_if_three_or_more_nonlast : str
        String separating each item *except* the last two of ``texts`` if
        ``len(texts) >= 3``.
    delimiter_if_three_or_more_last : str
        String separating the last two items of ``texts`` if
        ``len(texts) >= 3``.

    Returns
    ----------
    str
        Resulting string as described above.

    Examples
    ----------
    >>> join_as(
    ...     ('Fulgrim', 'Perturabo', 'Angron', 'Mortarion'),
    ...     delimiter_if_two=' and ',
    ...     delimiter_if_three_or_more_nonlast=', ',
    ...     delimiter_if_three_or_more_last=', and '
    ... )
    'Fulgrim, Perturabo, Angron, and Mortarion'
    '''

    # Number of passed strings.
    texts_count = len(texts)

    # If no strings are passed, return the empty string.
    if texts_count == 0:
        return ''
    # If one string is passed, return this string as is.
    elif texts_count == 1:
        return texts[0]
    # If two strings are passed, return these strings joined appropriately.
    elif texts_count == 2:
        return '{}{}{}'.format(
            texts[0], delimiter_if_two, texts[1])
    # If three or more strings are passed, return these strings joined
    # appropriately.
    else:
        # All such strings except the last two, joined appropriately.
        texts_nonlast = join_on(
            *texts[0:-2], delimiter=delimiter_if_three_or_more_nonlast)

        # The last two such strings, joined appropriately.
        texts_last = '{}{}{}'.format(
            texts[-2], delimiter_if_three_or_more_last, texts[-1])

        # Return these two substrings, joined appropriately.
        return '{}{}{}'.format(
            texts_nonlast, delimiter_if_three_or_more_nonlast, texts_last)

# ....................{ JOINERS ~ as : conjunction        }....................
@type_check
def join_as_conjunction(*texts: str) -> str:
    '''
    Conjunctively join all passed strings in a human-readable manner.

    Specifically:

    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, and ``.
    '''

    return join_as(
        *texts,
        delimiter_if_two=' and ',
        delimiter_if_three_or_more_nonlast=', ',
        delimiter_if_three_or_more_last=', and '
    )


@type_check
def join_as_conjunction_double_quoted(*texts: str) -> str:
    '''
    Conjunctively double-quote and join all passed strings in a human-readable
    manner.

    Specifically:

    * All passed strings are double-quoted.
    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, and ``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strs

    # Tuple of all passed strings double-quoted. Since the "*" operator applied
    # to this tuple below requires a sequence rather than generator, this is
    # the most space-efficient available sequence (i.e., frozen tuple).
    texts_quoted = tuple(strs.double_quote(text) for text in texts)

    # Conjunctively join these strings.
    return join_as_conjunction(*texts_quoted)


@type_check
def join_iterable_as_conjunction_double_quoted(iterable: IterableTypes) -> str:
    '''
    Conjunctively double-quote and join all string items of the passed iterable
    in a human-readable manner, implicitly converting each item of this
    iterable to a string as needed.

    See Also
    ----------
    :func:`join_as_conjunction_double_quoted`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iterables

    # Iterable of the same type as that of the passed iterable, implicitly
    # converting each item of this iterable to a string as needed.
    iterable_strs = iterables.to_iterable(iterable=iterable, item_cls=str)

    # Return these strings conjunctively double-quoted and joined.
    return join_as_conjunction_double_quoted(*iterable_strs)

# ....................{ JOINERS ~ as : conjunction        }....................
@type_check
def join_as_disjunction(*texts: str) -> str:
    '''
    Disjunctively join all passed strings in a human-readable manner.

    Specifically:

    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, or ``.
    '''

    return join_as(
        *texts,
        delimiter_if_two=' or ',
        delimiter_if_three_or_more_nonlast=', ',
        delimiter_if_three_or_more_last=', or '
    )


@type_check
def join_as_disconjunction_double_quoted(*texts: str) -> str:
    '''
    Disjunctively double-quote and join all passed strings in a human-readable
    manner.

    Specifically:

    * All passed strings are double-quoted.
    * All passed strings excluding the last two are joined with ``, ``.
    * The last two passed strings are joined with ``, or ``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strs

    # Tuple of all passed strings double-quoted. See
    # join_as_conjunction_double_quoted().
    texts_quoted = tuple(strs.double_quote(text) for text in texts)

    # Disjunctively join these strings.
    return join_as_disjunction(*texts_quoted)
