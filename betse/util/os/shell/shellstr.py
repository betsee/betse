#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level shell-specific string facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, IntOrNoneTypes

# ....................{ GETTERS                            }....................
#FIXME: Consider integrating this method with the main codebase, ideally by
#calling this function before running each simulation phase. This banner is
#sufficiently aesthetic that I'd certainly appreciate seeing it everywhere.
@type_check
def get_banner(
    # Mandatory arguments.
    title: str,

    # Optional arguments.
    padding: str = '=',
    length: IntOrNoneTypes = None,
) -> str:
    '''
    Single-line banner embedding the passed human-readable title.

    Parameters
    ----------
    title : str
        Human-readable substring to be embedded in this banner, ideally less
        than or equal to 40 characters in length (i.e., half the default
        terminal width and thus banner length of 80 characters).
    padding : optional[str]
        Single character with which to pad this title. To visually center this
        title in this banner, this character is repeated both before and after
        this title as necessary. Defaults to an ASCII punctuation character.
    length : IntOrNoneTypes
        Preferred length of the returned banner. If the length of this title
        does *not* already exceed this preferred length, the returned banner is
        guaranteed to have this exact length. Defaults to ``None``, in which
        case this length internally defaults to:
        * If the current shell is an interactive terminal, the **width** (i.e.,
          maximum number of columns currently displayed by) this terminal.
        * Else, the conventional default of 80.

    Returns
    ----------
    str
        Single-line banner embedding this title.

    Examples
    ----------
        >>> from betse.util.os.shell import shellstr
        >>> shellstr.get_banner(
        ...     title='The Love Song of J. Alfred Prufrock', padding='~')
        ~~~~~~~~~~~~~~~~~~~~~~ The Love Song of J. Alfred Prufrock ~~~~~~~~~~~~~~~~~~~~~
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import chars

    # If this padding is *NOT* a single character, raise an exception.
    chars.die_unless_char(padding)

    #FIXME: Generalize this hard-coded maximum length as follows:
    #
    #* If the current shell is interactive, set this length to this shell's
    #  current column length. Sadly, doing so requires platform-specific logic.
    #* Else, default this length to the current reasonable default of 80.

    # If the caller passed no preferred banner length...
    if length is None:
        # Default this length to the conventional terminal width.
        length = 80

    # Minimal-length banner containing this text and the smallest
    # possible prefixing and suffixing banner characters.
    banner_min = '{padding} {title} {padding}'.format(
        title=title, padding=padding)

    # If the length of this minimal-length banner exceeds the maximum
    # length, print this text with no such banner.
    if len(banner_min) > length:
        return title
    # Else, print this text within a banner.
    else:
        #FIXME: Again, generalize this hard-coded length in a manner dependent
        #on the current terminal width.

        # Length of the string preceding this text. To produce a uniform
        # rather than ragged left margin for aesthetic uniformity, this
        # length is a constant.
        BANNER_PREFIX_LEN = 30

        # String preceding this text.
        banner_prefix = padding * BANNER_PREFIX_LEN

        # Banner to be printed excluding all suffixing banner characters.
        banner_sans_suffix = '{} {} '.format(banner_prefix, title)

        # Length of the string following this text. To dynamically fill all
        # remaining line space, this length contextually depends on the
        # lengths of all other strings. By the above conditional, this
        # length is guaranteed to be non-zero and need *NOT* be tested.
        banner_suffix_len = length - len(banner_sans_suffix)

        # Nonetheless, test this length for sanity.
        assert banner_suffix_len > 0, (
            'Banner suffix length {} not positive.'.format(
                banner_suffix_len))

        # String following this text.
        banner_suffix = padding * banner_suffix_len

        # Single-line banner to be returned.
        banner = '\n{}{}'.format(banner_sans_suffix, banner_suffix)

        # Return this banner.
        return banner

# ....................{ QUOTERS                            }....................
@type_check
def shell_quote(text: str) -> str:
    '''
    Shell-quote the passed string in a platform-specific manner.

    If the current platform is:

    * *Not* Windows (e.g., Linux, OS X), the returned string is guaranteed to be
      suitable for passing as an arbitrary positional argument to external
      commands.
    * Windows, the returned string is suitable for passing *only* to external
      commands parsing arguments in the same manner as the Microsoft C runtime.
      While *all* applications on POSIX-compliant systems are required to parse
      arguments in the same manner (i.e., according to Bourne shell lexing), no
      such standard applies to Windows applications. Shell quoting is therefore
      fragile under Windows -- like pretty much everything.
    '''

    # Avoid circular import dependencies.
    from betse.util.os import oses

    # If the current OS is Windows, do *NOT* perform POSIX-compatible quoting.
    # Windows is POSIX-incompatible and hence does *NOT* parse command-line
    # arguments according to POSIX standards. In particular, Windows does *NOT*
    # treat single-quoted arguments as single arguments but rather as multiple
    # shell words delimited by the raw literal `'`. This is circumventable by
    # calling an officially undocumented Windows-specific function. (Awesome.)
    if oses.is_windows():
        import subprocess
        return subprocess.list2cmdline([text])
    # Else, perform POSIX-compatible quoting.
    else:
        import shlex
        return shlex.quote(text)
