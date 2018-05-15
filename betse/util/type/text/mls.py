#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

r'''
Low-level \*ML (e.g., HTML, SGML, XML) facilities.
'''

# ....................{ IMPORTS                            }....................
import html
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.types import type_check, RegexCompiledType

# ....................{ GETTERS                            }....................
@func_cached
def get_tag_regex() -> RegexCompiledType:
    r'''
    Compiled regular expression matching a syntactically but *not* necessarily
    semantically valid \*ML (e.g., HTML, SGML, XML) tag.

    For efficiency in downstream clients (e.g., BETSEE) frequently calling the
    :func:`is_ml` function and hence requiring this expression be compiled, this
    expression is intentionally pre-compiled rather than uncompiled and thus
    returned as a cached getter.

    See Also
    ----------
    https://kevin.deldycke.com/2008/07/python-ultimate-regular-expression-to-catch-html-tags
        Blog article strongly inspiring this regular expression.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # Create, return, and cache this expression.
    return regexes.compile_regex(
        # Tag prefix.
        r'<'
        # Tag closure prefix (optional).
        r'\/?'
        # Tag name.
        r'\w+'
        # Zero or more tag attributes.
        r'(?:'
            # Attribute name.
            r'\s+\w+'
            # Attribute value (optional).
            r'(?:'
                # Equals sign.
                r'\s*=\s*'
                # Either:
                r'(?:'
                    # Double-quoted string.
                    r'".*?"|'
                    # Single-quoted string.
                    r"'.*?'|"
                    # Unquoted word.
                    r"[^'"
                    r'">\s]+'
                r')'
            r')?'
        r')*'
        # Whitespace (optional).
        r'\s*'
        # Tag closure suffix (optional).
        r'\/?'
        # Tag suffix.
        r'>'
    )

# ....................{ TESTERS                            }....................
@type_check
def is_ml(text: str) -> bool:
    r'''
    ``True`` only if the passed string contains one or more \*ML tags and hence
    superficially appears to contain \*ML-specific syntax.

    Parameters
    ----------
    text : str
        String to be tested for the presence of \*ML-specific syntax.

    Returns
    ----------
    bool
        ``True`` only if this string contains one or more \*ML tags.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import regexes

    # Test this string against this expression.
    return regexes.is_match(text=text, regex=get_tag_regex())

# ....................{ ESCAPERS                           }....................
@type_check
def escape_ml(text: str) -> str:
    r'''
    Passed string with all \*ML-specific syntax escaped and hence guaranteed to
    contain *no* \*ML-specific syntax.

    Specifically, this function converts all:

    * ``<`` characters into ``&lt;`` substrings.
    * ``>`` characters into ``&gt;`` substrings.
    * ``&`` characters into ``&amp;`` substrings.
    * ``'`` characters into ``&#x27;`` substrings.
    * ``"`` characters into ``&quot;`` substrings.

    Caveats
    ----------
    This function renders arbitrary strings syntactically but *not* necessarily
    semantically suitable for use as valid \*ML (e.g., HTML, SGML, XML).
    Notably, returned strings may still contain unescaped JavaScript and hence
    cross-site scripting (XSS) vulnerabilities. These strings should *never* be
    directly displayed as HTML to end users. For example, the following call to
    this function preserves a JavaScript fragment unescaped:

        >>> from betse.util.type.text import mls
        >>> mls.escape_ml('<a href="javascript:alert()">')
        '&lt;a href=&quot;javascript:alert()&quot;&gt;'

    Parameters
    ----------
    text : str
        String to escape all \*ML-specific syntax of.

    Returns
    ----------
    str
        String with all \*ML-specific syntax escaped.
    '''

    # Praise be the stdlib.
    return html.escape(text)

# ....................{ TAGIFIERS                          }....................
@type_check
def tagify_newlines(text: str) -> str:
    '''
    Convert each newline in the passed string into a ``<br/>`` tag.

    Parameters
    ----------
    text : str
        String to convert all newlines of.

    Returns
    ----------
    str
        String with all newlines converted to corresponding tags.
    '''

    # Sometimes, it really is just that easy.
    return text.replace('\n', '<br/>')
