#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level logging formatter subclasses.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific modules
# must *ALWAYS* be able to safely import from this module at any level, these
# circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# from betse.util.type.types import type_check
from logging import Formatter

# ....................{ CLASSES                            }....................
#FIXME: Unfortunately, this fundamentally fails to work. The reason why? The
#"TextWrapper" class inserts spurious newlines *EVEN WHEN YOU EXPLICITLY TELL
#IT NOT TO*. This is crazy, but noted in the documentation:
#
#    "If replace_whitespace is False, newlines may appear in the middle of a
#     line and cause strange output. For this reason, text should be split into
#     paragraphs (using str.splitlines() or similar) which are wrapped
#     separately."
#
#Until this is resolved, the only remaining means of wrapping log messages will
#be to define new top-level module functions suffixed by "_wrapped" ensuring
#that the appropriate formatter is used (e.g., a new log_info_wrapped()
#function). For now, let's just avoid the topic entirely. It's all a bit
#cumbersome and we're rather weary of it.

class LogFormatterWrap(Formatter):
    '''
    Log formatter wrapping all lines in handled log records to a sane line
    length (e.g., 80 characters).

    Attributes
    ----------
    _text_wrapper : TextWrapper
        Object with which to wrap log messages, cached for efficiency.
    '''

    pass
    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self._text_wrapper = TextWrapper(
    #         drop_whitespace = False,
    #         replace_whitespace = False,
    #     )

    # def format(self, log_record: LogRecord) -> str:
    #     # Avoid circular import dependencies.
    #     from betse.util.type import strs
    #
    #     # Get such message by (in order):
    #     #
    #     # * Formatting such message according to our superclass.
    #     # * Wrapping such formatted message.
    #     return strs.wrap(
    #         text = super().format(log_record),
    #         text_wrapper = self._text_wrapper,
    #     )
