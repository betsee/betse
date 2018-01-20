#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level command-line argument facilities.
'''

#FIXME: It'd be great to contribute this back to the official "argparse" module
#by contributing a patch. See:
#    http://bugs.python.org/issue12806

# ....................{ IMPORTS                            }....................
import re
from argparse import HelpFormatter
# from betse.util.io.log import logs
from betse.util.type.text import strs
from betse.util.type.types import type_check, GeneratorType, SequenceTypes

# ....................{ CONSTANTS                          }....................
_BLANK_LINE_REGEX = re.compile(r'^\s*$')
'''
Regular expression matching only blank lines.
'''


_INDENTATION_REGEX = re.compile(r'^( *)')
'''
Regular expression capturing all spaces prefixing the subject string.
'''

# ....................{ FORMATTERS                         }....................
class SemicolonAwareHelpFormatter(HelpFormatter):
    '''
    Formatter wrapping all lines *not* prefixed by ``;`` in
    :class:`ArgumentParser` help text into paragraphs.

    Newlines in all lines prefixed by ``;`` will be preserved as is, preventing
    such lines from being wrapped. Newlines in all other lines will be
    implicitly converted to single spaces. (See examples below.)

    See Also
    ----------
    http://bugs.python.org/issue12806
    http://bugs.python.org/file28091/paraformatter.py
        argparse's lack of official paragraph formatter is long recognized, as
        documented by the above bug thread. A year and a half after the first
        post in such thread, rurpy2 contributed a module "paraformatter", which
        this module remains strongly inspired by.

    Examples
    ----------
        >>> from argparse import ArgumentParser
        >>> from betse.util.path.command import HelpFormatterParagraph
        >>> arg_parser = ArgumentParser(
        ...     formatter_class = HelpFormatterParagraph)
        >>> arg_subparsers = arg_parser.add_subparsers()
        >>> arg_subparsers.add_parser(
        ...     name = 'road',
        ...     description = """
        ...         The Road goes ever on and on
        ...         Down from the door where it began.
        ...         Now far ahead the Road has gone,
        ...         ;And I must follow, if I can,
        ...         ;Pursuing it with eager feet,
        ...         ;Until it joins some larger way,
        ...         ;Where many paths and errands meet.
        ...
        ...            The Road goes ever on and on
        ...            Down from the door where it began.
        ...            Now far ahead the Road has gone,
        ...            And I must follow, if I can,
        ...         ;  Pursuing it with weary feet,
        ...         ;  Until it joins some larger way,
        ...         ;  Where many paths and errands meet.
        ...         ;  And whither then? I cannot say.
        ...     """)
        >>> arg_parser.parse_args('road --help')
        The Road goes ever on and on Down from the door where it began. Now far
        ahead the Road has gone,
        And I must follow, if I can,
        Pursuing it with eager feet,
        Until it joins some larger way,
        Where many paths and errands meet.

        The Road goes ever on and on Down from the door where it began. Now far
        ahead the Road has gone, And I must follow, if I can,
        Pursuing it with weary feet,
        Until it joins some larger way,
        Where many paths and errands meet.
        And whither then? I cannot say.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _fill_text(self, text: str, width: int, indent: str) -> str:
        '''
        Passed string with all non-indented lines wrapped into paragraphs.
        '''

        return '\n'.join(self._split_lines_indented(text, width, indent))


    def _split_lines(self, text: str, width: int) -> list:
        '''
        List of lines produced by wrapping all non-indented lines in the passed
        text into paragraphs and splitting such string on newlines.
        '''

        return self._split_lines_indented(text, width)

    # ..................{ SUBCLASS                           }..................
    @type_check
    def _split_lines_indented(
        self, text: str, width: int, indent: str = '') -> list:
        '''
        List of lines produced by wrapping all non-indented lines in the passed
        text into paragraphs and splitting this string on newlines.
        '''

        # Initialize this text wrapper with the passed arguments.
        strs.text_wrapper.width = width
        strs.text_wrapper.initial_indent = indent
        strs.text_wrapper.subsequent_indent = indent

        # List of lines wrapped from this text.
        lines = []

        # Iteratively construct this list.
        for paragraph_wrapped_lines in self._paragraphs_wrapped_lines(text):
            lines.extend(paragraph_wrapped_lines)

        # Return this list.
        return lines


    @type_check
    def _paragraphs_wrapped_lines(self, text: str) -> GeneratorType:
        '''
        Generator yielding a list of all lines comprising the next paragraph
        wrapped from the passed text.
        '''

        # Number of spaces indenting the first line of this text.
        main_indent_len = _get_text_indent_len(text)

        # List of lines accumulating the current paragraph. For efficiency,
        # this would ideally be a StringIO-based string buffer rather than a
        # list; unfortunately, the former fails to provide a sane method for
        # testing whether this buffer is empty (i.e., has been written to),
        # whereas the latter implicitly provides such functionality.
        paragraph_lines = []

        # For each line split from this text...
        for line in text.splitlines():
            # If this line is a blank line, pretend this line consisted of a ";"
            # instead to force this line to stop wrapping the currently
            # accumulated paragraph if any.
            if _BLANK_LINE_REGEX.match(line):
                line = ';'
            assert len(line), 'Line empty.'

            # If this line is prefixed by ";", prevent this line from being
            # wrapped by yielding such line as a newly unwrapped paragraph.
            #
            # Since the prior conditional guarantees this line to now be
            # nonempty, the first character of this line is safely indexable.
            if line[0] == ';':
                # print('Non-wrappable line detected: ' + line)

                # Strip this ";".
                line = line[1:]

                # If a previously accumulated paragraph exists, yield this
                # paragraph *BEFORE* yielding this line.
                if paragraph_lines:
                    yield self._wrap_lines(paragraph_lines, main_indent_len)

                    # Clear such paragraph.
                    paragraph_lines = []

                # Yield such line as a discrete paragraph.
                yield [line]
            # Else, wrap this line. Specifically, accumulate this line for
            # subsequent wrapping and yielding.
            else:
                # print('Wrappable line detected: ' + line)
                paragraph_lines.append(line)

        # If a previously accumulated paragraph exists, yield this as the last.
        if paragraph_lines:
            yield self._wrap_lines(paragraph_lines, main_indent_len)


    @type_check
    def _wrap_lines(self, lines: SequenceTypes, main_indent_len: int) -> list:
        '''
        Join the passed list of lines on newline, reduce two or more consecutive
        whitespace characters in the passed string to one space, strip leading
        and trailing whitespace, wrap the result, and return a wrapped list of
        lines.

        Since the base class performs the same reduction and stripping of
        whitespace, so do we.
        '''

        # Difference in spaces between the indentation of the first line of the
        # entire text passed above and the indentation of the first line of the
        # currently passed paragraph.
        paragraph_indent_len = _get_text_indent_len(lines[0]) - main_indent_len

        # String of spaces of this length.
        paragraph_indent = ' ' * paragraph_indent_len

        # Indent each wrapped line by this difference.
        strs.text_wrapper.initial_indent    = paragraph_indent
        strs.text_wrapper.subsequent_indent = paragraph_indent

        # Strip and wrap these lines.
        return strs.text_wrapper.wrap(
            self._whitespace_matcher.sub(
                ' ', strs.join_on_newline(lines)).strip())

# ....................{ GETTERS                            }....................
@type_check
def _get_text_indent_len(text: str) -> int:
    '''
    Number of spaces indenting the first line of the passed text.
    '''

    return len(_INDENTATION_REGEX.match(text).group(1))
