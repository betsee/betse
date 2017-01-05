#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: It'd be great to contribute this back to the official "argparse" module
#by contributing a patch. See:
#    http://bugs.python.org/issue12806

'''
Low-level command-line argument facilities.
'''

# ....................{ IMPORTS                            }....................
from argparse import HelpFormatter
from betse.util.type import strs
import re

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
class HelpFormatterParagraph(HelpFormatter):
    '''
    Formatter wrapping all lines *not* prefixed by `;` in `ArgumentParser`
    help text into paragraphs.

    Newlines in all lines prefixed by `;` will be preserved as is, preventing
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
        Get the passed string with all non-indented lines wrapped into
        paragraphs.
        '''
        return '\n'.join(self._split_lines_indented(text, width, indent))

    def _split_lines(self, text: str, width: int) -> list:
        '''
        Get a list of lines produced by wrapping all non-indented lines in the
        passed text into paragraphs and splitting such string on newlines.
        '''
        return self._split_lines_indented(text, width)

    # ..................{ PRIVATE                            }..................
    def _split_lines_indented(self, text, width, indent = '') -> list:
        '''
        Get a list of lines produced by wrapping all non-indented lines in the
        passed text into paragraphs and splitting such string on newlines.
        '''
        # Initialize such text wrapper with the passed arguments.
        strs.text_wrapper.width = width
        strs.text_wrapper.initial_indent = indent
        strs.text_wrapper.subsequent_indent = indent

        # List of lines wrapped from such text.
        lines = []

        # Iteratively construct such list.
        for paragraph_wrapped_lines in self._paragraphs_wrapped_lines(text):
            lines.extend(paragraph_wrapped_lines)

        return lines

    def _paragraphs_wrapped_lines(self, text: str) -> tuple:
        '''
        Generator yielding a list of all lines comprising the next paragraph
        wrapped from the passed text.
        '''
        assert isinstance(text, str), '"{}" not a string.'.format(text)

        # Number of spaces indenting the first line of such text.
        main_indent_len = _get_text_indent_len(text)

        # List of lines accumulating the current paragraph. For efficiency,
        # this would ideally be a StringIO-based string buffer rather than a
        # list; unfortunately, the former fails to provide a sane method for
        # testing whether such buffer is empty (i.e., has been written to),
        # whereas the latter implicitly provides such functionality.
        paragraph_lines = []

        # For each line split from such text...
        for line in text.splitlines():
            # If such line is a blank line, pretend such line consisted of a ";"
            # instead to force such line to stop wrapping the currently
            # accumulated paragraph if any.
            if _BLANK_LINE_REGEX.match(line):
                line = ';'
            assert len(line), 'Line empty.'

            # If such line is prefixed by ";", prevent such line from being
            # wrapped by yielding such line as a newly unwrapped paragraph.
            #
            # Since the prior conditional guarantees such line to now be
            # nonempty, the first character of such line is safely indexable.
            if line[0] == ';':
                # Strip such ";".
                line = line[1:]

                # If a previously accumulated paragraph exists, yield such
                # paragraph *BEFORE* yielding such line.
                if paragraph_lines:
                    yield self._wrap_lines(paragraph_lines, main_indent_len)

                    # Clear such paragraph.
                    paragraph_lines = []

                # Yield such line as a discrete paragraph.
                yield [line]
            # Else, wrap such line. Specifically, accumulate such line for
            # subsequent wrapping and yielding.
            else:
                paragraph_lines.append(line)

        # If a previously accumulated paragraph exists, yield this as the last.
        if paragraph_lines:
            yield self._wrap_lines(paragraph_lines, main_indent_len)

    def _wrap_lines(self, lines: list, main_indent_len: int) -> list:
        '''
        Join the passed list of lines on newline, reduce two or more consecutive
        whitespace characters in the passed string to one space, strip leading
        and trailing whitespace, wrap the result, and return a wrapped list of
        lines.

        Since the base class performs the same reduction and stripping of
        whitespace, so do we.
        '''
        assert isinstance(lines, list), '"{}" not a list.'.format(lines)
        assert len(lines), 'List of lines empty.'
        assert isinstance(main_indent_len, int),\
            '"{}" not an integer.'.format(main_indent_len)

        # Difference in spaces between the indentation of the first line of the
        # entire text passed above and the indentation of the first line of the
        # currently passed paragraph.
        paragraph_indent_len = _get_text_indent_len(lines[0]) - main_indent_len

        # String of spaces of such length.
        paragraph_indent = ' ' * paragraph_indent_len

        # Indent each wrapped line by such difference.
        strs.text_wrapper.initial_indent    = paragraph_indent
        strs.text_wrapper.subsequent_indent = paragraph_indent

        # Strip and wrap such lines.
        return strs.text_wrapper.wrap(
            self._whitespace_matcher.sub(
                ' ', strs.join_on_newline(lines)).strip())

# ....................{ GETTERS                            }....................
def _get_text_indent_len(text: str) -> int:
    '''
    Get the number of spaces indenting the first line of the passed text.
    '''
    return len(_INDENTATION_REGEX.match(text).group(1))

# --------------------( WASTELANDS                         )--------------------
        # # True if the current line is to be wrapped.
        # is_line_wrapping = False
        #
        # # For each line split from such text...
        # for line in text.splitlines():
        #     # Wrap lines by default.
        #     is_line_wrapping = True
        #
        #     # If such line is prefixed by ";", strip such ";" and prevent
        #     # such line from being wrapped.
        #     if line[0] == ';':
        #         line = line[1:]
        #         is_line_wrapping = False
        #
        #     # If wrapping such line, accumulate such line for subsequent
        #     # wrapping and yielding.
        #     if is_line_wrapping:

#FUXME: Refactor according to the new syntax.
                    #FUXME: Both here and below, we need to prefix each
                    #indented line by the indentation difference: e.g.,
                    # _get_text_indent_len(line) - main_indent_len

                # # If such line is prefixed by ":" followed by one or more
                # # spaces, strip such ":" but *NOT* spaces and prevent such line
                # # from being wrapped.
                # if _COLON_INDENTATION_REGEX.match(line):
                #     line = line[1:]
                #     is_line_wrapping = False
                # # If such line is a blank line or has indentation exceeding that
                # # of the first line, strip such indentation and prevent such
                # # line from being wrapped.
                # elif _BLANK_LINE_REGEX.match(line) or\
                #      _get_text_indent_len(line) > main_indent_len:
                #     line = line.strip()
                #     is_line_wrapping = False
                #
                # # If wrapping such line, accumulate such line for subsequent
                # # wrapping and yielding.
                # if is_line_wrapping:
                #     paragraph_lines.append(line)
                # # Else, yield such line as an unwrapped paragraph.
                # else:
                #     # If a previously accumulated paragraph exists, yield such
                #     # paragraph *BEFORE* yielding such line.
                #     if paragraph_lines:
                #         yield self._wrap_lines(paragraph_lines)
                #
                #         # Clear such paragraph.
                #         paragraph_lines = []
                #
                #     # Yield such line as a discrete paragraph.
                #     yield [line]

# _COLON_INDENTATION_REGEX = re.compile(r'^(: *)')
# '''
# Regular expression capturing a `:` followed by one or more spaces prefixing the
# subject string.
# '''

# Lines prefixed by `:` followed by zero or more spaces will
    # *not* be wrapped, thereby preserving both indentation and newlines in such
    # lines.
    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     strs.text_wrapper = TextWrapper()

    # Attributes
    # ----------
    # _text_wrapper : TextWrapper
    #     Object with which to wrap text.

                    # line = [main_indent_len:]
        #     # For each line split from such text...
        #     for line in text.splitlines():
        #         # If such line is a blank line or has indentation exceeding that
        #         # of the first line, yield only such line.
        #         if _BLANK_LINE_REGEX.match(line) or\
        #            _get_text_indent_len(line) > main_indent_len:
        #             # If a previously accumulated paragraph exists, yield such
        #             # paragraph *BEFORE* yielding such line.
        #             if paragraph_lines:
        #                 yield '\n'.join(paragraph_lines), True
        #
        #                 # Clear such list of lines.
        #                 paragraph_lines = []
        #
        #             # Yield such line as a discrete paragraph.
        #             yield line, False
        #         # Else, append such line to such list of lines.
        #         else:
        #             paragraph_lines.append(line)
        #
        #     # If a previously accumulated paragraph exists, yield such
        #     # paragraph (as above).
        #     if paragraph_lines:
        #         yield '\n'.join(paragraph_lines), True
        #
        # # List of lines wrapped from such text.
        # lines = []
        #
        # # For each paragraph parsed from such text...
        # for paragraph, is_wrapping in _paragraphs():
        #     # If wrapping such paragraph, do so.
        #     if is_wrapping:
        #         # Reduce two or more consecutive whitespace characters to one
        #         # space, strip leading whitespace, wrap the result, and append
        #         # the resulting lines. Our base class performs the same
        #         # reduction and stripping of whitespace; so, so do we.
        #         lines.extend(
        #             strs.text_wrapper.wrap(
        #                 self._whitespace_matcher.sub(' ', paragraph).strip()))
        #     # Else, such paragraph *MUST* be a single line and hence is
        #     # appendable as is after stripping prefixing indentation.
        #     else:
        #         lines.append(paragraph[main_indent_len:])

# from io import StringIO
            # String buffer accumulating the current paragraph.
            # paragraph = StringIO()
            # True if such paragraph exists (i.e., if such buffer has been
            # written to at least once but *NOT* yet yielded).
            # is_paragraph = False
                        # Technically, StringIO() objects are also clearable
                        # without recreating such ojbects as follows:
                        #
                        #     paragraph.truncate(0)
                        #     paragraph.seek(0)
                        #
                        # That said, the top answer at the following
                        # stackoverflow thread shows the current approach to be
                        # both simpler and (surprisingly) more efficient:
                        #
                        #     https://stackoverflow.com/questions/4330812/how-do-i-clear-a-stringio-object
                        # paragraph = StringIO()
                        # is_paragraph = False
