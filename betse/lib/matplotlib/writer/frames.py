#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific classes writing animations as frames.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMatplotlibException
from betse.util.path import dirs, paths
from matplotlib import verbose
from matplotlib.animation import writers, MovieWriter

# ....................{ CLASSES                            }....................
# Subclass the base "MovieWriter" superclass rather than the child
# "FileMovieWriter" class (itself subclassing the former). Subclassing
# "FileMovieWriter" introduces considerable complexity for little gain.
@writers.register('frame')
class FileFrameWriter(MovieWriter):
    '''
    Object serializing each frame of an animation to a unique image file.

    For drop-in use as an animation writer (e.g., to the `Anim.save()`
    method), this class masquerades as a `MovieWriter` subclass by the unique
    name of `frame`. Nonetheless, no movie is written; only frames are written.

    Attributes
    -----------
    _frame_number : int
        0-based index of the next frame to be written.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @classmethod
    def isAvailable(cls):
        '''
        Notify users that this class does _not_ depend on external commands and
        hence is _always_ available.
        '''
        return True

    # ..................{ PUBLIC                             }..................
    def __init__(self, *args, **kwargs) -> None:
        '''
        Construct a new `FileFrameWriter` object.

        This constructor coerces the passed `extra_args` parameter if any to the
        empty list. Since no external process is forked to write animation
        frames, all arguments to be passed to that process are ignorable.
        '''

        # Ignoring "extra_args" is essential. Failing to do so results in
        # exceptions in the superclass constructor, which when "extra_args" is
        # None attempts to access the non-existent "args_key" class attribute.
        kwargs['extra_args'] = list()
        super().__init__(*args, **kwargs)


    def setup(self, *args, **kwargs) -> None:
        '''
        Prepare to write animation frames.

        This method is implicitly called by the superclass `saving()` method
        implicitly called by the `Anim.save()` method. Note that,
        unfortunately, the design of both methods prohibits this method from
        accepting subclass-specific parameters.

        Parameters
        -----------
        outfile : str
            `str.format()`-formatted template describing the filenames of the
            resulting frame images. This template is subject to the following
            `str.format()` formatting:
            * The first `{`- and `}`-delimited substring (e.g., `{:07d}`) will
              be replaced by the 0-based index of the current frame. If this
              substring does _not_ exist, an exception is raised.

        See the superclass `setup()` method for details on all other parameters.

        See Also
        -----------
        https://www.python.org/dev/peps/pep-3101
            For details on format specifiers.
        '''

        super().setup(*args, **kwargs)

        # Frame number of the next frame to be written.
        self._frame_number = 0

        # If this filename template is malformed, raise an exception.
        if not ('{' in self.outfile and '}' in self.outfile):
            raise BetseMatplotlibException(
                'Frame filename template "{}" contains no "{{"- and "}}"-'
                'delimited format specifier.'.format(self.outfile))

        # Output filetype. Override the superclass' awkward choice of "rgba" as
        # output filetype default.
        self.frame_format = paths.get_filetype(self.outfile)

        # List of all output filetypes supported by this class.
        #
        # Since this class serializes frames by calling the savefig() function
        # *AND* since that function supports all image filetypes supported by
        # the current backend, this syllogistically follows. (Q.E.D.)
        out_filetypes_supported = self.fig.canvas.get_supported_filetypes()

        # If this filetype is unsupported, raise an exception.
        if self.frame_format not in out_filetypes_supported:
            raise BetseMatplotlibException(
                'Frame filetype "{}" unsupported by the '
                'current Matplotlib backend (i.e., not in "{}").'.format(
                    self.frame_format, str(out_filetypes_supported)))

        # Parent directory of all output files.
        out_dirname = paths.get_dirname(self.outfile)

        # Create this directory if needed.
        dirs.make_parent_unless_dir(out_dirname)


    def grab_frame(self, **kwargs) -> None:
        '''
        Write the next frame for the current figure to the image file defined by
        the current filename template.

        The high-level `Animation.save()` method unpacks and passes all
        key-value pairs of the optional `savefig_kwargs` dictionary argument to
        this method as keyword arguments. This method then passes all passed
        keyword arguments to the `Figure.savefig()` method.
        '''

        # Leverage similar code as our superclass with exception of:
        #
        # * *NOT* catching and logging exceptions. The superclass implementation
        #   catches and logs exceptions by deferring to the output of an
        #   external process no longer forked by this subclass.
        # * Passing a filename rather than open file handle to the
        #   self.fig.savefig() method. While we could reimplement the superclass
        #   _frame_sink() method to do so, this seems simpler and more reliable.

        # Filename of the current frame to be written.
        frame_filename = self.outfile.format(self._frame_number)

        # Log this attempt via Matplotlib's logging API.
        verbose.report(
            'FileFrameWriter.grab_frame: saving frame %d to fname=%s' % (
                self._frame_number, frame_filename),
            level='debug')

        # Increment the number of the next frame to be written *AFTER* logging.
        self._frame_number += 1

        # Write the current frame.
        self.fig.savefig(
            filename=frame_filename,
            format=self.frame_format,
            dpi=self.dpi,
            **kwargs
        )

    # ..................{ IGNORE                             }..................
    def cleanup(self):
        '''
        Prevent the superclass `cleanup()` method from attempting to capture
        output from an external process no longer forked by this subclass.
        '''
        pass


    def _run(self) -> None:
        '''
        Prevent the superclass `_run()` method from attempting to fork an
        external process running a non-existent video encoding command.
        '''
        pass

# --------------------( WASTELANDS                         )--------------------
# from betse.util.type import types
# from betse.lib.matplotlib import mpl
        # # List of all output filetypes supported by this class.
        # #
        # # Since this class serializes frames by calling the savefig() function
        # # *AND* since that function supports all image filetypes supported by
        # # the current backend, the conclusion syllogistically follows. (Q.E.D.)
        # #
        # # This list is intentionally *NOT* declared as the standard Matplotlib
        # # animation writer class attribute "supported_formats". Doing so would
        # # require calling the following function and hence importing the
        # # "matplotlib.pyplot" module at the time of this module's importation,
        # # substantially complicating Matplotlib use. (Think backends.)
        # out_filetypes_supported = mpl.backend_figure_filetypes()
        #
        # # If this filetype is unsupported, raise an exception.
        # if self.frame_format not in out_filetypes_supported:
        #     raise BetseMatplotlibException(
        #         'Frame filetype "{}" unsupported by the '
        #         'current Matplotlib backend (i.e., not in "{}").'.format(
        #             self.frame_format, str(out_filetypes_supported)))

        # frame_prefix : str
        #     **Ignored.**
        # Prevent the superclass from deleting output files. This may or may not
        # actually be needed, but... hopefully never hurts.
        # self.clear_temp = False

    # def finish(self):
    #     '''
    #     Prevent the superclass `finish()` method from attempting to fork an
    #     external process running a non-existent video encoding command.
    #     '''
    #     pass

    # supported_formats = matplotlibs.get_backend_figure_filetypes()
    # '''
    # List of all image filetypes supported by this class.
    #
    # Since this class serializes frames by calling the `savefig()` function _and_
    # since that function supports all image filetypes supported by the current
    # backend, the conclusion syllogistically follows. **Q.E.D., yo.**
    # '''

        #FUXME: Is this actually needed?
        # Type of all output files.
        # self.frame_format = paths.get_filetype(self.outfile)

        #FUXME: Is this actually needed?
        # Parent directory of all output files. The oddball attribute name
        # satisfies equally oddball superclass requirements. Ugh, you!
        # self.temp_prefix = paths.get_dirname(self.outfile)

    # def __init__(self, *args, **kwargs):
    #     FileMovieWriter.__init__(self, *args, **kwargs)

        # Create the parent directory of output files if needed.
        # dirs.make_unless_parent_dir(self.outfile)

# from collections import namedtuple
    # def _run(self):
    #     '''
    #     Notify our superclass that this method succeeded by dynamically creating
    #     and instantiating a fake class providing the required attribute
    #     subsequently tested by the superclass `finish()` method.
    #
    #     That method expects this method to add a private `_proc` attribute to
    #     the current object whose public `returncode` attribute provides the exit
    #     status of the presumably run video encoding command. Since serializing
    #     frames requires no such command, this method coerces that attribute to
    #     be the standard exit status for success (i.e., `0`).
    #
    #     Welcome to Matplotlib Hell. We hope you enjoy your hellish stay.
    #     '''
    #     # Interestingly, note that the standard "errno" module provides no
    #     # integer constant for success. In theory, even Windows should use 0 to
    #     # denote success. (Let's hope Bill didn't screw that pooch too.)
    #     FakeProc = namedtuple('FakeProc', 'returncode')
    #     self._proc = FakeProc(0)
