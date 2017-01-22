#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific classes writing animations as frames.
'''

#FIXME: Consider contributing most or all of this submodule back to matplotlib.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMatplotlibException
from betse.util.path import dirs, paths
from matplotlib import verbose
from matplotlib.animation import writers, MovieWriter

# ....................{ CLASSES                            }....................
#FIXME: Rename to "NoopMovieWriter".
@writers.register('noop')
class NoopWriter(MovieWriter):
    '''
    Matplotlib animation writer ignoring rather than writing animation frames.

    For drop-in use as an animation writer (e.g., to the `Animation.save()`
    method), this class masquerades as a `MovieWriter` subclass by the unique
    name of `noop`. Nonetheless, no movies or files are written. This class is
    principally intended for use by callers:

    * Manually handling writing (e.g., by manually instantiating and calling
      methods of `MovieWriter` subclasses).
    * Calling the `Animation.save()` method, which obstructs manual handling of
      writing by:
      * Mandating that a writer be passed.
      * Automatically writing with that writer.

    Such callers may continue to manually handle writing by passing instances of
    this subclass to the `Animation.save()` method.
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
        Initialize this writer.

        Specifically, the passed `extra_args` parameter if any is coerced to the
        empty list. Since this subclass forks no external process to write
        animation frames, all arguments to be passed to that process are
        ignorable.
        '''

        # Ignoring "extra_args" is essential. Failing to do so results in
        # exceptions in the superclass constructor, which when "extra_args" is
        # None attempts to access the non-existent "args_key" class attribute.
        kwargs['extra_args'] = list()
        super().__init__(*args, **kwargs)

    # ..................{ IGNORE                             }..................
    def grab_frame(self, **kwargs) -> None:
        '''
        Prevent the superclass `grab_frame()` method from writing this frame to
        the current sink for this writer.

        This writer fails to redefine the `_frame_sink()` method and hence
        defaults to the sink provided by the superclass: stdin. Since stdin is
        typically an unsafe sink, this method circumvents the need to define a
        safe sink by reducing this method to a noop.
        '''
        pass


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


#FIXME: Rename to "ImageMovieWriter".

# This writer shares sufficiently many similarities with the "noop" writer to
# warrant inheriting that writer.
@writers.register('image')
class ImageWriter(NoopWriter):
    '''
    Matplotlib animation writer writing animation frames to still image files
    rather than animated videos.

    For drop-in use as an animation writer (e.g., to the `Animation.save()`
    method), this class masquerades as a `MovieWriter` subclass by the unique
    name of `image`. Nonetheless, no movie is written; only frames are written.

    Attributes
    -----------
    _frame_number : int
        0-based index of the next frame to be written.
    '''


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
        self.frame_format = paths.get_filetype_undotted_or_none(self.outfile)

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
            'ImageWriter.grab_frame: saving frame %d to fname=%s' % (
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
