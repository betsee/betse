#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific animation classes.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionFile
from betse.lib.matplotlib import matplotlibs
from betse.util.path import dirs, paths
# from betse.util.type import types
from matplotlib.animation import writers, FileMovieWriter

# ....................{ CLASSES                            }....................
@writers.register('frame')
class FileFrameWriter(FileMovieWriter):
    '''
    Object serializing each frame of an animation to a unique image file.

    For drop-in use as an animation writer (e.g., to the `Animation.save()`
    method), this class masquerades as a `MovieWriter` subclass by the unique
    name of `frame`. Nonetheless, no movie is written; only frames are written.
    '''


    supported_formats = matplotlibs.get_backend_figure_filetypes()
    '''
    List of all image filetypes supported by this class.

    Since this class serializes frames by calling the `savefig()` function _and_
    since that function supports all image filetypes supported by the current
    backend, the conclusion syllogistically follows. **Q.E.D., yo.**
    '''


    def setup(self, *args, **kwargs) -> None:
        '''
        Prepare to serialize animation frames.

        This method is implicitly called by the superclass `saving()` method
        implicitly called by the `Animation.save()` method. Note that,
        unfortunately, the design of both methods prohibits this method from
        accepting subclass-specific parameters.

        Parameters
        -----------
        frame_prefix : str
            **Ignored.**
        outfile : str
            `str.format()`-formatted template describing the filenames of the
            resulting frame images. This template is subject to the following
            `str.format()` formatting:
            * `{frame_number}` is replaced by the 0-based index of the current
              frame. If this substring is _not_ found, an exception will be
              thrown.

        See the superclass `setup()` method for details on all other parameters.
        '''
        super().setup(*args, **kwargs)

        # If this filename template is malformed, raise an exception.
        if '{frame_number}' not in self.outfile:
            raise BetseExceptionFile(
                'Frame filename template "{}" contains no '
                'substring "{{frame_number}}".'.format(self.outfile))

        # Filetype of all output files.
        out_filetype = paths.get_filetype(self.outfile)

        # If this filetype is unsupported, raise an exception.
        if out_filetype not in self.supported_formats:
            raise BetseExceptionFile(
                'Frame filetype "{}" unsupported by the '
                'current Matplotlib backend (i.e., not in "{}").'.format(
                    out_filetype, str(self.supported_formats)))

        # Prevent the superclass from deleting output files. This may or may not
        # actually be needed, but... hopefully never hurts.
        self.clear_temp = False

        # Parent directory of all output files.
        out_dirname = paths.get_dirname(self.outfile)

        # Create this directory if needed.
        dirs.make_unless_parent_dir(out_dirname)


    #FIXME: Implement me. This should be dramatically simpler than the
    #superclass implementation. Nonetheless, see that for inspiration.
    def grab_frame(self, **savefig_kwargs):
        pass


    def finish(self):
        '''
        Prevent the superclass `finish()` method from attempting to fork an
        external process running a non-existent video encoding command.
        '''
        pass

# --------------------( WASTELANDS                         )--------------------
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
