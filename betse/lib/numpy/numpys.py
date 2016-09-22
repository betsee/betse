#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Numpy, a mandatory runtime dependency.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.util.py import modules
from betse.util.os import libs, oses
from betse.util.type import iterables, regexes, strs
from betse.util.type.types import NoneType
from collections import OrderedDict

# Note that the optionally installed "numpy.distutils" subpackage, if available,
# contains the same "__config__" submodule. Since this subpackage is only
# optional, however, this submodule is imported from the top-level mandatory
# "numpy" package instead.
from numpy import __config__ as numpy_config

# ....................{ GLOBALS ~ opt_info                 }....................
# Fully initialized by the _init_globals() function below.
_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX = None
'''
Uncompiled regular expression heuristically matching the basenames of shared
libraries providing optimized BLAS shared libraries in the `libraries` list of
the global :data:`numpy.__config__.blas_opt_info` dictionary.

This expression does _not_ match the strict superset of optimized BLAS shared
libraries that are also optimized, as doing so in a cross-platform manner is
infeasible. Debian-based Linux distributions, for example, remove the prefixes
uniquely identifying the threaded variants of both ATLAS and OpenBLAS from the
basenames of their shared libraries (e.g., "libatlas.so" rather than
"libtatlas.so").

Optimized BLAS shared libraries matched by this expression include:

* **AMD Core Math Library (ACML), unconditionally GPU- and CPU- optimized in
  both OpenMP-based and non-OpenMP-based variants regardless of underlying
  compiler (e.g., GNU Fortran, Open64). Note, however, that ACML does _not_ ship
  with a CBLAS interface and hence is non-trivial to link Numpy against. While
  unlikely that any end users will ever do so, it nonetheless remains feasible
  and hence somewhat supported.
* **Automatically Tuned Linear Algebra Software (ATLAS),** both single- and
  multithreaded CBLAS- and Fortran-based variants for both ATLAS < 3.10 and
  ATLAS >= 3.10, which ships shared libraries under different basenames than
  ATLAS < 3.10. (Life complicates life.)
* **Intel Math Kernel Library (MKL),** unconditionally multithreaded in both
  OpenMP-based and non-OpenMP-based variants regardless of underlying compiler
  (e.g., dynamic, GCC, Intel).  Note that **Intel Vector Mathematical Library
  (VML)** is intentionally ignored.  Although also unconditionally
  multithreaded, VML does _not_ implement the BLAS API. Numpy currently contains
  no VML-specific handling, apart from (somewhat uselessly) detecting VML
  installation on reporting system diagnostics.
* **OpenBLAS,** both single- and multithreaded 32- and 64-bit variants.

This expression is typically only required once at application startup and hence
is conditionally compiled in a just-in-time (JIT) manner by the
:func:`_is_blas_optimized_python_general` function rather than
unconditionally compiled here.
'''
# print('blas regex: ' + _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX)


#FIXME: Actually use this, please.

# Fully initialized by the _init_globals() function below.
_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRNAME_REGEX = None
'''
Uncompiled regular expression heuristically matching the dirnames of shared
libraries providing optimized BLAS shared libraries in the `libraries` list of
the global :data:`numpy.__config__.blas_opt_info` dictionary.

See Also
----------
:data:`_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX`
    Further details.
'''


_OPTIMIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_OS_X = {
    # Accelerate. Although Accelerate is only conditionally multithreaded,
    # multithreading is enabled by default and hence a safe assumption.
    '-Wl,Accelerate',

    # vecLib. Similar assumptions as with Accelerate apply.
    '-Wl,vecLib',
}
'''
Set of all strings in the `extra_link_args` list of the global
:data:`numpy.__config__.blas_opt_info` dictionary heuristically corresponding to
optimized BLAS implementations under OS X.

Unlike all other BLAS implementations, Numpy does _not_ declare unique
dictionary globals describing these implementations when linked against. Ergo,
this lower-level solution.
'''

# ....................{ GLOBALS ~ filename                 }....................
_OPTIMIZED_BLAS_FILENAME_REGEX = r'^({}).*$'.format(
    r'|'.join((
        # AMD Core Math Library (ACML).
        r'acml',

        # Automatically Tuned Linear Algebra Software (ATLAS) >= 3.10.
        r'tatlas',

        # Automatically Tuned Linear Algebra Software (ATLAS) < 3.10.
        r'pt(c|f77)blas',

        # Intel Math Kernel Library (MKL). Thanks to the profusion of possible
        # library basenames, this regular expression fragment simplistically
        # assumes *ALL* library basenames prefixed by "mkl" to unconditionally
        # connote MKL. What could go wrong?
        r'mkl_',

        # OpenBLAS.
        r'openblas(_[^_]+)?_threads',

        #FIXME: Research. No idea if this basename substring is even accurate.

        # BLAS-like Library Instantiation Software (BLIS). Unconditionally
        # multithreaded. Technically, Numpy has yet to add official support for
        # BLIS. Since numerous contributors nonetheless perceive BLIS to be the
        # eventual successor of BLAS *AND* since Numpy currently hosts an open
        # pull request to explicitly add BLIS support under the sensible
        # subclass name "blis_info" (see Numpy PR #7294), explicitly listing
        # BLIS here should assist in future-proofing our multithreading
        # detection.
        #r'blis',
    ))
)
'''

This expression is typically only required once at application startup and hence
is conditionally compiled in a just-in-time (JIT) manner by the
:func:`_is_blas_optimized_linkage` function rather than unconditionally
compiled here.
'''

# ....................{ INITIALIZERS                       }....................
# For simplicity, this function is called below on the first importation of this
# submodule rather than explicitly called by callers.
def init() -> None:
    '''
    Initialize Numpy.

    Specifically:

    * If the currently installed version of Numpy was linked against an
      unoptimized BLAS implementation and is thus itself unoptimized, log a
      non-fatal warning.
    '''

    #FIXME: Excise this after is_blas_optimized() is known to reliably work
    #in a cross-platform manner.
    return

    _init_globals()

    # If Numpy linked against an unoptimized BLAS, log a non-fatal warning.
    if not is_blas_optimized():
        logs.log_warning(
            'Numpy unoptimized. '
            'Consider installing an optimized CBLAS implementation '
            '(e.g., OpenBLAS, ATLAS, ACML, MKL) and '
            'reinstalling Numpy against this implementation.'
        )


def _init_globals() -> None:
    '''
    Initialize all uninitialized global variables of this submodule.
    '''

    # Permit these globals to be redefined.
    global\
        _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX,\
        _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRNAME_REGEX

    # Redefine this global.
    _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX = (
        r'^({})(?:[_-].*)?$'.format(r'|'.join((
            # AMD Core Math Library (ACML).
            r'acml',

            # Automatically Tuned Linear Algebra Software (ATLAS) >= 3.10.
            r'[st]?atlas',

            # Automatically Tuned Linear Algebra Software (ATLAS) < 3.10.
            r'(?:pt)?f77blas',

            # Automatically Tuned Linear Algebra Software (ATLAS) < 3.10.  Although
            # some platforms (e.g., Ubuntu) distribute the CBLAS implementation of
            # ATLAS with the ambiguous basename of "cblas" rather than "ptcblas",
            # the former may also refer to the unoptimized reference CBLAS
            # implementation and is hence ignored.
            r'ptcblas',

            # Intel Math Kernel Library (MKL). Thanks to the profusion of possible
            # library basenames, this regular expression fragment simplistically
            # assumes *ALL* library basenames prefixed by "mkl" to unconditionally
            # connote MKL. What could go wrong?
            r'mkl',

            # OpenBLAS.
            r'openblas',

            #FIXME: Research. No idea if this basename substring is even accurate.

            # BLAS-like Library Instantiation Software (BLIS). Unconditionally
            # multithreaded. Technically, Numpy has yet to add official support for
            # BLIS. Since numerous contributors nonetheless perceive BLIS to be the
            # eventual successor of BLAS *AND* since Numpy currently hosts an open
            # pull request to explicitly add BLIS support under the sensible
            # subclass name "blis_info" (see Numpy PR #7294), explicitly listing
            # BLIS here should assist in future-proofing our multithreading
            # detection.
            #r'blis',
    ))))

    # Regular expression fragment matching the boundary of a dirname at which a
    # substring matching an optimized BLAS name may either begin or end.
    dirname_boundary_regex = r'[{}_.-])'.format(regexes.DIR_SEPARATOR_REGEX)

    # Redefine this global.
    _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRNAME_REGEX = (
        r'^(?:{boundary})?({name})(?:{boundary}.*)?$'.format(
            boundary=dirname_boundary_regex,
            name=r'|'.join((
                # AMD Core Math Library (ACML).
                r'acml',

                # Automatically Tuned Linear Algebra Software (ATLAS).
                r'atlas',

                # Intel Math Kernel Library (MKL).
                r'mkl',

                # OpenBLAS.
                r'openblas',
    ))))

# ....................{ TESTERS                            }....................
#FIXME: Revise docstring, which is pretty much completely wrong now.

def is_blas_optimized() -> bool:
    '''
    `True` only if the currently installed version of Numpy is linked against an
    optimized BLAS (Basic Linear Algebra Subprograms) implementation, ideally
    but _not_ necessarily parallelized across multiple processors.

    Optimized BLAS implementations are _strongly_ recommended over unoptimized
    BLAS implementations. The `numpy.dot()` operator, which is implicitly
    optimized when Numpy is linked against a optimized BLAS implementation, is
    frequently called by BETSE in its critical path.

    Note that testing for parallelized optimized BLAS implementations, while
    more specific and hence preferable, is infeasible for common edge-cases
    (e.g., Debian-based Linux distributions). For further details, see the
    :data:`_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX` string global.

    Heuristic
    ----------
    Numpy does _not_ provide a convenient API for readily querying this boolean.
    Numpy does, however, provide an admittedly inconvenient API for aggregating
    this boolean together from various sources: the
    :mod:`numpy.__config__` submodule. The
    :func:`numpy.distutils.misc_util.generate_config_py` function
    programmatically fabricates the contents of the :mod:`numpy.__config__`
    submodule at Numpy installation time.

    Specifically, for each subclass of the
    :class:`numpy.distutils.system_info.system_info` base class defined by the
    :mod:`numpy.distutils.system_info` submodule (e.g.,
    :class:`numpy.distutils.system_info.atlas_info`) whose corresponding shared
    library (e.g., ATLAS_) or feature (e.g., ATLAS_ multithreading) is available
    on the current system at Numpy installation time, a dictionary global of the
    same name as that subclass whose keys are the names of metadata types and
    values are metadata is programmatically added to the
    :mod:`numpy.__config__` submodule. Hence, this function returns
    `True` only if that submodule declares a global specific to a optimized
    BLAS implementation.

    .. _ATLAS: http://math-atlas.sourceforge.net
    '''

    # For each private tester implementing a heuristic for this public test (in
    # order of decreasing generality, portability, and reliability)...
    for tester_heuristic in (
        _is_blas_optimized_python_general,
        #FIXME: Uncomment these heuristics after retesting.
        # _is_blas_optimized_python_os_x,
        # _is_blas_optimized_linkage,
    ):
        # Attempt to...
        try:
            # Call this tester, capturing the result for subsequent handling.
            tester_result = tester_heuristic()

            # If this tester definitively identified Numpy as either
            # optimized or non-optimized, return this result.
            if tester_result is not None:
                return tester_result
            # Else, continue to the next tester.
        # If an error occurs, log that error *WITHOUT* raising an exception.
        # Detecting Numpy optimization is non-essential and hence hardly worth
        # halting the application over.
        except Exception as exception:
            logs.log_exception(exception)

    # Else, all heuristics failed to definitively identify Numpy to be either
    # optimized or non-optimized. For safety, assume the latter.
    return False


#FIXME: Document us up.
def _is_blas_optimized_python_general() -> (bool, NoneType):
    '''
    The `libraries` list of the global
    :data:`numpy.__config__.blas_opt_info` dictionary.
    '''

    # Global BLAS linkage dictionary for this Numpy installation if any or
    # "None" otherwise. Technically, this dictionary should *ALWAYS* be defined.
    # Reality probably occasionally begs to disagree, however.
    blas_lib = getattr(numpy_config, 'blas_opt_info', None)

    # If this dictionary is undefined, log a non-fatal warning and return False.
    # While unfortunate, this is *NOT* worth raising a fatal exception over.
    if blas_lib is None:
        logs.log_warning(
            'Numpy installation misconfigured: '
            '"numpy.__config__.blas_opt_info" dictionary not found.')
        return False

    # List of all substrings of BLAS library basenames this version of Numpy is
    # linked against if any or "None" otherwise.
    blas_basename_substrs = blas_lib.get('libraries', None)

    # If this list is either undefined or empty, log a non-fatal warning and
    # return False.  While unfortunate, this is *NOT* worth raising a fatal
    # exception over. (Unlike the optional "extra_link_args" metadata tested
    # for below, the "libraries" metadata is mandatory.)
    if not blas_basename_substrs:
        logs.log_warning(
            'Numpy installation misconfigured: '
            "\"numpy.__config__.blas_opt_info['libraries']\" "
            'dictionary key not found or empty.')
        return False
    # Else, this list is non-empty.

    # First element of this list. For simplicity, this function assumes the
    # BLAS library identified by this element currently exists. While
    # iteratively testing all listed BLAS libraries for existence would be
    # feasible, doing so is platform-specific and hence non-trivially fragile.
    blas_basename_substr = blas_basename_substrs[0]

    # If the BLAS library identified by this element is optimized, return
    # True. Since this element may be suffixed by non-identifying metadata
    # (e.g., version), a regular expression is leveraged.
    if regexes.is_match(
        text=blas_basename_substr,
        regex=_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX,
    ):
        return True

    # Else, instruct our caller to continue to the next heuristic.
    return None


#FIXME: Document us up.
def _is_blas_optimized_python_os_x() -> (bool, NoneType):

    # If the current platform is OS X, fallback to testing whether Numpy was
    # linked against a optimized BLAS implementation specific to OS X:
    # namely, "Accelerate" or "vecLib". Unlike all other BLAS implementations,
    # these implementations are linked against with explicit linker flags --
    # requiring further logic. For further confirmation that the
    # "numpy.__config__.blas_opt_info" dictionary gives this metadata when
    # linked against these implementations, see:
    #
    # * https://trac.macports.org/ticket/22200
    # * https://github.com/BVLC/caffe/issues/2677
    #
    # When life sells you cat food, you eat cat food.
    if oses.is_os_x():
        # List of all implementation-specific link arguments with which Numpy
        # linked against the current BLAS implementation if any or "None". Note
        # that the "blas_opt_info" dictionary global is guaranteed to exist due
        # to the previously called _is_blas_optimized_python_general()
        # function.
        blas_link_args_list = numpy_config.blas_opt_info.get(
            'extra_link_args', None)

        # If at least one such argument exists...
        if blas_link_args_list:
            # Set of these arguments, converted from this list for efficiency.
            blas_link_args = set(blas_link_args_list)

            # Subset of this set specific to multithreaded BLAS implementations.
            blas_link_args_multithreaded = (
                blas_link_args &
                _OPTIMIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_OS_X
            )

            # Return True only if this subset is nonempty.
            return len(blas_link_args_multithreaded) > 0

    # Else, instruct our caller to continue to the next heuristic.
    return None


#FIXME: Document us up.
def _is_blas_optimized_linkage() -> (bool, NoneType):

    # First element of the list of all substrings of BLAS library basenames this
    # version of Numpy is linked against. Note that this list is guaranteed to
    # be non-empty due to the previously called
    # _is_blas_optimized_python_general() function.
    blas_basename_substr = numpy_config.blas_opt_info['libraries']

    # If this appears to be either the reference BLAS or CBLAS implementations
    # *AND* this platform is POSIX-compliant and hence supports symbolic links,
    # fallback to testing whether this library is in fact a symbolic link to a
    # optimized BLAS implementation.
    #
    # Unfortunately, the "numpy.__config__" API fails to specify the absolute
    # paths of the libraries it links against. Since there exists no reliable
    # means of reverse engineering these paths from this API, these paths must
    # be obtained by another means: specifically, by querying the standard
    # "numpy.core.multiarray" C extension installed under all supported Numpy
    # for the absolute paths of all external shared libraries to which this
    # extension links -- exactly one of which is guaranteed to be the absolute
    # path of what appears to be a reference BLAS or CBLAS implementation.
    # if oses.is_posix():
    if (
        blas_basename_substr == 'blas' or
        blas_basename_substr == 'cblas'

        #FIXME: Generalize to OS X as well once the
        #libs.iter_linked_lib_filenames() function supports OS X.
    ) and oses.is_linux():
        # Absolute path of this shared library-based Numpy C extension.
        numpy_lib_filename = modules.get_filename('numpy.core.multiarray')

        # For the basename and absolute path of each shared library linked to
        # by this Numpy shared library...
        for (numpy_linked_lib_basename, numpy_linked_lib_filename) in (
            libs.iter_linked_lib_filenames(numpy_lib_filename)):

            #FIXME: Implement the following heuristic:
            #
            #* If "numpy_linked_lib_basename" sans filetype (e.g., via
            #  paths.get_pathname_sans_filetype()) is suffixed by "blas" *AND*
            #  "numpy_linked_lib_filename" is a symbolic link, then:
            #  * If the transitive target pathname of this symbolic link (e.g.,
            #    via paths.canonicalize()) matches a regular expression
            #    resembling (...but possibly differing from?)
            #    "_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX", return
            #    True.
            #* Else, continue to the next linked lib.

            # If this library is *NOT* a symbolic link to another library,
            # continue to the next library.
            if not files.is_symlink(numpy_linked_lib_filename):
                continue
            # Else, this library is actually a symbolic link to another library.

            # Basename excluding suffixing filetype of this library.
            numpy_linked_lib_rootname = (
                paths.get_pathname_sans_filetype(numpy_linked_lib_basename))

            # If this is neither the BLAS nor CBLAS reference library, continue
            # to the next library.
            if not numpy_linked_lib_rootname.endswith('blas'):
                continue
            # Else, this is either the BLAS or CBLAS reference library.

    # Else, instruct our caller to continue to the next heuristic.
    return None

# ....................{ GETTERS                            }....................
def get_metadatas() -> tuple:
    '''
    Tuple of 2-tuples `(metedata_name, metadata_value`), describing all
    currently installed third-party dependencies against which Numpy was linked
    (e.g., BLAS, LAPACK).
    '''

    #FIXME: Add LAPACK linkage metadata as well.
    return (
        ('numpy (blas)', get_blas_metadata()),
    )


def get_blas_metadata() -> OrderedDict:
    '''
    Ordered dictionary synopsizing the current Numpy installation with respect
    to BLAS linkage.
    '''

    # This dictionary.
    metadata = OrderedDict((
        ('optimized', is_blas_optimized()),
    ))

    # Set of all keys of the dictionary global synopsizing this metadata,
    # sorted in ascending lexicographic order for readability.
    blas_opt_info_keys = iterables.sort_lexicographic_ascending(
        numpy_config.blas_opt_info.keys())

    # For each such key...
    for blas_opt_info_key in blas_opt_info_keys:
        # The value of this key, unconditionally converted into a string and
        # then trimmed to a reasonable string length. The values of numerous
        # keys (e.g., "libraries", "sources") commonly exceed this length,
        # hampering readability for little to no gain. Excise them all.
        metadata[blas_opt_info_key] = strs.trim(
            obj=numpy_config.blas_opt_info[blas_opt_info_key],
            max_len=256,
        )

    # Return this dictionary.
    return metadata
