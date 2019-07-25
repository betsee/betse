#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Numpy, a mandatory runtime dependency.
'''

#FIXME: Add detection support for NVBLAS, the Nvidia GPU-specific equivalent of
#AMD's ACML. Naturally, further research is required.

#FIXME: Consider replacing bottleneck Numpy routines with routines imported
#from the following third-party Numpy-like frameworks:
#
#* "bottleneck", providing optimized routines accepting Numpy arrays --
#  implemented in Cython and hence faster than comparible Numpy routines.
#* "numexpr", a Theano-like framework accepting Numpy arrays -- performing
#  CPU-centric parallelization of expensive array operations. Whereas Theano
#  permits such operations to be conveniently expressed in pure-Python, however,
#  numexpr inconveniently requires such operations be expressed as... wait for
#  it, raw strings. So, that sucks. Nonetheless, worth a possible look.
#* "blaze", a purported second-gen Numpy replacement. We harbour sincere doubts,
#  but everything deserves its millisecond to shine in the light. Ah. We see.
#  Blaze is considerably more heavy-weight than Numpy, and largely serves a
#  completely different marketshare: supercomputing. That's well beyond our
#  means, at the moment. Numpy it is!
#
#In short, "bottleneck" is probably the only framework listed above of interest.

#FIXME: Consider optimizing frequently used matrix and vector computations with
#Theano, a general-purpose Python mathematical optimization framework. One
#particularly compelling use case for Theano is to portably distribute
#computational work across multiple GPUs. In general, Theano can be used to
#reduce arbitrarily complex symbolic expressions expressed in pure Python to
#dynamically compiled machine code on-the-fly. Fairly amazing, all around. For
#the high-level synopsis, see:
#
#    http://deeplearning.net/software/theano/introduction.html
#FIXME: Theano and Torch (a similar framewark) appear to now be subsumed by
#TensorFlow, a Google-backed framework originally implemented in support of
#machine learning workflows at Google (e.g., DeepMind), but sufficiently
#generalized as to support a wide variety of computational needs -- like ours.

# ....................{ IMPORTS                           }....................
import numpy
from betse.util.io.log import logs
from betse.util.os import dlls
from betse.util.os.brand import linux, macos, posix
from betse.util.path import dirs, files, pathnames
from betse.util.py import pys
from betse.util.py.module import pymodname, pymodule
from betse.util.type.decorator.decmemo import func_cached
from betse.util.type.iterable import itersort
from betse.util.type.iterable.mapping.mapcls import OrderedArgsDict
from betse.util.type.numeric import versions
from betse.util.type.text import regexes
from betse.util.type.types import BoolOrNoneTypes, ModuleType
from numpy import __config__ as numpy_config

# ....................{ GLOBALS                           }....................
VERSION = numpy.__version__
'''
Human-readable :mod:`numpy` version string (e.g., ``1.14.5``).
'''

# ....................{ GLOBALS ~ opt_info                }....................
# Fully initialized by the _init_globals() function below.
_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX = None
'''
Uncompiled regular expression heuristically matching the basenames of optimized
BLAS shared libraries in the ``libraries`` list of the global
:data:`numpy.__config__.blas_opt_info` dictionary.

This expression does *not* match the strict superset of optimized BLAS shared
libraries that are also optimized, as doing so in a cross-platform manner is
infeasible. Debian-based Linux distributions, for example, remove the prefixes
uniquely identifying the threaded variants of both ATLAS and OpenBLAS from the
basenames of their shared libraries (e.g., ``libatlas.so`` rather than
``libtatlas.so``).

Optimized BLAS shared libraries matched by this expression include:

* **AMD Core Math Library (ACML), unconditionally GPU- and CPU- optimized in
  both OpenMP-based and non-OpenMP-based variants regardless of underlying
  compiler (e.g., GNU Fortran, Open64). Note, however, that ACML does *not*
  ship with a CBLAS interface and hence is non-trivial to link Numpy against.
  While unlikely that any end users will ever do so, it nonetheless remains
  feasible and hence somewhat supported.
* **Automatically Tuned Linear Algebra Software (ATLAS),** both single- and
  multithreaded CBLAS- and Fortran-based variants for both ATLAS < 3.10 and
  ATLAS >= 3.10, which ships shared libraries under different basenames than
  ATLAS < 3.10. (Life complicates life.)
* **Intel Math Kernel Library (MKL),** unconditionally multithreaded in both
  OpenMP-based and non-OpenMP-based variants regardless of underlying compiler
  (e.g., dynamic, GCC, Intel). Note that **Intel Vector Mathematical Library
  (VML)** is intentionally ignored. Although also unconditionally
  multithreaded, VML does *not* implement the BLAS API. Numpy currently
  contains no VML-specific handling, apart from (somewhat uselessly) detecting
  VML installation on reporting system diagnostics.
* **OpenBLAS,** both single- and multithreaded 32- and 64-bit variants.

This expression is typically only required once at application startup and
hence is conditionally compiled in a just-in-time (JIT) manner by the
:func:`_is_blas_optimized_python_general` function rather than unconditionally
compiled here.
'''
# print('blas regex: ' + _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX)


# Fully initialized by the _init_globals() function below.
_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRS_REGEX = None
'''
Uncompiled regular expression heuristically matching the dirnames of optimized
BLAS shared libraries in the `libraries` list of the global
:data:`numpy.__config__.blas_opt_info` dictionary.

See Also
----------
:data:`_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX`
    Further details.
'''


_OPTIMIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_MACOS = {
    # Accelerate. Although Accelerate is only conditionally multithreaded,
    # multithreading is enabled by default and hence a safe assumption.
    '-Wl,Accelerate',

    # vecLib. Similar assumptions as with Accelerate apply.
    '-Wl,vecLib',
}
'''
Set of all strings in the `extra_link_args` list of the global
:data:`numpy.__config__.blas_opt_info` dictionary heuristically corresponding
to macOS-specific optimized BLAS shared libraries.

Unlike all other such libraries, Numpy does _not_ declare unique dictionary
globals describing macOS-specific BLAS shared libraries when linked against.
Hence, this lower-level solution.
'''

# ....................{ GLOBALS ~ linked lib              }....................
# Fully initialized by the _init_globals() function below.
_OPTIMIZED_BLAS_LINKED_LIB_BASENAME_REGEX = None
'''
Uncompiled regular expression heuristically matching the basenames of optimized
BLAS shared libraries dynamically linked to by Numpy.

See Also
----------
:data:`_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX`
    Further details.
'''


# Fully initialized by the _init_globals() function below.
_OPTIMIZED_BLAS_LINKED_LIB_DIRNAME_REGEX = None
'''
Uncompiled regular expression heuristically matching the dirnames of optimized
BLAS shared libraries dynamically linked to by Numpy.

See Also
----------
:data:`_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX`
    Further details.
'''

# ....................{ INITIALIZERS                      }....................
# For simplicity, this function is called below on the first importation of
# this submodule rather than explicitly called by callers.
def init() -> None:
    '''
    Initialize this submodule.

    Specifically (in order):

    #. Initialize all uninitialized global variables of this submodule.
    #. If the currently installed version of Numpy was linked against an
       unoptimized BLAS implementation and is thus itself unoptimized, log a
       non-fatal warning.
    '''

    # Log this initialization.
    logs.log_debug('Initializing NumPy...')

    # Initialize all uninitialized global variables of this submodule.
    _init_globals()

    # If Numpy linked against an unoptimized BLAS, log a non-fatal warning.
    if not is_blas_optimized():
        logs.log_warning(
            'Numpy unoptimized; scaling down to single-core operation. '
            'Consider installing an optimized multithreaded '
            'CBLAS implementation (e.g., OpenBLAS, ATLAS, ACML, MKL) and '
            'reinstalling Numpy to use this implementation.'
        )


def _init_globals() -> None:
    '''
    Initialize all uninitialized global variables of this submodule.
    '''

    # Permit these globals to be redefined.
    global\
        _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX,\
        _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRS_REGEX,\
        _OPTIMIZED_BLAS_LINKED_LIB_BASENAME_REGEX,\
        _OPTIMIZED_BLAS_LINKED_LIB_DIRNAME_REGEX

    # Regular expression fragment matching uniquely identifying substrings of
    # basenames of optimized BLAS shared libraries.
    blas_lib_basename_regex = r'|'.join((
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
    ))

    # Regular expression fragment matching uniquely identifying substrings of
    # dirnames of optimized BLAS shared libraries.
    blas_lib_dirname_regex = r'|'.join((
        # AMD Core Math Library (ACML).
        r'acml',

        # Automatically Tuned Linear Algebra Software (ATLAS).
        r'atlas',

        # Intel Math Kernel Library (MKL).
        r'mkl',

        # OpenBLAS.
        r'openblas',
    ))

    # Regular expression fragment matching the boundary of a dirname at which a
    # substring matching an optimized BLAS substring may either begin or end.
    dirname_boundary_regex = r'[{}_.-]'.format(dirs.SEPARATOR_REGEX)

    #FIXME: What about macOS? Are shared libraries prefixed by "lib" under that
    #platform as well?

    # Redefine this global.
    _OPTIMIZED_BLAS_LINKED_LIB_BASENAME_REGEX = (
        r'^lib({})(?:[_.-].*)?$'.format(blas_lib_basename_regex))

    # Redefine this global. Since Numpy has already stripped all
    # platform-specific prefixes (e.g., "lib") and suffixes (e.g., ".so") from
    # this basename, only this substring followed by an arbitrary bounded
    # suffix need be matched.
    _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX = (
        r'^({})(?:[_-].*)?$'.format(blas_lib_basename_regex))

    # Redefine these globals to the same regular expression.
    _OPTIMIZED_BLAS_LINKED_LIB_DIRNAME_REGEX = (
        r'^.*?{boundary}({blas_dirname})(?:{boundary}.*)?$'.format(
            boundary=dirname_boundary_regex,
            blas_dirname=blas_lib_dirname_regex,))
    _OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRS_REGEX = (
        _OPTIMIZED_BLAS_LINKED_LIB_DIRNAME_REGEX)

# ....................{ TESTERS                           }....................
@func_cached
def is_blas_optimized() -> bool:
    '''
    ``True`` only if the currently installed version of Numpy is linked against
    an optimized BLAS (Basic Linear Algebra Subprograms) implementation,
    ideally but *not* necessarily parallelized across multiple processors.

    Optimized BLAS implementations are *strongly* recommended over unoptimized
    BLAS implementations. The :func:`numpy.dot` operator, which is implicitly
    optimized when Numpy is linked against a optimized BLAS implementation, is
    frequently called by BETSE in its critical path.

    Note that testing for parallelized optimized BLAS implementations, while
    more specific and hence preferable, is infeasible for common edge-cases
    (e.g., Debian-based Linux distributions). For further details, see the
    :data:`_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX` string global.
    '''

    # For each private tester implementing a heuristic for this public test (in
    # order of decreasing generality, portability, and reliability)...
    for tester_heuristic in (
        # Detect conda-managed Numpy first, as doing so reduces to a single
        # well-defined filesystem access and hence is guaranteed to be both the
        # most optimal and portable solution.
        _is_blas_optimized_conda,
        _is_blas_optimized_opt_info_libraries,
        _is_blas_optimized_opt_info_library_dirs,
        _is_blas_optimized_opt_info_macos,
        _is_blas_optimized_posix_symlink,
    ):
        # Attempt to...
        try:
            # Log the current heuristic being attempted.
            logs.log_debug(
                'Detecting BLAS by heuristic %s()...',
                tester_heuristic.__name__)

            # Call this tester, capturing the result for subsequent handling.
            tester_result = tester_heuristic()

            # If this tester definitively identified Numpy as either
            # optimized or non-optimized...
            if tester_result is not None:
                # Log this result.
                logs.log_debug('BLAS optimization detected: %r', tester_result)

                # Return this result.
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

# ....................{ TESTERS ~ private                 }....................
def _is_blas_optimized_conda() -> BoolOrNoneTypes:
    '''
    ``True`` only if the active Python interpreter is managed by ``conda`` *or*
    ``None`` otherwise (i.e., if this interpreter is managed by any other
    means, typically a system-wide package manager).

    If the active Python interpreter is managed by ``conda``, then the current
    version of Numpy was necessarily installed from one of the following two
    prominent Anaconda channels:

    * ``anaconda``, the default proprietary Anaconda channel. In this case,
      Anaconda guarantees Numpy to be linked against Intel Math Kernel Library
      (MKL) and hence optimized.
    * ``conda-forge``, the most popular third-party open-source Anaconda
      channel. Since this and downstream applications are currently only
      available from this channel, Numpy is typically installed from this
      channel during this application's runtime. In this case, conda-forge
      guarantees Numpy to be linked against OpenBLAS and hence optimized.

    In either case, installing Numpy via ``conda`` effectively guarantees
    optimization in all sensible use cases.
    '''

    # Ultimate freedom is a working one-liner.
    #
    # Note that we intentionally return "None" rather than "False" in the event
    # that this interpreter is *NOT* managed by conda. (Returning "False" would
    # erroneously halt the detection process here.)
    return pys.is_conda() or None

# ....................{ TESTERS ~ private : opt_info      }....................
def _is_blas_optimized_opt_info_libraries() -> BoolOrNoneTypes:
    '''
    ``True`` only if the first item of the ``libraries`` list of the global
    :data:`numpy.__config__.blas_opt_info` dictionary heuristically
    corresponds to that of an optimized BLAS implementation, ``False`` if a
    non-fatal error condition arises (e.g., due this list or dictionary being
    undefined), *or* ``None`` otherwise.

    This function returns ``None`` when unable to deterministically decide this
    boolean, in which case a subsequent heuristic will attempt to do so.

    Numpy does *not* define a public API exposing this boolean to callers.
    Numpy only defines a private API defining a medley of metadata from which
    this boolean is indirectly derivable: the :mod:`numpy.__config__`
    submodule. The :func:`numpy.distutils.misc_util.generate_config_py`
    function programmatically fabricates the contents of the
    :mod:`numpy.__config__` submodule at Numpy installation time. Ergo, this
    function introspectively inspects these contents for uniquely identifying
    metadata in a portable manner.
    '''

    # Global BLAS linkage dictionary for this Numpy installation if any or
    # "None" otherwise. Technically, this dictionary should *ALWAYS* be
    # defined.  Reality probably occasionally begs to disagree, however.
    blas_opt_info = getattr(numpy_config, 'blas_opt_info', None)

    # If this dictionary is undefined, log a non-fatal warning and return
    # False. While sad, this is *NOT* worth raising an exception over.
    if blas_opt_info is None:
        logs.log_warning(
            'Numpy installation misconfigured: '
            '"numpy.__config__.blas_opt_info" dictionary not found.')
        return False

    # List of the uniquely identifying substrings of all BLAS library basenames
    # this version of Numpy is linked against in a high-level manner if any or
    # "None" otherwise.
    #
    # Note that this list is *NOT* guaranteed to exist. When this version of
    # Numpy is linked against a BLAS library in a low-level manner (e.g., via
    # "'extra_link_args': ['-Wl,-framework', '-Wl,Accelerate']" on macOS), this
    # list should *NOT* exist. In most other cases, this list should exist. To
    # avoid edge cases, this list is ignored if absent.
    blas_basename_substrs = blas_opt_info.get('libraries', None)

    # If this list is either undefined or empty, silently noop.
    if not blas_basename_substrs:
        return None
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
        regex=_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_REGEX,
    ):
        return True

    # Else, instruct our caller to continue to the next heuristic.
    return None


def _is_blas_optimized_opt_info_library_dirs() -> BoolOrNoneTypes:
    '''
    ``True`` only if the first element of the `library_dirs` list of the
    global :data:`numpy.__config__.blas_opt_info` dictionary heuristically
    corresponds to that of an optimized BLAS implementation, ``False`` if a
    non-fatal error condition arises (e.g., due this list or dictionary being
    undefined), *or* ``None`` otherwise.

    This function returns ``None`` when unable to deterministically decide this
    boolean, in which case a subsequent heuristic will attempt to do so.
    '''

    # List of the dirnames of all BLAS libraries this version of Numpy is
    # linked against in a high-level manner if any or "None" otherwise.
    #
    # Note that the "blas_opt_info" dictionary global is guaranteed to exist
    # due to the previously called _is_blas_optimized_opt_info_basename()
    # function.
    #
    # Note that this list is *NOT* guaranteed to exist. When this version of
    # Numpy is linked against a BLAS library in a low-level manner (e.g., via
    # "'extra_link_args': ['-Wl,-framework', '-Wl,Accelerate']" on macOS), this
    # list should *NOT* exist. In most other cases, this list should exist. To
    # avoid edge cases, this list is ignored if absent.
    blas_dirnames = numpy_config.blas_opt_info.get('library_dirs', None)

    # If this list is either undefined or empty, silently noop.
    if not blas_dirnames:
        return None
    # Else, this list is non-empty.

    # First element of this list. For simplicity, this function assumes the
    # BLAS library identified by this element currently exists. While
    # iteratively testing all listed BLAS libraries for existence would be
    # feasible, doing so is platform-specific and hence non-trivially fragile.
    blas_dirname = blas_dirnames[0]

    # If the BLAS library identified by this element is optimized, return
    # True. Since this element is an arbitrary pathname typically containing
    # non-identifying metadata, a regular expression is leveraged.
    if regexes.is_match(
        text=blas_dirname,
        regex=_OPTIMIZED_BLAS_OPT_INFO_LIBRARY_DIRS_REGEX,
    ):
        return True

    # Else, instruct our caller to continue to the next heuristic.
    return None


def _is_blas_optimized_opt_info_macos() -> BoolOrNoneTypes:
    '''
    ``True`` only if the current platform is macOS *and* the
    ``extra_link_args`` list of the global
    :data:`numpy.__config__.blas_opt_info` dictionary both exists *and*
    heuristically corresponds to that of an optimized BLAS implementation
    specific to macOS (e.g., Accelerate, vecLib), ``False`` if a non-fatal error
    condition arises (e.g., due this list or dictionary being undefined), *or*
    ``None`` otherwise.

    This function returns ``None`` when unable to deterministically decide this
    boolean, in which case a subsequent heuristic will attempt to do so.

    Unlike all other BLAS implementations, macOS-specific BLAS implementations
    are linked against with explicit linker flags rather than pathnames. For
    further confirmation that the :attr:`numpy.__config__.blas_opt_info`
    dictionary defines these flags when linked to these implementations, see:

    * https://trac.macports.org/ticket/22200
    * https://github.com/BVLC/caffe/issues/2677

    When life buys you cat food, you eat cat food.
    '''

    # If the current platform is *NOT* macOS, continue to the next heuristic.
    if not macos.is_macos():
        return None
    # Else, the current platform is macOS.

    # List of all implementation-specific link arguments with which Numpy
    # linked against the current BLAS implementation if any or "None".
    #
    # Note that the "blas_opt_info" dictionary global is guaranteed to exist
    # due to the previously called _is_blas_optimized_opt_info_basename()
    # function.
    blas_link_args_list = numpy_config.blas_opt_info.get(
        'extra_link_args', None)

    # If no such argument exists, continue to the next heuristic. Since this
    # list is strictly optional, no errors or warnings are logged.
    if not blas_link_args_list:
        return None

    # Set of these arguments, converted from this list for efficiency.
    blas_link_args = set(blas_link_args_list)
    # logs.log_info('blas_link_args: {}'.format(blas_link_args))

    # Subset of this set specific to multithreaded BLAS implementations.
    blas_link_args_multithreaded = (
        blas_link_args & _OPTIMIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_MACOS)

    # If this subset is nonempty, return True.
    if len(blas_link_args_multithreaded) > 0:
        return True

    # Else, instruct our caller to continue to the next heuristic.
    return None

# ....................{ TESTERS ~ private : linkage       }....................
def _is_blas_optimized_posix_symlink() -> BoolOrNoneTypes:
    '''
    ``True`` only if the current platform is POSIX-compliant and hence
    supports symbolic links *and* the first item of the ``libraries`` list of
    the global :data:`numpy.__config__.blas_opt_info` dictionary is a symbolic
    link masquerading as either the unoptimized reference BLAS implementation
    but in fact linking to an optimized BLAS implementation.

    This function returns ``None`` when unable to deterministically decide this
    boolean, in which case a subsequent heuristic will attempt to do so.
    '''

    # If the current platform is POSIX-incompatible and hence does *NOT*
    # support symbolic links, continue to the next heuristic.
    if not posix.is_posix():
        return None

    #FIXME: Generalize to macOS as well once the
    #libs.iter_linked_filenames() function supports macOS.

    # If the current platform is *NOT* Linux, continue to the next heuristic.
    #
    # The libs.iter_linked_filenames() function called below currently only
    # supports Linux.
    if not linux.is_linux():
        return None

    # First element of the list of uniquely identifying substrings of all BLAS
    # library basenames this version of Numpy is linked against.
    #
    # Note that this list is guaranteed to both exist and be non-empty due to
    # the previously called _is_blas_optimized_opt_info_basename() function.
    blas_basename_substr = numpy_config.blas_opt_info['libraries'][0]

    # If this element appears to be neither the reference BLAS or CBLAS
    # implementations (e.g., "blas", "cblas", "refblas", "refcblas"), continue
    # to the next heuristic.
    if not blas_basename_substr.endswith('blas'):
        return None

    # Arbitrary Numpy C extension.
    #
    # Unfortunately, the "numpy.__config__" API fails to specify the absolute
    # paths of the libraries it links against. Since there exists no reliable
    # means of reverse engineering these paths from this API, these paths must
    # be obtained by another means: specifically, by querying the standard
    # "numpy.core.multiarray" C extension installed under all supported Numpy
    # for the absolute paths of all external shared libraries to which this
    # extension links -- exactly one of which is guaranteed to be the absolute
    # path of what appears to be a reference BLAS or CBLAS implementation.
    numpy_lib = get_c_extension()

    # Absolute filename of this C extension.
    numpy_lib_filename = pymodule.get_filename(module=numpy_lib)

    # For the basename and absolute filename of each shared library linked to
    # by this Numpy shared library...
    for (numpy_linked_lib_basename, numpy_linked_lib_filename) in (
        dlls.iter_linked_filenames(numpy_lib_filename)):
        # Basename excluding all suffixing filetypes of this library.
        numpy_linked_lib_rootname = pathnames.get_pathname_sans_filetypes(
            numpy_linked_lib_basename)
        # logs.log_info('rootname: %s; basename: %s; filename: %s', numpy_linked_lib_rootname, numpy_linked_lib_basename, numpy_linked_lib_filename)

        # If this appears to be neither the BLAS nor CBLAS reference library,
        # continue to the next library.
        if not numpy_linked_lib_rootname.endswith('blas'):
            continue
        # Else, this is either the BLAS or CBLAS reference library.

        # Absolute filename of the target library to which this library links
        # if this library is a symbolic link *OR* of this library as is
        # otherwise (i.e., if this is a library rather than symbolic link).
        numpy_linked_lib_target_filename = pathnames.canonicalize(
            numpy_linked_lib_filename)
        # logs.log_info('target filename: %s', numpy_linked_lib_target_filename)

        # If either the basename or dirname of this path corresponds to that of
        # an optimized BLAS library, return True.
        if regexes.is_match(
            text=pathnames.get_basename(numpy_linked_lib_target_filename),
            regex=_OPTIMIZED_BLAS_LINKED_LIB_BASENAME_REGEX,
        ) or regexes.is_match(
            text=pathnames.get_dirname(numpy_linked_lib_target_filename),
            regex=_OPTIMIZED_BLAS_LINKED_LIB_DIRNAME_REGEX,
        ):
            return True

        # Else, Numpy links against an unoptimized BLAS implementation. Halt!
        break

    # Else, instruct our caller to continue to the next heuristic.
    return None

# ....................{ GETTERS                           }....................
@func_cached
def get_c_extension() -> ModuleType:
    '''
    Arbitrary Numpy-specific submodule guaranteed to be implemented as a C
    extension.

    Application startup typically tests platform-specific libraries linked
    against this C extension to attempt to dynamically detect whether this
    version of Numpy is multithreaded or not.

    :func:`get_c_extension_name_qualified`
        Further details.
    '''

    # One-liners for the greater glory of BETSE.
    return pymodname.import_module(get_c_extension_name_qualified())


@func_cached
def get_c_extension_name_qualified() -> str:
    '''
    Fully-qualified name of an arbitrary Numpy-specific submodule guaranteed to
    be implemented as a C extension.

    Specifically, this function returns either:

    * If this is Numpy >= 1.16.0, the newly unified
      :mod:`numpy.core._multiarray_umath` C extension. To quote a comment
      heading the pure-Python :mod:`numpy.core.multiarray` submodule in recent
      versions of Numpy:

          Create the ``numpy.core.multiarray`` namespace for backward
          compatibility. In v1.16 the ``multiarray`` and ``umath`` c-extension
          modules were merged into a single ``_multiarray_umath`` extension
          module. So we replicate the old namespace by importing from the
          extension module.

    * Else, the obsoleted :mod:`numpy.core.multiarray` C extension.
    '''

    # Return either...
    return (
        # If Numpy >= 1.16.0, this newly unified C extension.
        'numpy.core._multiarray_umath'
        if versions.is_greater_than_or_equal_to(VERSION, '1.16.0') else
        # Else, Numpy < 1.16.0. In this case, this obsoleted C extension.
        'numpy.core.multiarray'
    )

# ....................{ GETTERS ~ metadata                }....................
def get_metadatas() -> tuple:
    '''
    Tuple of 2-tuples ``(metedata_name, metadata_value)``, describing all
    currently installed third-party dependencies against which Numpy was linked
    (e.g., BLAS, LAPACK).
    '''

    #FIXME: Add LAPACK linkage metadata as well.
    return (
        ('numpy (blas)', get_blas_metadata()),
    )


def get_blas_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current Numpy installation with respect
    to BLAS linkage.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text.string import strs

    # This dictionary.
    metadata = OrderedArgsDict('optimized', is_blas_optimized())

    # Set of all keys of the dictionary global synopsizing this metadata,
    # sorted in ascending lexicographic order for readability.
    blas_opt_info_keys = itersort.sort_ascending(
        tuple(numpy_config.blas_opt_info.keys()))

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
