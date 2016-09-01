#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Numpy, a mandatory runtime dependency.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.py import modules
from betse.util.os import libs, oses
from betse.util.type import iterables, regexes, strs
from collections import OrderedDict
from numpy.distutils import __config__ as numpy_config

# ....................{ GLOBALS                            }....................
_PARALLELIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX = r'^({}).*$'.format(
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
Uncompiled regular expression matching the basename of a parallelized BLAS
shared library.

Parallelized BLAS shared libraries matched by this expression include:

* **AMD Core Math Library (ACML), unconditionally GPU- and CPU- parallelized in
  both OpenMP-based and non-OpenMP-based variants regardless of underlying
  compiler (e.g., GNU Fortran, Open64). Note, however, that ACML does _not_ ship
  with a CBLAS interface and hence is non-trivial to link Numpy against. While
  unlikely that any end users will ever do so, it nonetheless remains feasible
  and hence somewhat supported.
* **Automatically Tuned Linear Algebra Software (ATLAS),** multithreaded
  CBLAS-based and Fortran-based variants for both ATLAS < 3.10 and ATLAS >=
  3.10, which ships shared libraries under different basenames than ATLAS <
  3.10. (Life complicates life.)
* **Intel Math Kernel Library (MKL),** unconditionally multithreaded in both
  OpenMP-based and non-OpenMP-based variants regardless of underlying compiler
  (e.g., dynamic, GCC, Intel).  Note that **Intel Vector Mathematical Library
  (VML)** is intentionally ignored.  Although also unconditionally
  multithreaded, VML does _not_ implement the BLAS API. Numpy currently contains
  no VML-specific handling, apart from (somewhat uselessly) detecting VML
  installation on reporting system diagnostics.
* **OpenBLAS,** both multithreaded 32- and 64-bit variants. All single-threaded
  variants of OpenBLAS are ignored.

This expression is typically only required once at application startup and hence
is conditionally compiled in a just-in-time (JIT) manner by the
:func:`is_blas_parallelized` function rather than unconditionally compiled here.

See Also
----------
:data:`_PARALLELIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_OS_X`
    Set matching OS X-specific parallelized BLAS implementations -- which, for
    obscure reasons pertaining to Numpy internals, are _not_ matchable via this
    regular expression..
'''
# print('blas regex: ' + _PARALLELIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX)


_PARALLELIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_OS_X = {
    # Accelerate. Although Accelerate is only conditionally multithreaded,
    # multithreading is enabled by default and hence a safe assumption.
    '-Wl,Accelerate',

    # vecLib. Similar assumptions as with Accelerate apply.
    '-Wl,vecLib',
}
'''
Set of all uniquely identifying elements of the list value of the
`extra_link_args` key of the :data:`numpy.distutils.__config__.blas_opt_info`
dictionary specific to parallelized BLAS implementations under OS X.

Unlike all other BLAS implementations, Numpy does _not_ declare unique
dictionary globals describing these implementations when linked against. Ergo,
this lower-level solution.
'''

# ....................{ INITIALIZERS                       }....................
# For simplicity, this function is called below on the first importation of this
# submodule rather than explicitly called by callers.
def init() -> None:
    '''
    Initialize Numpy.

    Specifically:

    * If the currently installed version of Numpy was linked against an
      unparallelized BLAS implementation and is thus itself unparallelized, log
      a non-fatal warning.
    '''

    # If Numpy linked against an unparallelized BLAS, log a non-fatal warning.
    if not is_blas_parallelized():
        return

        #FIXME: Sadly, the is_blas_parallelized() is insufficiently granular to
        #support such behaviour at the moment. Improve that method substantially
        #before reenabling this warning.
        logs.log_warning(
            'Numpy not parallelized. '
            'Consider installing a parallelized BLAS implementation '
            '(e.g., OpenBLAS, ATLAS, ACML, MKL) and '
            'reinstalling Numpy against this implementation.'
        )

# ....................{ TESTERS                            }....................
#FIXME: Revise docstring, which is pretty much completely wrong now.

def is_blas_parallelized() -> bool:
    '''
    `True` only if the currently installed version of Numpy is linked against a
    BLAS (Basic Linear Algebra Subprograms) implementation implicitly
    parallelized across multiple processors -- either CPU- or GPU-based.

    Parallelized BLAS implementations are _strongly_ recommended over
    unparallelized BLAS implementations. The `numpy.dot()` operator, which is
    implicitly parallelized when Numpy is linked against a parallelized BLAS
    implementation, is frequently called by BETSE in its critical path.

    Heuristic
    ----------
    Numpy does _not_ provide a convenient API for readily querying this boolean.
    Numpy does, however, provide an admittedly inconvenient API for aggregating
    this boolean together from various sources: the
    :mod:`numpy.distutils.__config__` submodule. The
    :func:`numpy.distutils.misc_util.generate_config_py` function
    programmatically fabricates the contents of the
    :mod:`numpy.distutils.__config__` submodule at Numpy installation time.

    Specifically, for each subclass of the
    :class:`numpy.distutils.system_info.system_info` base class defined by the
    :mod:`numpy.distutils.system_info` submodule (e.g.,
    :class:`numpy.distutils.system_info.atlas_info`) whose corresponding shared
    library (e.g., ATLAS_) or feature (e.g., ATLAS_ multithreading) is available
    on the current system at Numpy installation time, a dictionary global of the
    same name as that subclass whose keys are the names of metadata types and
    values are metadata is programmatically added to the
    :mod:`numpy.distutils.__config__` submodule. Hence, this function returns
    `True` only if that submodule declares a global specific to a parallelized
    BLAS implementation.

    .. _ATLAS: http://math-atlas.sourceforge.net
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

    # If the BLAS library identified by this element is parallelized, return
    # True. Since this element may be suffixed by non-identifying metadata
    # (e.g., version), a regular expression is leveraged.
    if regexes.is_match(
        text=blas_basename_substr,
        regex=_PARALLELIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX,
    ):
        return True

    # If the current platform is OS X, fallback to testing whether Numpy was
    # linked against a parallelized BLAS implementation specific to OS X:
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
        # linked against the current BLAS implementation if any or "None"
        # otherwise.
        blas_link_args_list = blas_lib.get('extra_link_args', None)

        # If at least one such argument exists...
        if blas_link_args_list:
            # Set of these arguments, converted from this list for efficiency.
            blas_link_args = set(blas_link_args_list)

            # Subset of this set specific to multithreaded BLAS implementations.
            blas_link_args_multithreaded = (
                blas_link_args &
                _PARALLELIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_OS_X
            )

            # Return True only if this subset is nonempty.
            return len(blas_link_args_multithreaded) > 0

    # If this BLAS library appears to be either the reference BLAS or CBLAS
    # implementations *AND* this platform is POSIX-compliant and hence supports
    # symbolic links, fallback to testing whether this library is in fact a
    # symbolic link to a parallelized BLAS implementation.
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
    # if (
    #     blas_basename_substr == 'blas' or
    #     blas_basename_substr == 'cblas'
    # ) and oses.is_linux():
    #     # If the standard "numpy.core.multiarray" submodule is importable...
    #     if modules.is_module('numpy.core.multiarray'):
    #         # Do so.
    #         import numpy.core.multiarray as numpy_lib
    #
    #         # Absolute path of this submodule.
    #         numpy_lib_filename = modules.get_filename(numpy_lib)
    #
    #         #FIXME: Test whether or not this path is that of a shared library
    #         #first (e.g., is suffixed by ".so" under Linux).
    #
    #         # Absolute paths of all shared libraries required by this library.
    #         numpy_lib_libs = libs.get_dependency_filenames(numpy_lib_filename)

    # Else, all hope is lost.
    return False

# ....................{ GETTERS                            }....................
def get_metadatas() -> tuple:
    '''
    Tuple of 2-tuples `(metedata_name, metadata_value`), describing all
    currently installed third-party dependencies against which Numpy was linked
    (e.g., BLAS, LAPACK).
    '''

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
        ('parallelized', is_blas_parallelized()),
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
