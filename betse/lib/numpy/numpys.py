#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Numpy, a mandatory runtime dependency.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.os import oses
from betse.util.type import iterables, regexes, strs
from collections import OrderedDict
from numpy.distutils import __config__ as numpy_config

# ....................{ GLOBALS                            }....................
#FIXME: Replace with a new "_PARALLELIZED_BLAS_OPT_INFO_LIBRARIES" set. The
#current approach appears to be fundamentally incorrect, sadly. It appears
#feasible to have, for example, an "atlas_info" simply reporting "NOT
#AVAILABLE". Ergo, this set as currently defined must go.
#FIXME: Actually, for robustness, just replace this with a single regular
#expression. Attempting to match all possible basenames with a finite number of
#set entries is technically infeasible, as basenames may be versioned.

_PARALLELIZED_BLAS_OPT_INFO_LIBRARY_BASENAME_REGEX = (
    # Basename beginning, matching exactly one of:
    r'^('

        # OpenBLAS. Match only the multithreaded variant of the 32- or 64-bit
        # OpenBLAS implementation. The single-threaded variant is ignored.
        r'openblas(_int64)?_threads'
        # r'|'

    # Basename end, consuming all remaining optional metadata (e.g., version).
    r').*?$'
)
'''
Uncompiled regular expression matching the basename of a parallelized BLAS
shared library.

This expression is conditionally compiled in a just-in-time manner rather than
unconditionally compiled. This expression is typically only required once at
application startup.
'''


#FIXME: Excise after adding the corresponding support to the above regex.
# _CONFIG_BLAS_MULTITHREADED_GLOBAL_NAMES = {
    #FIXME: Research this.
    # AMD Core Math Library (ACML). Unconditionally GPU- and CPU-parallelized.
    #'????',

    # Automatically Tuned Linear Algebra Software (ATLAS). Frustratingly, ATLAS
    # is shipped in both single- and multithreaded variants. Furthermore, ATLAS
    # >= 3.10 ships shared libraries under different basenames than under ATLAS
    # < 3.10. With respect to BLAS multithreading, three discrete states exist:
    #
    # * Single-threaded ATLAS regardless of version, in which case neither the
    #   "atlas_blas_threads_info" nor
    #   "atlas_3_10_blas_threads_info" globals are defined.
    # * Multi-threaded ATLAS < 3.10, in which case only the
    #   "atlas_blas_threads_info" global is defined.
    # * Multi-threaded ATLAS >= 3.10, in which case only the
    #   "atlas_3_10_blas_threads_info" global is defined.

    #FIXME: Research this.
    # 'atlas_blas_threads_info',
    # 'atlas_3_10_blas_threads_info',

    #FIXME: Research this.
    # BLAS-like Library Instantiation Software (BLIS). Unconditionally
    # multithreaded. Technically, Numpy has yet to add official support for
    # BLIS. Since numerous contributors nonetheless perceive BLIS to be the
    # eventual successor of BLAS *AND* since Numpy currently hosts an open pull
    # request to explicitly add BLIS support under the sensible subclass name
    # "blis_info" (see Numpy PR #7294), explicitly listing BLIS here should
    # assist in future-proofing our multithreading detection.
    # 'blis_info',

    #FIXME: Research this.
    # Intel Math Kernel Library (MKL). Unconditionally multithreaded.
    # 'blas_mkl_info',

    # OpenBLAS. Conditionally multithreaded. Ignore the single-threaded variant.
#     'openblas_threads',
# }
'''
Set of the names of all possible dictionary globals declared by the
:mod:`numpy.distutils.__config__` submodule specific to parallelized BLAS
implementations.

Each element of this set is guaranteed to be the unqualified name of a subclass
of the :class:`numpy.distutils.system_info.system_info` base class.

See Also
----------
:meth:`numpy.distutils.system_info.blas_opt_info.calc_info`
    Method whose body demonstrates the canonical heurestic employed by Numpy to
    decide which BLAS implementation to link against at installation time.
'''


#FIXME: Rename to "_PARALLELIZED_BLAS_OPT_INFO_EXTRA_LINK_ARGS_OS_X".
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
#FIXME: Add support for detecting ACML (AMD Core Math Library), an OpenCL-based
#BLAS implementation parallelized over GPU shader units. While Numpy has no
#explicit support for detecting this library, Numpy appears to be trivially
#linkable against ACML by simply replacing the standard BLAS reference library
#with a symbolic link to ACML and reinstalling Numpy. Ergo, detecting ACML may
#reduce to:
#
#* Under POSIX-compliant platforms, determining whether:
#  * The currently linked BLAS library is reported as being the standard BLAS
#    reference library *AND*
#  * The shared BLAS library at this path (perhaps obtained via the
#    "numpy_config.blas_opt_info['libraries']" list) is a symbolic link *AND*
#  * This symbolic link refers to a pathname indicative of ACML.

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

    # If the current platform is OS X, test for whether Numpy was linked against
    # a parallelized BLAS implementation only available under OS X: namely,
    # either "Accelerate" or "vecLib". Unlike all other BLAS implementations,
    # these implementations are linked against with explicit linker flags --
    # requiring further logic. When life sells you cat food, you eat cat food.
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
