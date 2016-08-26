#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Numpy, a mandatory runtime dependency.
'''

# ....................{ IMPORTS                            }....................
from betse.util.py import modules
from numpy.distutils import __config_ as numpy_config

# ....................{ GLOBALS                            }....................
#FIXME: Sadly, a dictionary lookup does *NOT* suffice under OS X. This platform
#provides two platform-specific multithreaded BLAS implementations: "Accelerate"
#and "vecLib". Numpy does *NOT* currently provide separate dictionary globals
#describing these implementations if they are linked against -- unlike all other
#BLAS implementations. That said, there is still a reliable means of deciding
#whether these implementations have been linked against or not:
#
#* Under OS X, if:
#  * The "numpy.distutils.__config__" submodule declares no dictionary globals
#    whose names are elements of this set *AND*
#  * The "numpy.distutils.__config__.blas_opt_info" dictionary global contains
#    a "extra_link_args" key (which is *NOT* guaranteed to be the case) *AND*
#  * The value of that key is a list containing either of the following two
#    literal strings:
#    * '-Wl,Accelerate'
#    * '-Wl,vecLib'
#* Then Numpy is linked against a multithreaded BLAS implementation.
#
#A bit crazy, admittedly, but it's the only feasible approach.

_CONFIG_BLAS_MULTITHREADED_GLOBAL_NAMES = {
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
    'atlas_blas_threads_info',
    'atlas_3_10_blas_threads_info',

    # BLAS-like Library Instantiation Software (BLIS). Unconditionally
    # multithreaded. Technically, Numpy has yet to add official support for
    # BLIS. Since numerous contributors nonetheless perceive BLIS to be the
    # eventual successor of BLAS *AND* since Numpy currently hosts an open pull
    # request to explicitly add BLIS support under the sensible subclass name
    # "blis_info" (see Numpy PR #7294), explicitly listing BLIS here should
    # assist in future-proofing our multithreading detection.
    'blis_info',

    # Intel Math Kernel Library (MKL). Unconditionally multithreaded.
    'blas_mkl_info',

    # OpenBLAS. Unconditionally multithreaded.
    'openblas_info',
}
'''
Set of the names of all possible dictionary globals declared by the
:mod:`numpy.distutils.__config__` submodule specific to multithreaded BLAS
implementations.

Each element of this set is guaranteed to be the unqualified name of a subclass
of the :class:`numpy.distutils.system_info.system_info` base class.

See Also
----------
:meth:`numpy.distutils.system_info.blas_opt_info.calc_info`
    Method whose body demonstrates the canonical heurestic employed by Numpy to
    decide which BLAS implementation to link against at installation time.
'''

# ....................{ TESTERS                            }....................
#FIXME: Call this in the libs.init() method. It's hardly critical, so it can
#safely wait to the customary calling point.
def is_multithreaded() -> bool:
    '''
    `True` only if the currently installed version of Numpy is linked against a
    multithreaded BLAS (Basic Linear Algebra Subprograms) implementation.

    Multithreaded BLAS implementations are _strongly_ recommended over
    single-threaded BLAS implementations. BETSE frequently calls the `np.dot()`,
    operator, which is implicitly multithreaded when linked against a
    multithreaded BLAS implementation, in its critical path.

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
    `True` only if that submodule declares a global specific to a multithreaded
    BLAS library.

    .. _ATLAS: http://math-atlas.sourceforge.net
    '''

    # Set of the names of all dictionary globals concerning Numpy configuration.
    config_global_names = modules.get_global_names(numpy_config)

    # Subset of this set specific to multithreaded BLAS implementations.
    config_blas_multithreaded_global_names = (
        config_global_names & _CONFIG_BLAS_MULTITHREADED_GLOBAL_NAMES)

    # Return True only if this subset is nonempty.
    return len(config_blas_multithreaded_global_names) > 0
