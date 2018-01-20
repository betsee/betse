#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level logging handler subclasses.
'''

#FIXME: Transparently compress rotated logfiles with the standard ".gz"-style
#compression format. (Surely, StackOverflow has already solved this.)

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid circular import dependencies, avoid importing from *ANY*
# application-specific modules at the top-level -- excluding those explicitly
# known *NOT* to import from this module. Since all application-specific modules
# must *ALWAYS* be able to safely import from this module at any level, these
# circularities are best avoided here rather than elsewhere.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import errno, time
from betse.util.io import stderrs
from logging.handlers import RotatingFileHandler

# ....................{ SUBCLASSES                         }....................
class LogHandlerFileRotateSafe(RotatingFileHandler):
    '''
    Process-safe rotating file handler.

    The standard :class:`RotatingFileHandler` class is thread- but *not*
    process-safe. Concurrent attempts to log to the same physical file from
    multiple processes can and typically will produce fatal race conditions
    producing raised exceptions from one or more of these processes. On logfile
    rotation, each process will aggressively contend with each other process for
    write access to the same physical file to be rotated.

    This :class:`RotatingFileHandler` subclass is both thread- *and*
    process-safe, obviating these concerns.

    Caveats
    ----------
    This handler does *not* implement platform-specific file locking to
    constrain logfile access during rotation. While feasible (e.g., via the
    third-party ConcurrentLogHandler_ package), doing so is complicated by
    numerous non-trivial tradeoffs -- including:

    * **Efficiency.** File locking is notoriously inefficient. File locking
      logfile access, which necessitates acquiring and release a file lock for
      each log record written to disk by this handler, is doubly so. Ideally,
      logging should have negligible to no performance implications. When this
      ceases to be the case, logging ceases to be universally viable.
    * **Security.** File locking is notoriously insecure. Both BSD-style
      ``flock()`` locks *and* POSIX-style ``fcntl()`` locks on files having
      permissions more permissive than `0600` (e.g., group- or world-readable)
      expose applications to permanent deadlocks. Malicious users with group- or
      world-readable access to files to be locked can permanently halt the
      execution of applications locking those files by:

      #. Preemptively acquiring those locks first.
      #. Permanently preserving those locks (i.e., never releasing those locks).

    * **Thread safety.** Specifically:

      * POSIX-style ``fcntl()`` locks lock on process IDs (PIDs) and hence are
        implicitly non-thread-safe.
      * BSD-style ``flock()`` locks lock on file descriptors (FDs) and hence are
        implicitly thread-safe.

    * **Portability.** While the prior point implies BSD-style ``flock()`` locks
      to be preferable to POSIX-style ``fcntl()`` locks for purposes of
      thread-safety, the latter remain mildly more portable than the former.
      Naturally, neither are supported under Microsoft platforms, which provides
      entirely different synchronization platforms suffering completely
      different tradeoffs (e.g., inability to rename files concurrently opened
      by multiple processes). Python's stdlib implements low-level wrappers
      encapsulating all three platform-specific APIs but *no* high-level wrapper
      unifying these fundamentally dissimilar approaches under a common API.
      Optional third-party packages implementing such high-level wrappers exist
      but introduce additional tradeoffs (in addition to those documented
      above), including:

      * **Maintenance.** Most packages are poorly maintained at best, with
        stable releases few and far between.
      * **Documentation.** Most packages are poorly documented at best, with
        little to no human-readable documentation published online.
      * **Uncertainty.** While multiple packages exist, no package appears to
        predominate the others with respect to popularity or usage.

    These deficiencies and more are well-documented, long-standing, and unlikely
    to be resolved on any platform at any point in the near future. Until
    universally resolved on *all* supported platforms, the existence of these
    deficiencies implies file locking to be broken by design and unusable for
    applications in the "real world."

    Non-filesystem locking (e.g., :mod:`multiprocessing` module synchronization
    primitives) is strongly preferable. Unfortunately, such locking requires
    centralized coordination managed at the high-level Python level. No such
    coordination exists between independent BETSE processes run by external
    users at the low-level operating system level.

    Both file- and non-file locking are inapplicable within this context. Hence,
    this handler cannot reasonably constrain logfile access during rotation.
    Instead, on detecting exceptions produced by race conditions between
    multiple processes competing for access when attempting to emit log records,
    this handler temporarily halts the current process for a non-observable
    amount of the CPU timeslice (e.g., 50ms) and repeats the attempt.

    While non-ideal, no sane solutions exist. File locking is insane.

    .. _ConcurrentLogHandler:
       https://pypi.python.org/pypi/ConcurrentLogHandler

    See Also
    ----------
    http://0pointer.de/blog/projects/locking.html
        *On the Brokenness of File Locking.* Classic treatise on the
        long-standing deficiencies of both BSD- and POSIX-style file locking.
    https://loonytek.com/2015/01/15/advisory-file-locking-differences-between-posix-and-bsd-locks
        *Advisory File Locking – My take on POSIX and BSD locks.* Detailed
        analysis concluding with similar deficiencies and lack of solutions.
    '''

    # ..................{ EMITTERS                           }..................
    def emit(self, record) -> None:

        # Attempt to emit this record to this logfile and conditionally rotate
        # this logfile in the default non-process-safe manner.
        try:
            super().emit(record)
        # If an exception indicative of a race condition between multiple
        # processes competing for logfile rotation is raised, attempt to re-emit
        # this record to this logfile in a process-safe manner. Since this
        # necessarily entails inefficiency, we do so *ONLY* as needed.
        #
        # These exceptions are indicated by the following real-world tracebacks:
        #
        #     # Traceback implicating the "FileNotFoundError" exception.
        #     Traceback (most recent call last):
        #       File "/opt/shared/python/3.6.0/lib/python3.6/logging/handlers.py", line 72, in emit
        #         self.doRollover()
        #       File "/opt/shared/python/3.6.0/lib/python3.6/logging/handlers.py", line 169, in doRollover
        #         os.rename(sfn, dfn)
        #     FileNotFoundError: [Errno 2] No such file or directory: '/cluster/home/sburck01/.betse/betse.log.4' -> '/cluster/home/sburck01/.betse/betse.log.5'
        #     Call stack:
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/bin/betse", line 63, in <module>
        #         sys.exit(entry_module.main())
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/cli/__main__.py", line 45, in main
        #         return CLIMain().run(arg_list)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/cli/cliabc.py", line 156, in run
        #         matplotlib_backend_name=self._args.matplotlib_backend_name)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/libs.py", line 258, in reinit
        #         init(*args, **kwargs)
        #       File "<string>", line 15, in func_type_checked
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/libs.py", line 291, in init
        #         mpl_config.init(backend_name=matplotlib_backend_name)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/matplotlib/matplotlibs.py", line 226, in init
        #         self._init_backend(backend_name=backend_name)
        #
        #     # Traceback implicating the "OSError" exception.
        #     Traceback (most recent call last):
        #       File "/opt/shared/python/3.6.0/lib/python3.6/logging/handlers.py", line 71, in emit
        #         if self.shouldRollover(record):
        #       File "/opt/shared/python/3.6.0/lib/python3.6/logging/handlers.py", line 188, in shouldRollover
        #         self.stream.seek(0, 2)  #due to non-posix-compliant Windows feature
        #     OSError: [Errno 116] Stale file handle
        #     Call stack:
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/bin/betse", line 63, in <module>
        #         sys.exit(entry_module.main())
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/cli/__main__.py", line 45, in main
        #         return CLIMain().run(arg_list)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/cli/cliabc.py", line 156, in run
        #         matplotlib_backend_name=self._args.matplotlib_backend_name)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/libs.py", line 258, in reinit
        #         init(*args, **kwargs)
        #       File "<string>", line 15, in func_type_checked
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/libs.py", line 291, in init
        #         mpl_config.init(backend_name=matplotlib_backend_name)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/matplotlib/matplotlibs.py", line 226, in init
        #         self._init_backend(backend_name=backend_name)
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/matplotlib/matplotlibs.py", line 434, in _init_backend
        #         if self.is_backend_usable(_BACKEND_NAME_HEADLESS):
        #       File "<string>", line 15, in func_type_checked
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/lib/matplotlib/matplotlibs.py", line 760, in is_backend_usable
        #         logs.log_debug('Testing matplotlib backend "%s"...', backend_name)
        #       File "<string>", line 15, in func_type_checked
        #       File "/cluster/kappa/90-days-archive/levin/levinlab/sburck01/venv-betse/lib/python3.6/site-packages/betse-0.4.2-py3.6.egg/betse/util/io/log/logs.py", line 100, in log_debug
        #         logging.debug(message, *args, **kwargs)
        #
        # If this is a "FileNotFoundError" exception, unconditionally retry.
        except FileNotFoundError:
            self._emit_safely(record)
        # If this is an "OSError" exception...
        except OSError as exception:
            # If this is the same error as in the traceback above, retry.
            if exception.errno == errno.ESTALE:
                self._emit_safely(record)
            # Else, re-raise this exception.
            else:
                raise
        # Else, permit this exception to continue unwinding the call stack.

    # ..................{ PRIVATE                            }..................
    # Note that, while the emit() method defined above *COULD* be reimplemented
    # to perform the iteration performed by this method and this method then
    # removed, doing so would inefficiently incur the cost of such iteration on
    # ever logging call -- which is clearly unacceptable.
    #
    # Instead, we bite the DRY bullet and simply repeat ourselves below.
    def _emit_safely(record):
        '''
        Attempt to repeatedly emit the passed record in a process-safe manner
        *after* the parent :meth:`emit` call fails to do.
        '''

        # Arbitrary maximum number of times to attempt to re-emit this record.
        MAX_ATTEMPTS = 5

        # Number of seconds to temporarily halt this process, equal to 50ms.
        SLEEP_INTERVAL = 0.05

        # For each such attempt...
        for _ in range(MAX_ATTEMPTS):
            # Notify the user of this race condition. Due to the circumstances,
            # logging this message is right out.
            stderrs.output(
                'Logging race between multiple processes detected...')

            # Temporarily halt the current process in the hope that this race
            # condition will be resolved before this process is awakened.
            time.sleep(SLEEP_INTERVAL)

            # Attempt to re-emit this record to this logfile and conditionally
            # rotate this logfile in the default non-process-safe manner,
            # immediately returning on success.
            try:
                super().emit(record)
                return
            # Catch the same exceptions caught by the emit() method above.
            #
            # If this is a "FileNotFoundError" exception, unconditionally retry.
            except FileNotFoundError:
                continue
            # If this is an "OSError" exception...
            except OSError as exception:
                # If this is the same error as in the traceback above, retry.
                if exception.errno == errno.ESTALE:
                    continue
                # Else, re-raise this exception.
                else:
                    raise
            # Else, permit this exception to continue unwinding the call stack.

        # Attempt to re-emit this record to this logfile one last time. If this
        # attempt fails, a non-human-readable exception will be raised. This is
        # probably a good thing. Assuming this exception is reported back to us,
        # the above logic may be improved by receiving this exception.
        super().emit(record)
