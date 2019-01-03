#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **callback** (i.e., caller-defined object external to this
application whose methods are called by other methods internal to this
application as a means of facilitating inter-package communication)
functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseCallbackException
from betse.util.io.log import logs
from betse.util.type.types import type_check, StrOrNoneTypes

# ....................{ SUPERCLASSES                      }....................
class CallbacksBC(object):
    '''
    Concrete base class of all **callbacks** (i.e., caller-defined object
    external to this application whose methods are called by other methods
    internal to this application as a means of facilitating inter-package
    communication) subclasses.

    Each of the public methods defined by this superclass is a callback
    intended to be called *only* by the **source callable** (i.e., function or
    method calling each such callback), which receives this object from the
    **sink callable** (i.e., function or method instantiating this object and
    passing this object to the source callable).

    Attributes
    ----------
    _progress_min : int
        Minimum value passed to the most recent call to the
        :meth:`progress_ranged` callback by the source callable. While commonly
        0, the sink callable should *not* assume this to be the case.
    _progress_max : int
        Maximum value passed to the most recent call to the
        :meth:`progress_ranged` callback by the source callable.
    _progress_next : int
        Next value to be implicitly passed to the next call to the
        :meth:`progressed` callback by the :meth:`progressed_next` callback.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self) -> None:
        '''
        Initialize this callbacks object.
        '''

        # Nullify all instance variables for safety.
        self._clear_progress()

    # ..................{ CALLBACKS ~ progress              }..................
    @type_check
    def progress_ranged(
        self, progress_max: int, progress_min: int = 0) -> None:
        '''
        Callback passed the range of progress values subsequently passed to the
        :meth:`progressed` callback by the source callable.

        Design
        ----------
        The source callable should explicitly call the following callbacks
        *before* performing significant work (in order):

        #. This callback exactly once. Doing so enables the sink callable to
           configure external objects (e.g., progress bars) requiring this
           range.
        #. Either:

           * The :meth:`progressed_next` callback. (Recommended.)
           * The :meth:`progressed` callback with the minimum progress value
             passed to this callback. Since doing so is equivalent to calling
             the :meth:`progressed_next` callback, calling the latter comes
             recommended.

        In either case, the call to the second such callback initializes the
        same external objects (e.g., progress bars) to the minimum progress of
        this range. By design, this :meth:`progress_ranged` callback does *not*
        implicitly call either the :meth:`progressed` or
        :meth:`progressed_next` callbacks on your behalf. Doing so would invite
        chicken-and-egg subtleties with subclasses overriding this method to
        configure external objects.

        Subclasses are advised to reimplement this callback and call this
        superclass implementation in their own reimplementation.

        Parameters
        ----------
        progress_min : optional[int]
            Minimum value subsequently passed to the :meth:`progressed`
            callback by the source callable. Defaults to 0, ensuring the
            ``progress_max`` parameter to be exactly equal to the number of
            times that the :meth:`progressed_next` callback is called by the
            source callable. Nonetheless, the sink callable should *not* assume
            this value to always be 0.
        progress_max : int
            Maximum value subsequently passed to the :meth:`progressed`
            callback by the source callable. If the ``progress_min`` parameter
            is 0, this is equivalent to the number of times that the
            :meth:`progressed_next` callback is called by the source callable.

        Raises
        ----------
        :class:`BetseCallbackException`
            If either:

            * The source callable is still performing work (i.e., the
              :meth:`progressed` callback has yet to be called with a progress
              value equal to ``progress_max``).
            * ``progress_min`` is greater than or equal to ``progress_max``.
        '''

        # If this range is invalid, raise an exception.
        if progress_min > progress_max:
            raise BetseCallbackException(
                'Minimum progress {} exceeds maximum progress {}.'.format(
                    progress_min, progress_max))
        elif progress_min == progress_max:
            raise BetseCallbackException(
                'Minimum progress {} equals maximum progress {}.'.format(
                    progress_min, progress_max))

        # If the source callable is still performing work, raise an exception.
        self._die_if_progressing()

        # Classify all passed parameters for subsequent reference by the
        # progressed() callback.
        self._progress_min = progress_min
        self._progress_max = progress_max

        # Value to be subsequently passed to the progressed() callback by the
        # first call to the progressed_next() callback. See the
        # progressed_next() callback for a justification of the "+ 1" here.
        self._progress_next = progress_min + 1
        # logs.log_debug('next progress: %d', self._progress_next)


    @type_check
    def progress_stated(self, status: str) -> None:
        '''
        Callback passed the current state of progress for work performed by the
        source callable as an arbitrary human-readable string.

        The superclass implementation of this callback logs this string at the
        :attr:`LogLevel.INFO` level, under the generally safe assumption that
        this string is of reasonable relevance to developers *and* users alike.

        Design
        ----------
        The source callable should call this callback repeatedly while
        performing any significant work -- ideally either immediately before or
        after calling the :meth:`progressed`, :meth:`progressed_last`, or
        :meth:`progressed_next` callbacks. Doing so enables the sink callable
        to incrementally update external objects (e.g., status bars) in
        synchronicity with numeric-based progress updates.

        The source callable may technically call this callback at any arbitrary
        time. Unlike the :meth:`progress_ranged` callback and the
        aforementioned callbacks, this callback is purely subjective and hence
        *not* constrained by hard contractual obligations. For usability by the
        sink callable, however, this callable should typically *only* be called
        at approximately the same time as the aforementioned callbacks are
        called by the source callable; doing so increases the likelihood that
        the sink callable is prepared to handle this callback properly.

        Subclasses are advised to reimplement this callback and call this
        superclass implementation in their own reimplementation.

        Parameters
        ----------
        status : str
            Current state of progress for work performed by the source
            callable as an arbitrary human-readable string. To enable the sink
            callable to safely embed this string in size-constrained regions,
            this string should be terse and ideally at most a single line.

        See Also
        ----------
        :meth:`progressed_next`
        :meth:`progressed_last`
            Higher-level callbacks wrapping this lower-level callback in a
            safer manner complying with the above recommendations -- notably,
            that this callback be called in unison with the :meth:`progressed`
            callback by the source callable.
        '''

        # Log this string as an informational message.
        logs.log_info(status)

    # ..................{ CALLBACKS ~ progressed            }..................
    @type_check
    def progressed(self, progress: int) -> None:
        '''
        Callback passed the current state of progress for work performed by the
        source callable.

        If this is the last possible progress value (i.e., if the passed
        ``progress`` parameter is equal to the ``progress_max`` parameter
        previously passed to the :meth:`progress_ranged` callback), this method
        additionally records the source callable to no longer be performing
        work by calling the :meth:`_clear_progress` method. This preserves
        contractual guarantees, including the guarantee that the next call to
        either this or the :meth:`progressed_next` callback will raise an
        exception.

        Design
        ----------
        The source callable should call this callback repeatedly while
        performing any significant work. Doing so enables the sink callable to
        incrementally update external objects (e.g., progress bars).

        Subclasses are advised to reimplement this callback and call this
        superclass implementation in their own reimplementation.

        Parameters
        ----------
        progress : int
            Current state of progress for work performed by the source
            callable. If this integer is *not* in the range previously passed
            to the :meth:`progress_ranged` callback, an exception is raised. If
            this is the first call to this callback by the source callable,
            this integer should ideally be one larger than the minimum progress
            value previously passed to the :meth:`progress_ranged` callback.
            See the :meth:`progressed_next` callback for further discussion.

        Raises
        ----------
        :class:`BetseCallbackException`
            If either:

            * The :meth:`progress_ranged` callback has yet to be called.
            * This progress is *not* in the range previously passed to the
              :meth:`progress_ranged` callback by this source callable.

        See Also
        ----------
        :meth:`progressed_next`
        :meth:`progressed_last`
            Higher-level callbacks automatically passing the next progress
            value according to the the range previously passed to the
            :meth:`progress_ranged` callback.
        '''

        # If progress_ranged() has yet to be called, raise an exception.
        self._die_unless_progressing()

        # If this progress is *NOT* in the range previously passed to
        # progress_ranged(), raise an exception.
        if not (self._progress_min <= progress <= self._progress_max):
            raise BetseCallbackException(
                'Progress {} not in range [{}, {}].'.format(
                    progress, self._progress_min, self._progress_max))

        # If this is the last possible progress value, record the source
        # callable to no longer be performing work.
        if progress == self._progress_max:
            self._clear_progress()

    # ..................{ CALLBACKS ~ progress : next       }..................
    def progressed_next(self, status: StrOrNoneTypes = None) -> None:
        '''
        Higher-level callback wrapping the lower-level :meth:`progressed` and
        :meth:`progress_stated` callbacks by automatically passing the next
        progress value to tha former callback according to the range previously
        passed to the :meth:`progress_ranged` callback *and* passing the passed
        human-readable string if any to the :meth:`progress_stated` callback.

        Note that callers recalling this callback a predetermined number of
        times are advised to replace the last call to this callback
        with a call to the higher-level :meth:`progressed_last` callback; doing
        so guarantees that the last progress value passed to the sink callable
        is the maximum progress value previously passed to the
        :meth:`progress_ranged` callback, preserving an essential contractual
        obligation of this API.

        Guarantees
        ----------
        The first call to this callback is guaranteed to emit one greater than
        the minimum progress value previously passed to the
        :meth:`progress_ranged` callback. While arguably unintuitive, this
        behaviour corresponds with expectations in the sink callable.

        Assuming that the minimum progress value previously passed to the
        :meth:`progress_ranged` callback was 0 as is the common case, this
        behaviour has the additional advantage of guaranteeing that the maximum
        progress value previously passed to that callback is exactly equal to
        the number of times that this :meth:`progressed_next` callback is
        called by the source callable.

        Specifically, since the sink callable typically handles the
        :meth:`progress_ranged` callback by configuring a progress bar to the
        passed range, the sink callable is assumed to have already "consumed"
        the first progress value for all intents and purposes; instructing the
        sink callable to "reconsume" the same progress value would effectively
        hide the actual first iteration of work performed by the source
        callable from the end user perspective.

        Design
        ----------
        For generality, the sink callable should typically *only* handle the
        lower-level :meth:`progressed` and :meth:`progress_stated` callbacks.
        This callback is merely syntactic sugar conveniently simplifying the
        implementation of the source callback calling this callback.

        Parameters
        ----------
        status : StrOrNoneTypes
            Current state of progress for work performed by the source
            callable as an arbitrary human-readable string. To enable the sink
            callable to safely embed this string in size-constrained regions,
            this string should be terse and ideally at most a single line.
            Defaults to ``None``, in which case the :meth:`progress_stated`
            callback is *not* called by this higher-level callback.

        Raises
        ----------
        :class:`BetseCallbackException`
            If the :meth:`progress_ranged` callback has yet to be called.

        See Also
        ----------
        :meth:`progressed`
        :meth:`progress_stated`
            Lower-level callbacks wrapped by this higher-level callback.
        '''

        # If progress_ranged() has yet to be called, raise an exception.
        self._die_unless_progressing()

        # Notify the sink callback of the current state of progress.
        self.progressed(progress=self._progress_next)

        # If the source callback has work remaining to perform, increment this
        # state *AFTER* notifying this callback.
        if self._progress_next is not None:
            self._progress_next += 1

        # If passed status, notify the sink callback of this status.
        if status is not None:
            self.progress_stated(status)


    def progressed_last(self, *args, **kwargs) -> None:
        '''
        Higher-level callback wrapping the lower-level :meth:`progressed_next`
        callback by internally calling that callback and raising an exception
        if the current progress value is *not* the maximum progress value
        previously passed to the :meth:`progress_ranged` callback.

        This callback is intended to be called by callers re-calling the
        :meth:`progressed_next` callback a predetermined number of times as a
        safer replacement for the last such call, ensuring that the last
        progress value passed to the sink callable is indeed the maximum.

        Parameters
        ----------
        All positional and keyword arguments are passed as is to the
        lower-level :meth:`progressed_next` callback wrapped by this callback.

        Raises
        ----------
        :class:`BetseCallbackException`
            If either:

            * The :meth:`progress_ranged` callback has yet to be called.
            * The current progress value is *not* the maximum progress value.

        See Also
        ----------
        :meth:`progressed_next`
            Lower-level callback wrapped by this higher-level callback.
        '''

        # Notify the sink callback of the current state of progress.
        self.progressed_next(*args, **kwargs)

        # If the source callable is still performing work, raise an exception.
        self._die_if_progressing()

    # ..................{ EXCEPTIONS                        }..................
    def _die_if_progressing(self) -> None:
        '''
        Raise an exception if the source callable is still performing work.

        Raises
        ----------
        :class:`BetseCallbackException`
            If the source callable is currently performing work.

        See Also
        ----------
        :meth:`_is_progressing`
            Further details.
        '''

        if self._is_progressing:
            raise BetseCallbackException(
                'progress_ranged() callback called before '
                'progressed() callback called with maximum progress value.')


    def _die_unless_progressing(self) -> None:
        '''
        Raise an exception unless the source callable is still performing work.

        Raises
        ----------
        :class:`BetseCallbackException`
            If the source callable is *not* currently performing work.

        See Also
        ----------
        :meth:`_is_progressing`
            Further details.
        '''

        if not self._is_progressing:
            raise BetseCallbackException(
                'progressed()-style callback called before '
                'progress_ranged() callback '
                '(e.g., due to previously passing a maximum progress value to '
                'progressed(), resetting progress state and requiring calling '
                'progress_ranged() again).'
            )

    # ..................{ TESTERS                           }..................
    @property
    def _is_progressing(self) -> bool:
        '''
        ``True`` only if the source callable is still performing work.

        Equivalently, this method returns ``False`` only:

        * If the source callable is *not* currently performing work.
        * If the state of work performed by this callable has been cleared
          (i.e., to ``None`` values) rather than set to non-``None`` integers.
        * Unless the :meth:`progress_ranged` callback has been called more
          recently than the :meth:`_clear_progress` method.
        '''

        return (
            self._progress_min  is not None and
            self._progress_max  is not None and
            self._progress_next is not None
        )

    # ..................{ CLEARERS                          }..................
    def _clear_progress(self) -> None:
        '''
        Record the source callable to no longer be performing work by clearing
        (i.e., resetting) the state of work performed by this callable to
        ``None`` values.
        '''

        # Nullify all instance variables for safety.
        self._progress_min = None
        self._progress_max = None
        self._progress_next = None
