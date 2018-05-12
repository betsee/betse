#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **callback** (i.e., caller-defined object external to this
application whose methods are called by other methods internal to this
application as a means of facilitating inter-package communication)
functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta  #, abstractmethod
from betse.exceptions import BetseCallbackException
from betse.util.io.log import logs
from betse.util.type.types import type_check  #, NoneType

# ....................{ SUPERCLASSES                       }....................
class CallbacksABC(metaclass=ABCMeta):
    '''
    Abstract base class of all **callbacks** (i.e., caller-defined object
    external to this application whose methods are called by other methods
    internal to this application as a means of facilitating inter-package
    communication) subclasses.

    Design
    ----------
    Each of the public methods defined by this superclass is a callback
    intended to be called *only* by the **source callable** (i.e., function or
    method calling each such callback). The source callback receives this
    object from the **sink callable** (i.e., function or method instantiating
    this object and passing this object to the source callable).

    Attributes
    ----------
    _progress_min : int
        Minimum value subsequently passed to the :meth:`progressed` callback by
        the source callable. While typically 0, the sink callable should *not*
        necessarily assume this to be the case.
    _progress_max : int
        Maximum value subsequently passed to the :meth:`progressed` callback by
        the source callable.
    _progress_next : int
        Value to be subsequently passed to the :meth:`progressed` callback by
        the :meth:`progressed_next` callback.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this callbacks object.
        '''

        # Nullify all instance variables for safety.
        self._progress_min = None
        self._progress_max = None
        self._progress_next = None

    # ..................{ CALLBACKS ~ progress               }..................
    @type_check
    def progress_ranged(self, progress_min: int, progress_max: int) -> None:
        '''
        Callback passed the range of progress values subsequently passed to the
        :meth:`progressed` callback by the source callable.

        Design
        ----------
        The source callable should call this callback exactly once per call
        *before* either performing any significant work or calling the
        :meth:`progressed` callback. Doing so enables the sink callable to
        initialize external objects (e.g., progress bars) requiring this range.

        Subclasses are required to reimplement this callback and call this
        superclass implementation in their own reimplementation.

        Parameters
        ----------
        progress_min : int
            Minimum value subsequently passed to the :meth:`progressed`
            callback by the source callable. While typically 0, the sink
            callable should *not* necessarily assume this to be the case.
        progress_max : int
            Maximum value subsequently passed to the :meth:`progressed`
            callback by the source callable.

        Raises
        ----------
        :class:`BetseCallbackException`
            If ``progress_min`` is greater than or equal to ``progress_max``.
        '''

        # If this range is invalid, raise an exception.
        if self._progress_min > self._progress_max:
            raise BetseCallbackException(
                'Minimum progress {} exceeds maximum progress {}.'.format(
                    progress_min, progress_max))
        elif self._progress_min == self._progress_max:
            raise BetseCallbackException(
                'Minimum progress {} equals maximum progress {}.'.format(
                    progress_min, progress_max))

        # Classify all passed parameters for subsequent reference by the
        # progressed() callback.
        self._progress_min = progress_min
        self._progress_max = progress_max

        # Value to be subsequently passed to the progressed() callback by the
        # first call to the progressed_next() callback. See the
        # progressed_next() callback for a justification of the "+ 1" here.
        self._progress_next = progress_min + 1
        # logs.log_debug('next progress: %d', self._progress_next)


    #FIXME: If "progress" == "self._progress_max", we should probably nullify
    #everything for sanity. The show only goes one way, folks.
    @type_check
    def progressed(self, progress: int) -> None:
        '''
        Callback passed the current state of progress for work performed by the
        source callable.

        Design
        ----------
        The source callable should call this callback repeatedly while
        performing any significant work. Doing so enables the sink callable to
        incrementally update external objects (e.g., progress bars).

        Subclasses are required to reimplement this callback and call this
        superclass implementation in their own reimplementation.

        Parameters
        ----------
        progress : int
            Current state of progress for work performed by the source
            callable. If this integer is *not* in the range previously passed
            to the :attr:`progress_ranged` callback, an exception is raised.

        Raises
        ----------
        :class:`BetseCallbackException`
            If either:

            * The :attr:`progress_ranged` callback has yet to be called by this
              source callable.
            * This progress is *not* in the range previously passed to the
              :attr:`progress_ranged` callback by this source callable.

        See Also
        ----------
        :meth:`progressed_next`
            Higher-level callback automatically passing the next progress value
            according to the the range previously passed to the
            :attr:`progress_ranged` callback.
        '''

        # If progress_ranged() has yet to be called, raise an exception.
        self._die_unless_progress_ranged()

        # If this progress is *NOT* in the range previously passed to
        # progress_ranged(), raise an exception.
        if not (self._progress_min <= progress <= self._progress_max):
            raise BetseCallbackException(
                'Progress {} not in range [{}, {}].'.format(
                    progress, self._progress_min, self._progress_max))


    def progressed_next(self) -> None:
        '''
        Higher-level callback wrapping the lower-level :meth:`progressed`
        callback by automatically passing the next progress value to that
        callback according to the the range previously passed to the
        :attr:`progress_ranged` callback.

        Guarantees
        ----------
        The first call to this callback is guaranteed to emit one greater than
        the minimum progress value previously passed to the
        :attr:`progress_ranged` callback. While arguably unintuitive, this
        behaviour corresponds with expectations in the sink callable.

        Assuming that the minimum progress value previously passed to the
        :attr:`progress_ranged` callback is 0 as is the common case, this
        behaviour has the additional advantage of guaranteeing that the
        maximum progress value previously passed to that callback is exactly
        equal to the number of times that this :meth:`progressed_next` callback
        is called by the source callable.

        Specifically, since the sink callable typically handles the
        :attr:`progress_ranged` callback by initializing a progress bar to the
        passed range, the sink callable is assumed to have already "consumed"
        the first progress value for all intents and purposes; instructing the
        sink callable to "reconsume" the same progress value would effectively
        hide the actual first iteration of work performed by the source
        callable from the end user perspective.

        Design
        ----------
        For generality, the sink callable should typically *only* handle the
        :meth:`progressed` callback. This callback is merely syntactic sugar
        conveniently simplifying the implementation of the source callback
        calling this callback.

        Raises
        ----------
        :class:`BetseCallbackException`
            If the :attr:`progress_ranged` callback has yet to be called by
            this source callable.
        '''

        # If progress_ranged() has yet to be called, raise an exception.
        self._die_unless_progress_ranged()

        # Notify the sink callback of the current state of progress.
        self.progressed(progress=self._progress_next)

        # If the source callback has work remaining to perform, increment this
        # state *AFTER* notifying this callback.
        if self._progress_next is not None:
            self._progress_next += 1

    # ..................{ EXCEPTIONS                         }..................
    def _die_unless_progress_ranged(self) -> None:
        '''
        Raise an exception unless the :attr:`progress_ranged` callback has been
        called by this source callable.

        Raises
        ----------
        :class:`BetseCallbackException`
            If the :attr:`progress_ranged` callback has yet to be called by
            this source callable.
        '''

        if (self._progress_min is None or
            self._progress_max is None or
            self._progress_next is None):
            raise BetseCallbackException(
                'progressed() callback called before '
                'progress_ranged() callback.')
