#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all matplotlib-based layer subclasses  spatially
plotting onto the cell cluster.
'''

#FIXME: The current approach to implementing animation overlays is totally BALLS
#fundamentally flawed. We currently attempt to provide a crude form of plot
#composition (i.e., merging two or more types of plots together into a single
#plot) by adding new booleans to the "AnimCellsABC" base class (e.g.,
#"is_current_overlayable") -- a fundamentally unwieldy and ultimately
#unworkable approach. By definition, you cannot provide true composability frow
#within a single class hierarchy. Instead, we need to split the specific
#process of generating different types of artists (e.g., mesh plots, stream
#plots) from the general process of animating and saving frames and plots as
#follows:
#
#* Refactor all concrete subclasses of "AnimCellsABC" into one or more
#  subclasses of "LayerCellsABC" instead, which may then be instantiated and
#  composed together into a new "layers" list passed to
#  CellsLayerABC.__init__(). For example:
#  * Split the existing "AnimGapJuncTimeSeries" subclass into:
#    * A new "CellsLayerGapJunc" subclass plotting *ONLY* the gap junction
#      open state as a "LineCollection" overlay. This layer subclass would
#      probably only be used for this specific purpose.
#    * A new "CellsLayerTimeSeries" subclass plotting *ONLY* an arbitrary
#      time series for the cell cluster as a mesh plot underlay. Clearly, this
#      layer subclass would be extensively reused elsewhere as well.
#* Replace all current overlay functionality in "AnimCellsABC" with "layers".
#* Refactor the configuration file from the current hard-coded non-composable
#  approach to a dynamic list-based approach permitting zero or more
#  user-defined animations, each consisting of one or more stock BETSE-defined
#  layers, to be defined. Users would then be able to construct arbitrarily
#  simple or complex animations as required.
#
#So, yes. It's quite a bit of work. But it's absolutely essential as well,
#particularly for implementing a general-purpose BETSE GUI.

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
from betse.exceptions import BetseSimVisualLayerException
from betse.util.io.log import logs
from betse.util.py import pyref
from betse.util.type import iterables, types
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.types import (
    type_check, IterableTypes, SequenceOrNoneTypes,)

# ....................{ SUPERCLASS                         }....................
class LayerCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting a single feature of
    the cell cluster for a parent visual.

    Each subclass of this class plots the spatial distribution of a single
    modelled variable (e.g., membrane voltage) for one or more simulation time
    steps. Each instance of the higher-level
    :class:`betse.science.visual.visabc.VisualCellsABC` abstract base class
    contains one or more instances of subclasses of this lower-level class.

    Separating low-level layer logic from high-level visual logic (e.g., frame
    iteration, video compression) enables composition between otherwise
    unrelated types. Thanks to layers, two or more types of visuals may be
    trivially composed into a unique third type of visual with *no* modification
    to existing layers or visuals.

    Attributes
    ----------
    _is_layered : bool
        ``True`` only if the :meth:`layer` method has been called at least once
        for this layer instance.
    _phase : SimPhase
        Current simulation phase *or* ``None`` if the :meth:`prep` method has
        yet to be called. Note that this attribute is also accessible via the
        :meth:`_visual.phase` property and is thus technically redundant. For
        both convenience and orthogonality with similar ``_phase`` attributes
        elsewhere in the codebase, this attribute is nonetheless defined and
        should be used in place of the :meth:`_visual.phase` property.
    _visual : VisualCellsABC
        Plot or animation to layer onto *or* ``None`` if the :meth:`prep` method
        has yet to be called.
    _zorder : int
        **Z-order** (i.e., positive integer ordering artist drawing, such that
        artists with larger z-orders are drawn above artists with smaller
        z-orders) of all artists plotted by this layer. While this zorder is
        *not* strictly enforced by this abstract base class, layer subclasses
        are encouraged to voluntarily pass the ``zorder=self._zorder`` option to
        all matplotlib axis-specific artist creation methods.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this layer.

        This method intentionally accepts *no* parameters except constants
        parametrizing this layer's behaviour. In particular, this method accepts
        references to neither the parent
        :class:`betse.science.visual.visabc.VisualCellsABC` instance containing
        this layer instance *nor* any other instances also contained by that
        parent instance (e.g., Matplotlib figure or axes objects). Why? Because
        plotters are instantiated by callers *before* their parent
        :class:`VisualCellsABC` instances are instantiated.

        See Also
        ----------
        :meth:`layer`
            Further details on class design.
        '''

        # Default instance attributes.
        self._is_layered = False
        self._phase = None
        self._visual = None
        self._zorder = None


    @type_check
    def prep(
        self,
        visual: 'betse.science.visual.visabc.VisualCellsABC',
        zorder: int,
    ) -> None:
        '''
        Prepare this layer to be layered onto the passed plot or animation.

        Parameters
        ----------
        visual : VisualCellsABC
            Plot or animation to layer onto.
        zorder : int
            **Z-order** (i.e., positive integer ordering artist drawing, such
            that artists with larger z-orders are drawn above artists with
            smaller z-orders) of all artists plotted by this layer.
        '''

        # Classify this visual with a weak rather than strong (the default)
        # reference, thereby avoiding circular references and the resulting
        # complications thereof (e.g., increased memory overhead).  Since the
        # parent plot or animation necessarily lives significantly longer than
        # this layer, no complications arise. Ergo, this attribute *ALWAYS*
        # yields this object (rather than non-deterministically yielding "None"
        # if this object is unexpectedly garbage-collected).
        self._visual = pyref.proxy_weak(visual)

        # Classify this zorder unmodified.
        self._zorder = zorder

        # Alias the current simulation phase to a convenience variable.
        self._phase = self._visual.phase

        # Ensure the next call to the layer() method calls the _layer_first()
        # rather than _layer_next() method, ensuring layers to be safely
        # reusable between multiple parent visuals.
        self._is_layered = False

    # ..................{ LAYERS                             }..................
    def layer(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the current time step and each cell of the current
        cluster onto the figure axes of the current plot or animation.
        '''

        # If this method has yet to be called...
        if not self._is_layered:
            # Perform logic specific to this call.
            self._layer_first()

            # Prevent subsequent calls to this method from repeating this logic.
            self._is_layered = True
        # Else, this method has been called at least once.
        else:
            # Perform logic specific to all subsequent calls.
            self._layer_next()


    # Layer subclasses are recommended but *NOT* required to reimplement this
    # empty method. Since simplistic layers plotting artists whose appearance is
    # constant across all time steps only require the abstract _layer_first()
    # method to be implemented, this method is concrete rather than abstract.
    def _layer_next(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the next simulation time step onto the figure axes
        of the current plot or animation.
        '''

        pass

    # ..................{ SUBCLASS                           }..................
    # Subclasses are required to implement the following abstract methods.

    @abstractmethod
    def _layer_first(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the first simulation time step onto the figure
        axes of the current plot or animation.
        '''

        pass

# ....................{ SUBCLASSES ~ colorful              }....................
class LayerCellsColorfulABC(LayerCellsABC):
    '''
    Abstract base class of all classes spatially plotting a single modelled
    variable of the cell cluster (e.g., cell membrane voltage) whose values are
    mappable as colors onto the colorbar for a parent plot or animation.

    Attributes
    ----------
    _color_mappables : IterableTypes
        Iterable of all mappables internally cached and returned by the
        :meth:`color_mappables` property.
    _color_max : NumericSimpleTypes
        Maximum color value to be displayed by the parent visual's colorbar.
    _color_min : NumericSimpleTypes
        Minimum color value to be displayed by the parent visual's colorbar.
        Ignored if :attr:`_visual.conf.is_color_autoscaled` is ``True``.
        Minimum color value to be displayed by the colorbar. Ignored if
        :attr:`_conf.is_color_autoscaled` is ``True``.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Default all instance attributes.
        self._color_mappables = None
        self._color_min = None
        self._color_max = None

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Generalize this class to support layers recomputing color mappables
    #each time step. To do so:
    #
    #* Define a new abstract _layer_next_color_mappables() method.
    #* Refactor all subclasses to override that rather than the _layer_next()
    #  method. (Ugh. Apologies on that one.)
    #* Refactor all subclass implementations of _layer_next_color_mappables()
    #  to explicitly return None, for safety.
    #* Override the _layer_next() method here as follows:
    #
    #    def _layer_next(self) -> None:
    #
    #        # Iterable of mappables layered by the subclass for the current time
    #        # step if any or None otherwise.
    #        color_mappables = self._layer_next_color_mappables()
    #
    #        if color_mappables is not None:
    #            self._color_mappables = color_mappables
    #
    #Great. We're not quite done, however. The question then becomes: how do we
    #propagate this change of color mappables to the parent visual, which will
    #be required to reset its colorbar accordingly? *UGH.*

    def _layer_first(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g.,
        transmembrane voltage) for the first sampled time step onto the current
        visual's figure axes.

        This method internally calls the subclass
        :meth:`_layer_first_color_mappables` method and caches the
        :attr:`_color_mappables` iterable returned by that call.

        Caveats
        ----------
        **Subclasses should not override this method,** which preserves various
        invariants required for colorbar mapping. Instead, subclasses should
        *only* override the :meth:`_layer_first_color_mappables`,
        :meth:`_layer_next` method, and :meth:`_layer_next_color_mappables`
        methods as needed.
        '''

        # Color mappables layered by the subclass for the first time step.
        self._color_mappables = self._layer_first_color_mappables()

        # If no mappables were layered, raise an exception.
        if not self._color_mappables:
            raise BetseSimVisualLayerException(
                '_layer_first_color_mappables() returned no color mappables.')

        # Add a colorbar to this visual associated with these mappables.
        self._make_colorbar()


    def _layer_next(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g.,
        transmembrane voltage) for the next sampled time step onto the current
        visual's figure axes.

        This method internally calls the subclass
        :meth:`_layer_next_color_mappables` method and caches the
        :attr:`_color_mappables` iterable returned by that call.

        Caveats
        ----------
        **Subclasses should typically override this method,** which
        inefficiently recreates all color mappables each time step, which most
        layer subclasses do *not* require. To efficiently preserve the initial
        color mappables created and cached by the prior call to the
        :meth:`_layer_first` method, most layer subclasses should instead
        override the default implementation of this method rather than
        implementing the optional :meth:`_layer_next_color_mappables` method.
        '''

        # Color mappables layered by the subclass for the next time step.
        self._color_mappables = self._layer_next_color_mappables()

        # If no mappables were layered, raise an exception.
        if not self._color_mappables:
            raise BetseSimVisualLayerException(
                '_layer_next_color_mappables() returned no color mappables.')

        # Rescale these mappables to the previously established minimum and
        # maximum color values. Since matplotlib provides no convenient means
        # for efficiently replacing the existing colorbar added by the prior
        # call to the _layer_first() method *AND* since replacing this colorbar
        # would be generally undesirable, this colorbar is preserved as is.
        self._scale_color_mappables()

    # ..................{ SUBCLASS ~ mandatory               }..................
    # The following abstract methods *MUST* be implemented by subclasses.

    @abstractproperty
    def color_data(self) -> SequenceOrNoneTypes:
        '''
        Sequence of arbitrary dimensions of all possible color values for all
        time steps plotted by this layer *or* ``None`` if calculating these
        values is impractical (e.g., due to space or time constraints).

        If colorbar autoscaling is:

        * Disabled, this sequence is ignored.
        * Enabled *and* this sequence is:
          * ``None``, this sequence is ignored. In this case, the subclass is
            responsible for colorbar autoscaling.
          * Non-``None``, the colorbar is clipped to the minimum and maximum
            scalar values unravelled from this sequence.
        '''

        pass


    @abstractmethod
    def _layer_first_color_mappables(self) -> IterableTypes:
        '''
        Layer the spatial distribution of a single modelled variable (e.g.,
        transmembrane voltage) for the first sampled time step onto the current
        visual's figure axes, returning all matplotlib color mappables whose
        artists are to be mapped onto (i.e., coloured according to) this
        figure's colorbar.

        Returns
        ----------
        IterableTypes
            Iterable of all matplotlib color mappables, cached into the
            :attr:`_color_mappables` attribute by the parent
            :meth:`_layer_first` method. Note that, by matplotlib design, only
            the first mappable in this iterable is arbitrarily associated with
            this colorbar; all other mappables are ignored for this purpose.
        '''

        pass

    # ..................{ SUBCLASS ~ optional                }..................
    # The following non-abstract methods *MAY* be implemented by subclasses.

    def _layer_next_color_mappables(self) -> IterableTypes:
        '''
        Layer the spatial distribution of a single modelled variable (e.g.,
        transmembrane voltage) for the next sampled time step onto the current
        visual's figure axes, returning all matplotlib color mappables whose
        artists are to be mapped onto (i.e., coloured according to) this
        figure's colorbar.

        Design
        ----------
        All layer subclasses *must* implement either:

        * The :meth:`layer_next` method, in which case this layer is assumed to
          preserve the original color mappables returned by the initial call to
          the :meth:`_layer_first_color_mappables` method for each time step.
          This is the common (and most efficient) case.
        * This method, in which case this layer is assumed to behave as follows
          for each time step:

          * Destroy all color mappables created by the prior call to this or the
            :meth:`_layer_first_color_mappables` method.
          * Create and return new color mappables from the call to this method.

        Since the latter is substantially less efficient than the former, *only*
        layers recreating color mappables must implement this method; all other
        layers should implement the :meth:`layer_next` method.

        Returns
        ----------
        IterableTypes
            Iterable of all matplotlib color mappables, cached into the
            :attr:`_color_mappables` attribute by the parent
            :meth:`_layer_next` method.
        '''

        pass

    # ..................{ COLORS                             }..................
    def _make_colorbar(self) -> None:
        '''
        Add a colorbar to the parent visual's figure.

        This method automatically configures the range of colors displayed by
        this colorbar *and* associates all color mappables defined by this
        layer subclass with this colorbar.
        '''

        # Set the minimum and maximum colorbar values.
        self._set_color_range()

        # Scale all color mappables defined by this subclass by these values.
        self._scale_color_mappables()

        # First color mappable defined by this subclass.
        color_mappable_first = iterables.get_item_first(self._color_mappables)

        # Create a colorbar associated with this color mappable.
        self._visual.make_colorbar(color_mappable_first)


    def _set_color_range(self) -> None:
        '''
        Set the minimum and maximum color values to be displayed by the colorbar
        for the parent visual's figure.

        Specifically, these values are set to:

        * If colorbar autoscaling is enabled by this visual's configuration,
          the minimum and maximum values unravelled from the
          :meth:`color_data` Numpy array defined by this layer subclass.
        * Else, the minimum and maximum values hardcoded into this visual's
          configuration.
        '''

        # If colorbar autoscaling is enabled by this visual's configuration...
        if self._visual.conf.is_color_autoscaled:
            # One-dimensional Numpy array of all colors flattened from this
            # possibly multi-dimensional Numpy array of these values.
            color_data_flat = self.color_data.ravel()

            # Set the minimum and maximum colors to the minimum and maximum
            # values in this array.
            self._color_min = color_data_flat.min()
            self._color_max = color_data_flat.max()
            # self._color_min = np.ma.min(color_data_flat)
            # self._color_max = np.ma.max(color_data_flat)
        # Else, colorbar autoscaling is disabled. In this case, set the minimum
        # and maximum colors to those hardcoded into this configuration.
        else:
            self._color_min = self._visual.conf.color_min
            self._color_max = self._visual.conf.color_max

        # If these values are identical, coerce them to differ. Failing to do
        # so produces spurious visual artifacts in both the axes and colorbar.
        # if self._color_min == self._color_max:
        #     self._color_min = self._color_min - 1
        #     self._color_max = self._color_max + 1

        # Ensure sanity.
        assert types.is_numeric(self._color_min), (
            types.assert_not_numeric(self._color_min))
        assert types.is_numeric(self._color_max), (
            types.assert_not_numeric(self._color_max))


    def _scale_color_mappables(self) -> None:
        '''
        Scale all color mappables defined by this layer subclass to the current
        minimum and maximum color values.
        '''

        # Log this attempt.
        logs.log_debug(
            'Rescaling "%s" colors to [%d, %d]...',
            self._visual.name, self._color_min, self._color_max)

        # For each color mappable, clip that mappable to the minimum and
        # maximum values discovered above. Note this also has the beneficial
        # side-effect of establishing the colorbar's range.
        for color_mappable in self._color_mappables:
            # Ensure sanity.
            assert types.is_matplotlib_mappable(color_mappable), (
                types.assert_not_matplotlib_mappable(color_mappable))

            # Clip this mappable.
            color_mappable.set_clim(self._color_min, self._color_max)
