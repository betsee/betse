#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Post-simulation cell cluster plot pipelines** (i.e., pipelines plotting
simulated data of all cells in the cell cluster).
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.config.export.visual.confvisabc import (
    SimConfVisualCellsListItem)
from betse.science.export import expmath
from betse.science.phase.pipe.piperun import piperunner
from betse.science.phase.require import phasereqs
from betse.science.visual.plot import plotutil
from betse.science.visual.plot.pipe.plotpipeabc import PlotPipeABC
from betse.util.io.log import logs
from betse.util.type import iterables
from betse.util.type.types import IterableTypes
from collections import OrderedDict
from matplotlib import pyplot as pyplot
from matplotlib.collections import LineCollection, PolyCollection
# from scipy.ndimage.filters import gaussian_filter

# ....................{ SUBCLASSES                         }....................
#FIXME: Rename to "SimPipeExportPlotsCells" for disambiguity.
class PlotCellsPipe(PlotPipeABC):
    '''
    **Post-simulation cell cluster plot pipeline** (i.e., object iteratively
    displaying and/or saving all plots depicting all cells of the cell cluster,
    produced after initialization and simulation as specified by the current
    simulation configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, label_singular='cell cluster plot', **kwargs)

    # ..................{ SUPERCLASS                         }..................
    @property
    def _runners_conf(self) -> IterableTypes:
        return self._phase.p.plot.plots_cells_after_sim

    # ..................{ EXPORTERS ~ channel                }..................
    # @piperunner(
    #     categories=('Ion Channel', 'Density Factor',),
    #     requirements={phasereqs.ELECTROOSMOSIS,},
    # )
    # def export_channel_density(self, conf: SimConfVisualCellsListItem) -> None:
    #     '''
    #     Plot all cell membrane ion channel density factors for the cell cluster
    #     at the last time step.
    #     '''
    #
    #     # Prepare to export the current plot.
    #     self._export_prep()
    #
    #     plotutil.plotMemData(
    #         self._phase.cells,
    #         self._phase.p,
    #         zdata=self._phase.sim.rho_channel,
    #         clrmap=self._phase.p.default_cm,
    #     )
    #     pyplot.xlabel('Spatial Dimension [um]')
    #     pyplot.ylabel('Spatial Dimension [um]')
    #     pyplot.title('Membrane ion channel density factor')
    #
    #     # Export this plot to disk and/or display.
    #     self._export(basename='final_channels_2D')

    # ..................{ EXPORTERS ~ cluster                }..................
    # This exporter is solver- and feature-agnostic.

    @piperunner(categories=('Cell Cluster', 'Mask',))
    def export_cluster_mask(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot the **cell cluster image mask** (i.e., user-defined image whose
        pure-black pixels exactly correspond to the shape of the cell cluster).

        This plot is irrespective of time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        fig = pyplot.figure()
        pyplot.imshow(
            self._phase.cells.maskM,
            origin='lower',
            extent=self._phase.cache.upscaled.extent,
            cmap=self._phase.p.background_cm,
        )
        pyplot.colorbar()
        pyplot.title('Cluster Masking Matrix')

        # Export this plot to disk and/or display.
        self._export(basename='cluster_mask')

    # ..................{ EXPORTERS ~ cluster : tissue       }..................
    # This exporter is solver- and feature-agnostic.
    @piperunner(categories=('Cell Cluster', 'Tissue and Cut Profiles',))
    def export_tissue_cuts(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot a **tissue and cut profile tessellation** (i.e., tiled mosaic of
        all cells spatially subdivided into tissue and cut profile regions such
        that all cells in the same region share the same arbitrary color) for
        the cell cluster.

        This plot is irrespective of time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        # Localize frequently accessed fields for convenience.
        p = self._phase.p
        dyna = self._phase.dyna
        cells = self._phase.cells
        colormap = self._phase.p.background_cm

        col_dic = {}
        cb_ticks = []
        cb_tick_labels = []

        # Ordered dictionary mapping from the name of each tissue and cut
        # profile to a one-dimensional Numpy array of the 0-based indices of all
        # cells owned by this profile.
        #
        # Note that order is absolutely significant. The first tissue profile is
        # a pseudo-tissue defining the cell cluster itself. If order were *NOT*
        # preserved here, this tissue would be assigned an arbitrary z-order, in
        # which case all tissues assigned smaller z-orders would be entirely
        # obscured by that pseudo-tissue covering the cell cluster.
        profile_name_to_cells_index = OrderedDict()

        fig = pyplot.figure()
        ax = pyplot.subplot(111)

        # For the name and object encapsulating each tissue profile...
        for tissue_name, _ in dyna.tissue_name_to_profile.items():
            # One-dimensional Numpy array of the indices of all tissue cells.
            tissue_cells_index = dyna.cell_target_inds[tissue_name]

            # If this tissue contains no cells, skip to the next tissue.
            if not len(tissue_cells_index):
                logs.log_warning('Tissue "%s" contains no cells.', tissue_name)
                continue

            # Else, tissue contains one or more cells. Map this tissue to these
            # indices (for subsequent lookup).
            profile_name_to_cells_index[tissue_name] = tissue_cells_index

        #FIXME: The "p.plot_cutlines" boolean is only ever leveraged here and
        #thus is arguably extraneous. Consider refactoring as follows:
        #
        #* Remove the "p.plot_cutlines" boolean and corresponding YAML-formetted
        #  default configuration option.
        #* Split the existing "tissue_cuts" plot type in the "cell cluster
        #  pipeline" into the following two types:
        #  * "tissue", unconditionally plotting *ONLY* tissue profiles.
        #  * "tissue_cuts", unconditionally plotting both tissue and cut
        #    profiles.
        #* Define a new private _export_profiles() method as follows:
        #      @type_check
        #      _export_profiles(
        #          self,
        #          conf: SimConfVisualCellsListItem,
        #          is_tissue: bool,
        #          is_cuts: bool
        #       ) -> None:
        #* Implement export_tissue() to call _export_profiles() as:
        #    self._export_profiles(
        #        self,
        #        conf=conf,
        #        is_tissue=True,
        #        is_cuts=False,
        #    )
        #* Implement export_tissue_cuts() similarly.

        # If plotting cut profiles as well *AND* the cutting event is enabled...
        if p.plot_cutlines and dyna.event_cut is not None:
            # For the name and object encapsulating each cut profile...
            for cut_name, cut_profile in dyna.cut_name_to_profile.items():
                # Map this cut to the indices of all cells cut by this profile.
                profile_name_to_cells_index[cut_name] = (
                    cut_profile.picker.pick_cells(cells=cells, p=p))

        # Minimum and maximum 0-based integers uniquely identifying the first
        # and last tissue and cut profile (respoctively), localized for ordering
        # purposes in the colorbar legend.
        profile_zorder = 0
        profile_zorder_max = len(profile_name_to_cells_index)

        # For the name and one-dimensional Numpy array of the 0-based indices of
        # all cells in each tissue and/or cut profile...
        for profile_name, profile_cells_index in (
            profile_name_to_cells_index.items()):
            # logs.log_debug('Plotting tissue "%s"...', profile_name)
            profile_zorder += 1

            profile_points = expmath.upscale_coordinates(
                cells.cell_verts[profile_cells_index])

            z = np.zeros(len(profile_points))
            z[:] = profile_zorder

            col_dic[profile_name] = PolyCollection(
                profile_points, array=z, cmap=colormap, edgecolors='none')
            col_dic[profile_name].set_clim(0, profile_zorder_max)

            # col_dic[profile_name].set_alpha(0.8)
            col_dic[profile_name].set_zorder(profile_zorder)
            ax.add_collection(col_dic[profile_name])

            # Add this profile name to the colour legend.
            cb_ticks.append(profile_zorder)
            cb_tick_labels.append(profile_name)

        # logs.log_debug('Plotting colorbar ticks: %r', cb_ticks)
        # logs.log_debug('Plotting colorbar tick labels: %r', cb_tick_labels)

        ax_cb = None
        if dyna.tissue_name_to_profile:
            # Name of the first tissue profile.
            tissue_first_name = iterables.get_item_first(
                dyna.tissue_name_to_profile.keys())

            # Color mappable associated with this tissue profile, guaranteed in
            # this case to be a "PolyCollection" instance.
            tissue_first_mappable = col_dic[tissue_first_name]

            ax_cb = fig.colorbar(tissue_first_mappable, ax=ax, ticks=cb_ticks)
            ax_cb.ax.set_yticklabels(cb_tick_labels)

        if p.enumerate_cells:
            for i, cll in enumerate(cells.cell_centres):
                ax.text(
                    p.um*cll[0], p.um*cll[1], i,
                    ha='center', va='center', zorder=20)

        ax.set_xlabel('Spatial Distance [um]')
        ax.set_ylabel('Spatial Distance [um]')
        ax.set_title('Cell Cluster')

        ax.axis('equal')
        ax.axis(self._phase.cache.upscaled.extent)

        # Export this plot to disk and/or display.
        self._export(basename='cluster_mosaic')

    # ..................{ EXPORTERS ~ current                }..................
    @piperunner(
        categories=('Current Density', 'Intracellular',),
        requirements=phasereqs.ELECTRIC_CURRENT,
    )
    def export_currents_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular current densities for the cell cluster at the
        last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figI, axI, cbI = plotutil.plotStreamField(
            100*self._phase.sim.J_cell_x,
            100*self._phase.sim.J_cell_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=False,
            title='Intracellular Current Density',
            cb_title='Current Density [uA/cm2]',
            show_cells=False,
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        axI.set_xlabel('Spatial distance [um]')
        axI.set_ylabel('Spatial distance [um]')
        cbI.set_label('Current Density [uA/cm2]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Current_gj')


    @piperunner(
        categories=('Current Density', 'Extracellular',),
        requirements=phasereqs.ELECTRIC_CURRENT_EXTRA,
    )
    def export_currents_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular current densities for the environment at the last
        time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figI2, axI2, cbI2 = plotutil.plotStreamField(
            100*self._phase.sim.J_env_x,
            100*self._phase.sim.J_env_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=True,
            title='Extracellular Current Density',
            cb_title='Current Density [uA/cm2]',
            show_cells=False,
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        axI2.set_xlabel('Spatial distance [um]')
        axI2.set_ylabel('Spatial distance [um]')
        cbI2.set_label('Extracellular Current Density [uA/cm2]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Current_extracellular')

    # ..................{ EXPORTERS ~ deform                 }..................
    @piperunner(
        categories=('Deformation', 'Total',),
        requirements=phasereqs.DEFORM,
    )
    def export_deform_total(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all **total cellular displacements** (i.e., summations of all
        cellular deformations due to galvanotropic and osmotic pressure body
        forces) for the cell cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        plotutil.plotStreamField(
            self._phase.p.um*self._phase.sim.dx_cell_time[-1],
            self._phase.p.um*self._phase.sim.dy_cell_time[-1],
            self._phase.cells, self._phase.p,
            plot_ecm=False,
            title='Final Displacement of Cell Collective',
            cb_title='Displacement [um]',
            show_cells=self._phase.p.showCells,
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='final_displacement_2D')

    # ..................{ EXPORTERS ~ diffusion              }..................
    @piperunner(
        categories=('Diffusion Weights', 'Extracellular',),
        requirements=phasereqs.ECM,
    )
    def export_diffusion_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all logarithm-scaled extracellular diffusion weights for the
        environment at the last time step.
        '''

        # Prepare to export the electric plot.
        self._export_prep()

        fig = pyplot.figure()
        ax99 = pyplot.subplot(111)
        pyplot.imshow(
            np.log10(self._phase.sim.D_env_weight.reshape(
                self._phase.cells.X.shape)),
            origin='lower',
            extent=self._phase.cache.upscaled.extent,
            cmap=self._phase.p.background_cm,
        )
        pyplot.colorbar()

        cell_edges_flat = expmath.upscale_coordinates(
            self._phase.cells.mem_edges_flat)
        coll = LineCollection(cell_edges_flat, colors='k')
        coll.set_alpha(1.0)
        ax99.add_collection(coll)

        pyplot.title('Logarithm of Environmental Diffusion Weight Matrix')

        # Export this plot to disk and/or display.
        self._export(basename='env_diffusion_weights')

    # ..................{ EXPORTERS ~ electric               }..................
    @piperunner(
        categories=('Electric Field', 'Intracellular',),
        requirements=phasereqs.ELECTRIC_FIELD,
    )
    def export_electric_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular electric field lines for the cell cluster at the
        last time step.
        '''

        # Prepare to export the electric plot.
        self._export_prep()

        plotutil.plotVectField(
            self._phase.sim.E_cell_x,
            self._phase.sim.E_cell_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=False,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='Final_Electric_Field_GJ')


    @piperunner(
        categories=('Electric Field', 'Extracellular',),
        requirements=phasereqs.ELECTRIC_FIELD_EXTRA,
    )
    def export_electric_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular electric field lines for the environment at the
        last time step.
        '''

        # Prepare to export the electric plot.
        self._export_prep()

        plotutil.plotVectField(
            self._phase.sim.E_env_x,
            self._phase.sim.E_env_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=True,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='Final_Electric_Field_ECM')

    # ..................{ EXPORTERS ~ fluid                  }..................
    @piperunner(
        categories=('Fluid Flow', 'Intracellular',),
        requirements=phasereqs.FLUID,
    )
    def export_fluid_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular fluid flow field lines for the cell cluster at
        the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        plotutil.plotStreamField(
            1e6*self._phase.sim.u_cells_x,
            1e6*self._phase.sim.u_cells_y,
            self._phase.cells, self._phase.p,
            plot_ecm=False,
            title='Final Fluid Velocity in Cell Collective',
            cb_title='Velocity [um/s]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='final_vel_2D_gj')


    @piperunner(
        categories=('Fluid Flow', 'Extracellular',),
        requirements=phasereqs.FLUID_EXTRA,
    )
    def export_fluid_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular fluid flow field lines for the environment at the
        last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        plotutil.plotStreamField(
            1e9*self._phase.sim.u_env_x,
            1e9*self._phase.sim.u_env_y,
            self._phase.cells, self._phase.p,
            plot_ecm=True,
            title='Final Fluid Velocity in Environment',
            cb_title='Velocity [nm/s]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='final_vel_2D_env')

    # ..................{ EXPORTERS ~ ion : calcium          }..................
    @piperunner(
        categories=('Ion Concentration', 'Calcium', 'Intracellular',),
        requirements=phasereqs.ION_CALCIUM,
    )
    def export_ion_calcium_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular calcium (i.e., Ca2+) ion concentrations for the
        cell cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figCa, axCa, cbCa = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=self._phase.sim.cc_time[-1][self._phase.sim.iCa]*1e6,
            number_cells=self._phase.p.enumerate_cells,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        # Export this plot to disk and/or display.
        self._export(basename='final_Ca_2D')


    @piperunner(
        categories=('Ion Concentration', 'Calcium', 'Extracellular',),
        requirements=phasereqs.ION_CALCIUM_EXTRA,
    )
    def export_ion_calcium_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular calcium (i.e., Ca2+) ion concentrations for the
        environment at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        cc_Ca = self._phase.sim.cc_env[self._phase.sim.iCa].reshape(
            self._phase.cells.X.shape)

        pyplot.figure()
        pyplot.imshow(
            cc_Ca,
            origin='lower',
            extent=self._phase.cache.upscaled.extent,
            cmap=self._phase.p.default_cm,
        )
        pyplot.colorbar()
        pyplot.title('Environmental Calcium [mmol/L]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_environmental_calcium')

    # ..................{ EXPORTERS ~ gap junction           }..................
    # These exporters are solver- and feature-agnostic.

    @piperunner(categories=('Gap Junction', 'Connectivity Network',))
    def export_gj_connectivity(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot the **gap junction connectivity network** (i.e., graph of all gap
        junctions connecting cell membranes) for the cell cluster.

        This plot is irrespective of time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        fig = pyplot.figure()
        ax_x = pyplot.subplot(111)

        if self._phase.p.showCells:
            base_points = expmath.upscale_coordinates(
                self._phase.cells.cell_verts)
            col_cells = PolyCollection(
                base_points, facecolors='k', edgecolors='none')
            col_cells.set_alpha(0.3)
            ax_x.add_collection(col_cells)

        connects = expmath.upscale_coordinates(self._phase.cells.nn_edges)
        collection = LineCollection(connects, linewidths=1.0, color='b')
        ax_x.add_collection(collection)
        pyplot.axis('equal')
        pyplot.axis(self._phase.cache.upscaled.extent)

        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Cell Connectivity Network')

        # Export this plot to disk and/or display.
        self._export(basename='gj_connectivity_network')


    @piperunner(categories=('Gap Junction', 'Relative Permeability',))
    def export_gj_permeability(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all **gap junction connectivity states** (i.e., relative
        permeabilities of all gap junctions connecting cell membranes) for the
        cell cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        fig_x = pyplot.figure()
        ax_x = pyplot.subplot(111)

        con_segs = self._phase.cells.nn_edges
        connects = self._phase.p.um*np.asarray(con_segs)
        collection = LineCollection(
            connects,
            array=self._phase.sim.gjopen,
            cmap=self._phase.p.background_cm,
            linewidths=2.0,
        )
        # collection.set_clim(0, 1)

        ax_x.add_collection(collection)
        cb = fig_x.colorbar(collection)
        pyplot.axis('equal')
        pyplot.axis(self._phase.cache.upscaled.extent)

        cb.set_label('Relative Permeability')
        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Final Gap Junction Relative Permeability')

        # Export this plot to disk and/or display.
        self._export(basename='final_gjState')

    # ..................{ EXPORTERS ~ microtubule            }..................
    @piperunner(
        categories=('Microtubule', 'Coherence',),
        requirements=phasereqs.MICROTUBULE,
    )
    def export_microtubule(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot the coherence of all cellular microtubules for the cell cluster at
        the last time step.
        '''

        # Prepare to export the microtubules plot.
        self._export_prep()

        umtx, umty = self._phase.sim.mtubes.mtubes_to_cell(self._phase.cells, self._phase.p)

        pyplot.figure()
        ax = pyplot.subplot(111)

        plotutil.plotVectField(
            umtx,
            umty,
            self._phase.cells,
            self._phase.p,
            plot_ecm=False,
            title='Final Microtubule Alignment Field',
            cb_title='Aligned MT Fraction',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # pyplot.figure()
        # ax = pyplot.subplot(111)
        # umtx, umty = self._phase.sim.mtubes.mtubes_to_cell(self._phase.cells, self._phase.p)
        # plotutil.cell_quiver(umtx, umty, ax, self._phase.cells, self._phase.p)

        # plotutil.mem_quiver(
        #     self._phase.sim.mtubes.mtubes_x,
        #     self._phase.sim.mtubes.mtubes_y,
        #     ax,
        #     self._phase.cells,
        #     self._phase.p,
        # )

        ax.set_xlabel('X-Distance [um]')
        ax.set_ylabel('Y-Distance [um]')
        ax.set_title('Net microtubule alignment in cells')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Microtubules')

    # ..................{ EXPORTERS ~ pump                   }..................
    # @piperunner(
    #     categories=('Ion Pump', 'Density Factor',),
    #     requirements={phasereqs.ELECTROOSMOSIS,},
    # )
    # def export_pump_density(self, conf: SimConfVisualCellsListItem) -> None:
    #     '''
    #     Plot all cell membrane ion pump density factors for the cell cluster at
    #     the last time step.
    #     '''
    #
    #     # Prepare to export the current plot.
    #     self._export_prep()
    #
    #     plotutil.plotMemData(
    #         self._phase.cells, self._phase.p,
    #         zdata=self._phase.sim.rho_pump,
    #         clrmap=self._phase.p.default_cm,
    #     )
    #     pyplot.xlabel('Spatial Dimension [um]')
    #     pyplot.ylabel('Spatial Dimension [um]')
    #     pyplot.title('Membrane ion pump density factor')
    #
    #     # Export this plot to disk and/or display.
    #     self._export(basename='final_pumps_2D')

    @piperunner(
        categories=('Ion Pump', 'Pump Rate', 'Na-K-ATPase',),
        requirements=phasereqs.PUMP_NAKATPASE,
    )
    def export_pump_nakatpase(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all cell membrane Na-K-ATPase pump rates for the cell cluster at
        the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        pumpData = self._phase.sim.rate_NaKATP*1e9

        figPump, axPump, cbPump = plotutil.plotPrettyPolyData(
            pumpData,
            self._phase.sim, self._phase.cells, self._phase.p,
            number_cells=self._phase.p.enumerate_cells,
            clrmap=self._phase.p.default_cm,
        )

        axPump.set_title('Final Na/K-ATPase Pump Rate')
        axPump.set_xlabel('Spatial distance [um]')
        axPump.set_ylabel('Spatial distance [um]')
        cbPump.set_label('Pump Na+ Flux (nmol/m2*s)')

        # Export this plot to disk and/or display.
        self._export(basename='final_NaKPump_2D')

    # ..................{ EXPORTERS ~ pressure               }..................
    @piperunner(
        categories=('Pressure', 'Total',),
        requirements=phasereqs.PRESSURE_TOTAL,
    )
    def export_pressure_total(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all **cellular pressure totals** (i.e., summations of all cellular
        mechanical and osmotic pressures) for the cell cluster at the last time
        step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figP, axP, cbP = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=self._phase.sim.P_cells,
            number_cells=self._phase.p.enumerate_cells,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        axP.set_title('Final Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure [Pa]')

        # Export this plot to disk and/or display.
        self._export(basename='final_P_2D_gj')

    # ..................{ EXPORTERS ~ voltage                }..................
    @piperunner(
        categories=('Voltage', 'Extracellular',),
        requirements=phasereqs.VOLTAGE_EXTRA,
    )
    def export_voltage_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular voltages for the environment at the last time
        step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        vv = self._phase.sim.v_env.reshape(self._phase.cells.X.shape)

        pyplot.figure()
        pyplot.imshow(
            1e3*vv,
            origin='lower',
            extent=self._phase.cache.upscaled.extent,
            cmap=self._phase.p.default_cm,
        )
        pyplot.colorbar()
        pyplot.title('Environmental Voltage [mV]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_environmental_V')


    @piperunner(
        categories=('Voltage', 'Polarity',),
        requirements=phasereqs.VOLTAGE_POLARITY,
    )
    def export_voltage_polarity(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all cellular voltage polarities for the cell cluster at the last
        time step.
        '''

        # Prepare to export the polarization plot.
        self._export_prep()

        pyplot.figure()

        # Plot a background Vmem mesh.
        fig, ax, cb = plotutil.plotPrettyPolyData(
            1000*self._phase.sim.vm_time[-1],
            self._phase.sim, self._phase.cells, self._phase.p,
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        # Calculate the Vmem polarity vectors.
        polm = self._phase.sim.vm - (
            self._phase.sim.vm_ave_time[-1][self._phase.cells.mem_to_cells])
        polx = polm*self._phase.cells.mem_vects_flat[:,2]
        poly = polm*self._phase.cells.mem_vects_flat[:,3]

        pcx = np.dot(
            self._phase.cells.M_sum_mems,
            polx*self._phase.cells.mem_sa) / self._phase.cells.cell_sa
        pcy = np.dot(
            self._phase.cells.M_sum_mems,
            poly*self._phase.cells.mem_sa) / self._phase.cells.cell_sa

        plotutil.cell_quiver(pcx, pcy, ax, self._phase.cells, self._phase.p)

        ax.set_xlabel('X-Distance [um]')
        ax.set_ylabel('Y-Distance [um]')
        ax.set_title('Cell Vmem polarity')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Polarity')

    # ..................{ EXPORTERS ~ voltage : vmem         }..................
    @piperunner(categories=('Voltage', 'Transmembrane', 'Actual',))
    def export_voltage_membrane(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all transmembrane voltages (Vmem) for the cell cluster at the last
        time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figV, axV, cbV = plotutil.plotPrettyPolyData(
            1000*self._phase.sim.vm_time[-1],
            self._phase.sim, self._phase.cells, self._phase.p,
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        figV.suptitle('Final Vmem', fontsize=14, fontweight='bold')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')

        # Export this plot to disk and/or display.
        self._export(basename='final_Vmem_2D')


    @piperunner(categories=('Voltage', 'Transmembrane', 'Average',))
    def export_voltage_membrane_average(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot the averages of all transmembrane voltages (Vmem) for the cell
        cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figVa, axVa, cbVa = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=1000*self._phase.sim.vm_ave,
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        # axVa.quiver(
        #     p.um*cells.cell_centres[:,0],
        #     p.um*cells.cell_centres[:,1], sim.pol_cell_x, sim.pol_cell_y)

        figVa.suptitle('Final Average Vmem', fontsize=14, fontweight='bold')
        axVa.set_xlabel('Spatial distance [um]')
        axVa.set_ylabel('Spatial distance [um]')
        cbVa.set_label('Voltage [mV]')

        # Export this plot to disk and/or display.
        self._export(basename='final_AverageVmem_2D')


    @piperunner(
        categories=('Voltage', 'Transmembrane', 'GHK',),
        requirements=phasereqs.VOLTAGE_MEMBRANE_GHK,
    )
    def export_voltage_membrane_ghk(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all transmembrane voltages (Vmem) calculated by the
        Goldman-Hodgkin-Katz (GHK) equation for the cell cluster at the last
        time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figV_ghk, axV_ghk, cbV_ghk = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=1000*self._phase.sim.vm_GHK_time[-1],
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        figV_ghk.suptitle(
            'Final Vmem using Goldman Equation', fontsize=14, fontweight='bold')
        axV_ghk.set_xlabel('Spatial distance [um]')
        axV_ghk.set_ylabel('Spatial distance [um]')
        cbV_ghk.set_label('Voltage [mV]')

        # Export this plot to disk and/or display.
        self._export(basename='final_Vmem_GHK_2D')
