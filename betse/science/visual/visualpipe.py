#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation plots and animations.
'''

# ....................{ IMPORTS                            }....................
from betse.science.simulate.simphase import SimPhaseABC
from betse.science.visual.anim import animpipe
from betse.science.visual.plot import plotpipe
from betse.util.io.log import logs
from betse.util.type.types import type_check
from matplotlib import pyplot as plt

# ....................{ PIPELINES                          }....................
#FIXME: Refactor the "plot_type" parameter to be somewhat less crazy. This
#parameter is passed by the "simrunner" submodule. Ideally, this parameter
#should be a typesafe enum rather than a non-typesafe string or, preferably,
#simply go away entirely. Related commentary follows on how to best achieve the
#latter goal.
#FIXME: I don't quite grok our usage of "sim.run_sim". This undocumented
#attribute appears to be internally set by the Simulator.run_phase_sans_ecm()
#method. That makes sense; however, what's the parallel "p.run_sim" attribute
#for, then?  Interestingly, the "SimRunner" class sets "p.run_sim" as follows:
#
#* To "False" if an initialization is being performed.
#* To "True" if a simulation is being performed.
#
#This doesn't seem quite ideal, however. Ideally, there would exist one and only
#one attribute whose value is an instance of a multi-state "SimPhaseKind" class
#rather than two binary boolean attributes. Possible enum values might include:
#
#* "SimPhaseKind.seed" when seeding a new cluster.
#* "SimPhaseKind.init" when initializing a seeded cluster.
#* "SimPhaseKind.sim" when simulating an initialized cluster.
#
#This attribute would probably exist in the "Simulator" class -- say, as
#"sim.phase". In light of that, consider the following refactoring:
#
#* Define a new "SimPhaseKind" class in the "sim" module with the above
#  attributes.
#* Define a new "Simulator.phase" attribute initialized to None.
#* Replace all existing uses of the "p.run_sim" and "sim.run_sim" booleans with
#  "sim.phase" instead. Note that only the:
#  * "SimRunner" class sets "p.run_sim".
#  * "Simulator" class sets "sim.run_sim".
#
#Note also the:
#
#* "plot_type" parameter passed to the pipeline_plots() function by the
#  "SimRunner" class.
#* The seemingly duplicate "p.plot_type" attribute internally set by the
#  pipeline_plots() function, which is frankly crazy.
#
#Both parameters should probably receive similar treatment and be replaced
#entirely by use of the new "sim.phase" attribute.
#
#Wonder temptress at the speed of light and the sound of love!

#FIXME: Refactor to accept only a single "SimPhaseABC" instance.
@type_check
def pipeline(phase: SimPhaseABC, plot_type: str) -> None:
    '''
    Display and/or save all currently enabled plots and animations for the
    passed simulation phase.

    Parameters
    ----------------------------
    phase: SimPhaseABC
        Current simulation phase.
    plot_type : str
        String constant corresponding to the current simulation phase. Valid
        values include:
        * `init`, for plotting simulation initialization results.
        * `sim`, for plotting simulation run results.
    '''

    #FIXME: This is terrible. I don't even.
    phase.p.plot_type = plot_type

    # Display and/or save all animations.
    animpipe.pipeline(phase)

    #FIXME: Display plots *BEFORE* animations. Doing so safely will probably
    #necessitate resolving why we're calling plt.show() below, however.

    # Display and/or save all plots.
    plotpipe.pipeline(phase)

    #FIXME: What is this? What requires showing? Are we finalizing some
    #previously displayed visual artifact? We suspect this to be safely
    #jettisoned deadweight, but... let's verify that, please. If this *IS*
    #required, it's probably only required for plots -- in which case this logic
    #should be appended onto the plotpipe.pipeline() function.

    # If displaying plots and animations, display... something? I guess?
    if phase.p.turn_all_plots_off is False:
        plt.show()
    # Else, log the directory to which results were exported.
    else:
        logs.log_info('Results exported to: %s', phase.save_dirname)
