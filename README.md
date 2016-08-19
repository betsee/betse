[![build status](https://gitlab.com/betse/betse/badges/master/build.svg)](https://gitlab.com/betse/betse/commits/master)

BETSE
===========

**BETSE** (**B**io **E**lectric **T**issue **S**imulation **E**ngine) is an
open-source cross-platform [finite
volume](https://en.wikipedia.org/wiki/Finite_volume_method) simulator for 2D
computational multiphysics problems in the life sciences, including
[electrodiffusion](https://en.wikipedia.org/wiki/Nernst%E2%80%93Planck_equation),
[electro-osmosis](https://en.wikipedia.org/wiki/Electro-osmosis),
[galvanotaxis](https://en.wiktionary.org/wiki/galvanotaxis), [voltage-gated ion
channels](https://en.wikipedia.org/wiki/Voltage-gated_ion_channel), [gene
regulatory networks](https://en.wikipedia.org/wiki/Gene_regulatory_network),
and [biochemical reaction
networks](http://www.nature.com/subjects/biochemical-reaction-networks) (e.g.,
metabolism). BETSE is associated with the [Paul Allen Discovery Center at Tufts
University](http://www.alleninstitute.org/what-we-do/frontiers-group/discovery-centers/allen-discovery-center-tufts-university/)
and is supported by a Paul Allen Discovery Center award from the [Paul G. Allen
Frontiers Group](https://www.alleninstitute.org/what-we-do/frontiers-group).

BETSE is [portably implemented](betse) in pure [Python
3](https://en.wikipedia.org/wiki/History_of_Python), [continuously
stress-tested](#testing) with [GitLab-CI](https://about.gitlab.com/gitlab-ci)
**+** [py.test](http://pytest.org), and [permissively distributed](#license)
under the [BSD 2-clause license](https://opensource.org/licenses/BSD-2-Clause).
While a high-level graphical user interface (GUI) supporting all popular
platforms is planned, BETSE currently _only_ provides a low-level command
line interface (CLI).

## Installation

BETSE is installable under **Linux**, **OS X**, and **Windows** as follows:

1. Install **[Git](https://git-scm.com/downloads).**
1. Install the **Python 3.x** (e.g., 3.5) variant of
   **[Anaconda](https://www.continuum.io/downloads).** <sup>Do _not_ install the
   Python 2.7 variant of Anaconda. BETSE requires Python 3.x.</sup>
1. Run the following commands from within a command-line terminal:
   1. Download the live version of BETSE.

            git clone https://gitlab.com/betse/betse.git

   1. Install BETSE.

            cd betse && sudo python3 setup.py install

   1. **_(Optional)_** Test BETSE by running a sample simulation.

            cd /tmp && betse try

For a comprehensive list of all system requirements, software dependencies, and
platform-specific installation instructions, see [**BETSE
Installation**](doc/md/INSTALL.md).

## Introduction

BETSE simulates biorealistic electrochemical phenomena in gap
junction-networked 2D cellular collectives. To predict bioelectric patterns
and their spatio-temporal dynamics, BETSE:

* Models [ion channel](https://en.wikipedia.org/wiki/Ion_channel) and
  [gap junction](https://en.wikipedia.org/wiki/Gap_junction) activity.
* Tracks changes in ion concentration and net ionic charge.
* Calculates endogenous voltages and currents.
* Exports simulation results to a variety of output formats, including:
  * Publication-quality plots and graphs.
  * Internet-friendly compressed video.
  * Post-processable tabular data (e.g., [comma-separated values
    (CSV)](https://en.wikipedia.org/wiki/Comma-separated_values)).
* Imports bitmask images defining the shapes of:
  * Cell clusters.
  * Cell cluster regions localizing ion channel activity, typically signifying
    different types of adjacent tissues (e.g., organs).

BETSE offers a rich biological toolset for aggregating simple concepts into
complex simulations.

### Ions

Simulations may enable arbitrary combinations of the seven principal ions
underpinning bioelectrical signals:
[Na<sup>+</sup>](https://en.wikipedia.org/wiki/Sodium_in_biology),
[K<sup>+</sup>](https://en.wikipedia.org/wiki/Potassium_in_biology),
[Cl<sup>-</sup>](https://en.wikipedia.org/wiki/Chloride),
[Ca<sup>2+</sup>](https://en.wikipedia.org/wiki/Calcium_in_biology),
[H<sup>+</sup>](https://en.wikipedia.org/wiki/Hydron_(chemistry)),
[HCO<sub>3</sub><sup>-</sup>](https://en.wikipedia.org/wiki/Bicarbonate_transporter_protein),
and [anionic proteins
(P<sup>-</sup>)](https://en.wikipedia.org/wiki/Gibbs%E2%80%93Donnan_effect).

### Ion Channels

Individual cells in simulations may enable arbitrary combinations of
[voltage-gated ion
channels](https://en.wikipedia.org/wiki/Voltage-gated_ion_channel) implementing
the [Hodgkin-Huxley (HH)
formalism](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) with
experimentally-derived variables sourced from
[Channelpedia](http://channelpedia.epfl.ch). Currently implemented channel
types include [Nav1.2](http://channelpedia.epfl.ch/ionchannels/121),
[Nav1.3](http://channelpedia.epfl.ch/ionchannels/122),
[Nav1.6](http://channelpedia.epfl.ch/ionchannels/125),
[Kv1.1](http://channelpedia.epfl.ch/ionchannels/1),
[Kv1.2](http://channelpedia.epfl.ch/ionchannels/2),
[Kv1.5](http://channelpedia.epfl.ch/ionchannels/5),
[Kv3.3](http://channelpedia.epfl.ch/ionchannels/13),
[Kv3.4](http://channelpedia.epfl.ch/ionchannels/14),
[Kir2.1](http://channelpedia.epfl.ch/ionchannels/42), [L-type
Ca](http://channelpedia.epfl.ch/ionchannels/212), [T-type
Ca](https://en.wikipedia.org/wiki/T-type_calcium_channel), [P/Q-type
Ca](http://channelpedia.epfl.ch/ionchannels/78),
[HCN1](http://channelpedia.epfl.ch/ionchannels/61),
[HCN2](http://channelpedia.epfl.ch/ionchannels/62), and
[HCN4](http://channelpedia.epfl.ch/ionchannels/64). Built-in and custom
[leak](https://en.wikipedia.org/wiki/Leak_channel) and [ligand-gated
channels](https://en.wikipedia.org/wiki/Ligand-gated_ion_channel) are
available, including [calcium-gated K<sup>+</sup>
channels](https://en.wikipedia.org/wiki/Voltage-dependent_calcium_channel).

For added control over bioelectric cell dynamics, notable ion pumps and
exchangers such as
[Na<sup>+</sup>/K<sup>+</sup>-ATPase](https://en.wikipedia.org/wiki/Na%2B/K%2B-ATPase),
[H<sup>+</sup>/K<sup>+</sup>-ATPase](https://en.wikipedia.org/wiki/Hydrogen_potassium_ATPase),
[V-ATPase](https://en.wikipedia.org/wiki/V-ATPase), and
[Ca<sup>2+</sup>-ATPase](https://en.wikipedia.org/wiki/Calcium_ATPase) are
also available.

### Extracellular Space

Cells form interconnected intracellular networks via [voltage-sensitive gap
junction connections]((https://en.wikipedia.org/wiki/Gap_junction)) embedded
within an [extracellular space](https://en.wikipedia.org/wiki/Extracellular)
maintained by [tight junctions](https://en.wikipedia.org/wiki/Tight_junction)
at the cell cluster periphery. Simulation of the extracellular environment
enables exploration of [local field
potentials](https://en.wikipedia.org/wiki/Local_field_potential), the
[transepithelial
potential](https://en.wikipedia.org/wiki/Transepithelial_potential_difference),
and [ephaptic coupling](https://en.wikipedia.org/wiki/Ephaptic_coupling)
between cells.

### Networks

Built-in and custom [gene regulatory
networks]((https://en.wikipedia.org/wiki/Gene_regulatory_network)) and
[biochemical reaction
networks]((http://www.nature.com/subjects/biochemical-reaction-networks))
(emphasizing metabolism) are fully supported. To unite these powerful control
systems with bioelectrical signaling, the activity-modulated interaction
between [gene products](https://en.wikipedia.org/wiki/Gene_product) and other
biochemicals is fully integrated with [ion
channels](https://en.wikipedia.org/wiki/Ion_channel),
[pumps](https://en.wikipedia.org/wiki/Active_transport), and [gap
junctions](https://en.wikipedia.org/wiki/Gap_junction).

## Validation

BETSE is peer-reviewed software continuing to receive evidence-based scrutiny.
Simulation output is reproducibly synchronized with experimental observations
on [membrane
permeability](https://en.wikipedia.org/wiki/Cell_membrane#Permeability),
[resting potential](https://en.wikipedia.org/wiki/Resting_potential), ion
concentration, and adjunct biophysical quantities. Predictable outcomes have
been demonstrated for well-known cases, including:

* Prediction of the correct [transmembrane
  voltage](https://en.wikipedia.org/wiki/Membrane_potential) changes on
  perturbations to single cell membrane states and environmental ion
  concentrations.
* Development of realistic
  [transepithelial potential differences
  (TEPD)](https://en.wikipedia.org/wiki/Transepithelial_potential_difference).
* Development of realistic bioelectric signals on large-scale cellular wounds.

For details, see our recently published [introductory paper](#reference).

## License

BETSE is open-source software [licensed](LICENSE) under the permissive [BSD
2-clause license](https://opensource.org/licenses/BSD-2-Clause).

## Reference

When leveraging BETSE in your own publications, consider citing our
[introductory
paper](http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract):

> [Pietak, Alexis](https://www.researchgate.net/profile/Alexis_Pietak) and [Levin, Michael](https://ase.tufts.edu/biology/labs/levin)
> (2016). [**Exploring Instructive Physiological
> Signaling with the Bioelectric Tissue Simulation Engine
> (BETSE)**](http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract).
> [_Frontiers in Bioengineering and
> Biotechnology_](http://journal.frontiersin.org/journal/bioengineering-and-biotechnology)
> 4, 55. `doi:10.3389/fbioe.2016.00055`

## Authors

BETSE comes courtesy the contributions of a cadre of [authors](AUTHORS.md) –
without whom this engine would be computationally impoverished, intellectually
diminished, and decrepit beyond all unusable compare. **Thanks, all.**

## See Also

For prospective users:

* [**Installation**](doc/md/INSTALL.md), detailing BETSE's installation with
  exhaustive platform-specific instructions.
* [**Usage**](doc/md/USAGE.md), detailing BETSE's command-line interface (CLI)
  with human-readable explanation and examples.

For prospective contributors:

* [**Development**](doc/md/DEVELOP.md), detailing development of the BETSE
  codebase – philosophy, workflow, and otherwise.
* [**Testing**](doc/md/TEST.md), detailing testing of the BETSE codebase –
  [continuous integration
  (CI)](https://en.wikipedia.org/wiki/Continuous_integration), manual testing,
  and otherwise.
* [**Freezing**](doc/md/FREEZE.md), detailing conversion of the BETSE
  codebase into redistributable platform-specific executable binaries.
