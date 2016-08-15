|build-status|

.. |build-status| image:: https://gitlab.com/betse/betse/badges/master/build.svg
   :target: https://gitlab.com/betse/betse/commits/master

=====
BETSE
=====

**BETSE** (**B**\ io\ **E**\ lectric **T**\ issue **S**\ imulation **E**\ ngine)
is an open-source cross-platform `finite volume`_ simulator for 2D computational
multiphysics problems in the life sciences – including electrodiffusion_,
electro-osmosis_, galvanotaxis_, `voltage-gated ion channels`_, `gene regulatory
networks`_, and `biochemical reaction networks`_ (e.g., metabolism). BETSE is
associated with the `Paul Allen Discovery Center`_ at `Tufts University`_ and
supported by a `Paul Allen Discovery Center award`_ from the `Paul G. Allen
Frontiers Group`_.

BETSE is `portably implemented <codebase_>`__ in pure `Python 3`_, `continuously
stress-tested <testing_>`__ with GitLab-CI_ **+** py.test_, and `permissively
distributed <License_>`__ under the `BSD 2-clause license`_. While a high-level
graphical user interface (GUI) supporting all popular platforms is planned,
BETSE currently *only* provides a low-level command line interface (CLI).

.. # Table of contents, excluding the document heading above.
.. contents:: :local:

Installation
============

BETSE is installable under **Linux**, **OS X**, and **Windows** as
follows:

#. Install Git_.
#. Install the **Python 3.x** (e.g., 3.5) variant of Anaconda_. :sup:`Do not
   install the Python 2.7 variant of Anaconda_. BETSE requires Python 3.x.`
#. Run the following commands from within a command-line terminal:

   #. Download the live version of BETSE.

      .. code:: bash
   
         git clone https://gitlab.com/betse/betse.git

   #. Install BETSE.

      .. code:: bash
   
         cd betse && sudo python3 setup.py install

   #. *(Optional)* Test BETSE by running a sample simulation.

      .. code:: bash
   
         cd /tmp && betse try

For a comprehensive list of system requirements and software dependencies by
platform, see our `complete installation instructions <doc/md/INSTALL.md>`__.

Introduction
============

BETSE simulates biorealistic electrochemical phenomena in `gap junction`_\
-networked 2D cellular collectives. To predict `bioelectric patterns
<bioelectricity_>`__ and their spatio-temporal dynamics, BETSE:

-  Models `ion channel`_ and `gap junction`_ activity.
-  Tracks changes in ion concentration and net ionic charge.
-  Calculates endogenous voltages and currents.
-  Exports simulation results to a variety of output formats, including:

   -  Publication-quality plots and graphs.
   -  Internet-friendly compressed video.
   -  Post-processable tabular data (e.g., `comma-separated values
      (CSV) <https://en.wikipedia.org/wiki/Comma-separated_values>`__).

-  Imports bitmask images defining the shapes of:

   -  Cell clusters.
   -  Cell cluster regions localizing `ion channel`_ activity, typically
      signifying disparate types of adjacent tissue.

To assemble simple concepts into complex simulations, BETSE provides a rich
biological toolset consisting of:

Ions
----

Simulations may enable arbitrary combinations of the principal ions implicated
in bioelectrical signaling – including:

- Calcium_ (*Ca*\ :sup:`2+`).
- Chloride_ (*Cl*\ :sup:`-`).
- Hydron_ (*H*\ :sup:`+`).
- Potassium_ (*K*\ :sup:`+`).
- Sodium_ (*Na*\ :sup:`+`).
- `Anionic protein`_ (*P*\ :sup:`-`).
- `Bicarbonate transporter protein`_ (*HCO*\ :sup:`-`\ :sub:`3`).

Ion Channels
------------

Individual cells in simulations may enable arbitrary combinations of
`voltage-gated ion channels`_, each implementing the `Hodgkin-Huxley (HH)
formalism`_ with experimentally-derived parameters sourced from reputable
`knowledge-based systems`_ (e.g., Channelpedia_). Explicitly supported channel
types include:

- HCN1_, HCN2_, and HCN4_.
- `L-type Ca`_, `T-type Ca`_, and |P/Q-type Ca|_.
- Kir2.1_.
- Kv1.1_, Kv1.2_, Kv1.5_. Kv3.3_, and Kv3.4_.
- Nav1.2_, Nav1.3_, and Nav1.6_.
- Built-in and custom `leak <leak channels_>`__ and `ligand-gated channels`_,
  including:

  - |Calcium-gated K+ channels|_.

Ion Pumps and Exchangers
------------------------

For fine-grained control over cell dynamics, notable ion pumps and exchangers
may also be selectively enabled – including:

- |Ca2+-ATPase|_.
- |H+/K+-ATPase|_.
- |Na+/K+-ATPase|_.
- V-ATPase_.

Extracellular Space
-------------------

Cells form interconnected intracellular networks via `voltage-sensitive
gap junction
connections <(https://en.wikipedia.org/wiki/Gap_junction)>`__ embedded
within an `extracellular
space <https://en.wikipedia.org/wiki/Extracellular>`__ maintained by
`tight junctions <https://en.wikipedia.org/wiki/Tight_junction>`__ at
the cell cluster periphery. Simulation of the extracellular environment
enables exploration of `local field
potentials <https://en.wikipedia.org/wiki/Local_field_potential>`__, the
`transepithelial
potential <https://en.wikipedia.org/wiki/Transepithelial_potential_difference>`__,
and `ephaptic
coupling <https://en.wikipedia.org/wiki/Ephaptic_coupling>`__ between
cells.

Networks
--------

Built-in and custom `gene regulatory
networks <(https://en.wikipedia.org/wiki/Gene_regulatory_network)>`__
and `biochemical reaction
networks <(http://www.nature.com/subjects/biochemical-reaction-networks)>`__
(emphasizing metabolism) are fully supported. To unite these powerful
control systems with bioelectrical signaling, the activity-modulated
interaction between `gene
products <https://en.wikipedia.org/wiki/Gene_product>`__ and other
biochemicals is fully integrated with `ion
channels <https://en.wikipedia.org/wiki/Ion_channel>`__,
`pumps <https://en.wikipedia.org/wiki/Active_transport>`__, and `gap
junctions <https://en.wikipedia.org/wiki/Gap_junction>`__.

Validation
==========

BETSE is peer-reviewed software continuing to receive evidence-based
scrutiny. Simulation output is reproducibly synchronized with
experimental observations on `membrane
permeability <https://en.wikipedia.org/wiki/Cell_membrane#Permeability>`__,
`resting potential <https://en.wikipedia.org/wiki/Resting_potential>`__,
ion concentration, and adjunct biophysical quantities. Predictable
outcomes have been demonstrated for well-known cases, including:

-  Prediction of the correct `transmembrane
   voltage <https://en.wikipedia.org/wiki/Membrane_potential>`__ changes
   on perturbations to single cell membrane states and environmental ion
   concentrations.
-  Development of realistic `transepithelial potential differences
   (TEPD) <https://en.wikipedia.org/wiki/Transepithelial_potential_difference>`__.
-  Development of realistic bioelectric signals on large-scale cellular
   wounds.

For details, see our recently published `introductory
paper <#reference>`__.

License
=======

BETSE is open-source software `licensed <LICENSE>`__ under the
permissive `BSD 2-clause
license <https://opensource.org/licenses/BSD-2-Clause>`__.

Reference
=========

If leveraging BETSE in your own work, consider citing our `introductory
paper`_:

    `Pietak, Alexis`_ and `Levin, Michael`_ (\ *2016*\ ). |article name|_
    |journal name|_ 4, 55. ``doi:10.3389/fbioe.2016.00055``

Authors
=======

BETSE comes courtesy a dedicated community of authors_ and contributors_ –
without whom this project would be computationally impoverished, intellectually
bankrupt, and unusable beyond reason. **Thanks, all.**

See Also
========

For prospective users:

-  `Installation <doc/md/INSTALL.md>`__, detailing BETSE's
   installation with exhaustive platform-specific instructions.
-  `Usage <doc/md/USAGE.md>`__, detailing BETSE's command-line
   interface (CLI) with human-readable explanation and examples.

For prospective contributors:

-  `Development <doc/md/DEVELOP.md>`__, detailing development of the
   BETSE codebase – philosophy, workflow, and otherwise.
-  `Testing <doc/md/TEST.md>`__, detailing testing of the BETSE
   codebase – `continuous integration
   (CI) <https://en.wikipedia.org/wiki/Continuous_integration>`__,
   manual testing, and otherwise.
-  `Freezing <doc/md/FREEZE.md>`__, detailing conversion of the
   BETSE codebase into redistributable platform-specific executable
   binaries.

.. # ------------------( LINKS ~ academia                   )------------------
.. _Pietak, Alexis:
   https://www.researchgate.net/profile/Alexis_Pietak
.. _Levin, Michael:
   https://ase.tufts.edu/biology/labs/levin
.. _Channelpedia:
   http://channelpedia.epfl.ch
.. _Paul Allen Discovery Center:
   http://www.alleninstitute.org/what-we-do/frontiers-group/discovery-centers/allen-discovery-center-tufts-university
.. _Paul Allen Discovery Center award:
   https://www.alleninstitute.org/what-we-do/frontiers-group/news-press/press-resources/press-releases/paul-g-allen-frontiers-group-announces-allen-discovery-center-tufts-university
.. _Paul G. Allen Frontiers Group:
   https://www.alleninstitute.org/what-we-do/frontiers-group
.. _Tufts University:
   https://www.tufts.edu

.. # ------------------( LINKS ~ citation                   )------------------
.. _introductory paper:
   http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract

.. |article name| replace::
   **Exploring Instructive Physiological Signaling with the Bioelectric Tissue
   Simulation Engine (BETSE).**
.. _article name:
   http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract

.. |journal name| replace::
   *Frontiers in Bioengineering and Biotechnology.*
.. _journal name:
   http://journal.frontiersin.org/journal/bioengineering-and-biotechnology

.. # ------------------( LINKS ~ project                    )------------------
.. _authors:
   AUTHORS.md
.. _contributors:
   https://gitlab.com/betse/betse/graphs/master
.. _testing:
   https://gitlab.com/betse/betse/pipelines
.. _codebase:
   https://gitlab.com/betse/betse/tree/master

.. # ------------------( LINKS ~ science                    )------------------
.. _bioelectricity:
   https://en.wikipedia.org/wiki/Bioelectromagnetics
.. _biochemical reaction networks:
   http://www.nature.com/subjects/biochemical-reaction-networks
.. _electrodiffusion:
   https://en.wikipedia.org/wiki/Nernst%E2%80%93Planck_equation
.. _electro-osmosis:
   https://en.wikipedia.org/wiki/Electro-osmosis
.. _finite volume:
   https://en.wikipedia.org/wiki/Finite_volume_method
.. _galvanotaxis:
   https://en.wiktionary.org/wiki/galvanotaxis
.. _gap junction:
   https://en.wikipedia.org/wiki/Gap_junction
.. _gene regulatory networks:
   https://en.wikipedia.org/wiki/Gene_regulatory_network
.. _Hodgkin-Huxley (HH) formalism:
   https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model
.. _knowledge-based systems:
   https://en.wikipedia.org/wiki/Knowledge-based_systems

.. # ------------------( LINKS ~ science : ions             )------------------
.. _calcium:   https://en.wikipedia.org/wiki/Calcium_in_biology
.. _chloride:  https://en.wikipedia.org/wiki/Chloride
.. _hydron:    https://en.wikipedia.org/wiki/Hydron_(chemistry)
.. _sodium:    https://en.wikipedia.org/wiki/Sodium_in_biology
.. _potassium: https://en.wikipedia.org/wiki/Potassium_in_biology
.. _anionic protein:
   https://en.wikipedia.org/wiki/Gibbs%E2%80%93Donnan_effect
.. _bicarbonate transporter protein:
   https://en.wikipedia.org/wiki/Bicarbonate_transporter_protein

.. # ------------------( LINKS ~ science : channels         )------------------
.. _ion channel:
   https://en.wikipedia.org/wiki/Ion_channel
.. _leak channels:
   https://en.wikipedia.org/wiki/Leak_channel
.. _ligand-gated channels:
   https://en.wikipedia.org/wiki/Ligand-gated_ion_channel
.. _voltage-gated ion channels:
   https://en.wikipedia.org/wiki/Voltage-gated_ion_channel

.. |calcium-gated K+ channels| replace::
   Calcium-gated K\ :sup:`+` channels
.. _calcium-gated K+ channels:
   https://en.wikipedia.org/wiki/Calcium-activated_potassium_channel

.. # ------------------( LINKS ~ science : channels : type  )------------------
.. _HCN1:   http://channelpedia.epfl.ch/ionchannels/61
.. _HCN2:   http://channelpedia.epfl.ch/ionchannels/62
.. _HCN4:   http://channelpedia.epfl.ch/ionchannels/64
.. _Kir2.1: http://channelpedia.epfl.ch/ionchannels/42
.. _Kv1.1:  http://channelpedia.epfl.ch/ionchannels/1
.. _Kv1.2:  http://channelpedia.epfl.ch/ionchannels/2
.. _Kv1.5:  http://channelpedia.epfl.ch/ionchannels/5
.. _Kv3.3:  http://channelpedia.epfl.ch/ionchannels/13
.. _Kv3.4:  http://channelpedia.epfl.ch/ionchannels/14
.. _Nav1.2: http://channelpedia.epfl.ch/ionchannels/121
.. _Nav1.3: http://channelpedia.epfl.ch/ionchannels/122
.. _Nav1.6: http://channelpedia.epfl.ch/ionchannels/125
.. _L-type Ca:   http://channelpedia.epfl.ch/ionchannels/212
.. _T-type Ca:   https://en.wikipedia.org/wiki/T-type_calcium_channel

.. |P/Q-type Ca| replace:: :sup:`P`\ /\ :sub:`Q`-type Ca
.. _P/Q-type Ca:
   http://channelpedia.epfl.ch/ionchannels/78

.. # ------------------( LINKS ~ science : pumps : type     )------------------
.. _V-ATPase: https://en.wikipedia.org/wiki/V-ATPase

.. |Ca2+-ATPase| replace:: Ca\ :sup:`2+`-ATPase
.. _Ca2+-ATPase: https://en.wikipedia.org/wiki/Calcium_ATPase

.. |H+/K+-ATPase| replace:: H\ :sup:`+`/K\ :sup:`+`-ATPase
.. _H+/K+-ATPase: https://en.wikipedia.org/wiki/Hydrogen_potassium_ATPase

.. |Na+/K+-ATPase| replace:: Na\ :sup:`+`/K\ :sup:`+`-ATPase
.. _Na+/K+-ATPase: https://en.wikipedia.org/wiki/Na%2B/K%2B-ATPase

.. # ------------------( LINKS ~ software                   )------------------
.. _Anaconda:
   https://www.continuum.io/downloads
.. _BSD 2-clause license:
   https://opensource.org/licenses/BSD-2-Clause
.. _Git:
   https://git-scm.com/downloads
.. _GitLab-CI:
   https://about.gitlab.com/gitlab-ci
.. _Python 3:
   https://www.python.org
.. _py.test:
   http://pytest.org
