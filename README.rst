.. # ------------------( BADGES                             )------------------
.. image::  https://gitlab.com/betse/betse/badges/master/build.svg
   :target: https://gitlab.com/betse/betse/pipelines
   :alt: Linux Build Status
.. image::  https://ci.appveyor.com/api/projects/status/mow7y8k3vpfu30c6/branch/master?svg=true
   :target: https://ci.appveyor.com/project/betse/betse/branch/master
   :alt: Windows Build Status

.. # ------------------( SYNOPSIS                           )------------------

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
stress-tested <testing_>`__ with GitLab-CI_ **×** Appveyor_ **+** py.test_, and
`permissively distributed <License_>`__ under the `BSD 2-clause license`_.

======
BETSEE
======

BETSEE_ (\ **BETSE E**\ nvironment) is the official open-source cross-platform
graphical user interface (GUI) for BETSE. BETSEE_ wraps the low-level
command-line interface (CLI) bundled with BETSE in a high-level interactive
modelling environment optimized for both new and advanced users alike.

Like BETSE, BETSEE_ is `portably implemented <BETSEE codebase_>`__ in `Python
3`_ and `permissively distributed <License_>`__ under the `BSD 2-clause
license`_. Unlike BETSE, BETSEE_ leverages the industry-standard PySide2_-based
`Qt 5 <Qt_>`_ application framework to deliver a modern scientific workflow.

.. # ------------------( TABLE OF CONTENTS                  )------------------
.. # Blank line. By default, Docutils appears to only separate the subsequent
.. # table of contents heading from the prior paragraph by less than a single
.. # blank line, hampering this table's readability and aesthetic comeliness.

|

.. # Table of contents, excluding the above document heading. While the
.. # official reStructuredText documentation suggests that a language-specific
.. # heading will automatically prepend this table, this does *NOT* appear to
.. # be the case. Instead, this heading must be explicitly declared.

.. contents:: **Contents**
   :local:

.. # ------------------( DESCRIPTION                        )------------------

Installation
============

BETSE currently supports **Linux**, **macOS**, and **Windows** out-of-the-box.

Simple
--------

For new users, BETSE is portably installable as follows:

- [\ *Windows*\ ] Emulate **Ubuntu Linux** via the `Windows Subsystem for Linux
  (WSL) <WSL_>`__. [#windows_not]_
- Install the **Python 3.x** [#python2_not]_ (e.g., 3.6) variant of Anaconda_.
  [#anaconda_not]_
- Open a **Bash terminal.** [#terminal]_
- Run the following commands.
  
  - Enable conda-forge_.
  
    .. code:: bash
  
       conda config --add channels conda-forge
  
  - Install **BETSE.** [#conda_package]_
  
    .. code:: bash
  
       conda install betse
  
  - [\ *Optional*\ ] Test **BETSE.** [#test_command]_

    .. code:: bash

       betse -v try

.. [#windows_not]
   The `Windows Subsystem for Linux (WSL) <WSL_>`__ and – hence BETSE itself –
   is *only* installable under **Windows 10.** Under older Windows versions,
   BETSE may be installed from a `virtual Linux guest <VirtualBox_>`__.

.. [#python2_not]
   Do *not* install the **Python 2.7** variant of Anaconda_. BETSE requires
   **Python 3.x.**

.. [#anaconda_not]
   If you prefer *not* to install Anaconda_, BETSE dependencies are also
   `manually installable <Advanced_>`__ via your platform-specific package
   manager (e.g., Homebrew_ on macOS, APT_ on Ubuntu Linux). Doing so is
   non-trivial and, where performed incorrectly, could produce a
   performance-crippled single-core installation of BETSE. Anaconda_ suffers no
   such issues and is guaranteed to produce a performance-optimized multicore
   installation of BETSE on *all* supported platforms. We strongly recommend
   Anaconda_ – even when you think you know better.

.. [#terminal]
   To open a `POSIX`_\ -compatible terminal under:

   - **Windows:**

     #. Install **Ubuntu Linux** via the `Windows Subsystem for Linux (WSL) <WSL_>`__.
     #. Open an *Ubuntu Linux terminal.*

   - **macOS:**

     #. Open the *Finder*.
     #. Open the *Applications* folder.
     #. Open the *Utilities* folder.
     #. Open *Terminal.app*.

   - **Ubuntu Linux:**

     #. Type ``Ctrl``\ +\ ``Alt``\ +\ ``t``.

.. [#conda_package]

   This command installs both the `most recent stable release of BETSE <conda
   package_>`__ *and* all mandatory and most optional dependencies of this
   release. Older stable releases are installable in a similar manner (e.g.,
   ``conda install betse=0.7.0`` for BETSE 0.7.0). All `Anaconda packages`_ are
   kindly hosted by the `non-profit conda-forge organization <conda-forge_>`__.

.. [#test_command]

   This test command:

   #. Enables verbosity with the ``-v`` option, simplifying issue reporting in
      the event of an unexpected error.
   #. Creates a ``sample_sim/`` subdirectory in the current directory, providing
      the default simulation for this release of BETSE. This includes *all*
      configuration files and resources referenced by these files.
   #. Runs all simulation phases (e.g., seed, initialization) of this simulation.

   When finished, you may safely either:
   
   - Remove this subdirectory.
   - Rename this subdirectory (e.g., to ``my_sim/``) to serve as a basis for
     subsequent simulations.
   - Preserve this subdirectory as is.

Advanced
--------

For developers and advanced users, *any* version of BETSE – including the live
repository and prior stable releases – is manually installable as follows:

#. Install **Python 3.x** and `all dependencies <dependencies_>`__ required by
   BETSE. Under:

   - **Linux,** install `these dependencies <dependencies_>`__ via your
     distribution-specific package manager (e.g., APT_ under Debian-based
     distributions). Do *not* use ``pip``.\ [#pip_not]_
   - **macOS,** either:

     - (\ *Recommended*\ ) Install the **Python 3.x** variant of Anaconda_.
     - Or both:

       #. Install a third-party package manager (e.g., Homebrew_, MacPorts_).
          Apple does *not* provide a package manager out-of-the-box.
       #. Install `these dependencies <dependencies_>`__ via that package
          manager. Do *not* use ``pip``.\ [#pip_not]_

   - **Windows,** install the **Python 3.x** variant of Anaconda_.\ [#windows]_

#. Open a **terminal.**
#. **Download** either:

   - **The unstable BETSE repository** as follows:

     #. Install Git_.
     #. Clone the ``master`` branch of this repository.

        .. code:: bash

           git clone https://gitlab.com/betse/betse.git

     #. Prepare for installation.

        .. code:: bash

           cd betse

   - **Any stable BETSE release,** including the most recent, as follows:

     #. Visit our `source tarball archive <tarballs_>`__.
     #. Click the download icon to the right of the desired release and select
        *Download tar.gz*.
     #. Extract the downloaded tarball into the current directory.

        .. code:: bash

           tar -xvzf betse-*.tar.gz

     #. (\ *Optional*\ ) Remove this tarball.

        .. code:: bash

           rm betse-*.tar.gz

     #. Prepare for installation.

        .. code:: bash

           cd betse-*

#. **Install BETSE** either:

   - (\ *Recommended*\ ) **Editably,** installing a cross-platform symbolic link
     to the current BETSE codebase. Modifications to this code are applied
     immediately *without* requiring reinstallation.

     .. code:: bash

        sudo python3 setup.py develop

   - **Non-editably,** installing a physical copy of the current BETSE codebase.
     Modifications to this code are ignored and thus require reinstallation.

     .. code:: bash

        sudo python3 setup.py install

#. (\ *Optional*\ ) **Test BETSE,** running all modelling phases of a sample
   simulation from a new directory.

   .. code:: bash

      cd /tmp && betse try


.. [#pip_not]
   Do *not* install scientific dependencies (e.g., NumPy_, SciPy_) with either
   ``pip`` or ``easy_install``; doing so typically degrades BETSE to single-core
   performance. To optimize BETSE across multiple cores, *always* install these
   dependencies with your platform-specific package manager (e.g., Homebrew_,
   APT_).

.. [#windows]
   Unlike Linux and macOS, Anaconda_ is (\ *effectively*\ ) required under
   Windows. Due to Microsoft's lack of support for `POSIX`_\ -compliant
   toolchains, *no* reasonable alternatives for installing multicore-aware
   scientific dependencies exist.

Usage
============

BETSE installs the ``betse`` command, a low-level command line interface (CLI)
designed for experienced users. The typical CLI-based workflow resembles:

#. Create a default simulation configuration in a new ``my_sim/`` subdirectory
   of the current directory.

   .. code:: bash

      betse config my_sim/sim_config.yml

#. [\ *Optional*\ ] Edit the newly created ``sim_config.yml`` file with your
   favourite text editor.
#. Create a pseudo-random cell cluster from this simulation configuration.

   .. code:: bash

      betse seed my_sim/sim_config.yml

#. [\ *Optional*\ ] Plot this cell cluster.

   .. code:: bash

      betse plot seed my_sim/sim_config.yml

#. Initialize (i.e., calculate steady-state concentrations for) this cell cluster.

   .. code:: bash

      betse init my_sim/sim_config.yml

#. [\ *Optional*\ ] Plot this initialized cell cluster.

   .. code:: bash

      betse plot init my_sim/sim_config.yml

#. Simulate this initialized cell cluster.

   .. code:: bash

      betse sim my_sim/sim_config.yml

#. [\ *Optional*\ ] Plot this simulated cell cluster.

   .. code:: bash

      betse plot sim my_sim/sim_config.yml

See the following external documents for detailed usage instructions – complete
with explanatory examples, sample plots, and ample screenshots:

- Official `BETSE 0.4 documentation`_. (\ *PDF format; 72 pages.*\ )
- Official `BETSE 0.3 documentation`_. (\ *PDF format; 77 pages.*\ )

Alternately, our sister project BETSEE_ installs the ``betsee`` command, a
high-level graphical user interface (GUI) designed for *all* users – regardless
of prior experience or familiarity with the command line.

Introduction
============

BETSE simulates biorealistic electrochemical phenomena in `gap junction`_\
-networked 2D cellular collectives. To predict `bioelectric patterns
<bioelectricity_>`__ and their spatio-temporal dynamics, BETSE:

- Models `ion channel`_ and `gap junction`_ activity.
- Tracks changes in ion concentration and net ionic charge.
- Calculates endogenous voltages and currents.
- Accepts simulation parameters, variables, and options as human-readable,
  well-commented configuration files in YAML_ format.
- Exports simulation results to a variety of output formats, including:

  - Publication-quality:

    - Plots, charts, and animations driven by Matplotlib_, the industry
      standard for open-source plot visualization.
    - `Directed graphs`_ (i.e., networks) driven by Graphviz_, the industry
      standard for open-source graph visualization.

  - Internet-friendly compressed video driven by any of various popular
    open-source video encoders, including FFmpeg_, Libav_, and MEncoder_.
  - Post-processable tabular data (e.g., `comma-separated values (CSV)
    <comma-separated values_>`__).

- Imports bitmask images defining the shapes of:

  - Cell clusters.
  - Cell cluster regions localizing `ion channel`_ activity, typically
    signifying disparate types of adjacent tissue.

To assemble simple concepts into complex simulations, BETSE supplies a richly
configurable, highly scalable biological toolset consisting of:

Ions
----

Simulations may enable arbitrary combinations of the principal ions implicated
in bioelectrical signaling – including:

- Sodium_ (*Na*\ :sup:`+`).
- Potassium_ (*K*\ :sup:`+`).
- Chloride_ (*Cl*\ :sup:`-`).
- Calcium_ (*Ca*\ :sup:`2+`).
- Hydrogen_ (*H*\ :sup:`+`).
- `Anionic proteins`_ (*P*\ :sup:`-`).
- Bicarbonate_ (*HCO*\ :sup:`-`\ :sub:`3`).

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
- `Leak <leak channels_>`__ and `ligand-gated channels`_, including:

  - |Calcium-gated K+ channels|_.

Custom ion channels parametrized by user-selected constants may be trivially
defined in the same manner (e.g., via a YAML_\ -formatted configuration file).

Ion Pumps and Exchangers
------------------------

For fine-grained control over cell dynamics, notable ion pumps and exchangers
may also be selectively enabled – including:

- |Ca2+-ATPase|_.
- |H+/K+-ATPase|_.
- |Na+/K+-ATPase|_.
- V-ATPase_.

Custom ion pumps and exchangers parametrized by user-selected constants may be
trivially defined in the same manner (e.g., via a YAML_\ -formatted
configuration file).

Extracellular Space
-------------------

Cells form interconnected intracellular networks via voltage-sensitive `gap
junction connections <gap junction_>`__ embedded within an `extracellular
environment`_, maintained by `tight junctions`_ at the cell cluster periphery.
Simulation of this environment enables exploration of `local field
potentials`_, `transepithelial potential`_, and `ephaptic coupling`_ between
cells.

Biological Networks
-------------------

Simulation of `gene regulatory <gene regulatory networks_>`__ and `biochemical
reaction networks`_ at both the cellular and mitochondrial level supports deep
spatial analysis of otherwise intractable biological processes. Metabolism,
disease, aging, and other `genetic <genetics_>`__ and `epigenetic
<epigenetics_>`__ phenomena commonly associated with quasi-`Big Data`_ are all
valid targets for exhaustive study with BETSE.

To integrate these potent control systems with bioelectrical signaling, the
`activity <enzyme activity_>`__-modulated interaction between `gene products`_
and similar biochemicals is fully integrated with `ion channels <ion
channel_>`__, `ion pumps`_, and `gap junctions`_.

Validation
==========

BETSE is peer-reviewed software receiving continual evidence-based scrutiny.
Simulation output is reproducibly synchronized with experimental observations on
`membrane permeability`_, `resting potential`_, ion concentration, and similar
real-world biophysical quantities. Predictable outcomes have been demonstrated
for such well-known cases as:

-  `Transmembrane voltage changes <transmembrane voltage_>`__ on perturbations
   to single cell membrane states and environmental ion concentrations.
-  `Transepithelial potential differences (TEPD) <transepithelial
   potential_>`__.
-  Bioelectrical signals at large-scale cellular wound sites.

Citation
========

BETSE is formally described in our `introductory paper <2016 article_>`__.
Third-party papers, theses, and other texts leveraging BETSE should ideally
cite the following:

    `Alexis Pietak`_ and `Michael Levin`_, 2016. |2016 article name|_
    |2016 article supplement|_ [#supplement]_ |2016 journal name|_ *4*\ (55).
    :sup:`https://doi.org/10.3389/fbioe.2016.00055`

Subsequent papers expanding the BETSE architecture with additional theory,
experimental results, and comparative metrics include:

    * `Alexis Pietak`_ and `Michael Levin`_, 2017. |2017 article name|_
      |2017 article supplement|_ [#supplement]_ |2017 journal name|_ *14*\ (134),
      p.20170425. :sup:`https://doi.org/10.1098/rsif.2017.0425`
    * Vaibhav P. Pai, `Alexis Pietak`_, Valerie Willocq, Bin Ye, Nian-Qing Shi,
      and `Michael Levin`_, 2018. |2018 hcn2 article name|_ |2018 hcn2 article
      supplement|_ [#supplement]_ |2018 hcn2 journal name|_ *9*\ (1), p.998.
      :sup:`https://doi.org/10.1038/s41467-018-03334-5`
    * Javier Cervera, `Alexis Pietak`_, `Michael Levin`_, and Salvador Mafe,
      2018. |2018 coupling article name|_ |2018 coupling journal name|_ *128*,
      pp.45—61. :sup:`https://doi.org/10.1016/j.bioelechem.2018.04.013`
    * `Alexis Pietak`_ and `Michael Levin`_, 2018. |2018 review article name|_
      |2018 review journal name|_
      :sup:`https://doi.org/10.1016/j.pbiomolbio.2018.03.008`

.. # FIXME: Add an image thumbnail for the first article above displaying the
.. # cover image selected by the prior journal for that edition's cover article.

.. # Note that, for unknown reasons, this footnote *MUST* be refenced above and
.. # defined here rather than in the supplement replacements defined below.

.. [#supplement]
   This article's supplement extends the cursory theory presented by this
   article with a rigorous treatment of the mathematics, formalisms, and
   abstractions required to fully reproduce this work. If theoretical questions
   remain after completing the main article, please consult this supplement.

Contact
=======

To contact `Dr. Pietak`_, the principal developer of the BETSE codebase and
first or second author of all above papers, cordially direct correspondence to
either:

* The `anonymous contact form <Organic Mechanics Contact_>`__ at `Organic
  Mechanics`_, the personal archives of all material published to date by
  `Dr. Pietak`_ – including papers, presentations, textbooks, and additional
  theoretical work.
* The personal e-mail account of `Dr. Pietak`_: [#e-mail]_

  * *Username:* **alexis** ``{dot}`` **pietak**
  * *Hostname:* **gmail** ``{dot}`` **com**

To report a software issue (e.g., bug, crash, or other unexpected behaviour)
*or* request a new feature in BETSE, consider `submitting a new issue <issue
submission_>`__ to our `issue tracker`_. Thanks in advance; it's only through
generous user contributions that your user experience can continue to improve.

.. [#e-mail]
   To protect Dr. Pietak's e-mail address against `automated harvesting <e-mail
   harvesting_>`__, this address has been intentionally obfuscated. Reconstruct
   the original address by:
   
   * Replacing the ``{dot}`` substring with the ``.`` charecter.
   * Concatenating the username and hostname strings with the ``@`` character.
     
   For this reason, consider directing correspondence to the `anonymous contact
   form <Organic Mechanics Contact_>`__ at `Organic Mechanics`_ instead.

Authors
=======

BETSE comes courtesy a dedicated community of `authors <author list_>`__ and
contributors_ – without whom this project would be computationally impoverished,
biologically misaligned, and simply unusable.

**Thanks, all.**

See Also
========

For prospective users:

-  `Installation <dependencies_>`__, detailing BETSE's installation with
   exhaustive platform-specific instructions.

For prospective contributors:

-  `Development <doc/md/DEVELOP.md>`__, detailing development of the BETSE
   codebase – philosophy, workflow, and otherwise.
-  `Testing <doc/md/TEST.md>`__, detailing testing of the BETSE codebase –
   `continuous integration`_, manual testing, and otherwise.
-  `Freezing <doc/md/FREEZE.md>`__, detailing conversion of the BETSE codebase
   into redistributable platform-specific executable binaries.

License
=======

BETSE is open-source software `released <license_>`__ under the permissive `BSD
2-clause license`_.

The logo prominently displayed on this `project page <project_>`__ is a flat
`Noun Project`_ icon entitled `"Cow," <Cows collection_>`__ `kindly released
<Noun Project license_>`__ under the permissive `BSD-compatible <license
compatibility_>`__ `CC BY 3.0 license`_ by `Maxim Kulikov`_.

.. # ------------------( LINKS ~ betse                      )------------------
.. _author list:
   doc/md/AUTHORS.md
.. _codebase:
   https://gitlab.com/betse/betse/tree/master
.. _conda package:
   https://anaconda.org/conda-forge/betse
.. _contributors:
   https://gitlab.com/betse/betse/graphs/master
.. _dependencies:
   doc/md/INSTALL.md
.. _issue submission:
   https://gitlab.com/betse/betse/issues/new?issue%5Bassignee_id%5D=&issue%5Bmilestone_id%5D=
.. _issue tracker:
   https://gitlab.com/betse/betse/issues
.. _license:
   LICENSE
.. _project:
   https://gitlab.com/betse/betse
.. _testing:
   https://gitlab.com/betse/betse/pipelines
.. _tarballs:
   https://gitlab.com/betse/betse/tags

.. # ------------------( LINKS ~ betse : docs               )------------------
.. _BETSE 0.4 documentation:
   https://www.dropbox.com/s/n8qfms2oks9cvv2/BETSE04_Documentation_Dec1st2016.pdf?dl=0
.. _BETSE 0.3 documentation:
   https://www.dropbox.com/s/fsxhjpipbiog0ru/BETSE_Documentation_Nov1st2015.pdf?dl=0

.. # ------------------( LINKS ~ betsee                     )------------------
.. _BETSEE:
   https://gitlab.com/betse/betsee
.. _BETSEE codebase:
   https://gitlab.com/betse/betsee/tree/master

.. # ------------------( LINKS ~ academia                   )------------------
.. _Michael Levin:
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

.. # ------------------( LINKS ~ academia : ally            )------------------
.. _Alexis Pietak:
.. _Pietak, Alexis:
.. _Dr. Pietak:
   https://www.researchgate.net/profile/Alexis_Pietak
.. _Organic Mechanics:
   https://www.omecha.org
.. _Organic Mechanics Contact:
   https://www.omecha.org/contact

.. # ------------------( LINKS ~ paper : 2016               )------------------
.. _2016 article:
   http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract

.. |2016 article name| replace::
   **Exploring instructive physiological signaling with the bioelectric tissue
   simulation engine (BETSE).**
.. _2016 article name:
   http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract

.. |2016 article supplement| replace::
   **(**\ Supplement\ **).**
.. _2016 article supplement:
   https://www.frontiersin.org/articles/file/downloadfile/203679_supplementary-materials_datasheets_1_pdf/octet-stream/Data%20Sheet%201.PDF/1/203679

.. |2016 journal name| replace::
   *Frontiers in Bioengineering and Biotechnology,*
.. _2016 journal name:
   http://journal.frontiersin.org/journal/bioengineering-and-biotechnology

.. # ------------------( LINKS ~ paper ~ 2017               )------------------
.. |2017 article name| replace::
   **Bioelectric gene and reaction networks: computational modelling of genetic, biochemical and bioelectrical dynamics in pattern regulation.**
.. _2017 article name:
   http://rsif.royalsocietypublishing.org/content/14/134/20170425

.. |2017 article supplement| replace::
   **(**\ Supplement\ **).**
.. _2017 article supplement:
   https://figshare.com/collections/Supplementary_material_from_Bioelectric_gene_and_reaction_networks_computational_modelling_of_genetic_biochemical_and_bioelectrical_dynamics_in_pattern_regulation_/3878404

.. |2017 journal name| replace::
   *Journal of The Royal Society Interface,*
.. _2017 journal name:
   http://rsif.royalsocietypublishing.org

.. # ------------------( LINKS ~ paper ~ 2018 : hcn2        )------------------
.. |2018 hcn2 article name| replace::
   **HCN2 rescues brain defects by enforcing endogenous voltage pre-patterns.**
.. _2018 hcn2 article name:
   https://www.nature.com/articles/s41467-018-03334-5

.. |2018 hcn2 article supplement| replace::
   **(**\ Supplement\ **).**
.. _2018 hcn2 article supplement:
   https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03334-5/MediaObjects/41467_2018_3334_MOESM1_ESM.pdf

.. |2018 hcn2 journal name| replace::
   *Nature Communications.*
.. _2018 hcn2 journal name:
   https://www.nature.com

.. # ------------------( LINKS ~ paper ~ 2018 : coupling   )------------------
.. |2018 coupling article name| replace::
   **Bioelectrical coupling in multicellular domains regulated by gap junctions: A conceptual approach.**
.. _2018 coupling article name:
   https://www.sciencedirect.com/science/article/pii/S156753941830063X?via%3Dihub

.. |2018 coupling journal name| replace::
   *Bioelectrochemistry.*
.. _2018 coupling journal name:
   https://www.sciencedirect.com/journal/bioelectrochemistry

.. # ------------------( LINKS ~ paper ~ 2018 : review     )------------------
.. |2018 review article name| replace::
   **Bioelectrical control of positional information in development and regeneration: A review of conceptual and computational advances.**
.. _2018 review article name:
   https://www.sciencedirect.com/science/article/pii/S0079610718300415

.. |2018 review journal name| replace::
   *Progress in Biophysics and Molecular Biology.*
.. _2018 review journal name:
   https://www.sciencedirect.com/journal/progress-in-biophysics-and-molecular-biology

.. # ------------------( LINKS ~ science                    )------------------
.. _bioelectricity:
   https://en.wikipedia.org/wiki/Bioelectromagnetics
.. _biochemical reaction networks:
   http://www.nature.com/subjects/biochemical-reaction-networks
.. _electrodiffusion:
   https://en.wikipedia.org/wiki/Nernst%E2%80%93Planck_equation
.. _electro-osmosis:
   https://en.wikipedia.org/wiki/Electro-osmosis
.. _enzyme activity:
   https://en.wikipedia.org/wiki/Enzyme_assay
.. _ephaptic coupling:
   https://en.wikipedia.org/wiki/Ephaptic_coupling
.. _epigenetics:
   https://en.wikipedia.org/wiki/Epigenetics
.. _extracellular environment:
   https://en.wikipedia.org/wiki/Extracellular
.. _finite volume:
   https://en.wikipedia.org/wiki/Finite_volume_method
.. _galvanotaxis:
   https://en.wiktionary.org/wiki/galvanotaxis
.. _gap junction:
.. _gap junctions:
   https://en.wikipedia.org/wiki/Gap_junction
.. _gene products:
   https://en.wikipedia.org/wiki/Gene_product
.. _gene regulatory networks:
   https://en.wikipedia.org/wiki/Gene_regulatory_network
.. _genetics:
   https://en.wikipedia.org/wiki/Genetics
.. _genetic algorithms:
   https://en.wikipedia.org/wiki/Genetic_algorithm
.. _Hodgkin-Huxley (HH) formalism:
   https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model
.. _local field potentials:
   https://en.wikipedia.org/wiki/Local_field_potential
.. _membrane permeability:
   https://en.wikipedia.org/wiki/Cell_membrane
.. _resting potential:
   https://en.wikipedia.org/wiki/Resting_potential
.. _tight junctions:
   https://en.wikipedia.org/wiki/Tight_junction
.. _transmembrane voltage:
   https://en.wikipedia.org/wiki/Membrane_potential
.. _transepithelial potential:
   https://en.wikipedia.org/wiki/Transepithelial_potential_difference

.. # ------------------( LINKS ~ science : ions             )------------------
.. _anionic proteins:
   https://en.wikipedia.org/wiki/Ion#anion
.. _bicarbonate: https://en.wikipedia.org/wiki/Bicarbonate
.. _calcium:     https://en.wikipedia.org/wiki/Calcium_in_biology
.. _chloride:    https://en.wikipedia.org/wiki/Chloride
.. _hydrogen:    https://en.wikipedia.org/wiki/Hydron_(chemistry)
.. _sodium:      https://en.wikipedia.org/wiki/Sodium_in_biology
.. _potassium:   https://en.wikipedia.org/wiki/Potassium_in_biology

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
.. _ion pumps:
   https://en.wikipedia.org/wiki/Active_transport

.. # ------------------( LINKS ~ science : pumps : type     )------------------
.. _V-ATPase: https://en.wikipedia.org/wiki/V-ATPase

.. |Ca2+-ATPase| replace:: Ca\ :sup:`2+`-ATPase
.. _Ca2+-ATPase: https://en.wikipedia.org/wiki/Calcium_ATPase

.. |H+/K+-ATPase| replace:: H\ :sup:`+`/K\ :sup:`+`-ATPase
.. _H+/K+-ATPase: https://en.wikipedia.org/wiki/Hydrogen_potassium_ATPase

.. |Na+/K+-ATPase| replace:: Na\ :sup:`+`/K\ :sup:`+`-ATPase
.. _Na+/K+-ATPase: https://en.wikipedia.org/wiki/Na%2B/K%2B-ATPase

.. # ------------------( LINKS ~ science : computer         )------------------
.. _Big Data:
   https://en.wikipedia.org/wiki/Big_data
.. _comma-separated values:
   https://en.wikipedia.org/wiki/Comma-separated_values
.. _continuous integration:
   https://en.wikipedia.org/wiki/Continuous_integration
.. _directed graphs:
   https://en.wikipedia.org/wiki/Directed_graph
.. _e-mail harvesting:
   https://en.wikipedia.org/wiki/Email_address_harvesting
.. _genenic algorithms:
   https://en.wikipedia.org/wiki/Genetic_algorithm
.. _knowledge-based systems:
   https://en.wikipedia.org/wiki/Knowledge-based_systems

.. # ------------------( LINKS ~ os : linux                 )------------------
.. _APT:
   https://en.wikipedia.org/wiki/Advanced_Packaging_Tool
.. _POSIX:
   https://en.wikipedia.org/wiki/POSIX
.. _Ubuntu:
.. _Ubuntu Linux:
   https://www.ubuntu.com
.. _Ubuntu Linux 16.04 (Xenial Xerus):
   http://releases.ubuntu.com/16.04

.. # ------------------( LINKS ~ os : macos                 )------------------
.. _Homebrew:
   http://brew.sh
.. _MacPorts:
   https://www.macports.org

.. # ------------------( LINKS ~ os : windows               )------------------
.. _WSL:
   https://msdn.microsoft.com/en-us/commandline/wsl/install-win10

.. # ------------------( LINKS ~ software                   )------------------
.. _Appveyor:
   https://ci.appveyor.com/project/betse/betse/branch/master
.. _dill:
   https://pypi.python.org/pypi/dill
.. _FFmpeg:
   https://ffmpeg.org
.. _Git:
   https://git-scm.com/downloads
.. _GitLab-CI:
   https://about.gitlab.com/gitlab-ci
.. _Graphviz:
   http://www.graphviz.org
.. _imageio:
   https://imageio.github.io
.. _Libav:
   https://libav.org
.. _Matplotlib:
   http://matplotlib.org
.. _NumPy:
   http://www.numpy.org
.. _MEncoder:
   https://en.wikipedia.org/wiki/MEncoder
.. _Python 3:
   https://www.python.org
.. _py.test:
   http://pytest.org
.. _SciPy:
   http://www.scipy.org
.. _VirtualBox:
   https://www.virtualbox.org
.. _YAML:
   http://yaml.org

.. # ------------------( LINKS ~ software : conda           )------------------
.. _Anaconda:
   https://www.anaconda.com/download
.. _Anaconda packages:
   https://anaconda.org
.. _conda-forge:
   https://conda-forge.org

.. # ------------------( LINKS ~ software : pyside2         )------------------
.. _PySide2:
   https://wiki.qt.io/PySide2
.. _PySide2 5.6:
   https://code.qt.io/cgit/pyside/pyside.git/log/?h=5.6
.. _PySide2 installation:
   https://wiki.qt.io/PySide2_GettingStarted
.. _PySide2 PPA:
   https://launchpad.net/~thopiekar/+archive/ubuntu/pyside-git
.. _Qt:
   https://www.qt.io
.. _Qt 5.6:
   https://wiki.qt.io/Qt_5.6_Release

.. # ------------------( LINKS ~ software : icons           )------------------
.. _Cows collection:
   https://thenounproject.com/maxim221/collection/cows
.. _Maxim Kulikov:
   https://thenounproject.com/maxim221
.. _Noun Project:
   https://thenounproject.com
.. _Noun Project license:
   https://thenounproject.com/legal

.. # ------------------( LINKS ~ software : licenses        )------------------
.. _license compatibility:
   https://en.wikipedia.org/wiki/License_compatibility#Compatibility_of_FOSS_licenses
.. _BSD 2-clause license:
   https://opensource.org/licenses/BSD-2-Clause
.. _CC BY 3.0 license:
   https://creativecommons.org/licenses/by/3.0
