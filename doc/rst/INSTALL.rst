.. # ------------------( SYNOPSIS                           )------------------

.. # FIXME: Merge ``doc/md/DEVELOP.md`` into this file; then, remove that file.
.. # FIXME: Merge ``doc/md/INSTALL.rst`` into this file; then, remove both that
.. # and the ``doc/md/INSTALL.md`` file.

============
Installation
============

BETSE requires:

-  Either **Microsoft Windows,** **Apple macOS,** or a **Linux distribution.**
   [#os_not]_
-  At least **Python 3.5** (e.g., 3.5, 3.8). [#python_not]_
-  At least **4GB RAM;** ideally, at least **16GB RAM.** [#thirtytwobit_not]_

.. [#os_not]
   All other platforms (e.g., Android, FreeBSD) are explicitly unsupported at
   this time.

.. [#python_not]
   All prior Python versions (e.g., Python 2.7, 3.4) are explicitly
   unsupported as well.

.. [#thirtytwobit_not]
   Ergo, BETSE effectively requires a **64-bit system.**  32-bit systems impose
   a so-called `"3GB barrier" <https://en.wikipedia.org/wiki/3_GB_barrier>`__
   preventing usage of more than 3—4GB of available RAM, which rarely suffices
   for even small-scale BETSE simulations.

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

Linux
=====

BETSE is installable under modern Linux distributions as follows:

#. Install **Python 3.x.** [#python2_not]_ [#python3_default]_
#. Open a **terminal.** [#linux_terminal]_
#. Copy-and-paste these commands into this terminal:

   #. Install **BETSE.**

      .. code-block:: console

         pip3 install betse

   #. [\ *Optional*\ ] Test **BETSE.**

      .. code-block:: console

         betse -v try

.. [#python2_not]
   Do *not* install **Python 2.7.** BETSE strictly requires **Python 3.x.**

.. [#python3_default]
   Most Linux distributions install Python 3.x by default. Since Python ≥ 3.4
   installs pip_ by default, the Python 3.x version of pip_ should already be
   installed by default – in theory, anyway.

.. [#linux_terminal]
   To open a terminal under:

   - **Ubuntu Linux:**

     #. Type ``Ctrl``\ +\ ``Alt``\ +\ ``t``.

macOS
=====

BETSE is installable under Apple macOS as follows:

#. Install the **Python 3.x** [#python2_not]_ (e.g., 3.7) variant of Anaconda_.
   [#macos_anaconda_not]_
#. Open the **Finder.**
#. Open the **Applications** folder.
#. Open the **Utilities** folder.
#. Open **Terminal.app.**
#. Copy-and-paste these commands into this terminal:

   #. **Enable** conda-forge_.

      .. code-block:: console

         conda config --add channels conda-forge

   #. Install **BETSE.** [#conda_package]_

      .. code-block:: console

         conda install betse

   #. [\ *Optional*\ ] Test **BETSE.**

      .. code-block:: console

         betse -v try

.. [#macos_anaconda_not]
   If you prefer *not* to install Anaconda_, BETSE is also `manually installable
   via third-party package managers for macOS <#macos-homebrew_> (e.g.,
   Homebrew_, MacPorts_). Doing so is non-trivial and, if performed improperly,
   could produce a performance-crippled single-core installation of BETSE.
   Anaconda_ suffers no such issues and is guaranteed to produce a
   performance-optimized multicore installation of BETSE. We strongly recommend
   Anaconda_ – even when you think you know better.

.. [#conda_package]
   This command installs both the `most recent stable release of BETSE <conda
   package_>`__ *and* all mandatory and most optional dependencies of this
   release. Older stable releases are installable in a similar manner (e.g.,
   ``conda install betse=0.7.0`` for BETSE 0.7.0). All `Anaconda packages`_ are
   kindly hosted by the `non-profit conda-forge organization <conda-forge_>`__.

Windows
=======

BETSE is installable under Microsoft Windows 10 as follows: [#windows_not]_

#. Emulate **Ubuntu Linux** via the `Windows Subsystem for Linux (WSL)
   <WSL_>`__ bundled with Windows 10.
#. Open an **Ubuntu Linux terminal.**
#. Copy-and-paste these commands into this terminal:

   #. **Update** the `Advanced Package Tool (APT) <APT_>`__ package cache.
 
     .. code-block:: console
 
        sudo apt update
   
   #. **Upgrade** all previously installed packages.
 
     .. code-block:: console
 
        sudo apt upgrade
 
   #. Install **Python 3.x.**\ [#python2_not]_
 
     .. code-block:: console
 
        sudo apt install python3-pip
 
   #. Install **BETSE.**
 
     .. code-block:: console
 
        pip3 install betse
 
   #. [\ *Optional*\ ] Test **BETSE.**
 
     .. code-block:: console
 
        betse -v try

.. [#windows_not]
   The `Windows Subsystem for Linux (WSL) <WSL_>`__ and (thus BETSE itself) is
   *only* installable under **Windows 10.** Under older Windows versions, BETSE
   may be installed from a `virtual Linux guest <VirtualBox_>`__.

Developers
==========

BETSE also supports Git_\ -based development as follows: [#dev_portable]_

#. Install **Python 3.x.** [#python2_not]_
#. Install Git_.
#. **Register** and **sign in** to a `GitLab account`_.
#. `Fork BETSE <BETSE fork_>`__ via this account.
#. Click the **Clone** button on this fork's front page.
#. Copy the URL given under the **Clone with SSH** heading.
#. Open a **terminal.**
#. Copy-and-paste these commands from this terminal:

   #. **Clone** the ``master`` branch of this fork into the current directory,
      replacing ``{fork_url}`` below with the URL copied above (e.g., ``git
      clone git@gitlab.com:muh_username/betse.git``).

      .. code-block:: console

         git clone {fork_url}

   #. Install BETSE **editably.** [#dev_editable]_

      .. code-block:: console

         cd betse && pip3 install -e .

   #. Create a new **feature branch,** replacing ``{branch_name}`` below with a
      string unique to this fork (e.g., ``git checkout -b muh_branch``).
      [#feature_branch] 

      .. code-block:: console

         git checkout -b {branch_name}
   
   #. **Change** this fork as desired.
   #. **Stage** and **commit** these changes.

      .. code-block:: console

         git commit -a

   #. **Push** these changes to the remote copy of your fork hosted at GitLab_.

      .. code-block:: console

         git push origin {branch_name}

#. **Browse** back to `BETSE's official project page <BETSE_>`__.
#. Click the **Create merge request** button.
#. Voilà! Instant open-source volunteerism.

.. [#dev_portable]
   For portability, these instructions assume the **Python 3.x** version of
   pip_ has already been installed in a platform-specific manner. While
   Anaconda_ may also be leveraged for Git_\ -based development, doing so
   exceeds the limited scope of these instructions.

.. [#dev_editable]
   An editable installation creates a symbolic link to this clone such that
   code changes are applied immediately *without* requiring reinstallation. A
   standard installation only installs a physical copy of this clone such that
   code changes are ignored until this clone is reinstalled.

.. [#feature_branch]
   Changes should *never* be committed directly to the ``master`` branch.
   Changes should *only* be committed to **feature branches** (i.e., branches
   of this clone containing all code changes needed to implement a new
   feature and/or improve an existing feature). The name of each feature branch
   *must* be unique to this clone but is otherwise arbitrary.

Docker
======

BETSE is also installable into a `Docker container`_, circumventing the need to
install BETSE directly into a host operating system. For simplicity, the
following instructions assume a modern Linux distribution: [#docker_not]_

#. `Install Docker <Docker install_>`__.
#. Instruct the **Xauthority security mechanism** to ignore hostnames, enabling
   `Docker containers`_ with different hostnames than that of the local host to
   access the current X11 socket. [#docker_thanks]_

   .. code-block:: console

      touch /tmp/.docker.xauth && xauth nlist :0 |
          sed -e 's/^..../ffff/' |
          xauth -f /tmp/.docker.xauth nmerge -

#. **Download** the latest version of the official `Anaconda 3 Docker image`_
   and **instantiate** this image as a new `Docker container`_ named ``betse``
   running an interactive Bash session mounting the X11 socket of the host's
   current X11 session.

   .. code-block:: console

      docker run -it\
          --name betse\
          -v /tmp/.X11-unix:/tmp/.X11-unix\
          -v /tmp/.docker.xauth:/tmp/.docker.xauth\
          -e DISPLAY=$DISPLAY\
          -e XAUTHORITY=/tmp/.docker.xauth\
          continuumio/anaconda3 bash

#. Copy-and-paste these commands into this container's terminal:

   #. [\ *Optional*\ ] Test the **X11 connection** by running ``xeyes``.

      .. code-block:: console

         apt-get update && apt-get install -y x11-apps && xeyes

   #. Download the live version of **BETSE** into the ``${HOME}`` directory of
      the current user (i.e., ``root``).

      .. code-block:: console

         cd ~ && git clone https://gitlab.com/betse/betse.git

   #. Install **BETSE** editably.

      .. code-block:: console

         cd betse && pip3 install -e .

   #. [\ *Optional*\ ] Test **BETSE** by running a sample simulation.

      .. code-block:: console

         cd /tmp && betse try && rm -rf sample_sim

   #. **Exit** this session.

      .. code-block:: console

         exit

To resume the previously instantiated container:

#. **Restart** this container.

   .. code-block:: console

      docker start betse

#. **Reenter** this container by running another interactive Bash session.

   .. code-block:: console

      docker attach betse

.. [#docker_not]
   While Docker_ is also installable under macOS or Windows, doing so exceeds
   the limited scope of these instructions.

.. [#docker_thanks]
   Thanks to `Jürgen Weigert
   <https://stackoverflow.com/users/3936284/j%c3%bcrgen-weigert>`__ for his
   observant `Stackoverflow answer <https://stackoverflow.com/a/25280523>`__
   inspiring this snippet.

Alternative Instructions
========================

.. warning::

   **The alternative installation instructions listed below are outdated.**
   Unlike the official instructions above, these unofficial instructions are
   intended more as inspiration for command-line aficionados than gospel.
   Copy-and-pasting these instructions into a terminal blindly *will* yield a
   broken installation of BETSE in the best case and a broken operating system,
   terminal, or Python environment in the worst case. *There be dragons here.*

BETSE is installable with only `one or two simple commands on all supported
platforms <BETSE install_>`__ – complete with multicore-aware hardware
optimizations. These commands strictly adhere to scientific standards for
Python packaging, including the standard pip_ Python package manager *and*
cross-platform Anaconda_ Python distribution.

For advanced users preferring to manually install dependencies with
platform-specific package managers (e.g., APT_), BETSE may also be manually
installed in a platform-specific manner. This approach has the obvious
advantage of cleanly integrating with existing packaging regimes but the
non-obvious disadvantage of typically installing a single-core version of BETSE
with *no* multicore-aware hardware optimizations. *That's bad.* 

Due to the difficulty of manually installing BETSE in a multicore-aware manner,
our `simple installation instructions <BETSE install_>`__ are *strongly*
recommended. For completeness, these instructions nonetheless detail the manual
approach for several popular package managers.

Linux
-----

BETSE is manually installable with *most* Linux package managers.

Debian
~~~~~~

Under `Debian <https://www.debian.org>`__-based Linux distributions
(e.g., `Linux Mint <https://www.linuxmint.com>`__,
`Ubuntu <https://www.ubuntu.com>`__), all mandatory dependencies are
installable in a system-wide manner as follows:

.. code-block:: console

   sudo apt-get install python3-dev python3-dill python3-matplotlib \
       python3-numpy python3-pil python3-pip python3-scipy python3-setuptools \
       python3-six python3-yaml tcl tk

Under some (especially older) `Debian <https://www.debian.org>`__-based Linux
distributions, the above instructions may not entirely suffice to satisfy all
installation-time or runtime requirements. Under these distributions,
dependencies may require some form of recompilation, relinking, or
reinstallation.

Matplotlib
++++++++++

BETSE requires a fairly recent version of Matplotlib_. If the newest version of
Matplotlib_ installed by your distribution is insufficient, the newest version
of Matplotlib_ is installable in a system-wide manner as follows:

.. code-block:: console

   sudo apt-get uninstall python3-matplotlib &&
   sudo apt-get install gcc gfortran libfreetype6-dev libpng-dev \
       libpython3-all-dev tcl-dev tk-dev &&
   sudo pip3 install matplotlib[all]

BLAS and LAPACK
+++++++++++++++

BETSE strongly recommends that optimized (rather than the unoptimized default)
implementations of the BLAS and LAPACK APIs for linear algebra be used. While
there exist numerous alternatives both open-source (e.g., CBLAS) and
proprietary (e.g., MKL), the following instructions assume use of either ATLAS
or OpenBLAS.

ATLAS
^^^^^

Automatically Tuned Linear Algebra Software (ATLAS) is the standard
baseline for all optimized BLAS and LAPACK implementations. ATLAS is
installable in a system-wide manner as follows:

.. code-block:: console

   sudo apt-get install build-essential libatlas-base-dev

Note that OpenBLAS and ATLAS *cannot* be installed at the same time.

OpenBLAS
^^^^^^^^

OpenBLAS is a more performant (*but arguably less stable*) optimized
BLAS and LAPACK implementation. While ATLAS is recommended for new
users, experienced users requiring improved performance may benefit from
installing OpenBLAS instead. OpenBLAS is installable in a system-wide
manner as follows:

.. code-block:: console

   sudo apt-get install build-essential libopenblas-dev

Note that OpenBLAS and ATLAS *cannot* be installed at the same time.

Gentoo
~~~~~~

Under `Gentoo <https://www.gentoo.org>`__-based Linux distributions
(e.g., `Chrome OS <https://en.wikipedia.org/wiki/Chrome_OS>`__,
`Sabayon <https://www.sabayon.org>`__), all mandatory and optional
dependencies are installable in a system-wide manner as follows:

#. **Install** the ``eselect repository`` module.

   .. code-block:: console

      emerge --ask app-eselect/eselect-repository &&
          mkdir -p /etc/portage/repos.conf

#. **Add** and **synchronize** the `raiagent overlay`_, religiously maintained
   by a `BETSE co-maintainer <https://github.com/leycec>`__.

   .. code-block:: console

      eselect repository enable raiagent && emerge --sync raiagent

#. Install both **BETSE** and **BETSEE,** our official PySide2_-based GUI.

   .. code-block:: console

      emerge --autounmask betsee

macOS (Homebrew)
----------------

Under Apple macOS, all mandatory dependencies are installable in a
system-wide manner with either:

-  [\ *Recommended*\ ] Homebrew_, an unofficial package manager for macOS.
   Homebrew_ provides robust support for features commonly required by BETSE
   developers, including the capacity to install older rather than merely the
   newest versions of packages.
-  MacPorts_, another unofficial package manager for macOS. MacPorts_ lacks
   robust support for features commonly required by BETSE developers, as
   described above. Since Homebrew_ and MacPorts_ install packages into
   different system directories (i.e., ``/usr/local`` for Homebrew_, ``/opt``
   for MacPorts_), the two *can* technically be used on the same system.
   However, this is generally discouraged. If you currently use and prefer
   MacPorts_, consider adopting the following instructions to use MacPorts_
   rather than Homebrew_.

For simplicity, the following instructions assume use of Homebrew_:

#. Register as an `Apple Developer <https://developer.apple.com>`__. While
   free, registration requires an existing Apple ID and hence ownership of an
   existing Apple product. *We don't make the awful rules. We only complain
   about them.*
#. **Upgrade your system** to the most recently released minor version for
   your currently installed major version of macOS. For example, if your
   system is macOS **10.8.3** (\ *Mountain Lion*\ ), upgrade to **10.8.5**
   (\ *Mountain Lion*\ ). Homebrew_ requires recent command-line tools (e.g.,
   ``clang``, ``gcc``), requiring requires recent XCode Command Line Tools
   (CLT), requiring a recent version of XCode, requiring a recent version of
   macOS. Provided your system meets the minimum requirements noted above, it
   should *not* be necessary to upgrade your system to a newer major version
   of macOS (e.g., from 10.8.5 to 10.9.5).
#. **Open a terminal window** (e.g., by running the pre-bundled
   ``Applications/Utilities/Terminal.app`` application).
#. If an older version of the XCode Command Line Tools (CLT) has already been
   installed, `manually uninstall
   <https://stackoverflow.com/questions/27438457/xcode-6-1-how-to-uninstall-command-line-tools>`__
   the CLT. While XCode itself is safely upgradable merely by installing a new
   version, the CLT generally is not. *You can thank Apple for that.*
#. Download and install the most recent version of `XCode
   <https://developer.apple.com/downloads>`__ available for your version of
   macOS. While free, this download requires an Apple Developer login.
#. [\ *Optional*\ ] After installing Xcode, perform the following *before*
   running Xcode:
#. **Instruct Gatekeeper to implicitly trust Xcode.** Gatekeeper is the macOS
   application security manager. By default, Gatekeeper uselessly verifies
   Xcode via a labouriously time-consuming and safely skippable bureaucratic
   process requiring in upwards of twenty minutes on lower-end laptops. Note
   that verification is *not* safely skippable for "dubious" applications
   downloaded from third-party sources.

   .. code-block:: console

      $ sudo xattr -d com.apple.quarantine /Applications/Xcode.app

#. **Run Xcode** (e.g., by double-clicking ``Applications/Xcode`` from the
   Finder). If you did *not* instruct Gatekeeper to implicitly trust this
   application as described above, grab a bag of greasy popcorn and `Blade
   Runner (The Final Cut)
   <https://en.wikipedia.org/wiki/Versions_of_Blade_Runner>`__. *You'll need
   both.*
#. **Agree to the Xcode license.** This *must* be done before attempting to run
   any Xcode-bundled commands from the terminal (e.g., ``clang``, ``gcc``,
   ``git``).
#. [\ *Optional*\ ] Close Xcode.
#. Download and install the exact same version of the `XCode Command Line Tools
   <https://developer.apple.com/downloads>`__ (CLT) as the installed version of
   XCode. Attempting to install an older or newer version of the CLT may
   superficially succeed but *will* result in obscure and difficult-to-debug
   issues on attempting to install dependencies with Homebrew_ or MacPorts_.
   There are various approaches to installing the correct version of the CLT –
   some inherently safer than others. Either:

   - [\ *Recommended*\ ] **Manually download and install the CLT:**

     #. Browse to the `Apple Developer Downloads
        <https://developer.apple.com/downloads>`__ site.
     #. Enter ``xcode`` into the search bar.
     #. Manually search the resulting hits for the installed version of XCode.
     #. Note the official date of this version's release (e.g., June 12, 2013
        for XCode 4.6.3).
     #. Manually search the resulting hits for the most recent version of the
        CLT *preceding* this date (e.g., April 11, 2013 for the CLT
        corresponding to XCode 4.6.3).
     #. Download and install this version.

   - [\ *Not cecommended*\ ] **Automatically download and install the CLT.**
     While error-prone and hence discouraged, automatically downloading and
     installing the CLT with Apple-based automation *is* technically feasible
     in common edge cases. If your system has been upgraded to both the most
     recently released minor version of the currently installed major version
     of macOS *and* the most recently released version of XCode for that
     version of macOS, the following command *should* suffice. If in doubt,
     prefer the manual approach listed above instead.

     .. code-block:: console

        xcode-select –install

#. Download and install Homebrew_. While these dependencies are also
   technically installable via MacPorts_, Homebrew_ provides significantly more
   robust support for features of interest to BETSE users. Critically, this
   includes the capacity to install alternative versions of dependencies rather
   than merely the newest.

   .. code-block:: console

      ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

#. Manually prepend the current ``${PATH}`` by the absoute paths of all
   directories to which Homebrew_ installs packages. To do so permanently,
   append the following line to the appropriate startup dotfile in your home
   directory for your preferred shell (e.g., ``.bashrc`` for Bash, the default
   macOS shell).

   .. code-block:: console

      export PATH="/usr/local/bin:/usr/local/sbin:${PATH}"

#. Activate this ``${PATH}`` change. Either:

   -  [\ *Recommended*\ ] Close the current terminal window and open a new
      terminal window.
   -  [\ *Not recommended*\ ] Manually source the modified dotfile: e.g.,
   
      .. code-block:: console
   
         source ~/.bashrc

1. [\ *Optional*\ ] Inspect your Homebrew_ installation for potential issues.
   The following command should report that ``"Your system is ready to brew."``
   If it does *not*, consider resolving all reported issues before continuing.

   .. code-block:: console

      brew doctor

2. **Install all dependencies.**

   .. code-block:: console

      brew tap homebrew/python &&
          brew install python3 &&
          pip3 install --upgrade pip setuptools wheel &&
          brew install matplotlib --with-python3 --without-python &&
          brew install numpy --with-python3 --without-python &&
          brew install pillow --with-python3 --without-python &&
          brew install scipy --with-python3 --without-python &&
          brew install libyaml &&
          pip3 install dill pyyaml

Note that Homebrew_ is a source-based package manager and hence relatively slow.
Expect the installation process to require anywhere from several hours to
several days, depending on hardware performance. We wish we were kidding.

Note also that these instructions link NumPy_ against the most optimized
multicore implementation of the BLAS and LAPACK APIs available under macOS as
of this writing: Apple's `Accelerate Framework
<https://developer.apple.com/reference/accelerate/1668466-blas>`__. No further
BLAS or LAPACK configuration is required or recommended.

Optional Dependencies
=====================

BETSE optionally leverages (but does *not* strictly require) the following
dependencies where available at runtime:

- `NetworkX <https://networkx.github.io>`__ >= 1.11, for optionally
  analyzing BETSE networks.
- `pprofile <https://github.com/vpelletier/pprofile>`__ >= 1.8, for
  optionally profiling BETSE in a line-granular manner.
- `ptpython <https://github.com/jonathanslenders/ptpython>`__ >= 0.29,
  for optionally wrapping the BETSE REPL with an improved interface. By
  default, the BETSE REPL leverages the stock Python REPL.
- `py.test <http://pytest.org>`__ >= 2.8.0, for optionally running unit
  tests.
- `PyDot <https://github.com/erocarrera/pydot>`__ >= 1.0.28 and
  `GraphViz <http://www.graphviz.org>`__ >= 2.38, for optionally
  visualizing BETSE networks.
- `PyInstaller <http://www.pyinstaller.org>`__ >= 3.0, for optionally
  freezing BETSE.
- `UPX <http://upx.sourceforge.net>`__ (any version), for optionally
  compressing frozen BETSE executables.

These dependencies are installable as follows.

NetworkX
--------

To optionally analyze networks (e.g., gene regulatory, biochemical reaction),
BETSE requires NetworkX, a pure-Python graph theoretic framework. This
dependency is installable in a system-wide manner as follows:

- Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

  .. code-block:: console

     $ sudo apt-get install python3-networkx

- Under all other supported platforms:

  .. code-block:: console

     $ pip3 install networkx

pprofile
--------

To optionally profile the BETSE codebase with line-granularity into `callgrind
<http://kcachegrind.sourceforge.net/>`__-compatible profile files, BETSE
requires ``pprofile``, an advanced pure-Python line profiler. This dependency
is installable in a system-wide manner as follows:

- Under all supported platforms:

  .. code-block:: console

     $ pip3 install pprofile

ptpython
--------

To optionally wrap the BETSE REPL with an improved interface providing syntax
highlighting, multiline editing, autocompletion, and presumably more, BETSE
requires ``ptpython``, an advanced pure-Python REPL. This dependency is
installable in a system-wide manner as follows:

- Under all supported platforms:

  .. code-block:: console

     $ pip3 install ptpython

py.test
-------

To optionally `run tests <#testing>`__, BETSE requires ``py.test``, a
pure-Python test harness. This dependency is installable in a system-wide
manner as follows:

- Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

  .. code-block:: console

     $ sudo apt-get install python3-pytest

- Under all other supported platforms:

  .. code-block:: console

     $ pip3 install pytest

Plugins
~~~~~~~

While optional, BETSE provides out-of-the-box support for the following
third-party ``py.test`` plugins:

- ``pytest-xdist``, parallelizing test runs across all available processors.
  ``py.test`` itself provides *no* built-in support for parallelization! Since
  BETSE's test suite is computationally expensive (if not prohibitive), this
  plugin is a hard prerequisite for sanity preservation.

Contributors are strongly encouraged to install these optional dependencies,
which BETSE's test suite will then implicitly detect and set accordingly. These
dependencies are installable in a system-wide manner as follows:

- Under all other supported platforms:

  .. code-block:: console

     $ pip3 install pytest-xdist

PyDot + GraphViz
----------------

To optionally visualize networks (e.g., gene regulatory, biochemical reaction),
BETSE requires both:

- PyDot, a high-level pure-Python GraphViz wrapper.
- GraphViz, a low-level C-based graph theoretic visualizer.

These dependencies are installable in a system-wide manner as follows:

- For PyDot:

  - Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):
  
    .. code-block:: console
  
       $ sudo apt-get install python3-pydot
  
  - Under all other supported platforms:
  
    .. code-block:: console
  
       $ pip3 install pydot

- For GraphViz:

  - Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):
  
    .. code-block:: console
  
       $ sudo apt-get install graphviz
  
  - Under Apple macOS:
  
    .. code-block:: console
  
       $ brew install graphviz

PyInstaller
-----------

To optionally `freezing BETSE <BETSE freeze_>`__, BETSE requires PyInstaller, a
non-pure-Python cross-platform command-line utility for freezing Python
applications. This dependency is installable in a system-wide manner as
follows:

- Under all supported platforms:

  .. code-block:: console

     $ pip3 install pyinstaller

UPX 
---

To optionally compress executables while `freezing BETSE <BETSE freeze_>`__,
BETSE requires the Ultimate Packer for eXecutables (UPX), a non-Python
cross-platform command-line utility for compressing arbitrary executables. This
dependency is installable in a system-wide manner as follows:

- Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

  .. code-block:: console

     $ sudo apt-get install upx-ucl

- Under Apple OS X:

  .. code-block:: console

     $ brew install upx

.. # ------------------( LINKS ~ betse                      )------------------
.. _BETSE:
   https://gitlab.com/betse/betse
.. _BETSE fork:
   https://gitlab.com/betse/betse/-/forks/new
.. _BETSE freeze:
   /doc/md/FREEZE.md
.. _BETSE install:
   /README.rst#installation

.. # ------------------( LINKS ~ host : gitlab              )------------------
.. _GitLab:
   https://gitlab.com
.. _GitLab account:
   https://gitlab.com/users/sign_in
.. _GitLab-CI:
   https://about.gitlab.com/gitlab-ci

.. # ------------------( LINKS ~ os : linux                 )------------------
.. _APT:
   https://en.wikipedia.org/wiki/Advanced_Packaging_Tool
.. _POSIX:
   https://en.wikipedia.org/wiki/POSIX

.. # ------------------( LINKS ~ os : gentoo                )------------------
.. _raiagent overlay:
   https://github.com/leycec/raiagen

.. # ------------------( LINKS ~ os : ubuntu                )------------------
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

.. # ------------------( LINKS ~ soft                       )------------------
.. _Atom:
   https://atom.io
.. _dill:
   https://pypi.python.org/pypi/dill
.. _FFmpeg:
   https://ffmpeg.org
.. _Git:
   https://git-scm.com/downloads
.. _Graphviz:
   http://www.graphviz.org
.. _Libav:
   https://libav.org
.. _MEncoder:
   https://en.wikipedia.org/wiki/MEncoder
.. _VirtualBox:
   https://www.virtualbox.org
.. _YAML:
   http://yaml.org

.. # ------------------( LINKS ~ soft : docker              )------------------
.. _Docker:
.. _Docker container:
.. _Docker containers:
   https://www.docker.com
.. _Docker install:
   https://docs.docker.com/install

.. # ------------------( LINKS ~ soft : py                  )------------------
.. _imageio:
   https://imageio.github.io
.. _Matplotlib:
   http://matplotlib.org
.. _NumPy:
   http://www.numpy.org
.. _Python 3:
   https://www.python.org
.. _pip:
   https://pip.pypa.io
.. _py.test:
   http://pytest.org
.. _SciPy:
   http://www.scipy.org

.. # ------------------( LINKS ~ soft : py : conda          )------------------
.. _Anaconda:
   https://www.anaconda.com/download
.. _Anaconda packages:
   https://anaconda.org
.. _Anaconda 3 Docker image:
   https://hub.docker.com/r/continuumio/anaconda3
.. _conda-forge:
   https://conda-forge.org

.. # ------------------( LINKS ~ soft : py : pyside2        )------------------
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

