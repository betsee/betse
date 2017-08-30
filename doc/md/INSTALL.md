Installation
===========

BETSE is installable with only [two simple steps](/README.rst) on all supported
platforms – complete with multicore-aware hardware optimizations. These steps
leverage scientific standards for Python packaging, including the cross-platform
**[Anaconda](https://www.continuum.io/downloads)** Python distribution *and*
**[`pip`](https://pypi.python.org/pypi/pip)** Python package manager.

For advanced users preferring to manually install dependencies with existing
package managers (e.g.,
[APT](https://en.wikipedia.org/wiki/Advanced_Packaging_Tool)) rather than
[Anaconda](https://www.continuum.io/downloads), BETSE may also be manually
installed in a platform-specific manner. This approach has the obvious advantage
of cleanly integrating with existing packaging regimes but the non-obvious
disadvantage of typically installing a single-core version of BETSE with *no*
multicore-aware hardware optimizations. <sup>_That's bad._</sup>

Due to the difficulty of manually installing BETSE in a multicore-aware manner,
the [simple installation instructions](/README.rst) are *strongly* recommended.
For completeness, this document nonetheless details the manual approach for
several popular package managers.

## Requirements

BETSE requires:

* Either **Windows,** **macOS,** or **Linux.** All other platforms (e.g.,
  Android, FreeBSD) are explicitly unsupported.
* At least **Python 3.4** (e.g., 3.4, 3.5, 3.6). All prior Python versions (e.g.,
  Python 2.7, 3.3) are explicitly unsupported.
* A **64-bit system.** This is principally due to the increasing obsolescence and
  hence irrelevance of 32-bit systems for scientific work. [Read: no clients or
  developers still use 32-bit systems.] To a lesser extent, this is due to the
  so-called ["3GB barrier"](https://en.wikipedia.org/wiki/3_GB_barrier) imposed
  by most existing 32-bit systems -- including *all* non-server 32-bit editions
  of Microsoft Windows. This barrier prevents usage of more than 3 to 4GB of
  available RAM, which rarely suffices for even small-scale tissue simulations.

## Recommendations

For optimal performance, BETSE also recommends:

* At least **4GB RAM;** ideally, at least **16GB RAM.**

## Dependencies

Like most Python applications, BETSE has both mandatory dependencies that must
be installed before BETSE itself is installed *and* optional dependencies safely
installable at any time.

Unlike most Python applications, these dependencies should *not* be installed
with [`pip`](https://pypi.python.org/pypi/pip), the standard Python package
manager. This is especially the case for scientific dependencies partially
implemented in C or Fortran (e.g., [NumPy](http://www.numpy.org),
[SciPy](http://www.scipy.org)). While technically feasible, installing these
dependencies with [`pip`](https://pypi.python.org/pypi/pip) typically results in
an unoptimized single-core installation of BETSE. <sup>_Again, that's bad._</sup>

### Mandatory

BETSE requires the following non-pure-Python packages – which themselves
require non-Python precompiled libraries (e.g., C, Fortran) and hence are best
installed manually via the system-wide package manager for your current
operating system:

* Python >= 3.4.
* [Dill](https://github.com/uqfoundation/dill) >= 0.2.3.
* [NumPy](http://www.numpy.org) >= 1.8.0.
* [Pillow](https://python-pillow.github.io) >= 2.3.0.
* [PyYaml](http://pyyaml.org) >= 3.10.
* [SciPy](http://www.scipy.org) >= 0.12.0.
* [matplotlib](http://matplotlib.org) >= 1.5.0.
* [setuptools](https://pythonhosted.org/setuptools) >= 3.3.
* [six](https://pythonhosted.org/six) >= 1.5.2.

To install these dependencies, the following instructions also install `pip3`,
the Python 3-specific release of the popular Python package manager `pip`.
<sup>Naturally, BETSE itself does _not_ require `pip3` at runtime.</sup>

#### Matplotlib Backend

Under Linux and Windows, BETSE also requires at least one of the following
third-party widget toolkits and corresponding matplotlib backends (*in
descending order of preference*):

* **[Tcl/Tk](https://www.tcl.tk)** *and* matplotlib compiled with Tcl/Tk support
  (i.e., the `TkAgg` backend). This backend comes highly recommended for optimal
  display of BETSE plots and animations under both Linux and Windows.
* **[PyQt5](https://www.riverbankcomputing.com/software/pyqt/download5)** *and*
  matplotlib compiled with Qt5 support (i.e., the `Qt5Agg` backend).
* **[PyQt4](https://www.riverbankcomputing.com/software/pyqt/download)** *and*
  matplotlib compiled with Qt4 support (i.e., the `Qt4Agg` backend).

Under macOS, BETSE instead prefers the `MacOSX` matplotlib backend. This backend
comes highly recommended for optimal display of BETSE plots and animations under
macOS, requiring no additional third-party widget toolkits.

## Linux

BETSE is manually installable with *most* Linux-centric package managers.

### Debian

Under [Debian](https://www.debian.org)-based Linux distributions (e.g., [Linux
Mint](https://www.linuxmint.com), [Ubuntu](https://www.ubuntu.com)), all
mandatory dependencies are installable in a system-wide manner as follows:

    $ sudo apt-get install python3-dev python3-dill python3-matplotlib \
      python3-numpy python3-pil python3-pip python3-scipy python3-setuptools \
      python3-six python3-yaml tcl tk

Under some (especially older) [Debian](https://www.debian.org)-based Linux
distributions, the above instructions may not entirely suffice to satisfy all
installation-time or runtime requirements. Under these distributions,
dependencies may require some form of recompilation, relinking, or
reinstallation.

##### Updated Matplotlib

BETSE requires a fairly recent version of matplotlib. If the newest version of
matplotlib installed by your distribution is insufficient, the newest version of
matplotlib is installable in a system-wide manner as follows:

    $ sudo apt-get uninstall python3-matplotlib &&
      sudo apt-get install gcc gfortran libfreetype6-dev libpng-dev \
        libpython3-all-dev tcl-dev tk-dev &&
      sudo pip3 install matplotlib[all]

##### Optimized BLAS and LAPACK

BETSE strongly recommends that optimized (rather than the unoptimized default)
implementations of the BLAS and LAPACK APIs for linear algebra be used. While
there exist numerous alternatives both open-source (e.g., CBLAS) and
proprietary (e.g., MKL), the following instructions assume use of either ATLAS
or OpenBLAS.

###### ATLAS

Automatically Tuned Linear Algebra Software (ATLAS) is the standard baseline
for all optimized BLAS and LAPACK implementations. ATLAS is installable in a
system-wide manner as follows:

    $ sudo apt-get install build-essential libatlas-base-dev &&
      sudo update-alternatives --set libblas.so.3 /usr/lib/atlas-base/atlas/libblas.so.3 &&
      sudo update-alternatives --set liblapack.so.3 /usr/lib/atlas-base/atlas/liblapack.so.3

Note that OpenBLAS and ATLAS *cannot* be installed at the same time.

###### OpenBLAS

OpenBLAS is a more performant (_but arguably less stable_) optimized BLAS and
LAPACK implementation. While ATLAS is recommended for new users, experienced
users requiring improved performance may benefit from installing OpenBLAS
instead. OpenBLAS is installable in a system-wide manner as follows:

    $ sudo apt-get install build-essential libopenblas-dev &&
      sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3 &&
      sudo update-alternatives --set liblapack.so.3 /usr/lib/lapack/liblapack.so.3

Note that OpenBLAS and ATLAS _cannot_ be installed at the same time.

### Gentoo

Under [Gentoo](https://www.gentoo.org)-based Linux distributions (e.g., [Chrome
OS](https://en.wikipedia.org/wiki/Chrome_OS),
[Sabayon](https://www.sabayon.org)), all mandatory and optional dependencies are
installable in a system-wide manner as follows:

1. **Install [`layman`](https://wiki.gentoo.org/wiki/Layman),** Gentoo's
   official overlay manager.

        $ sudo emerge layman
        $ sudo echo 'source /var/lib/layman/make.conf' >> /etc/portage/make.conf

* **Add the [`raiagent`](https://github.com/leycec/raiagent) overlay,**
  religiously maintained by a [BETSE co-maintainer](https://github.com/leycec).

        $ sudo layman -a raiagent

* **Synchronize overlays.**

        $ sudo layman -S

* Either:

  * **_(Recommended)_** Install the [optimized BLAS and LAPACK
    stack](https://wiki.gentoo.org/wiki/User_talk:Houseofsuns) published by the
    [`science`](https://github.com/gentoo-science/sci) overlay. While
    technically optional, failing to do so *will* reduce BETSE to unoptimized
    single-core behavior. To properly install this stack, see these
    [authoritative
    instructions](https://wiki.gentoo.org/wiki/User_talk:Houseofsuns).
  * Disable the `smp` USE flag enabled by default for BETSE. In this case, the
    default unoptimized BLAS and LAPACK stack will be linked against instead.

            $ sudo echo 'sci-biology/betse -smp' >> /etc/portage/package.use

* **Unmask BETSE.** Either:

  * **_(Recommended)_** Unmask the most recent stable release of BETSE.

            $ sudo echo '>=sci-biology/betse-0.4.1' >> /etc/portage/package.accept_keywords

  * Unmask the most recent unstable commit to the BETSE `git` repository.

            $ sudo echo '>=sci-biology/betse-0.4.1 **' >> /etc/portage/package.accept_keywords

* **Install BETSE.**

        $ sudo emerge betse

#### macOS

Under Apple macOS, all mandatory dependencies are installable in a system-wide
manner with either:

* **_(Recommended)_** **[Homebrew](http://brew.sh),** an unofficial 
  package manager for macOS. Homebrew provides robust support for features
  commonly required by BETSE developers, including the capacity to install older
  rather than merely the newest versions of packages.
* **[MacPorts](https://www.macports.org),** another unofficial package manager
  for macOS. MacPorts lacks robust support for features commonly required by
  BETSE developers, as described above. Since Homebrew and MacPorts install
  packages into different system directories (i.e., `/usr/local` for Homebrew
  and `/opt` for MacPorts), the two _can_ technically be used on the same
  system. However, this is generally discouraged. If you currently use and
  prefer MacPorts, we recommend adopting the following instructions to use
  MacPorts rather than Homebrew.

For simplicity, the following instructions assume use of Homebrew:

1. **Register as an [Apple Developer](https://developer.apple.com).** While
   free, registration requires an existing Apple ID and hence ownership of an
   existing Apple product. <sup>_We don't make the awful rules. We only complain
   about them._</sup>
1. **Upgrade your system** to the most recently released minor version for your
   currently installed major version of macOS. For example, if your system is
   macOS **10.8.3** (_Mountain Lion_), upgrade to **10.8.5** (_Mountain Lion_).
   Homebrew requires recent command-line tools (e.g., `clang`, `gcc`),
   requiring requires recent XCode Command Line Tools (CLT), requiring a recent
   version of XCode, requiring a recent version of macOS. Provided your system
   meets the minimum requirements noted above, it should _not_ be necessary to
   upgrade your system to a newer major version of macOS (e.g., from 10.8.5 to
   10.9.5).
1. **Open a terminal window** (e.g., by running the pre-bundled
   `Applications/Utilities/Terminal.app` application). All commands prefixed by
   `$` below _must_ be run from within a terminal window. Note that, by Unix
   convention, the `$` prefix only denotes the default Bash shell prompt and
   should _not_ actually be typed (e.g., type `xcode-select –install` rather
   than `$ xcode-select –install` when asked to do so below). Likewise, the
   `<return>` key should be typed after each such command.
1. If an older version of the XCode Command Line Tools (CLT) has already been
   installed, **[manually
   uninstall](https://stackoverflow.com/questions/27438457/xcode-6-1-how-to-uninstall-command-line-tools)
   the CLT.** While XCode itself is safely upgradable merely by installing a
   new version, the CLT generally is _not_. <sup>_You can thank Apple for
   that._</sup>
1. **Download and install the most recent version of
   [XCode](https://developer.apple.com/downloads)** available for your version
   of macOS. While free, this download requires an Apple Developer login.
1. **_(Optional)_** After installing Xcode, perform the following _before_
   running Xcode:
   1. **Instruct Gatekeeper to implicitly trust Xcode.** Gatekeeper is the macOS
      application security manager. By default, Gatekeeper uselessly verifies
      Xcode via a labouriously time-consuming and safely skippable bureaucratic
      process requiring in upwards of twenty minutes on lower-end laptops. Note
      that verification is _not_ safely skippable for "dubious" applications
      downloaded from third-party sources.

            $ sudo xattr -d com.apple.quarantine /Applications/Xcode.app

1. **Run Xcode** (e.g., by double-clicking `Applications/Xcode` from the
   Finder). If you did _not_ instruct Gatekeeper to implicitly trust this
   application as described above, grab a bag of greasy popcorn and [_Blade
   Runner (The Final
   Cut)_](https://en.wikipedia.org/wiki/Versions_of_Blade_Runner). You'll need
   both.
1. **Agree to the Xcode license.** This _must_ be done before attempting to run
   any Xcode-bundled commands from the terminal (e.g., `clang`, `gcc`, `git`).
1. **_(Optional)_ Close Xcode.**
1. **Download and install the exact same version of the [XCode Command Line
   Tools](https://developer.apple.com/downloads) (CLT)** as the installed
   version of XCode. Attempting to install an older or newer version of the CLT
   may superficially succeed but _will_ result in obscure and difficult-to-debug
   issues on attempting to install dependencies with Homebrew or MacPorts.
   There are various approaches to installing the correct version of the CLT –
   some inherently safer than others. Either:
   1. **_(Recommended)_** **Manually download and install the CLT:**
      1. Browse to the [Apple Developer
         Downloads](https://developer.apple.com/downloads) site.
      1. Enter `xcode` into the search bar.
      1. Manually search the resulting hits for the installed version of XCode.
      1. Note the official date of this version's release (e.g., June 12, 2013
         for XCode 4.6.3).
      1. Manually search the resulting hits for the most recent version of the
         CLT _preceding_ this date (e.g., April 11, 2013 for the CLT
         corresponding to XCode 4.6.3).
      1. Download and install this version.
   1. **_(Not recommended)_** **Automatically download and install the CLT.**
      While error-prone and hence discouraged, automatically downloading and
      installing the CLT with Apple-based automation _is_ technically feasible
      in common edge cases. If your system has been upgraded to both the most
      recently released minor version of the currently installed major version
      of macOS _and_ the most recently released version of XCode for that
      version of macOS, the following command _should_ suffice. If in doubt,
      prefer the manual approach listed above instead.

            $ xcode-select –install

1. **Download and install [Homebrew](http://brew.sh).** While these dependencies
   are also technically installable via [MacPorts](https://www.macports.org),
   Homebrew provides significantly more robust support for features of interest
   to BETSE users. Critically, this includes the capacity to install
   alternative versions of dependencies rather than merely the newest.

        $ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

1. **Manually prepend the current `${PATH}`** by the absoute paths of all
   directories to which Homebrew installs packages. To do so permanently, append
   the following line to the appropriate startup dotfile in your home directory
   for your preferred shell (e.g., `.bashrc` for Bash, the default macOS shell).

        export PATH="/usr/local/bin:/usr/local/sbin:${PATH}"

1. **Activate this `${PATH}` change.** Either:
   * **_(Recommended)_** Close the current terminal window and open a new terminal
     window.
   * **_(Not recommended)_** Manually source the modified dotfile: e.g.,

            $ source ~/.bashrc

1. **_(Optional)_** Inspect your Homebrew installation for potential issues.
   The following command should report that `"Your system is ready to brew."`
   If it does _not_, consider resolving all reported issues before continuing.

        $ brew doctor

1. **Install all dependencies.**

        $ brew tap homebrew/python &&
          brew install python3 &&
          pip3 install --upgrade pip setuptools wheel &&
          brew install matplotlib --with-python3 --without-python &&
          brew install numpy --with-python3 --without-python &&
          brew install pillow --with-python3 --without-python &&
          brew install scipy --with-python3 --without-python &&
          brew install libyaml &&
          pip3 install dill pyyaml

Note that Homebrew is a source-based package manager and hence relatively slow.
Expect the installation process to require anywhere from several hours to
several days, depending on hardware performance. We wish we were kidding.

Note also that these instructions link `numpy` against the most optimized
multicore implementation of the BLAS and LAPACK APIs available under macOS as of
this writing: Apple's **[Accelerate
Framework](https://developer.apple.com/reference/accelerate/1668466-blas).** No
further BLAS or LAPACK configuration is required or recommended.

#### Windows

> **Note:** these instructions are _woefully_ inadequate at present, encouraging
> installation of the [Cygwin](https://www.cygwin.com)-based
> **[Babun](https://babun.github.io)** wrapper rather than usage of the existing
> **[Bash on ubuntu on
> Windows](https://msdn.microsoft.com/en-us/commandline/wsl/about)** environment
> bundled with the Windows 10's Anniversary Update. Until these instructions are
> updated accordingly, Windows users are _strongly_ encouraged to follow the
> [simple installation instructions](/README.rst) instead.

Under Microsoft Windows, all mandatory dependencies are installable in a
system-wide manner via any number of POSIX compatibility layers. For simplicity,
the following instructions assume use of the
[Miniconda](http://conda.pydata.org/miniconda.html) Python distribution *and*
[Babun](http://babun.github.io) POSIX compatibility layer under 64-bit Windows:

1. Download and install **[Babun](https://babun.github.io)**, an open-source
   [Cygwin](https://www.cygwin.com) convenience wrapper complete with `pact`, a
   CLI-based package manager for Windows. Due to unreconcilable flaws in
   Windows' non-POSIX-compatible process model, Cygwin and hence Babun is
   incompatible with all Windows applications on the [Big List of Dodgy Apps
   (BLODA)](https://cygwin.com/faq/faq.html#faq.using.bloda).  Unfortunately,
   this includes most antivirus software. If Babun begins behaving erratically,
   consider temporarily disabling such software for the duration of Babun usage.
   (This is the fault of neither Babun nor Cygwin!)
1. Download and install the 64-bit Python 3 Windows version of **Miniconda**.
   (See the "Wine" subsection below for further details.)
1. Double click the desktop shortcut `babun` to open a new terminal window.
1. Prioritize Miniconda- over Babun-installed Python packages. By default, Babun
   prioritizes Babun- over Miniconda-installed Python packages. Since Babun
   packages only a subset of the dependencies required by BETSE, Miniconda's
   `conda` rather than Babun's `pact` package manager must be used to install
   such dependencies. To permit this, modify the `${PATH}` global exported at
   Babun startup by editing the `.zshrc` file in your home directory as follows:

        # Alter this...
        export PATH=$HOME/bin:/usr/local/bin:$PATH

        # ...to this.
        export MINICONDA_HOME="/cygdrive/c/Miniconda3"
        export PATH="${MINICONDA_HOME}:${MINICONDA_HOME}/Scripts:${HOME}/bin:${PATH}"

1. Apply such changes to the current shell session.

        $ source ~/.zshrc

1. Install Python dependencies via `conda`, Miniconda's package manager:

        $ conda install dill numpy matplotlib pyside pyyaml pywin32 scipy

1. Symbolically link the Python 3 executable `python.exe` installed by Miniconda
   to `python3`. For disambiguity, numerous core scripts including BETSE's
   `setup.py` installer run Python 3 as `python3` rather than `python`. For
   unknown reasons, the Python 3-specific version of Miniconda under Windows
   does *not* appear to provide a `python3` executable. To fix this:

        $ ln -s /c/Miniconda3/python.exe /c/Miniconda3/python3

### Optional

BETSE optionally leverages (but does _not_ strictly require) the following
dependencies where available at runtime:

* [NetworkX](https://networkx.github.io) >= 1.11, for optionally analyzing
  BETSE networks.
* [pprofile](https://github.com/vpelletier/pprofile) >= 1.8, for optionally
  profiling BETSE in a line-granular manner.
* [ptpython](https://github.com/jonathanslenders/ptpython) >= 0.29, for
  optionally wrapping the BETSE REPL with an improved interface. <sup>By
  default, the BETSE REPL leverages the stock Python REPL.</sup>
* [py.test](http://pytest.org) >= 2.8.0, for optionally running unit tests.
* [PyDot](https://github.com/erocarrera/pydot) >= 1.0.28 and
  [GraphViz](http://www.graphviz.org) >= 2.38, for optionally visualizing BETSE
  networks.
* [PyInstaller](http://www.pyinstaller.org) >= 3.0, for optionally freezing
  BETSE.
* [UPX](http://upx.sourceforge.net) (any version), for optionally compressing
  frozen BETSE executables.

These dependencies are installable as follows.

#### NetworkX

To optionally analyze networks (e.g., gene regulatory, biochemical reaction),
BETSE requires NetworkX, a pure-Python graph theoretic framework. This
dependency is installable in a system-wide manner as follows:

* Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

        $ sudo apt-get install python3-networkx

* Under all other supported platforms: <sup>Under Linux, additionally prefix
  this command by `sudo`.</sup>

        $ pip3 install networkx

#### `pprofile`

To optionally profile the BETSE codebase with line-granularity into
[callgrind](http://kcachegrind.sourceforge.net/)-compatible profile files,
BETSE requires `pprofile`, an advanced pure-Python line profiler. This
dependency is installable in a system-wide manner as follows:

* Under all supported platforms: <sup>Under Linux, additionally prefix this
  command by `sudo`.</sup>

        $ pip3 install pprofile

#### `ptpython`

To optionally wrap the BETSE REPL with an improved interface providing syntax
highlighting, multiline editing, autocompletion, and presumably more, BETSE
requires `ptpython`, an advanced pure-Python REPL. This dependency is
installable in a system-wide manner as follows:

* Under all supported platforms: <sup>Under Linux, additionally prefix this
  command by `sudo`.</sup>

        $ pip3 install ptpython

#### `py.test`

To optionally [run tests](#testing), BETSE requires `py.test`, a pure-Python
test harness. This dependency is installable in a system-wide manner as follows:

* Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

        $ sudo apt-get install python3-pytest

* Under all other supported platforms: <sup>Under Linux, additionally prefix
  this command by `sudo`.</sup>

        $ pip3 install pytest

##### `py.test` Plugins

While optional, BETSE provides out-of-the-box support for the following
third-party `py.test` plugins:

* `pytest-xdist`, parallelizing test runs across all available processors.
  `py.test` itself provides _no_ built-in support for parallelization! Since
  BETSE's test suite is computationally expensive (if not prohibitive), this
  plugin is a hard prerequisite for sanity preservation.

Contributors are strongly encouraged to install these optional dependencies,
which BETSE's test suite will then implicitly detect and set accordingly.
These dependencies are installable in a system-wide manner as follows:

* Under all other supported platforms: <sup>Under Linux, additionally prefix
  this command by `sudo`.</sup>

        $ pip3 install pytest-xdist

#### PyDot + GraphViz

To optionally visualize networks (e.g., gene regulatory, biochemical reaction),
BETSE requires:

* PyDot, a high-level pure-Python GraphViz wrapper.
* GraphViz, a low-level C-based graph theoretic visualizer.

These dependencies are installable in a system-wide manner as follows:

* For PyDot:
  * Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

            $ sudo apt-get install python3-pydot

  * Under all other supported platforms: <sup>Under Linux, additionally prefix
    this command by `sudo`.</sup>

            $ pip3 install pydot

* For GraphViz:
  * Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

            $ sudo apt-get install graphviz

  * Under Apple OS X:

            $ brew install graphviz

#### PyInstaller

To optionally [freeze BETSE](#freezing), BETSE requires PyInstaller, a
non-pure-Python cross-platform command-line utility for freezing Python
applications. This dependency is installable in a system-wide manner as
follows:

* Under all supported platforms: <sup>Under Linux, additionally prefix this
  command by `sudo`.</sup>

        $ sudo pip3 install pyinstaller

#### UPX (Ultimate Packer for eXecutables)

To optionally compress executables while [freezing BETSE](#freezing), BETSE
requires UPX, a non-Python cross-platform command-line utility for compressing
arbitrary executables. This dependency is installable in a system-wide manner
as follows:

* Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

        $ sudo apt-get install upx-ucl

* Under Apple OS X:

        $ brew install upx

## Installation

BETSE itself is installable into either:

* A system-wide directory accessible to all users of the current system.
* A venv (i.e., virtual environment) isolated to the current user.

The latter has the advantage of avoiding conflicts with already installed
system-wide Python and non-Python packages (e.g., in the event that BETSE
requires different versions of such packages), but the corresponding
disadvantage of requiring reinstallation of such packages and all transitive
dependencies of such packages. Since several dependencies are heavy-weight
(e.g., Qt4) and hence costly to reinstall, this is a notable disadvantage.

Note that the string `${BETSE\_DIR}` should be replaced everywhere below by the
absolute path of the top-level directory containing the source for BETSE.

### System-wide

BETSE is installable into a system-wide directory as follows:

* Compile BETSE.

        $ cd "${BETSE_DIR}"
        $ python3 setup.py build

* Install BETSE.

        $ sudo python3 setup.py easy_install --no-deps .

Curiously, although the `develop` command for `setuptools` provides a
`--no-deps` option, the `install` command does not. Hence, the `easy_install`
command is called above instead.

BETSE is subsequently uninstallable via `pip` as follows:

    $ sudo pip uninstall betse

### User-specific

BETSE is installable into a user-specific
[venv](https://docs.python.org/3/library/venv.html) by running the following
commands **from within that venv**:

    $ cd "${BETSE_DIR}"
    $ ./setup.py install

This command should *not* be run outside of a venv. Doing so will reinstall all
dependencies of BETSE already installed by the system-wide package manager
(e.g., `apt-get`). This may superficially appear to work but invites obscure and
difficult to debug conflicts at BETSE runtime between dependencies reinstalled
by `setuptools` and dependencies already installed by such package maneger.

BETSE is subsequently uninstallable via `pip` as follows:

    $ pip uninstall betse

### Docker

BETSE is also installable into a [Docker](https://www.docker.com)-hosted Linux
distribution contained within an existing Linux distribution, circumventing the
need to install BETSE directly into an existing system. For simplicity, the
following instructions assume usage of Docker's official Ubuntu image:

1. Install [**Docker**](https://docs.docker.com/engine/installation).
1. Instruct the Xauthority security mechanism to ignore hostnames, permitting
   Docker containers with different hostnames than that of the local host to
   access the current X11 socket. <sup>_Thanks to [Jürgen
   Weigert](https://stackoverflow.com/users/3936284/j%c3%bcrgen-weigert) for his
   observant [Stackoverflow answer](https://stackoverflow.com/a/25280523)
   inspiring this snippet.</sup>

        $ touch /tmp/.docker.xauth && xauth nlist :0 |
              sed -e 's/^..../ffff/' |
              xauth -f /tmp/.docker.xauth nmerge -

1. Download the latest version of [Continuum
   Analytics](https://www.continuum.io/downloads)' official [Anaconda 3 Docker
   image](https://hub.docker.com/r/continuumio/anaconda3) and instantiate this
   image as a new Docker container named `betse` running an interactive Bash
   session mounting the X11 socket of the host's current X11 session.

        $ docker run -it\
              --name betse\
              -v /tmp/.X11-unix:/tmp/.X11-unix\
              -v /tmp/.docker.xauth:/tmp/.docker.xauth\
              -e DISPLAY=$DISPLAY\
              -e XAUTHORITY=/tmp/.docker.xauth\
              continuumio/anaconda3 bash

1. Run the following commands from within this container:
   1. **_(Optional)_** Test the X11 connection by running `xeyes`.

            $ apt-get update && apt-get install -y x11-apps && xeyes

   1. Download the live version of BETSE into the `${HOME}` directory of the
      current user (i.e., `root`).

            $ cd ~ && git clone https://gitlab.com/betse/betse.git

   1. Install BETSE.

            $ cd betse && python3 setup.py install

   1. **_(Optional)_** Test BETSE by running a sample simulation.

            $ cd /tmp && betse try && rm -rf sample_sim

   1. Exit this session.

            $ exit

To resume the previously instantiated container:

1. Restart this container.

        $ docker start betse

1. Reenter this container by running another interactive Bash session.

        $ docker attach betse
