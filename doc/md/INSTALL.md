Installation
===========

BETSE is manually installable and uninstallable as follows.

## Requirements

BETSE currently runs _only_ on:

* **64-bit systems**. This is principally due to the increasing obsolescence and
  hence irrelevance of 32-bit systems for scientific work. [Read: no clients or
  developers still use 32-bit systems.] To a lesser extent, this is due to the
  so-called ["3GB barrier"](https://en.wikipedia.org/wiki/3_GB_barrier) imposed
  by most existing 32-bit systems -- including *all* non-server 32-bit editions
  of Microsoft Windows. This barrier prevents usage of more than 3 to 4GB of
  available RAM, which rarely suffices for even small-scale tissue simulations.
* **Python 3.3** or newer (e.g., 3.4, 3.5, 3.6).
* Operating systems matching either:
  * **Microsoft Windows XP** or newer (e.g., Vista, 7, 8, 10).
  * **Apple OS X 10.8.5** (Mountain Lion) or newer (e.g., 10.9, 10.10, 10.11).
  * **Linux distributions providing at least `glibc` 2.19** or newer. (The
    currently installed version of `glibc` is printable by running the following
    command at the command line: `ldd --version`.) This includes but is *not*
    limited to the following Linux distributions:
    * **Linux Mint 17.1** (Rebecca) or newer.
    * **Ubuntu 14.10** (Utopic Unicorn) or newer.

## Recommendations

BETSE currently recommends:

* **At least 4GB RAM**. 

## Dependencies

BETSE has both:

* Mandatory dependencies that must be installed _before_ BETSE itself is
  installed.
* Optional dependencies that may be installed at any time, at the user's
  discretion.

### Mandatory

BETSE requires the following non-pure-Python packages – which themselves
require non-Python precompiled libraries (e.g., C, Fortran) and hence are best
installed manually via the system-wide package manager for your current
operating system:

* Python >= 3.3.
* [Dill](https://github.com/uqfoundation/dill) >= 0.2.3.
* [Matplotlib](http://matplotlib.org) >= 1.4.0.
* [NumPy](http://www.numpy.org) >= 1.8.0.
* [Pillow](https://python-pillow.github.io) >= 2.3.0.
* [PyYaml](http://pyyaml.org) >= 3.10.
* [SciPy](http://www.scipy.org) >= 0.12.0.
* [setuptools](https://pythonhosted.org/setuptools) >= 3.3.
* [six](https://pythonhosted.org/six) >= 1.5.2.

To install these dependencies, the following instructions also install `pip3`,
the Python 3-specific release of the popular Python package manager `pip`.
<sup>Naturally, BETSE itself does _not_ require `pip3` at runtime.</sup>

### Linux

Under Linux, BETSE also requires:

* [Tcl/Tk](https://www.tcl.tk).
* Matplotlib compiled with Tcl/Tk support (i.e., the `tkagg` backend).

#### Linux Debian

Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu), these
dependencies are installable in a system-wide manner as follows:

    $ sudo apt-get install python3-dev python3-dill python3-matplotlib python3-numpy python3-pil python3-pip python3-pyqt5 python3-scipy python3-setuptools python3-six python3-yaml tcl tk

Under some (especially older) Debian-based Linux distributions, the above
instructions may not entirely suffice to satisfy all installation-time or
runtime requirements. Under these ditributions, dependencies may require some
form of recompilation, relinking, or reinstallation.

##### Updated Matplotlib

BETSE requires a fairly recent version of Matplotlib. If the newest version of
Matplotlib installed by your distribution is insufficient, the newest version of
Matplotlib is installable in a system-wide manner as follows:

    $ sudo apt-get uninstall python3-matplotlib &&
      sudo apt-get install gcc gfortran libfreetype6-dev libpng-dev libpython3-all-dev tcl-dev tk-dev &&
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

    $ sudo apt-get install build-essential libatlas-dev libatlas-base-dev libatlas3gf-base &&
      sudo update-alternatives --set libblas.so.3 /usr/lib/atlas-base/atlas/libblas.so.3 &&
      sudo update-alternatives --set liblapack.so.3 /usr/lib/atlas-base/atlas/liblapack.so.3

Note that OpenBLAS and ATLAS _cannot_ be installed at the same time.

###### OpenBLAS

OpenBLAS is a more performant (_but arguably less stable_) optimized BLAS and
LAPACK implementation. While ATLAS is recommended for new users, experienced
users requiring improved performance may benefit from installing OpenBLAS
instead. OpenBLAS is installable in a system-wide manner as follows:

    $ sudo apt-get install build-essential libopenblas-dev &&
      sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3 &&
      sudo update-alternatives --set liblapack.so.3 /usr/lib/lapack/liblapack.so.3

Note that OpenBLAS and ATLAS _cannot_ be installed at the same time.

#### Apple macOS

Under Apple macOS, these dependencies are installable in a system-wide manner
via either:

* **_(Recommended)_ [Homebrew](http://brew.sh),** an unofficial macOS package
  manager. Homebrew provides robust support for features commonly required by
  BETSE developers, including the capacity to install older rather than
  merely the newest versions of packages.
* **_(Not recommended)_ [MacPorts](https://www.macports.org),** an alternative
  unofficial macOS package manager. MacPorts lacks robust support for features
  commonly required by BETSE developers, as described above. Since
  Homebrew and MacPorts install packages into different system directories
  (i.e., `/usr/local` for Homebrew and `/opt` for MacPorts), the two _can_
  technically be used on the same system. However, this is generally
  discouraged. If you currently use and prefer MacPorts, we recommend adopting
  the following instructions to use MacPorts rather than Homebrew.

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
   1. **_(Recommended)_ Manually download and install the CLT:**
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
   1. **_(Not recommended)_ Automatically download and install the CLT.** While
      error-prone and hence discouraged, automatically downloading and
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
   for your preferred shell (e.g., `.bashrc` for Bash, the default OS X shell).

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

#### Microsoft Windows

Under Windows, BETSE additionally requires:

* [Qt4](https://www.qt.io).
* Matplotlib compiled with Qt4 support (i.e., the `qt4agg` backend).

These dependencies are installable under both Microsoft Windows _and_ Wine
prefixes emulating Windows on non-Windows systems (e.g., Linux, OS X).

##### Native

Under native Microsoft Windows, such dependencies are installable as follows.
For simplicity, the following instructions assume use of the
[Miniconda](http://conda.pydata.org/miniconda.html) Python distribution *and*
[Babun](http://babun.github.io) POSIX compatibility layer under 64-bit Windows:

1. Download and install **Babun**, an open-source Cygwin convenience wrapper
   complete with `pact`, a CLI-based package manager for Windows. Due to
   unreconcilable flaws in Windows' non-POSIX-compatible process model, Cygwin
   and hence Babun is incompatible with all Windows applications on the [Big List
   of Dodgy Apps (BLODA)](https://cygwin.com/faq/faq.html#faq.using.bloda).
   Unfortunately, this includes most antivirus software. If Babun begins behaving
   erratically, consider temporarily disabling such software for the duration of
   Babun usage. (This is the fault of neither Babun nor Cygwin!)
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

##### Wine

> **NOTE** While BETSE may indeed be installable _and_ runnable under Wine,
> there's little point in doing so, as the resulting PyInstaller-frozen binaries
> are likely to embed Wine-specific shared libraries unlikely to behave as
> expected under actual Windows systems. Excise this entire subsection, please.

Under non-Windows systems, these dependencies are installable in a system-wide
manner via Wine emulation. This emulation requires the following packages:

* Wine >= 1.7.41. Prior versions of Wine fail to implement Windows API functions
  transitively required by Anaconda (e.g., `GetSystemTime()`).

For simplicity, the following instructions assume use of the
[Miniconda](http://conda.pydata.org/miniconda.html) Python distribution *and*
[PlayOnLinux](https://www.playonlinux.com) Wine manager under 64-bit Ubuntu
Linux. Thanks to the cross-platform portability of both Wine and Anaconda (if
not PlayOnLinux, for obvious reasons), these instructions should trivially
generalize to alternate setups (e.g., 32-bit OS X) as well:

1. Download the 64-bit Python 3 Windows installer for Miniconda, an open-source,
   cross-platform, binary-based Python package manager, from
   `http://conda.pydata.org/miniconda.html`. Miniconda is a minimalist version of
   Anaconda, a popular full-stack SciPy distribution. Whereas Anaconda comes
   pre-bundled with numerous Python and non-Python dependencies, Miniconda
   requires manual installation of such dependencies. We prefer the latter.
   Either suffices, however.
1. If using Linux *and* your preferred Linux distribution provides a readily
   installable package for Wine Staging, consider installing Wine Staging.
1. Install PlayOnLinux, an open-source Wine manager simplifying Wine usage.

        $ sudo apt-get install playonlinux

1. Install the newest 64-bit version of Wine via PlayOnLinux.
  1. Run PlayOnLinux.
  1. Select the *Tools* -> *Manage Wine versions* menu item.
  1. Select the *Wine versions (amd64)* tab.
  1. Select the topmost list item under *Available Wine versions*. Note that this
     version of Wine *must* be greater than or equal to that stipulated above.
  1. Click the right arrow.
  1. Click the *Next* button until complete.
1. Create a new Wine prefix named `betse` via PlayOnLinux.
  1. Click the *Configure* toolbar button.
  1. Click the *New* button.
  1. Click the *Next* button.
  1. Select the *64 bits windows installation* list item and click the *Next*
     button.
  1. Select the list item corresponding to the newly installed version of Wine
     and click the *Next* button.
  1. Enter `betse` and click the *Next* button.
1. Open a terminal window.
1. Activate the newly installed version of Wine, where `${WINE\_VERSION}` should
   be replaced by the installed version number (e.g., `1.7.40`).

        $ export WINEDEBUG='-all'
        $ export WINEPREFIX="${HOME}/.PlayOnLinux/wineprefix/betse"
        $ export PATH="${HOME}/.PlayOnLinux/wine/linux-amd64/${WINE_VERSION}/bin:${PATH}"

1. Install Miniconda via Wine, where `${MINICONDA\_INSTALLER}` should be replaced
  by the path to the previously downloaded Miniconda installer (e.g.,
  `~/Downloads/Miniconda3-latest-Windows-x86\_64.exe`).

        $ wine "${MINICONDA_INSTALLER}"

1. Active Miniconda.

        $ export MINICONDA_HOME="${WINEPREFIX}/drive_c/Miniconda3"
        $ export PATH="${MINICONDA_HOME}:${MINICONDA_HOME}/Scripts:${PATH}"

1. Install dependencies via `conda`, Miniconda's package manager:

        $ wine conda install numpy matplotlib pyside pyyaml pywin32 scipy

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
* [PyDot](https://github.com/erocarrera/pydot) >= 1.0.29 and
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
