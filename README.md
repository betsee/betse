betse
===========

`betse` (BioElectric Tissue Simulation Engine) bio-realistic modelling
of dynamic electrochemical phenomena in gap junction networked cell collectives,
with a focus on spatio-temporal pattern formation.

## Requirements

`betse` currently runs *only* on:

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

`betse` currently recommends but does *not* require:

* **At least 8GB RAM**. Again, this is due to the memory intensiveness of
  even small-scale tissue simulations.

## Dependencies

`betse` has both mandatory dependencies that must be installed before `betse`
itself is installed _and_ optional dependencies that may be installed at the
user's discretion at any time. The uncomputable fun begins now!

### Mandatory

`betse` requires the following non-pure-Python packages – which themselves
require non-Python precompiled libraries (e.g., C, Fortran) and hence are best
installed manually via the system-wide package manager for your current
operating system:

* Python >= 3.3.
* [Matplotlib](http://matplotlib.org) >= 1.4.0.
* [NumPy](http://www.numpy.org) >= 1.8.0.
* [Pillow](https://python-pillow.github.io) >= 2.3.0.
* PySide >= 1.1.0.
* [PyYaml](http://pyyaml.org) >= 3.10.
* [SciPy](http://www.scipy.org) >= 0.12.0.
* [Yamale](https://github.com/23andMe/Yamale) >= 1.5.3.
* [setuptools](https://pythonhosted.org/setuptools) >= 3.3.
* [six](https://pythonhosted.org/six) >= 1.5.2.

To install these dependencies, the following instructions also install `pip3`,
the Python 3-specific release of the popular Python package manager `pip`.
That said, `betse` itself does _not_ require `pip3` at runtime.

### Linux

Under Linux, `betse` also requires:

* [Tcl/Tk](https://www.tcl.tk).
* Matplotlib compiled with Tcl/Tk support (i.e., the `tkagg` backend enabled).

#### Linux Debian

Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu), these
dependencies are installable in a system-wide manner as follows:

    $ sudo apt-get install python3-dev python3-matplotlib python3-numpy python3-pil python3-pip python3-pyqt5 python3-scipy python3-setuptools python3-six python3-yaml tcl tk &&
      sudo pip3 install yamale

Under some (especially older) Debian-based Linux distributions, the above
instructions may not entirely suffice to satisfy all installation-time or
runtime requirements. Under these ditributions, dependencies may require some
form of recompilation, relinking, or reinstallation.

##### Updated Matplotlib

`betse` requires a fairly recent version of Matplotlib. If the newest version of
Matplotlib installed by your distribution is insufficient, the newest version of
Matplotlib may be manually installed as follows:

    $ sudo apt-get uninstall python3-matplotlib &&
      sudo apt-get install gcc gfortran libfreetype6-dev libpng-dev libpython3-all-dev tcl-dev tk-dev &&
      sudo pip3 install matplotlib[all]

##### Optimized BLAS and LAPACK

`betse` strongly recommends that optimized (rather than the unoptimized default)
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

OpenBLAS is a more performant (but often less stable) optimized BLAS and LAPACK
implementation. While ATLAS is recommended for new users, users requiring
improved performance may benefit from installing OpenBLAS instead. OpenBLAS is
installable in a system-wide manner as follows:

    $ sudo apt-get install build-essential libopenblas-dev &&
      sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libopenblas.so.0 &&
      sudo update-alternatives --set liblapack.so.3 /usr/lib/lapack/liblapack.so.3 

Note that OpenBLAS and ATLAS _cannot_ be installed at the same time.

#### Apple OS X

Under Apple OS X, these dependencies are installable in a system-wide manner
via either:

* **(Recommended)** [Homebrew](http://brew.sh), an unofficial OS X package
  manager. Homebrew provides robust support for features commonly required by
  `betse` developers, including the capacity to install older rather than
  merely the newest versions of packages.
* **(Not recommended)** [MacPorts](https://www.macports.org), an alternative
  unofficial OS X package manager. MacPorts lacks robust support for features
  commonly required by `betse` developers, as described above. Since
  Homebrew and MacPorts install packages into different system directories
  (i.e., `/usr/local` for Homebrew and `/opt` for MacPorts), the two _can_
  technically be used on the same system. However, this is generally
  discouraged. If you currently use and prefer MacPorts, we recommend adopting
  the following instructions to use MacPorts rather than Homebrew.

For simplicity, the following instructions assume use of Homebrew:

1. Upgrade your system to the most recently released minor version for your
   currently installed major version of OS X. For example, if your system is
   OS X **10.8.3** (Mountain Lion), upgrade to **10.8.5** (Mountain Lion).
   Homebrew requires recent command-line tools (e.g., `clang`, `gcc`),
   requiring requires recent XCode Command Line Tools (CLT), requiring a recent
   version of XCode, requiring a recent version of OS X. Provided your system
   meets the minimum requirements noted above, it should _not_ be necessary to
   upgrade your system to a newer major version of OS X (e.g., from 10.8.5 to
   10.9.5).
1. Register as an [Apple Developer](https://developer.apple.com). While free,
   registration requires an existing Apple ID and hence ownership of an existing
   Apple product. (We don't make the awful rules. We only complain about them.)
1. If an older version of the XCode Command Line Tools (CLT) has already been
   installed, manually uninstall it _before_ proceeding. While XCode itself is
   safely upgradable merely by installing a new version, the CLT is generally
   _not_. (You can thank Apple for that.)
1. Download and install the most recent version of
   [XCode](https://developer.apple.com/downloads) available for your version of
   OS X. While free, this download requires an [Apple Developer] login. After
   installing Xcode, we recommend performing the following instructions _before_
   attempting to open the installed Xcode application.
1. Open a terminal window (e.g., by running the pre-bundled
   `Applications/Utilities/Terminal.app` application).
1. **(Optional)** Instruct Gatekeeper, the OS X application security manager,
   to implicitly trust without attempting to explicitly verify the installed
   Xcode application. Verification uselessly consumes non-trivial time (in
   upwards of ten minutes, on some systems) _and_ is safely skippable for the
   specific case of Xcode. Note that verification is _not_ safely skippable for
   arbitrary applications downloaded from non-Apple sites.

        $ sudo xattr -d com.apple.quarantine /Applications/Xcode.app

1. Open the installed Xcode application (e.g., by double-clicking
   `Applications/Xcode` from the Finder). If you did _not_ instruct Gatekeeper
   to implicitly trust this application as described above, we recommend a bag
   of greasy popcorn and the "Blade Runner" director's cut. You'll need both.
1. Agree to the Xcode license. This _must_ be done before attempting to run any
   Xcode-provided commands from the terminal (e.g., `clang`, `gcc`, `git`).
1. **(Optional)** Close Xcode.
1. Download and install the **exact same version** of the [XCode Command Line
   Tools](https://developer.apple.com/downloads) (CLT) as the installed version
   of XCode. Attempting to install an older or newer version of the CLT will
   typically superficially succeed but result in obscure and difficult-to-debug
   issues on attempting to install dependencies with Homebrew. Naturally, there
   are numerous approaches to installing the correct version of the CLT:
   1. **(Recommended)** We strongly recommend manually downloading and
      installing the CLT rather than relying on Apple-based automation to do so:
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
   1. **(Not recommended)** The CLT is also automatically downloadable and
      installable via Apple-based automation. If your system has been upgraded
      to both the most recently released minor version of your currently
      installed major version of OS X _and_ to the most recently released
      version of XCode for that version of OS X, the following command _should_
      suffice. When in doubt, prefer the manual approach above instead.

            $ xcode-select –install

1. Download and install [Homebrew](http://brew.sh). While these dependencies are
   also technically installable via [MacPorts](https://www.macports.org),
   Homebrew provides significantly more robust support for features of interest
   to `betse` users. Critically, this includes the capacity to install
   alternative versions of dependencies rather than merely the newest.

        $ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

1. Manually prepend the current `${PATH}` by the absoute paths of all
   directories to which Homebrew installs packages. To do so permanently, append
   the following line to the appropriate startup dotfile in your home directory
   for your preferred shell (e.g., `.bashrc` for Bash, the default OS X shell).

        export PATH="/usr/local/bin:/usr/local/sbin:${PATH}"

1. Activate this `${PATH}` change. Specifically, either:
   * **(Recommended)** Close the current terminal window and open a new terminal
     window.
   * **(Not recommended)** Manually source the modified dotfile: e.g.,

            $ source ~/.bashrc

1. **(Optional)** Inspect your Homebrew installation for potential issues.
   The following command should report that `"Your system is ready to brew."`
   If it does _not_, consider resolving all reported issues before continuing.

        $ brew doctor

1. Install Python 3, update Python package managers, and install all remaining
   dependencies.

        $ brew tap homebrew/python &&
          brew install python3 &&
          pip3 install --upgrade pip setuptools &&
          brew install matplotlib --with-python3 --without-python &&
          brew install numpy --with-python3 --without-python &&
          brew install pillow --with-python3 --without-python &&
          brew install pyside --with-python3 --without-python &&
          brew install scipy --with-python3 --without-python &&
          brew install libyaml &&
          pip3 install pyyaml yamale

Note that Homebrew is a source-based package manager and hence relatively slow.
Expect the installation process to require anywhere from several hours to
several days, depending on hardware performance. We wish we were kidding.

#### Microsoft Windows

Such dependencies are installable under both Microsoft Windows *and* Wine
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
   packages only a subset of the dependencies required by `betse`, Miniconda's
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

        $ conda install numpy matplotlib pyside pyyaml pywin32 scipy

1. Symbolically link the Python 3 executable `python.exe` installed by Miniconda
   to `python3`. For disambiguity, numerous core scripts including `betse`'s
   `setup.py` installer run Python 3 as `python3` rather than `python`. For
   unknown reasons, the Python 3-specific version of Miniconda under Windows
   does *not* appear to provide a `python3` executable. To fix this:

        $ ln -s /c/Miniconda3/python.exe /c/Miniconda3/python3

##### Wine

**FIXME:** While `betse` may indeed be installable _and_ runnable under Wine,
there's little point in doing so, as the resulting PyInstaller-frozen binaries
are likely to embed Wine-specific shared libraries unlikely to behave as
expected under actual Windows systems. Excise this entire subsection, please.

Under non-Windows systems, such dependencies are installable in a system-wide
manner via Wine emulation. Such emulation requires the following packages:

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

`betse` optionally leverages (but does _not_ strictly require) the following
dependencies where available at runtime:

* [py.test](http://pytest.org) >= 2.5.0, for optionally running unit tests.
* [PyInstaller](http://www.pyinstaller.org) >= 3.0, for optionally freezing
  `betse`.
* [UPX](http://upx.sourceforge.net) (any version), for optionally compressing
  frozen `betse` executables.

These dependencies are installable as follows.

#### `py.test`

To optionally [run tests](#testing), `betse` requires `py.test`, a pure-Python 
test harness. This dependency is installable in a system-wide manner as follows:

* Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

        $ sudo apt-get install python3-pytest

* Under Apple OS X:

        $ pip3 install pytest

##### `py.test` Plugins

While optional, `betse` provides out-of-the-box support for the following
third-party `py.test` plugins:

* `pytest-xdist`, parallelizing test runs across all available processors.
  `py.test` itself provides _no_ built-in support for parallelization! Since
  `betse`'s test suite is computationally expensive (if not prohibitive), this
  plugin is a hard prerequisite for sanity preservation.

Contributors are strongly encouraged to install these optional dependencies,
which `betse`'s test suite will then implicitly detect and set accordingly.
These dependencies are installable in a system-wide manner as follows:

    $ pip3 install pytest-xdist

#### PyInstaller

To optionally [freeze `betse`](#freezing), `betse` requires PyInstaller, a
non-pure-Python cross-platform command-line utility for freezing Python
applications. This dependency is installable in a system-wide manner as
follows:

* Under any Linux distribution:

        $ sudo pip3 install pyinstaller

* Under Apple OS X:

        $ pip3 install pyinstaller

#### UPX (Ultimate Packer for eXecutables)

To optionally compress executables while [freezing `betse`](#freezing), `betse`
requires UPX, a non-Python cross-platform command-line utility for compressing
arbitrary executables. This dependency is installable in a system-wide manner
as follows:

* Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

        $ sudo apt-get install upx-ucl

* Under Apple OS X:

        $ brew install upx

## Installation

`betse` itself is installable into either:

* A system-wide directory accessible to all users of the current system.
* A venv (i.e., virtual environment) isolated to the current user.

The latter has the advantage of avoiding conflicts with already installed
system-wide Python and non-Python packages (e.g., in the event that `betse`
requires different versions of such packages), but the corresponding
disadvantage of requiring reinstallation of such packages and all transitive
dependencies of such packages. Since several dependencies are heavy-weight
(e.g., Qt4) and hence costly to reinstall, this is a notable disadvantage.

Note that the string `${BETSE\_DIR}` should be replaced everywhere below by the
absolute path of the top-level directory containing the source for `betse`.

### System-wide

`betse` is installable into a system-wide directory as follows:

* Compile `betse`.

        $ cd "${BETSE_DIR}"
        $ python3 setup.py build

* Install `betse`.

        $ sudo python3 setup.py easy_install --no-deps .

Curiously, although the `develop` command for `setuptools` provides a
`--no-deps` option, the `install` command does not. Hence, the `easy_install`
command is called above instead.

`betse` is subsequently uninstallable via `pip` as follows:

    $ sudo pip uninstall betse

### User-specific

`betse` is installable into a user-specific venv by running the following
command **from within such venv**:

    $ cd "${BETSE_DIR}"
    $ ./setup.py install

This command should *not* be run outside of a venv. Doing so will reinstall all
dependencies of `betse` already installed by the system-wide package manager
(e.g., `apt-get`). This may superficially appear to work but invites obscure and
difficult to debug conflicts at `betse` runtime between dependencies reinstalled
by `setuptools` and dependencies already installed by such package maneger.

`betse` is subsequently uninstallable via `pip` as follows:

    $ pip uninstall betse

## Usage

`betse` is a front-facing application rather than backend framework. While
`betse`'s Python packages are importable by other packages, `betse` is typically
run by executing Python wrapper scripts installed to the current `${PATH}`.

### CLI

`betse` installs a low-level command line interface (CLI), runnable as follows:

    $ betse

Such CLI is also runnable by invoking the main CLI module under the active
Python 3 interpreter as follows:

    $ cd "${BETSE_DIR}"
    $ python3 -m betse.cli

This has the minor advantage of working regardless of whether `betse` has been
installed or not, but the corresponding disadvantage of requiring more typing.

### GUI

`betse` also installs a high-level graphical user interface (GUI) implemented in
the popular cross-platform windowing toolkit Qt4, runnable as follows:

    $ betse-qt &

Such GUI is also runnable by invoking the main GUI module under the active
Python 3 interpreter as follows:

    $ cd "${BETSE_DIR}"
    $ python3 -m betse.gui

## Development

For development purposes, `betse` is *editably installable* (i.e., as a symbolic
link rather than physical copy). As the name implies, editable installations are
modifiable at runtime and hence suitable for development. Thanks to the magic
of symbolic links, changes to the copy of `betse` from which an editable
installation was installed will be silently prograpagated back to such
installation.

### System-wide

`betse` is installable into a system-wide directory as follows:

* **(Optional).** Set the current umask to "002" as above.
    $ umask 002
* Editably install `betse`.
    $ cd "${BETSE_DIR}"
    $ sudo python3 setup.py symlink

The `symlink` command is a `betse`-specific `setuptools` command inspired by the
IPython `setuptools` command of the same name, generalizing the behaviour of the
default `develop` command to system-wide editable installations.

Why? Because the `develop` command is suitable *only* for user-specific editable
installations. While both `pip` and `setuptools` provide commands for performing
editable installations (e.g., `sudo pip3 install --no-deps --editable .` and
`sudo python3 setup.py develop --no-deps`, respectively), executable scripts
installed by such commands raise fatal exceptions on failing to find
`setuptools`-installed dependencies regardless of whether such dependencies have
already been installed in a system-wide manner. To quote [IPython developer
MinRK](http://mail.scipy.org/pipermail/ipython-dev/2014-February/013209.html):

    So much hate for setuptools right now.  I can't believe `--no-deps` skips
    dependency installation, but still adds a redundant check to entry points.

### User-specific

`betse` is editably installable into a user-specific venv via either `pip` or
`setuptools` **from within such venv.** While there appears to be no particular
advantage to using one over the other, it remains helpful to note that both
apply. In either case, external executables (e.g., `betse`, `betse-qt`) will
also be installed and usable in the expected manner.

#### pip

`betse` is editably installable into a user-specific venv via `pip` as follows:

    $ cd "${BETSE_DIR}"
    $ pip3 install --no-deps --editable .

Such installation is uninstallable as follows:

    $ pip3 uninstall betse

#### setuptools

`betse` is editably installable into a user-specific venv via `setuptools` as
follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py develop --no-deps

Such installation is uninstallable as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py develop --uninstall

## Testing

`betse` is testable via `py.test` as follows. Either:

* Run all available tests.

        $ cd "${BETSE_DIR}"
        $ ./setup.py test

* Run all tests matching a passed Python-evaluatable expression. For example, to
  run all test functions and classes whose names contain either `test_tartarus`
  _or_ `test_thalassa`:

        $ cd "${BETSE_DIR}"
        $ ./setup.py test -k 'test_tartarus or test_thalassa'

## Freezing

`betse` is **freezable** (i.e., convertable to platform-specific executable
binaries distributable to end users) via PyInstaller, a cross-platform open-
source utility for performing such freezing. While other such utilities exist
(e.g., the Windows-specific `py2exe`, the OS X-specific `py2app`), only
PyInstaller satisfies the basic requirements of both platform portability and
freedom (as in both beer and speech). `betse` offers setuptools-automated
support for *only* PyInstaller.

### Dependencies

`betse`'s setuptools-automated freezing support has additional dependencies
beyond those required for core `betse` installation and usage. In no particular
order, these are:

#### `betse` Installation

`betse` *must* be installed in a system-wide manner first -- either:

* Editably (e.g., via `sudo python3 setup.py symlink`).
* Non-editably (e.g., via `sudo python3 setup.py install`).

See above for further details.

#### PyInstaller

The `scipy` branch of [Cecil Curry](https://github.com/leycec)'s [unofficial
PyInstaller fork](https://github.com/leycec/pyinstaller) *must* be installed as
follows:

* Make a local directory to which the PyInstaller codebase will be downloaded
  and change to such directory: e.g.,

        $ mkdir ~/py
        $ cd ~/py

* Download the PyInstaller codebase.

        $ git clone --branch scipy https://github.com/leycec/pyinstaller.git

* Change to the downloaded directory.

        $ cd pyinstaller

* Install PyInstaller.

        $ sudo python3 setup.py install

This branch provides critical patches submitted to but *not* yet accepted into
PyInstaller's [official Python 3
branch](https://github.com/pyinstaller/pyinstaller/tree/python3). If this branch
is remotely updated (e.g., with new patches{, consider reinstalling PyInstaller
as follows:

* Change to the PyInstaller codebase: e.g.,

        $ cd ~/py/pyinstaller

* Update the PyInstaller codebase.
* Reinstall PyInstaller.

        $ sudo python3 setup.py install

#### UPX

The Ultimate Packer for Executables (UPX) should ideally also be installed,
though this is *not* strictly necessary. UPX is a cross-platform utility for
compressing executables. When installed, PyInstaller automatically compresses
all written executables.

##### Windows

UPX does *not* appear to be readily available under Microsoft Windows. Next!

##### Linux

Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu), UPX is
installable in a system-wide manner as follows:

    $ sudo apt-get install upx-ucl

### Usage

`betse` is freezable in either:

* One-file mode, in which case PyInstaller generates one executable file for
  each of `betse`'s wrapper scripts (e.g., `betse`, `betse-qt`).
* One-directory mode, in which case PyInstaller generates one directory
  containing one executable file (along with various ancillary libraries and
  archives) for each of `betse`'s wrapper scripts.

`betse` is freezable in one-file mode as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py freeze_file

`betse` is freezable in one-directory mode as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py freeze_dir

### Paths

Assuming one-file mode, executables will be frozen to the following paths:

* Under Linux:
  * `freeze/dist/betse` for `betse`'s CLI wrapper.
  * `freeze/dist/betse-qt` for `betse`'s GUI wrapper.
* Under OS X:
  * `freeze/dist/betse` for `betse`'s CLI wrapper.
  * `freeze/dist/betse-qt.app` for `betse`'s GUI wrapper.
* Under Windows:
  * `freeze/dist/betse.exe` for `betse`'s CLI wrapper.
  * `freeze/dist/betse-qt.exe` for `betse`'s GUI wrapper.

Such executables may be moved (and hence distributed to end users) as is but
should *not* be renamed. Doing so typically invalidates code signing (as well as
other embedded metadata).

### Caveats

While commonly regarded as the best utility for freezing Python applications,
PyInstaller is *not* without the occasional caveat. Most if not all such caveats
apply to alternative utilities (e.g., `cx_Freeze`, `py2app`, `py2exe`).

#### Switching Freeze Modes

For safety, attempting to switch freeze modes (i.e., to refreeze `betse` in
one-directory mode when already frozen in one-file mode *or* in one-file mode
when already frozen in one-directory mode) currently results in an error.

This is correctable by either manually removing the previously frozen file or
directory *or* passing option `--clean` to the desired `freeze_file` or
`freeze_dir` command, which automatically removes such path on your behalf. For
example, `betse` is switchable from one-directory to one-file mode and then back
again as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py freeze_dir
    $ ./setup.py freeze_file --clean
    $ ./setup.py freeze_dir --clean

#### Forwards Incompatibilities

Executables frozen on older versions of supported operating systems are
*typically* compatible with newer versions of the same systems. For example,
executables frozen on Ubuntu 12.04 (Precise Pangolin) are typically compatible
with Ubuntu 12.10, 14.04, 14.10, and newer. This is commonly referred to as
**backwards compatibility**.

The converse is *not* the case. That is, executables frozen on newer versions of
supported operating systems are guaranteeably incompatible with older versions
of the same systems. For example, executables frozen on OS X 10.10 (Yosemite)
are guaranteeably incompatible with OS X 10.9, 10.8, 10.7, and older. This is
commonly referred to as **forwards incompatibility**.

For this reason, executables should be frozen on the oldest possible versions
of supported operating systems themselves supporting all application
dependencies. The simplest means of doing so is to **host** (i.e., install and
run) such versions of such systems as virtual guests of the current system via a
**hypervisor** (i.e., an application creating and running virtual machines).

For this purpose, we recommend Oracle's open-source hypervisor VirtualBox rather
than VMware's closed-source hypervisor VMware Workstation. The former supports
all operating systems supported by `betse` out of the box, including OS X; the
latter supports Linux and Windows but *not* OS X out of the box. Moreover, the
former is free (as in both beer and speech); the latter is non-free (as in both
beer and speech).

#### No Cross-freezing or -compilation

PyInstaller supports neither **cross-freezing** (i.e., generation of
executables intended for execution on operating systems other than the current)
nor **cross-compilation** (i.e., generating of executables intended for
execution on CPU architectures other than the current) out of the box. Hence,
executables will be executable only under systems matching the current operating
system and architecture. As example, if the current system is 64-bit Linux,
executables generated under such system will be executable only under other
64-bit Linux systems.

In practice, the lack of cross-compilation support is ignorable. `betse` is
sufficiently memory-intensive that running under a 32-bit architecture would be
strongly inadvisable and hence unsupported.

The lack of cross-freezing support is *not* ignorable. While there exist well-
documented means of cross-freezing Windows executables on Linux systems (e.g.,
via the popular emulation layer Wine), there exist *no* means of cross-freezing
OS X executables on Linux. The conventional solution is to install OS X as a
virtual guest of an existing Linux host and perform freezing under such guest.
We recommend the open-source product VMware for this purpose.

## License

`betse` is open-source software licensed under the BSD-compatible [NetBSD
license](https://www.netbsd.org/about/redistribution.html). equivalent to the
[2-clause BSD license](http://opensource.org/licenses/BSD-2-Clause) excluding
the following proviso:

    The views and conclusions contained in the software and documentation are
    those of the authors and should not be interpreted as representing official
    policies, either expressed or implied, of the FreeBSD Project.

We regard this proviso as evasive legerdemain. If "the views and conclusions
contained in the software and documentation" become sufficiently desynchronized
from official project policies as to constitute a legal risk (e.g., due to
permissive defamation law such as exists in the United Kingdom), the remedy is
to bridge the gap -- not to decry its existence. That is, to amend the software,
its policies, or both until no desynchronization between the two remains.

The software _is_ the project. Strident claims to the contrary notwithstanding,
the software _is_ fundamentally representative of project policies. This is
always and irrevocably the case, and claims otherwise merit little standing.
