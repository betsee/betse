betse
===========

`betse` (BioElectric Tissue Simulation Engine) bio-realistic modelling
of dynamic electrochemical phenomena in gap junction networked cell collectives,
with a focus on spatio-temporal pattern formation.

## System Requirements

`betse` currently runs *only* on:

* **64-bit systems**. This is principally due to the increasing obsolescence and
  hence irrelevance of 32-bit systems for scientific work. [Read: no clients or
  developers use 32-bit systems.] To a lesser extent, this is due to the so-
  called ["3GB barrier"](https://en.wikipedia.org/wiki/3_GB_barrier) imposed by
  most existing 32-bit systems -- including *all* non-server 32-bit editions of
  Microsoft Windows. Such barrier prevents usage of more than 3 to 4GB of
  available RAM, which rarely suffices for even small-scale tissue simulations.
* Operating systems matching either:
  * **Microsoft Windows XP** or newer.
  * **Apple OS X 10.6** (Snow Leopard) or newer.
  * **Linux distributions providing at least `glibc` 2.19** or newer. (The
    currently installed version of `glibc` is printable by running the following
    command at the command line: `ldd --version`.) This includes but is *not*
    limited to the following Linux distributions:
    * **Linux Mint 17.1** (Rebecca) or newer.
    * **Ubuntu 14.10** (Utopic Unicorn) or newer.

`betse` currently recommends but does *not* require:

* **At least 16GB RAM**. Again, this is due to the memory intensiveness of
  even small-scale tissue simulations.

## Program Dependencies

`betse` has both mandatory dependencies that must be installed *before*
installing `betse` and optional dependencies that may be installed *after*
installing `betse`. Let's begin!

### Mandatory

`betse` requires the following non-pure-Python packages – which themselves
require non-Python libraries (e.g., C, Fortran) and hence are best installed
manually via the package manager specific to your current operating system:

* Python >= 3.3.
* Matplotlib >= 1.3.0.
* NumPy >= 1.8.0.
* PySide >= 1.1.0.
* PyYaml >= 3.10.
* SciPy >= 0.12.0.
* Voluptuous >= 0.8.7.
* setuptools >= 7.0.

#### Linux Debian

Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu), such
dependencies are installable in a system-wide manner as follows:

    >>> sudo apt-get install python3-dev python3-matplotlib python3-numpy python3-pyside python3-scipy python3-setuptools python3-yaml

#### Apple OS X

Under Apple OS X, such dependencies are installable in a system-wide manner as
follows:

. Register as an [Apple Developer](https://developer.apple.com). While free,
  such registration requires an existing Apple ID.
. Download [XCode](https://developer.apple.com/xcode). While free, such
  download requires an [Apple Developer] login.
. Install XCode, ensuring the "UNIX Development Support" checkbox is checked.
. Download and install [MacPorts](https://www.macports.org).
. Open a terminal window (e.g., by running the pre-bundled
  `Applications/Utilities`Terminal.app` application).
. Install dependencies:
    >>> sudo port install py34-matplotlib py34-numpy py34-pyside py34-scipy py34-setuptools py34-yaml
. Activate the version of Python required by betse:
    >>> sudo port select --set python python34
. Close the terminal, if you like:
    >>> exit
. Manually add the
  "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/bin" directory
  to the current ${PATH}. To do so permanently, edit the ".profile" file in your
  home directory and change the line beginning with "export PATH" to resemble
  the following:
    export PATH="/opt/local/bin:/opt/local/sbin:/opt/local/Library/Frameworks/Python.framework/Versions/3.4/bin:$PATH"

Note that MacPorts is a source-based package manager and hence extremely slow.
Expect the installation of dependencies to take several hours to several days.
(We're not kidding.)

### Optional

`betse` optionally benefits from the following mostly pure-Python packages –
which Python-specific package managers (e.g., `pip`, `setuptools`) install for
you and hence require no manual installation:

* PyInstaller `python3` branch, for optionally freezing `betse`. (See below.)
* nose >= 1.3.0, for optionally running unit tests. (See below.)

Assuming the common `git` and `pip3` commands to already be installed, such
dependencies are installable in a system-wide manner as follows:

    >>> sudo pip3 install nose
    >>> git clone --branch python3 https://github.com/pyinstaller/pyinstaller.git
    >>> cd pyinstaller
    >>> sudo python3 setup.py install

## Program Installation

`betse` is installable into either:

* A system-wide directory accessible to all users of such system.
* A venv (i.e., virtual environment) isolated to the current user.

The latter has the advantage of avoiding conflicts with already installed
system-wide Python and non-Python packages (e.g., in the event that `betse`
requires different versions of such packages), but the corresponding
disadvantage of requiring reinstallation of such packages and all transitive
dependencies of such packages. Since several dependencies are heavy-weight
(e.g., Qt4) and hence costly to reinstall, this is a notable disadvantage.

Note that the string `${BETSE\_DIR}` should be replaced everywhere below by the
absolute path of the directory containing this file.

### System-wide

`betse` is installable into a system-wide directory as follows:

* Compile `betse`.
    >>> cd "${BETSE_DIR}"
    >>> python3 setup.py build
* Install `betse`.
    >>> sudo python3 setup.py easy_install --no-deps .

Curiously, although the `develop` command for `setuptools` provides a
`--no-deps` option, the `install` command does not. Hence, the `easy\_install`
command is called above instead.

### User-specific

`betse` is installable into a user-specific venv by running the following
command **from within such venv**:

    >>> cd "${BETSE_DIR}"
    >>> ./setup.py install

This command should *not* be run outside of a venv. Doing so will reinstall all
dependencies of `betse` already installed by the system-wide package manager
(e.g., `apt-get`). This may superficially appear to work but invites obscure and
difficult to debug conflicts at `betse` runtime between dependencies reinstalled
by `setuptools` and dependencies already installed by such package maneger.

## Program Usage

`betse` is a front-facing application rather than backend framework. While
`betse`'s Python packages are importable by other packages, `betse` is typically
run by executing Python wrapper scripts installed to the current `${PATH}`.

### CLI

`betse` installs a low-level command line interface (CLI), runnable as follows:

    >>> betse

Such CLI is also runnable by invoking the main CLI module under the active
Python 3 interpreter as follows:

    >>> cd "${BETSE_DIR}"
    >>> python3 -m betse.cli

This has the minor advantage of working regardless of whether `betse` has been
installed or not, but the corresponding disadvantage of requiring more typing.

### GUI

`betse` also installs a high-level graphical user interface (GUI) implemented in
the popular cross-platform windowing toolkit Qt4, runnable as follows:

    >>> betse-qt &

Such GUI is also runnable by invoking the main GUI module under the active
Python 3 interpreter as follows:

    >>> cd "${BETSE_DIR}"
    >>> python3 -m betse.gui

## Program Development

For development purposes, `betse` is *editably installable* (i.e., as a symbolic
link rather than physical copy). As the name implies, editable installations are
modifiable at runtime and hence suitable for development. Thanks to the magic
of symbolic links, changes to the copy of `betse` from which an editable
installation was installed will be silently prograpagated back to such
installation.

### System-wide

`betse` is installable into a system-wide directory as follows:

* **(Optional).** Set the current umask to "002" as above.
    >>> umask 002
* Editably install `betse`.
    >>> cd "${BETSE_DIR}"
    >>> sudo python3 setup.py symlink

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

    >>> cd "${BETSE_DIR}"
    >>> pip3 install --no-deps --editable .

Such installation is uninstallable as follows:

    >>> pip3 uninstall betse

#### setuptools

`betse` is editably installable into a user-specific venv via `setuptools` as
follows:

    >>> cd "${BETSE_DIR}"
    >>> ./setup.py develop --no-deps

Such installation is uninstallable as follows:

    >>> cd "${BETSE_DIR}"
    >>> ./setup.py develop --uninstall

## Program Testing

`betse` is testable via `nose` as follows:

    >>> cd "${BETSE_DIR}"
    >>> ./setup.py test

If `nose` is already installed under the active Python 3 interpreter, improved
unit test output is available by either of the following two (effectively)
equivalent commands:

    >>> nosetests             # this works...
    >>> ./setup.py nosetest   # ...as does this.

## Program Freezing

`betse` is **freezable** (i.e., convertable to platform-specific executable
binaries distributable to end users) via the optional dependency PyInstaller in
one of two ways:

* One-file mode, in which case PyInstaller will generate one executable file for
  each of `betse`'s wrapper scripts (e.g., `betse`, `betse-qt`).
* One-directory mode, in which case PyInstaller will generate one directory
  containing one executable file (along with various ancillary libraries and
  archives) for each of `betse`'s wrapper scripts.

`betse` is freezable in one-file mode as follows:

    >>> cd "${BETSE_DIR}"
    >>> ./setup.py freeze_file

`betse` is freezable in one-directory mode as follows:

    >>> cd "${BETSE_DIR}"
    >>> ./setup.py freeze_dir

### Output

Assuming one-file mode, executables will be output as:

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

