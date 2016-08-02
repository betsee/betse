Freezing
===========

BETSE is freezable (i.e., convertable to platform-specific executable binaries
distributable to end users with _no_ installation requirements) via
[**PyInstaller**](http://www.pyinstaller.org), an open-source cross-platform
utility performing such freezing. While other such utilities exist (e.g., the
Windows-specific `py2exe`, the OS X-specific `py2app`), only PyInstaller
satisfies the basic requirement of both platform portability and freedom (as in
both beer and speech).

BETSE offers setuptools-automated support for _only_ PyInstaller.

## Dependencies

BETSE's setuptools-automated freezing support has additional dependencies
beyond those required for core BETSE installation and usage. In no particular
order, these are:

### BETSE Installation

BETSE _must_ be installed in a system-wide manner first -- either:

* Editably (e.g., via `sudo python3 setup.py symlink`).
* Non-editably (e.g., via `sudo python3 setup.py install`).

See [**Development**](DEVELOP.md) for further details.

### UPX

The [Ultimate Packer for Executables (UPX)](http://upx.sourceforge.net) should
ideally also be installed, though this is _not_ strictly necessary. UPX is a
cross-platform utility for compressing executables. When installed, PyInstaller
automatically compresses all written executables.

### Windows

UPX does _not_ appear to be readily available under Microsoft Windows. Next!

### Linux

Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu), UPX is
installable in a system-wide manner as follows:

    $ sudo apt-get install upx-ucl

## Usage

BETSE is freezable in either:

* One-file mode, in which case PyInstaller generates one executable file for
  each of BETSE's wrapper scripts (e.g., `betse`, `betse-qt`).
* One-directory mode, in which case PyInstaller generates one directory
  containing one executable file (along with various ancillary libraries and
  archives) for each of BETSE's wrapper scripts.

### One-file Mode

BETSE is freezable in one-file mode as follows. Either:

* Run `freeze_file`, the provided Bash shell script wrapper. For convenience,
  this script is runnable from any directory (including the top-level BETSE
  directory) _and_ accepts all arguments accepted by the `freeze_file`
  subcommand:

        $ ./freeze_file

* Run the `setuptools`-driven `freeze_file` subcommand. Due to `setuptools`
  constraints, this subcommand is runnable _only_ from the top-level BETSE
  directory:

        $ cd "${BETSE_DIR}"
        $ python3 setup.py freeze_file

### One-directory Mode

BETSE is freezable in one-directory mode as follows. Either:

* Run `freeze_dir`, the provided Bash shell script wrapper. For convenience,
  this script is runnable from any directory (including the top-level BETSE
  directory) _and_ accepts all arguments accepted by the `freeze_dir`
  subcommand:

        $ ./freeze_dir

* Run the `setuptools`-driven `freeze_dir` subcommand. Due to `setuptools`
  constraints, this subcommand is runnable _only_ from the top-level BETSE
  directory:

        $ cd "${BETSE_DIR}"
        $ python3 setup.py freeze_dir

## Paths

Assuming one-file mode, executables will be frozen to the following paths:

* Under Linux:
  * `freeze/dist/betse` for BETSE's CLI wrapper.
  * `freeze/dist/betse-qt` for BETSE's GUI wrapper.
* Under OS X:
  * `freeze/dist/betse` for BETSE's CLI wrapper.
  * `freeze/dist/betse-qt.app` for BETSE's GUI wrapper.
* Under Windows:
  * `freeze/dist/betse.exe` for BETSE's CLI wrapper.
  * `freeze/dist/betse-qt.exe` for BETSE's GUI wrapper.

These executables may be moved (and hence distributed to end users) as is but
should _not_ be renamed. Doing so typically invalidates code signing (as well as
other embedded metadata).

## Caveats

While commonly regarded as the best utility for freezing Python applications,
PyInstaller is _not_ without the occasional caveat. Most if not all such caveats
apply to alternative utilities (e.g., `cx_Freeze`, `py2app`, `py2exe`).

### Switching Freeze Modes

For safety, attempting to switch freeze modes (i.e., to refreeze BETSE in
one-directory mode when already frozen in one-file mode _or_ in one-file mode
when already frozen in one-directory mode) currently results in an error.

This is correctable by either manually removing the previously frozen file or
directory _or_ passing option `--clean` to the desired `freeze_file` or
`freeze_dir` command, which automatically removes such path on your behalf. For
example, BETSE is switchable from one-directory to one-file mode and then back
again as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py freeze_dir
    $ ./setup.py freeze_file --clean
    $ ./setup.py freeze_dir --clean

### Forwards Incompatibilities

Executables frozen on older versions of supported operating systems are
_typically_ compatible with newer versions of the same systems. For example,
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
dependencies. The simplest means of doing so is to host (i.e., install and
run) such versions of such systems as virtual guests of the current system via
either:

* A full-fledged [hypervisor](https://en.wikipedia.org/wiki/Hypervisor) (i.e.,
  an application creating and running virtual machines).
* A minimalistic [Docker](https://www.docker.com) container.

To that end, we recommend:

* Under Linux, [**Docker**](https://www.docker.com).
* Under all other platforms, Oracle's open-source
  [VirtualBox](https://www.virtualbox.org) hypervisor rather than VMware's
  closed-source [VMware Workstation](https://www.vmware.com) hypervisor. The
  former supports all operating systems supported by BETSE out-of-the-box,
  including OS X; the latter supports Linux and Windows but _not_ OS X
  out-of-the-box. Moreover, the former is free (as in both beer and speech); the
  latter is non-free (as in both beer and speech).

### No Cross-freezing or Cross-compiling

PyInstaller supports neither **cross-freezing** (i.e., generation of
executables intended for execution on operating systems other than the current)
nor **cross-compilation** (i.e., generating of executables intended for
execution on CPU architectures other than the current) out of the box. Hence,
executables will be executable only under systems matching the current operating
system and architecture. As example, if the current system is 64-bit Linux,
executables generated under such system will be executable only under other
64-bit Linux systems.

In practice, the lack of cross-compilation support is ignorable. BETSE is
sufficiently memory-intensive that running under a 32-bit architecture would be
strongly inadvisable and hence unsupported.

The lack of cross-freezing support is _not_ ignorable. While there exist well-
documented means of cross-freezing Windows executables on Linux systems (e.g.,
via the popular emulation layer Wine), there exist _no_ means of cross-freezing
OS X executables on Linux. The conventional solution is to install OS X as a
virtual guest of an existing Linux host and perform freezing under such guest.
We recommend the open-source product VirtualBox for this purpose.
