[![build status](https://gitlab.com/betse/betse/badges/master/build.svg)](https://gitlab.com/betse/betse/commits/master)

BETSE
===========

**BETSE** (**B**io **E**lectric **T**issue **S**imulation **E**ngine) is a
cross-platform, open-source [finite
volume](https://en.wikipedia.org/wiki/Finite_volume_method) simulator for 2D
computational multiphysics problems in the life sciences, including
[electrodiffusion](https://en.wikipedia.org/wiki/Nernst%E2%80%93Planck_equation),
[electro-osmosis](https://en.wikipedia.org/wiki/Electro-osmosis),
[galvanotaxis](https://en.wiktionary.org/wiki/galvanotaxis), [voltage-gated ion
channels](https://en.wikipedia.org/wiki/Voltage-gated_ion_channel), [gene
regulatory networks](https://en.wikipedia.org/wiki/Gene_regulatory_network),
and [biochemical reaction
networks](http://www.nature.com/subjects/biochemical-reaction-networks) (e.g.
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
line interface (CLI) supporting Linux and OS X.
<sup>_Windows is currently unsupported._</sup>

## Installation

BETSE is installable under **Linux** and **OS X** as follows:

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

For a comprehensive list of system requirements, software dependencies, and
platform-specific installation instructions, see [**BETSE
Installation**](doc/md/INSTALL.md).

## Description

BETSE simulates biorealistic electrochemical phenomena in gap
junction-networked, 2D cellular collectives. To predict bioelectric patterns
and their spatio-temporal dynamics, BETSE:

* Models [ion channel](https://en.wikipedia.org/wiki/Ion_channel) and
  [gap junction](https://en.wikipedia.org/wiki/Gap_junction) activity.
* Tracks changes in ion concentration and net ionic charge.
* Calculates endogenous voltages and currents.
* Imports bitmask images defining the shapes of:
  * Cell clusters.
  * Cell cluster regions localizing ion channel activity, typically signifying
    different types of adjacent tissues (e.g., organs).
* Exports simulation results to a variety of output formats, including:
  * Publication-quality plots and graphs.
  * Internet-friendly compressed video.
  * Post-processable tabular data (e.g., [comma-separated values
    (CSV)](https://en.wikipedia.org/wiki/Comma-separated_values)).

Simulations can optionally include the seven ions underpinning bioelectrical
signals: [Na<sup>+</sup>](https://en.wikipedia.org/wiki/Sodium_in_biology),
[K<sup>+</sup>](https://en.wikipedia.org/wiki/Potassium_in_biology),
[Cl<sup>-</sup>](https://en.wikipedia.org/wiki/Chloride),
[Ca<sup>2+</sup>](https://en.wikipedia.org/wiki/Calcium_in_biology),
[H<sup>+</sup>](https://en.wikipedia.org/wiki/Hydron_(chemistry)),
[HCO<sub>3</sub><sup>-</sup>](https://en.wikipedia.org/wiki/Bicarbonate_transporter_protein),
and [anionic proteins
(P<sup>-</sup>)](https://en.wikipedia.org/wiki/Gibbs%E2%80%93Donnan_effect).

Individual cells can include a variety of voltage-gated ion channels, which are
based on [Hodgkin-Huxley
formalism](https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model) with
experimentally-derived variables sourced from
[Channelpedia](http://channelpedia.epfl.ch). Currently implemented
voltage-gated channel types include Nav1.2, Nav1.3, Nav1.6, Kv1.2, Kv1.2, Kv.1.5,
Kv3.3, Kv3.4, Kir2.1, L-type Ca, T-type Ca, P/Q-type Ca, HCN1, HCN2, and HCN4.
A variety of built-in and custom leak and ligand gated channels are available,
including calcium-gated K+ channels. Further control of bioelectric cell
dynamics is enabled via a range of ion pumps and exchangers, including
Na/K-ATPase, H/K-ATPase, V-ATPase, and Ca-ATPase.

Cells form an interconnected intracellular network via voltage-sensitive
gap junction connections, and a unique extracellular environment can be
maintained by tight junctions at the cluster periphery. Simulation of the
extracellular environment allows for exploration of local field potentials, the
transepithelial potential, and ephaptic coupling between cells.

BETSE also has biosystems modeling capabilities to implement custom
gene regulatory networks, biochemical reaction networks (with an emphasis on
metabolism) and integrate the interaction (i.e. activity modulation) between
gene-products and other biochemicals with ion channels, pumps and gap junctions,
thereby allowing for the study of relationships between the powerful control
systems of both gene/biochemical regulatory networks and bioelectrical signals.

## Validation

Simulation performance is perpetually validated by matching experimentally
observed data on membrane permeability, ion concentration, resting potential,
and related biophysical quantities to simulation output. Expected outcomes have
been demonstrated for a range of well-known cases, including:

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

BETSE is open-source software licensed under the permissive [BSD 2-clause
license](https://opensource.org/licenses/BSD-2-Clause). See [`LICENSE`](LICENSE)
for exhaustive details.

## Reference

If utilizing BETSE in your own work, please cite the following
[open-access introductory
paper](http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract):

> Pietak Alexis and [Levin Michael](https://ase.tufts.edu/biology/labs/levin/)
> (2016). [**Exploring Instructive Physiological
> Signaling with the Bioelectric Tissue Simulation Engine
> (BETSE)**](http://journal.frontiersin.org/article/10.3389/fbioe.2016.00055/abstract).
> [_Frontiers in Bioengineering and
> Biotechnology_](http://journal.frontiersin.org/journal/bioengineering-and-biotechnology)
> 4, 55. `doi:10.3389/fbioe.2016.00055`

## Usage

BETSE is a front-facing application rather than backend framework. While
BETSE's Python packages are importable by other packages, BETSE is typically
run by executing Python wrapper scripts installed to the current `${PATH}`.

### CLI

BETSE installs a low-level command line interface (CLI), runnable as follows:

    $ betse

This CLI is also runnable by invoking the main CLI module under the active
Python 3 interpreter as follows:

    $ cd "${BETSE_DIR}"
    $ python3 -m betse.cli

This has the minor advantage of working regardless of whether BETSE has been
installed or not, but the corresponding disadvantage of requiring more typing.

### GUI

BETSE also installs a high-level graphical user interface (GUI) implemented in
the popular cross-platform windowing toolkit Qt4, runnable as follows:

    $ betse-qt &

This GUI is also runnable by invoking the main GUI module under the active
Python 3 interpreter as follows:

    $ cd "${BETSE_DIR}"
    $ python3 -m betse.gui

## Development

For development purposes, BETSE is *editably installable* (i.e., as a symbolic
link rather than physical copy). As the name implies, editable installations are
modifiable at runtime and hence suitable for development. Thanks to the magic
of symbolic links, changes to the copy of BETSE from which an editable
installation was installed will be silently prograpagated back to such
installation.

### System-wide

BETSE is installable into a system-wide directory as follows:

* **(Optional).** Set the current umask to "002" as above.
    $ umask 002
* Editably install BETSE.
    $ cd "${BETSE_DIR}"
    $ sudo python3 setup.py symlink

The `symlink` command is a BETSE-specific `setuptools` command inspired by the
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

BETSE is editably installable into a user-specific venv via either `pip` or
`setuptools` **from within such venv.** While there appears to be no particular
advantage to using one over the other, it remains helpful to note that both
apply. In either case, external executables (e.g., `betse`, `betse-qt`) will
also be installed and usable in the expected manner.

#### pip

BETSE is editably installable into a user-specific venv via `pip` as follows:

    $ cd "${BETSE_DIR}"
    $ pip3 install --no-deps --editable .

This installation is uninstallable as follows:

    $ pip3 uninstall betse

#### setuptools

BETSE is editably installable into a user-specific venv via `setuptools` as
follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py develop --no-deps

This installation is uninstallable as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py develop --uninstall

## Testing

BETSE is testable via `py.test` as follows. Either:

* Run all available tests. Either:
  * Run `test`, the provided Bash shell script wrapper. For convenience, this
    script is runnable from any directory (including this directory) _and_
    accepts all arguments accepted by the `test` subcommand (detailed below):

            $ ./test

  * Run the `setuptools`-driven `test` subcommand. Due to `setuptools`
    constraints, this subcommand is runnable _only_ from this directory:

            $ cd "${BETSE_DIR}"
            $ python3 setup.py test

* Run all tests matching a passed Python-evaluatable expression. For example, to
  run all test functions and classes whose names contain either `test_tartarus`
  _or_ `test_thalassa`:

        $ cd "${BETSE_DIR}"
        $ python3 setup.py test -k 'test_tartarus or test_thalassa'

## Freezing

BETSE is **freezable** (i.e., convertable to platform-specific executable
binaries distributable to end users) via PyInstaller, a cross-platform open-
source utility for performing such freezing. While other such utilities exist
(e.g., the Windows-specific `py2exe`, the OS X-specific `py2app`), only
PyInstaller satisfies the basic requirements of both platform portability and
freedom (as in both beer and speech). BETSE offers setuptools-automated
support for *only* PyInstaller.

### Dependencies

BETSE's setuptools-automated freezing support has additional dependencies
beyond those required for core BETSE installation and usage. In no particular
order, these are:

#### BETSE Installation

BETSE *must* be installed in a system-wide manner first -- either:

* Editably (e.g., via `sudo python3 setup.py symlink`).
* Non-editably (e.g., via `sudo python3 setup.py install`).

See above for further details.

#### PyInstaller

**NOTE: _This section is largely obsolete, and should either be significantly
revised or completely deleted._**

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
though this is _not_ strictly necessary. UPX is a cross-platform utility for
compressing executables. When installed, PyInstaller automatically compresses
all written executables.

##### Windows

UPX does _not_ appear to be readily available under Microsoft Windows. Next!

##### Linux

Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu), UPX is
installable in a system-wide manner as follows:

    $ sudo apt-get install upx-ucl

### Usage

BETSE is freezable in either:

* One-file mode, in which case PyInstaller generates one executable file for
  each of BETSE's wrapper scripts (e.g., `betse`, `betse-qt`).
* One-directory mode, in which case PyInstaller generates one directory
  containing one executable file (along with various ancillary libraries and
  archives) for each of BETSE's wrapper scripts.

#### One-file Mode

BETSE is freezable in one-file mode as follows. Either:

* Run `freeze_file`, the provided Bash shell script wrapper. For convenience,
  this script is runnable from any directory (including this directory) _and_
  accepts all arguments accepted by the `freeze_file` subcommand:

        $ ./freeze_file

* Run the `setuptools`-driven `freeze_file` subcommand. Due to `setuptools`
  constraints, this subcommand is runnable _only_ from this directory:

        $ cd "${BETSE_DIR}"
        $ python3 setup.py freeze_file

#### One-directory Mode

BETSE is freezable in one-directory mode as follows. Either:

* Run `freeze_dir`, the provided Bash shell script wrapper. For convenience,
  this script is runnable from any directory (including this directory) _and_
  accepts all arguments accepted by the `freeze_dir` subcommand:

        $ ./freeze_dir

* Run the `setuptools`-driven `freeze_dir` subcommand. Due to `setuptools`
  constraints, this subcommand is runnable _only_ from this directory:

        $ cd "${BETSE_DIR}"
        $ python3 setup.py freeze_dir

### Paths

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

### Caveats

While commonly regarded as the best utility for freezing Python applications,
PyInstaller is _not_ without the occasional caveat. Most if not all such caveats
apply to alternative utilities (e.g., `cx_Freeze`, `py2app`, `py2exe`).

#### Switching Freeze Modes

For safety, attempting to switch freeze modes (i.e., to refreeze BETSE in
one-directory mode when already frozen in one-file mode *or* in one-file mode
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

#### Forwards Incompatibilities

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
dependencies. The simplest means of doing so is to **host** (i.e., install and
run) such versions of such systems as virtual guests of the current system via a
**hypervisor** (i.e., an application creating and running virtual machines).

For this purpose, we recommend Oracle's open-source hypervisor VirtualBox rather
than VMware's closed-source hypervisor VMware Workstation. The former supports
all operating systems supported by BETSE out of the box, including OS X; the
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

In practice, the lack of cross-compilation support is ignorable. BETSE is
sufficiently memory-intensive that running under a 32-bit architecture would be
strongly inadvisable and hence unsupported.

The lack of cross-freezing support is *not* ignorable. While there exist well-
documented means of cross-freezing Windows executables on Linux systems (e.g.,
via the popular emulation layer Wine), there exist *no* means of cross-freezing
OS X executables on Linux. The conventional solution is to install OS X as a
virtual guest of an existing Linux host and perform freezing under such guest.
We recommend the open-source product VMware for this purpose.
