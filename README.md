betse
===========

`betse` (BioElectric Tissue Simulation Engine) simulates propagation of
electrical phenomena within biological tissue (e.g., ion channel-gated current
flow).

## Dependencies

`betse` requires the following non-pure-Python packages – which themselves
require non-Python libraries (e.g., C, Fortran) and hence are best installed
manually via the package manager specific to your current operating system:

* setuptools >= 7.0.
* NumPy >= 1.9.0.
* PySide >= 1.1.0.
* PyYaml >= 3.10.
* SciPy >= 0.12.0.

`betse` also requires the following (mostly) pure-Python packages – which `pip`
will install for you and hence require no manual intervention:

* Matplotlib >= 1.3.0.
* PyYaml >= 3.11.

## Installation

`betse` is installable via `setuptools` as follows:

    >>> ./setup.py install

`betse`'s pure-Python source tree will be installed into the `site-packages`
subdirectory of the the current Python 3 interpreter. Likewise, the following
Python scripts will be installed into a platform-specific directory accessible
in the current `${PATH}`:

* `betse`, a low-level command line interface (CLI).
* `betse-qt`, a high-level graphical user interface (GUI) implemented in the
  cross-platform windowing toolkit Qt4.

## Usage

`betse` is a front-facing application rather than backend framework. While
`betse`'s Python packages are importable by other packages, `betse` is typically
directly run via the aforementioned Python scripts -- either interactively from
the command line or desktop environment or non-interactively from shell scripts.

### CLI

`betse`'s CLI is runnable as follows (assuming either `./setup.py install` or
`./setup.py develop` have been run):

    betse

`betse`'s CLI is also runnable by invoking the main CLI module under the active
Python 3 interpreter. This has the minor benefit of requiring neither
`./setup.py install` nor `./setup.py develop` to have been run, but the
corresponding disadvantage of requiring somewhat more typing. Run the following
pair of commands, replace `${BETSE_PATH}` by `betse`'s top-level directory
(e.g., containing the `.git/` subdirectory):

    cd "${BETSE_PATH}"
    python3 -m betse

### GUI

`betse`'s GUI is runnable as follows (assuming either `./setup.py install` or
`./setup.py develop` have been run):

    betse-qt

## Development

For development purposes, `betse` is installable in symbolically linked fashion
via `setuptools` as follows:

    >>> ./setup.py develop

This also installs the above two executables in the expected manner.

## Testing

`betse` is testable via `nose` as follows:

    >>> ./setup.py test

If `nose` is already installed under the active Python 3 interpreter, improved
unit test output is available by either of the following two (effectively)
equivalent commands:

    >>> nosetests             # this works...
    >>> ./setup.py nosetest   # ...as does this.

## License

