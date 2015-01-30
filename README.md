betse
===========

`betse` (Bioelectric Tissue Simulation Environment) simulates propagation of
electrical phenomena within biological tissue (e.g., ion channel-gated current
flow).

## Dependencies

`betse` requires the following non-pure-Python packages – which themselves
require non-Python libraries (e.g., C, Fortran) and hence are best installed
manually via the package manager specific to your current operating system:

* NumPy >= 1.9.0.
* PySide >= 1.1.0.
* PyYaml >= 3.10.
* SciPy >= 0.12.0.

`betse` also requires the following pure-Python packages – which `pip` will
install for you and hence require no manual installation:

* Matplotlib >= 1.3.0.

## Installation

`betse` is installable via `setuptools` as follows:

    >>> ./setup.py install

`betse`'s pure-Python source tree will be installed into the `site-packages/`
subdirectory of the the active Python 3 interpreter. Likewise, the following two
executables will be installed into a platform-specific directory accessible from
the current `${PATH}`:

* `betse`, a low-level command line interface (CLI).
* `betse-qt`, a high-level graphical user interface (GUI) implemented in the
  cross-platform windowing toolkit Qt4.

## Usage

`betse`'s CLI is intended to be run from both interactive shells and non-
interactive shell scripts as follows, where `$` signifies an ignorable shell
prompt (for readability):

    # Assuming either "./setup.py install" or "./setup.py develop" have been
    # run, the following command suffices.
    $ betse

    # Regardless of whether such "setup.py" tasks have been run, the following
    # pair of commands also suffices, where ${BETSE_PATH} should be replaced by
    # the top-level betse directory (e.g., containing the ".git/" subdirectory).
    $ cd "${BETSE_PATH}"
    $ python -m betse

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

