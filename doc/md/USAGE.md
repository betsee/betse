Usage
===========

> **NOTE:** This synopsis is _painfully_ inadequate. For detailed usage
> instructions complete with explanatory examples, plots, and screenshots, see
> our [77-page
> PDF](https://www.dropbox.com/s/fsxhjpipbiog0ru/BETSE_Documentation_Nov1st2015.pdf?dl=0)
> instead.

BETSE is usable as follows.

## CLI

BETSE is currently _only_ available as a low-level command-line interface (CLI)
named `betse`, a Python wrapper script installed to the current `${PATH}` on
[BETSE installation](INSTALL.md).

## GUI

A high-level graphical user interface (GUI) named `betse-qt` implemented via
non-[GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License) Python
bindings (e.g., [PySide](https://wiki.qt.io/PySide),
[PySide2](https://wiki.qt.io/PySide2)) to the cross-platform
non-[GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License) [Qt]
(https://www.qt.io) windowing toolkit is planned **but currently
unimplemented.**

Consider submitting a [feature request](https://gitlab.com/betse/betse/issues)
if you would prefer to see this process prioritized.

## Python

BETSE is intended to be used as a front-facing interactive application rather
than as a backend non-interactive library.

The **BETSE API** comprises the top-level
[`betse` Python package](https://gitlab.com/betse/betse/tree/master/betse) and
all subpackages and submodules of that package installed with BETSE. While the
BETSE API _is_ externally importable by other Python packages, doing so is
currently unsupported. The BETSE API is intended to be imported and used _only_
by the first-party interfaces listed above.

**The BETSE API is _not_ intended to be imported or used by third parties** –
not because we don't like third parties,<sup>_We do!_</sup> but because we lack
sufficient resources (both grant funding _and_ developer time) to stabilize the
BETSE API for public consumption. Until resources arrive, breaking API changes
are likely to be a permanent fixture of the BETSE landscape.

No [backward](https://en.wikipedia.org/wiki/Backward_compatibility) or
[forward compatibility](https://en.wikipedia.org/wiki/Forward_compatibility)
guarantees are currently provided. Third-party Python scripts (e.g., [IPython
notebooks](http://jupyter.org)), frameworks, and applications attempting to
leverage BETSE for internal use are likely to be disappointed. 

Consider submitting a [feature request](https://gitlab.com/betse/betse/issues)
if you would prefer to see the BETSE API stabilized.
