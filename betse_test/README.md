Tests
===========

[`py.test`](http://pytest.org)-driven
[functional](https://en.wikipedia.org/wiki/Functional_testing) and [unit
tests](https://en.wikipedia.org/wiki/Unit_testing).

## Structure

For collective sanity, tests are rigorously structured as follows:

* All functional test reside in the `betse_test_func` subdirectory of this
  directory.
* All unit tests reside in the `betse_test_unit` subdirectory of this directory.
* This directory should directly contain _no_ files except this otherwise
  ignorable file. In particular, this directory should _not_ contain an
  `__init__.py` file.
* This directory should directly contain _no_ subdirectories except the
  aforementioned subdirectories.

### Nomenclature

The somewhat obscure choice of subdirectory names is intentional, ensuring:

* Functional tests may import other functional tests or fixtures via the topmost 
  `betse_test_func` package name, which `py.test` dynamically injects at test
  discovery time into Python's current `sys.path`.
* Unit tests may import other unit tests or fixtures via the topmost 
  `betse_test_unit` package name, which `py.test` dynamically injects at test
  discovery time into Python's current `sys.path`.
* No import collision between functional and unit tests and fixtures. Functional
  tests and fixtures reside in one Python namespace; unit tests and fixtures
  reside in another.
* No import collision between tests and fixtures and the remainder of the Python
  ecosystem (e.g., the first-party stdlib or third-party packages). Tests and
  fixtures reside in `betse_test_` namespaces. All other packages presumably do
  _not_.

### Importability

To permit any test or fixture to import any other test or fixture, this
directory and all subdirectories of this directory _must_ contain `__init__.py`
files. As no parent directories of this directory contain these files, `py.test`
[assigns](https://pytest.org/latest/goodpractices.html) this directory a topmost
package with the same basename in Python's import namespace. All submodules
(i.e., `.py`-suffixed files in subpackages) and subpackages (i.e., subdirectories
with `__init__.py` files) of this package may then be imported from any test or
fixture via this package name.

For example, any test or fixture may explicitly import both the test-specific
`betse_test/util/metafixture.py` module _and_ the test-agnostic
`betse/metadata.py` module defined in the main codebase as follows:

```
# This imports "betse_test/util/metafixture.py", thanks to magic.
from betse_test.util import metafixture

# This imports "betse/metadata.py", thanks to magic.
from betse import metadata
```

## See Also

For further details on test discovery, see the
[paragraph](https://pytest.org/latest/goodpractices.html) beginning _"If pytest
finds a “a/b/test_module.py” test file..."_.
