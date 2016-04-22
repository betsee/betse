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
