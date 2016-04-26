Unit Tests
===========

[`py.test`](http://pytest.org)-driven [unit
tests](https://en.wikipedia.org/wiki/Functional_testing) exercising BETSE's
programmatic Python API.

## Structure

For collective sanity, tests are rigorously structured as follows:

### Unit Tests

Each unit test should be:

* Defined as a function named `test_unit_{test_name}` (e.g.,
  `test_unit_mpl_import`, a unit test exercising the importability of BETSE's
  Matplotlib API), where `{test_name}` is this test's unique name.
* In a module with basename `test_{tests_name}.py` (e.g., `test_mpl.py`, a suite
  of unit tests exercising BETSE's Matplotlib API), where `{tests_name}` is any
  arbitrary non-empty string.

### Unit Test Fixtures

Each **unit test fixture** (i.e., callable depended upon by one or more unit
tests, typically preparing the external environment or filesystem for subsequent
test execution) should be:

* Defined as a function named `betse_{fixture_name}` (e.g., `betse_sim_context`,
  a fixture establishing a simulation configuration for subsequent reuse by
  several functional tests), where `{fixture_name}` is any arbitrary non-empty
  string. The `betse_` prefix future-proofs testing by preventing collisions
  between BETSE-specific fixtures and existing or future fixtures provided by
  either `py.test` itself or a third-party `py.test` plugin.
* In a module with any basename residing in the `fixture` subpackage.
* Imported by this package's `conftest` module, effectively
  implicitly this fixture in all unit test modules. While this fixture is
  explicitly importable in these modules, there are no tangible benefits to
  doing so and several tangible detriments -- including:
  * **Refactorability.** Since `py.test` silently adds this directory to
    Python's `sys.path` when testing, refactoring tools and IDEs have no means
    of automatically parsing and thus refactoring test and fixture imports.
