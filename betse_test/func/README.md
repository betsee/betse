Functional Tests
===========

[`py.test`](http://pytest.org)-driven [functional
tests](https://en.wikipedia.org/wiki/Functional_testing) exercising the
currently installed versions of BETSE's external CLI and GUI commands.

## Structure

For collective sanity, tests are rigorously structured as follows:

### Functional Tests

Each functional test should be:

* Defined as a function named `test_{ui_type}_{test_name}` (e.g.,
  `test_cli_betse_info`, a functional test exercising the `betse info`
  subcommand of BETSE's CLI), where:
  * `{ui_type}` is either:
    * `cli`, for functional tests exercising the BETSE CLI.
    * `gui`, for functional tests exercising the BETSE GUI.
  * `{test_name}` is this test's unique name.
* In a module with basename `test_{tests_name}.py` (e.g., `test_cli.py`, a suite
  of functional tests exercising BETSE's CLI), where `{tests_name}` is any
  arbitrary non-empty string, residing in either the `cli` or `gui` subpackages.

### Functional Test Fixtures

Each **functional test fixture** (i.e., callable depended upon by one or more
functional tests, typically preparing the external environment or filesystem for
subsequent test execution) should be:

* Defined as a function named `betse_{fixture_name}` (e.g., `betse_sim_context`,
  a fixture establishing a simulation configuration for subsequent reuse by
  several functional tests), where `{fixture_name}` is any arbitrary non-empty
  string. The `betse_` prefix future-proofs testing by preventing collisions
  between BETSE-specific fixtures and existing or future fixtures provided by
  either `py.test` itself or a third-party `py.test` plugin.
* In a module with any basename residing in the `fixture` subpackage.
* Imported by the topmost `conftest` module of this package, effectively
  implicitly this fixture in all functional test modules. While this fixture is
  explicitly importable in these modules, there are no tangible benefits to
  doing so and several tangible detriments -- including:
  * **Refactorability.** Since `py.test` silently adds this directory to
    Python's `sys.path` when testing, refactoring tools and IDEs have no means
    of automatically parsing and thus refactoring test and fixture imports.
