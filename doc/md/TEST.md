Testing
===========

BETSE is rigorously tested with a [comprehensive test
suite](https://gitlab.com/betse/betse/tree/master/betse_test) comprising both
[functional](https://en.wikipedia.org/wiki/Functional_testing) and [unit
tests](https://en.wikipedia.org/wiki/Unit_testing).

## Manual Testing

BETSE is manually testable by [**py.test**](http://pytest.org) at the command
line. Either:

* Run all available tests. Either:
  * Run `test`, the provided Bash shell script wrapper. For convenience, this
    script is runnable from any directory (including the top-level BETSE
    directory) _and_ accepts all arguments accepted by the `test` subcommand
    (detailed below):

            $ ./test

  * Run the `setuptools`-driven `test` subcommand. Due to `setuptools`
    constraints, this subcommand is runnable _only_ from the top-level BETSE
    directory:

            $ cd "${BETSE_DIR}"
            $ python3 setup.py test

* Run all tests matching a passed Python-evaluatable expression. For example, to
  run all test functions and classes whose names contain either `test_tartarus`
  _or_ `test_thalassa`:

        $ cd "${BETSE_DIR}"
        $ python3 setup.py test -k 'test_tartarus or test_thalassa'

## Continuous Integration (CI)

BETSE is continuously integrated by
[**GitLab-CI**](https://about.gitlab.com/gitlab-ci) on each commit pushed to
each branch of BETSE's [GitLab](https://gitlab.com)-hosted [Git
repository](https://gitlab.com/betse/betse/tree/master).
