.. # ------------------( SYNOPSIS                            )------------------

=========
Releasing
=========

Stable releases of BETSE are produced with a rigorous procedure reminiscent to
that of most open-source software: *namely,*

#. A stable commit is **tagged** with an annotated Git tag.
#. That commit is **packaged** into a compressed source tarball.
#. That tarball is **uploaded** to third-party services for public consumption.

While technically optional, this procedure reduces the likelihood of
installation and usage woes by downstream consumers (\ *e.g.,* end users,
package maintainers) and is thus effectively mandatory.

.. # ------------------( TABLE OF CONTENTS                   )------------------
.. # Blank line. By default, Docutils appears to only separate the subsequent
.. # table of contents heading from the prior paragraph by less than a single
.. # blank line, hampering this table's readability and aesthetic comeliness.

|

.. # Table of contents, excluding the above document heading. While the
.. # official reStructuredText documentation suggests that a language-specific
.. # heading will automatically prepend this table, this does *NOT* appear to
.. # be the case. Instead, this heading must be explicitly declared.

.. contents:: **Contents**
   :local:

.. # ------------------( DESCRIPTION                        )------------------

Procedure
============

BETSE is releasable to all supported platforms as follows:

#. **Install** wheel_, a third-party pure-Python package permitting this release
   to be packaged into a cross-platform pre-compiled binary distribution
   supported by both PyPI_ and `pip`. This mandatory dependency augments
   setuptools with the ``bdist_wheel`` subcommand invoked below.

   .. code:: bash

      $ sudo pip3 install wheel

#. (\ *Optional*\ ) **Install** twine_, a third-party pure-Python package
   simplifying remote communication with both `PyPI itself <PyPI_>`__ and
   `Test PyPI`_. While optional, securely registering and uploading PyPI
   distributions *without* twine_ is typically non-trivial, banal, and tedious.
   :sup:`Your mileage may vary.`
#. (\ *Optional*\ ) **Validate reStructuredText (reST) rendering.** The
   human-readable description for this release derives directly from `the
   top-level README.rst file <readme_>`__ for this project. Sadly, PyPI's reST
   renderer supports only a proper subset of the syntax supported by the reST
   standard – itself only a proper subset of the syntax supported by Sphinx. If
   `this file <readme_>`__ contains syntax unsupported by PyPI's reST renderer,
   PyPI erroneously preserves this file as plaintext rather than rendering this
   file as HTML. To avoid this:

   #. Install the ``collective.checkdocs`` Python package.

      .. code:: bash

         $ sudo pip3 install collective.checkdocs

   #. Validate the PyPI-specific compatilibility of `this file <readme_>`__.

      .. code:: bash

         $ python3 setup.py checkdocs

   #. After submitting this release to PyPI via ``twine`` below, manually browse
      to `the PyPI-hosted page <PyPI BETSE_>`__ for this project and verify by
      cursory inspection that this project's description is rendered as HTML.

#. (\ *Optional*\ ) **Bump release metadata.** Assuming the prior release
   followed these instructions, release metadata has already been bumped in
   preparation for the next (i.e., this) release. If another bump is required
   (e.g., to upgrade this release from a patch to a minor or even major update),
   this bump should be performed *before* tagging this release. For details, see
   see the eponymous *"Bump release metadata."* instructions below.
#. (\ *Optional*\ ) **List all existing tags.** For reference, listing all
   previously created tags *before* creating new tags is often advisable.

   .. code:: bash

      $ git tag

#. **Create an announcement commit.** This commit should have a message whose:

   * First line is of the format ``"BETSE {version} ({codename}) released."``,
     where:

     * ``{version}`` is the current value of the ``betse.metadata.__version__``
       global.
     * ``{codename}`` is the current value of the ``betse.metadata.CODENAME``
       global.

   * Remaining lines are a changelog synopsizing the most significant changes
     implemented by this release – ideally in the enumerated format given below.

   For example::

       BETSE 0.4.0 (Glad Galvani) released.

       Significant changes include:

       * Tissue profiles generalized.
       * Animation video encoding supported.
       * Simulation stability and efficiency improved.

#. **Tag this commit.** An annotated tag\ [#tags]_ should be created whose:

   * Name is ``v{version}``, where:

     * ``v`` is an arbitrary prefix preserving historical consistency with
       previous tag names in this repository.
     * ``{version}`` is the current value of the ``betse.metadata.__version__``
       global.

   * Message is the same commit message created above.

   .. code:: bash

      $ git tag -a v{version}

#. **Package both a source tarball and binary wheel.**

   .. code:: bash

      $ python3 setup.py sdist bdist_wheel

#. (\ *Optional*\ ) **List the contents of this source tarball,** where
   ``${version}`` is the purely numeric version of this release (e.g.,
   ``0.4.1``). Verify by inspection that no unwanted paths were packaged.

   .. code:: bash

      $ tar -tvzf dist/betse-${version}.tar.gz | less

#. (\ *Optional*\ ) **Test the local installation of this release.** If
   installation of this release differs from that of prior releases, testing
   *before* publishing this release to PyPI_ and elsewhere is advisable.

   #. **Test this source tarball locally.**

      #. **Create a new empty (venv)** (i.e., virtual environment).

         .. code:: bash

            $ python3 -m venv --clear /tmp/betse-sdist

      #. **Install this source tarball into this venv.**\ [#venv]_

         .. code:: bash

            $ /tmp/betse-sdist/bin/pip3 install dist/betse-${version}.tar.gz

      #. **Test this release from this venv.**

         .. code:: bash

            $ cd /tmp && /tmp/betse-sdist/bin/betse try

      #. **Remove this venv and return to the prior directory.**

         .. code:: bash

            $ rm -rf /tmp/betse-sdist && cd -

   #. **Test this binary wheel locally.**

      #. **Create a new empty venv.**

         .. code:: bash

            $ python3 -m venv --clear /tmp/betse-wheel

      #. **Install this binary wheel into this venv.**\ [#venv]_

         .. code:: bash

            $ /tmp/betse-wheel/bin/pip3 install \
              dist/betse-${version}-py3-none-any.whl

      #. **Test this release from this venv.**

         .. code:: bash

            $ cd /tmp && /tmp/betse-wheel/bin/betse try

      #. **Remove this venv and sample simulation and return to the prior
         directory.**

         .. code:: bash

            $ rm -rf /tmp/betse-wheel /tmp/sample_sim && cd -

#. **Bump release metadata.** In preparation for developing the next release:

   #. The ``betse.metadata.__version__`` global should be incremented according
      to the `best practices <Version Nomenclature_>`__ provided below.
   #. The ``betse.metadata.CODENAME`` global should be incremented according
      to the `best practices <Codename Nomenclature_>`__ provided below.

#. (\ *Optional*\ ) **Bump downstream metadata.** As example, if the current
   version of BETSEE_ strictly requires the current version of BETSE, the
   ``betsee.guimetadeps.BETSE_VERSION_REQUIRED_MIN`` global string variable of
   the former should be incremented to reflect the latter.

#. **Create another announcement commit.** This commit should have a message
   whose first line is of the format ``"BETSE {version} ({codename})
   started."``, where:

     * ``{version}`` is the new value of the ``betse.metadata.__version__``
       global.
     * ``{codename}`` is the new value of the ``betse.metadata.CODENAME``
       global.

   Since no changelog for this release yet exists, a single-line message
   suffices for this commit. For example::

       BETSE 0.4.1 (Gladder Galvani) started.

#. **Push this tagged commit.** After doing so, Gitlab will automatically
   publish source tarballs in various formats (e.g., ``.zip``, ``.tar.bz2``)
   containing the contents of this repository at this tagged commit in this
   project's `source tarball archive <tarballs_>`__. No further work is required
   to distribute source tarballs via Gitlab.

   .. code:: bash

      $ git push && git push --tags

#. (\ *Optional*\ ) **Test the remote installation of this release.**

   #. **Test this release on** `Test PyPI`_. Note that, as this server is a
      moving target, the `official instructions <Test PyPI instructions_>`__
      *always* supersede those listed for convenience below.

      #. **Create a** `Test PyPI user`_.
      #. **Create a** ``~/.pypirc`` **dotfile,** ideally by following the
         `official instructions <Test PyPI instructions_>`__ for doing so.
      #. **Register this project with** `Test PyPI`_.

         .. code:: bash

            $ python3 setup.py register -r testpypi

      #. **Browse to this project on** `Test PyPI`_. Verify by inspection all
         identifying metadata at the following URL:

         https://testpypi.python.org/pypi/betse

      #. **Upload this source tarball and binary wheel to**  `Test PyPI`_.

         .. code:: bash

            $ twine upload -r testpypi dist/betse-${version}*

      #. **Create a new empty venv.**

         .. code:: bash

            $ python3 -m venv --clear /tmp/betse-pypi

      #. **Install this release into this venv.**\ [#venv]_

         .. code:: bash

            $ /tmp/betse-pypi/bin/pip3 install \
              install -i https://testpypi.python.org/pypi betse

      #. **Test this release from this venv.**

         .. code:: bash

            $ cd /tmp && /tmp/betse-pypi/bin/betse try

      #. **Remove this venv and sample simulation and return to the prior
         directory.**

         .. code:: bash

            $ rm -rf /tmp/betse-pypi /tmp/sample_sim && cd -

#. **Publish this release to** `PyPI`_.

   #. **Create a** `PyPI user`_.
   #. **Validate the primary e-mail address associated with this account,**
      which `PyPI`_ requires as a hard prerequisite to performing the first
      upload (and hence creation) for this project.
   #. **Create a** ``~/.pypirc`` **dotfile,** ideally by following the
      `official instructions <Test PyPI instructions_>`__ for doing so.
   #. **Upload this source tarball and binary wheel to** `PyPI`_. If this is the
      first such upload for this project, a `PyPI`_-hosted project page will be
      implicitly created by this upload. `PyPI` neither requires, recommends,
      nor supports end user intervention in this process.

      .. code:: bash

         $ twine upload dist/betse-${version}*

   #. (\ *Optional*\ ) **Browse to this project on** `PyPI`_. Verify by
      inspection all identifying metadata at the following URL:

      https://pypi.python.org/pypi/betse

   #. (\ *Optional*\ ) **Test the installation of this release from** `PyPI`_.

      #. **Create a new empty venv.**

         .. code:: bash

            $ python3 -m venv --clear /tmp/betse-pypi

      #. **Install this release into this venv.**\ [#venv]_

         .. code:: bash

            $ /tmp/betse-pypi/bin/pip3 install betse

      #. **Test this release from this venv.**

         .. code:: bash

            $ cd /tmp && /tmp/betse-pypi/bin/betse try

      #. **Remove this venv and sample simulation and return to the prior
         directory.**

         .. code:: bash

            $ rm -rf /tmp/betse-pypi /tmp/sample_sim && cd -

#. (\ *Optional*\ ) **Update third-party packages.** As of this writing, these
   include (in no particular order):

   * Our official `Anaconda package`_, automatically produced for all supported
     platforms from the `conda recipe`_ hosted at the `conda-forge feedstock`_
     maintained by a co-maintainer of BETSE. Updating this package thus reduces
     to updating this recipe. To do so, avoid directly pushing to any branch
     (including ``master``) of the `feedstock repository`_, as doing so
     conflicts with `conda-forge`_ automation; instead (in order):

     #. Remotely create a `GitHub`_ account.
     #. Remotely login to this account.
     #. Remotely fork our `feedstock repository`_.
     #. Locally clone this forked feedstock repository.
     #. Locally create a new branch of this repository specific to this update.
     #. Locally update this recipe from this branch (typically, by editing the
        ``recipe/meta.yaml`` file). When doing so, note that:

        * The sha256 hash of the updated tarball *must* be manually embedded in
          this recipe. To obtain this hash remotely (in order):

          * Browse to `the PyPI-hosted page <PyPI BETSE_>`__ for this project.
          * Click the *Download Files* link.
          * Click the *SHA256* link to the right of the updated tarball.
          * Paste the resulting string as the value of the ``sha256`` Jinja2
            templated variable in this recipe.

     #. Locally commit these changes.
     #. Locally push these changes to the upstream fork.
     #. Remotely open a pull request (PR) from the upstream fork against the
        `original repository <feedstock repository_>`__.

     See also the `conda-forge FAQ`_ entry `"Using a fork vs a branch when
     updating a recipe." <conda-forge update recipe_>`__

   * Our official `Gentoo Linux ebuild`_, currently hosted at the `raiagent
     overlay`_ maintained by a co-maintainer of BETSE.

Thus begins the dawn of a new scientific epoch.

.. [#tags]
   Do *not* create a lightweight tag, which omits critical metadata (e.g.,
   author identity, descriptive message). *Always* create an annotated tag
   containing this metadata by explicitly passing the ``-a`` option to the
   ``git tag`` subcommand.
.. [#venv]
   Installing this release into a venv requires installing *all* mandatory
   dependencies of this release into this venv from either binary wheels or
   source tarballs. In either case, expect installation to consume non-trivial
   space and time. The cheese shop was not instantiated in a day.

Version Nomenclature
====================

This application should be **versioned** (i.e., assigned a new version)
according to the `Semantic Versioning`_ schema. Each version *must* consist of
three ``.``-delimited integers ``{major}.{minor}.{patch}``, where:

* ``{major}`` is the **major version,** incremented only when either:

  * **Breaking backward compatibility with existing simulation configurations.**
    The public API of this application is its configuration file format rather
    than the public subset of its codebase (e.g., public submodules or classes).
    No codebase change can be considered to break backward compatibility unless
    also changing the simulation configuration file format in a manner rendering
    existing files in the prior format unusable. Note that doing so is
    unequivocally bad and hence *much* discouraged.
  * **Implementing headline-worthy functionality** (e.g., a GUI). Technically,
    this condition breaks the `Semantic Versioning`_ schema, which stipulates
    that *only* changes breaking backward compatibility warrant major bumps.
    But this is the real world. In the real world, significant improvements
    are rewarded with significant version changes.

  In either case, the minor and patch versions both reset to 0.

* ``{minor}`` is the **minor version,** incremented only when implementing
  customary functionality in a manner preserving backward compatibility. In this
  case, only the patch version resets to 0.
* ``{patch}`` is the **patch version,** incremented only when correcting
  outstanding issues in a manner preserving backward compatibility.

When in doubt, bump only the minor version and reset only the patch version.

Codename Nomenclature
=====================

This application should be **code named** (i.e., assigned a new human-readable
code name) according to the following crude distortion of the `Ubuntu code name
schema`_. Each code name *must* consist of two capitalized English words
``{adjective} {bioelectrician}``, where:

* ``{adjective}`` is an arbitrary adjective whose first letter is the same as
  that of the first character of the subsequent ``{bioelectrician}``.
* ``{bioelectrician}`` is the last name of an arbitrary academic associated with
  the long-standing field of bioelectricity.

Unlike the `Ubuntu code name schema`_, the first letter of the code name for
each version need *not* succeed the first letter of the code name for the prior
version. For our insignificant purposes, preserving alphabetization across code
names is a fruitless and hence worthless goal.

.. # ------------------( LINKS ~ betse                       )------------------
.. _readme:
   https://gitlab.com/betse/betse/blob/master/README.rst
.. _tarballs:
   https://gitlab.com/betse/betse/tags
.. _PyPI BETSE:
   https://pypi.python.org/pypi/betse

.. # ------------------( LINKS ~ betse : gentoo              )------------------
.. _Gentoo Linux ebuild:
   https://github.com/leycec/raiagent/tree/master/sci-biology/betse
.. _raiagent overlay:
   https://github.com/leycec/raiagent

.. # ------------------( LINKS ~ betse : conda               )------------------
.. _Anaconda package:
   https://anaconda.org/conda-forge/betse
.. _conda recipe:
   https://github.com/leycec/betse-feedstock/blob/master/recipe/meta.yaml
.. _conda-forge feedstock:
.. _feedstock repository:
   https://github.com/leycec/betse-feedstock

.. # ------------------( LINKS ~ betsee                      )------------------
.. _BETSEE:
   https://gitlab.com/betse/betsee

.. # ------------------( LINKS ~ python                     )------------------
.. _Semantic Versioning:
   http://semver.org
.. _twine:
   https://pypi.python.org/pypi/twine
.. _wheel:
   https://wheel.readthedocs.io

.. # ------------------( LINKS ~ python : conda             )------------------
.. _conda-forge:
   https://conda-forge.org
.. _conda-forge FAQ:
   https://conda-forge.org/docs/conda-forge_gotchas.html
.. _conda-forge update recipe:
   https://conda-forge.org/docs/conda-forge_gotchas.html#using-a-fork-vs-a-branch-when-updating-a-recipe

.. # ------------------( LINKS ~ python : pypi              )------------------
.. _Test PyPI:
   https://testpypi.python.org/pypi
.. _Test PyPI instructions:
   https://wiki.python.org/moin/TestPyPI
.. _Test PyPI user:
   https://testpypi.python.org/pypi?%3Aaction=register_form
.. _PyPI:
   https://pypi.python.org/pypi
.. _PyPI user:
   https://pypi.python.org/pypi?%3Aaction=register_form

.. # ------------------( LINKS ~ software                   )------------------
.. _GitHub:
   https://github.com
.. _Ubuntu code name schema:
   https://wiki.ubuntu.com/DevelopmentCodeNames
