.. # ------------------( SYNOPSIS                            )------------------

=========
Releasing
=========

BETSE is typically released (i.e., versioned and tagged at a stable commit,
packaged into a source tarball, and uploaded for public consumption) according
to a rigorous procedure detailed below. While technically optional, this
procedure reduces the likelihood of installation and usage issues by downstream
consumers (e.g., end users, package maintainers) and hence is effectively
mandatory.

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

Distribution
============

BETSE is releasable to all supported platforms as follows:

#. (\ *Optional*\ ) **List all existing tags.** For reference, listing all
   previously created tags *before* creating new tags is often advisable. 

  .. code:: bash
   
     $ git tag

#. **Tag the current commit.** For completeness, create an annotated tag with
   name ``v{version}``, where ``v`` is an arbitrary prefix preserving historical
   consistency with prior tag names and  ``{version}`` is the current value of
   the public ``betse.metadata.__version__`` specifier. :sup:`Do not create a
   lightweight tag, which omits such critical metadata as author identity and
   descriptive message.`

   .. code:: bash

      $ git tag -a v{version}

#. **Push all tags,** including the most recently created tag. Gitlab then
   implicitly exposes source tarballs in various formats (e.g., `.zip`,
   `.tar.bz2`) containing the contents of this repository at this tagged commit
   within the `source tarball archive <tarballs_>`__ of this project. No further
   work is required to upload a source tarball to Gitlab.

   .. code:: bash

      $ git push --tags

#. **Package the codebase into a tarball.**
#. **Upload this tarball to PyPI_.**
#. **Bump the application version.** In preparation for subsequent development
   of the next version of this application:

   #. The public ``betse.metadata.__version__`` specifier should be incremented
      according to the `best practices <Versioning_>`__ provided below.
   #. This modification should be committed and pushed (possibly accompanying
      other modifications).

Versioning
============

This application should typically be **versioned** (i.e., assigned a new
version) in a manner conforming to the `Semantic Versioning`_ schema. Each
version *must* consist of three ``.``-delimited integers
``{major}.{minor}.{patch}``, where:

- ``{major}`` provides the **major version,** incremented only when either:
  - **Breaking backward compatibility with existing simulation configurations.**
    The public API of this application is its configuration file format rather
    than the public subset of its codebase (e.g., public submodules or classes).
    No codebase change can be considered to break backward compatibility unless
    also changing the simulation configuration file format in a manner rendering
    existing files in the prior format unusable. Note that doing so is
    unequivocally bad and hence *much* discouraged.
  - **Implementing headline-worthy functionality** (e.g., a GUI). Technically,
    this condition breaks the `Semantic Versioning`_ schema, which stipulates
    that *only* changes breaking backward compatibility warrant major bumps.
    But this is the real world. In the real world, significant improvements
    are rewarded with significant version changes.
  In either case, the minor and patch versions both reset to 0.
- ``{minor}`` provides the **minor version,** incremented only when implementing
  customary functionality in a manner preserving backward compatibility. In this
  case, only the patch version resets to 0.
- ``{patch}`` provides the **patch version,** incremented only when correcting
  outstanding issues in a manner preserving backward compatibility.

When in doubt, bump only the minor version and reset only the patch version.

.. # ------------------( LINKS                               )------------------
.. _PyPI:
   https://pypi.python.org/pypi
.. _Semantic Versioning:
   http://semver.org
