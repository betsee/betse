Development
===========

For development purposes, BETSE is _editably installable_ (i.e., as a symbolic
link rather than physical copy). As the name implies, editable installations are
modifiable at runtime and hence suitable for development. By the magic of
setuptools eggs and symbolic links, modifications to the copy of BETSE from
which an editable installation originated are propagated back to that
installation implicitly.

## System-wide

BETSE is installable into a system-wide directory as follows:

* **_(Optional)._** Set the current umask to `002`.

        $ umask 002

* Editably install BETSE.

        $ cd "${BETSE_DIR}"
        $ sudo python3 setup.py symlink

The `symlink` command is a BETSE-specific `setuptools` command inspired by the
IPython `setuptools` command of the same name, generalizing the behaviour of the
default `develop` command to system-wide editable installations.

Why? Because the `develop` command is suitable _only_ for user-specific
editable installations. While both `pip` and `setuptools` provide commands for
performing editable installations (e.g., `sudo pip3 install --no-deps
--editable .` and `sudo python3 setup.py develop --no-deps`, respectively),
executable scripts installed by these commands raise fatal exceptions on
failing to find `setuptools`-installed dependencies regardless of whether these
dependencies have already been installed in a system-wide manner. To quote
[IPython developer
MinRK](http://mail.scipy.org/pipermail/ipython-dev/2014-February/013209.html):

    So much hate for setuptools right now.  I can't believe `--no-deps` skips
    dependency installation, but still adds a redundant check to entry points.

## User-specific

BETSE is editably installable into a user-specific venv via either `pip` or
`setuptools` **from within such venv.** While there appears to be no particular
advantage to using one over the other, it remains helpful to note that both
apply. In either case, external executables (e.g., `betse`, `betse-qt`) will
also be installed and usable in the expected manner.

### pip

BETSE is editably installable into a user-specific venv via `pip` as follows:

    $ cd "${BETSE_DIR}"
    $ pip3 install --no-deps --editable .

This installation is uninstallable as follows:

    $ pip3 uninstall betse

### setuptools

BETSE is editably installable into a user-specific venv via `setuptools` as
follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py develop --no-deps

This installation is uninstallable as follows:

    $ cd "${BETSE_DIR}"
    $ ./setup.py develop --uninstall
