.. # ------------------( SYNOPSIS                           )------------------

.. # FIXME: Merge ``doc/md/DEVELOP.md`` into this file; then, remove that file.
.. # FIXME: Merge ``doc/md/INSTALL.rst`` into this file; then, remove both that
.. # and the ``doc/md/INSTALL.md`` file.

============
Installation
============

BETSE requires:

-  Either **Microsoft Windows,** **Apple macOS,** or a **Linux distribution.**
   All other platforms (e.g., Android, FreeBSD) are explicitly unsupported at
   this time.
-  At least **Python 3.5** (e.g., 3.5, 3.8). All prior Python versions (e.g.,
   Python 2.7, 3.4) are explicitly unsupported.
-  At least **4GB RAM;** ideally, at least **16GB RAM.**\ [#thirtytwobit_not]_

.. [#thirtytwobit_not]
   Ergo, BETSE effectively requires a **64-bit system.**  32-bit systems impose
   a so-called `"3GB barrier" <https://en.wikipedia.org/wiki/3_GB_barrier>`__
   preventing usage of more than 3—4GB of available RAM, which rarely suffices
   for even small-scale BETSE simulations. This constraint extends to *all*
   non-server 32-bit editions of Microsoft Windows.

Linux
=====

BETSE is installable under *most* Linux distributions as follows:

#. Open a **terminal.**\ [#linux_terminal]
#. Copy-and-paste the following commands into this terminal:

   #. Install **BETSE.**

      .. code-block:: console

         pip3 install betse

   #. [\ *Optional*\ ] Test **BETSE.**

      .. code-block:: console

         betse -v try

.. [#linux_terminal]
   To open a terminal under:

   - **Ubuntu Linux:**

     #. Type ``Ctrl``\ +\ ``Alt``\ +\ ``t``.

macOS
=====

BETSE is installable under Apple macOS as follows:

#. Install the **Python 3.x** [#python2_not]_ (e.g., 3.7) variant of
   Anaconda_.\ [#anaconda_not]_
#. Open the **Finder**.
#. Open the **Applications** folder.
#. Open the **Utilities** folder.
#. Open *Terminal.app*.
#. Copy-and-paste the following commands into this terminal:

   #. Enable conda-forge_.

      .. code-block:: console

         conda config --add channels conda-forge

   #. Install **BETSE.**\ [#conda_package]_

      .. code-block:: console

         conda install betse

   #. [\ *Optional*\ ] Test **BETSE.**

      .. code-block:: console

         betse -v try

.. [#python2_not]
   Do *not* install **Python 2.7.** BETSE strictly requires **Python 3.x.** 

.. # FIXME: Change the <Advanced_> link to a macOS-specific link.

.. [#anaconda_not]
   If you prefer *not* to install Anaconda_, BETSE is also `manually
   installable <Advanced_>`__ via a third-party package manager (e.g.,
   Homebrew_, MacPorts_). Doing so is non-trivial and, if performed improperly,
   could produce a performance-crippled single-core installation of BETSE.
   Anaconda_ suffers no such issues and is guaranteed to produce a
   performance-optimized multicore installation of BETSE. We strongly recommend
   Anaconda_ – even when you think you know better.

.. [#conda_package]
   This command installs both the `most recent stable release of BETSE <conda
   package_>`__ *and* all mandatory and most optional dependencies of this
   release. Older stable releases are installable in a similar manner (e.g.,
   ``conda install betse=0.7.0`` for BETSE 0.7.0). All `Anaconda packages`_ are
   kindly hosted by the `non-profit conda-forge organization <conda-forge_>`__.

Windows
=======

BETSE is installable under Microsoft Windows 10 as follows:\ [#windows_not]_

#. Emulate **Ubuntu Linux** via the `Windows Subsystem for Linux (WSL)
   <WSL_>`__ bundled with Windows 10.
#. Open an **Ubuntu Linux terminal.**
#. Copy-and-paste the following commands into this terminal:

   #. Update the `Advanced Package Tool (APT) <APT_>`__ package cache.
 
     .. code-block:: console
 
        sudo apt update
   
   #. Upgrade all previously installed packages.
 
     .. code-block:: console
 
        sudo apt upgrade
 
   #. Install **Python 3.x.**\ [#python2_not]_
 
     .. code-block:: console
 
        sudo apt install python3-pip
 
   #. Install **BETSE.**
 
     .. code-block:: console
 
        pip3 install betse
 
   #. [\ *Optional*\ ] Test **BETSE.**
 
     .. code-block:: console
 
        betse -v try

.. [#windows_not]
   The `Windows Subsystem for Linux (WSL) <WSL_>`__ and (thus BETSE itself) is
   *only* installable under **Windows 10.** Under older Windows versions, BETSE
   may be installed from a `virtual Linux guest <VirtualBox_>`__.

.. # FIXME: Simplify this subsection to be Git-specific.

Developers
----------

For developers and advanced users, *any* version of BETSE – including the live
repository and prior stable releases – is manually installable as follows:

.. # FIXME: Consider embedding the subset of ``doc/md/INSTALL.md`` pertaining
.. # to dependency installation here.

#. Install **Python 3.x** and all dependencies required by BETSE. Under:

   - **Linux,** install these dependencies via your distribution-specific
     package manager (e.g., APT_ under Debian-based distributions). Avoid other
     package managers (e.g., ``pip``, ``conda``) where feasible.\ [#pip_not]_
   - **macOS,** either:

     - (\ *Recommended*\ ) Install the **Python 3.x** variant of Anaconda_.
     - Or both:

       #. Install a third-party package manager (e.g., Homebrew_, MacPorts_).
          Apple does *not* provide a package manager out-of-the-box.
       #. Install these dependencies via that package manager. Avoid other
          package managers (e.g., ``pip``, ``conda``) where feasible.\
          [#pip_not]_

   - **Windows,** install the **Python 3.x** variant of Anaconda_.\ [#windows]_

#. Open a **terminal.**
#. **Download** either:

   - **The unstable BETSE repository** as follows:

     #. Install Git_.
     #. Clone the ``master`` branch of this repository.

        .. code-block:: console

           git clone https://gitlab.com/betse/betse.git

     #. Prepare for installation.

        .. code-block:: console

           cd betse

   - **Any stable BETSE release,** including the most recent, as follows:

     #. Visit our `source tarball archive <tarballs_>`__.
     #. Click the download icon to the right of the desired release and select
        *Download tar.gz*.
     #. Extract the downloaded tarball into the current directory.

        .. code-block:: console

           tar -xvzf betse-*.tar.gz

     #. (\ *Optional*\ ) Remove this tarball.

        .. code-block:: console

           rm betse-*.tar.gz

     #. Prepare for installation.

        .. code-block:: console

           cd betse-*

#. **Install BETSE** either:

   - (\ *Recommended*\ ) **Editably,** installing a cross-platform symbolic link
     to the current BETSE codebase. Modifications to this code are applied
     immediately *without* requiring reinstallation.

     .. code-block:: console

        sudo pip install --editable .

   - **Non-editably,** installing a physical copy of the current BETSE codebase.
     Modifications to this code are ignored and thus require reinstallation.

     .. code-block:: console

        sudo pip install .

#. (\ *Optional*\ ) **Test BETSE,** running all modelling phases of a sample
   simulation from a new directory.

   .. code-block:: console

      cd /tmp && betse try

.. # FIXME:  Actually, "pip" should now install OpenBLAS-optimized wheels for
.. # NumPy and SciPy. Ergo, the discussion below no longer applies.

.. [#pip_not]
   Do *not* install scientific dependencies (e.g., NumPy_, SciPy_) with either
   ``pip`` or ``easy_install``; doing so typically degrades BETSE to
   single-core performance. To optimize BETSE across multiple cores, *always*
   install these dependencies with your platform-specific package manager
   (e.g., Homebrew_, APT_).

.. [#windows]
   Unlike Linux and macOS, Anaconda_ is (\ *effectively*\ ) required under
   Windows. Due to Microsoft's lack of support for `POSIX`_\ -compliant
   toolchains, *no* reasonable alternatives for installing multicore-aware
   scientific dependencies exist.


.. # FIXME: Refactor everything below to provide alternative installation
.. # instructions to those delineated above.

Advanced
--------

BETSE is installable with only `two simple steps </README.rst>`__ on all
supported platforms – complete with multicore-aware hardware
optimizations. These steps leverage scientific standards for Python
packaging, including the cross-platform
**`Anaconda <https://www.continuum.io/downloads>`__** Python
distribution *and* **```pip`` <https://pypi.python.org/pypi/pip>`__**
Python package manager.

For advanced users preferring to manually install dependencies with
existing package managers (e.g.,
`APT <https://en.wikipedia.org/wiki/Advanced_Packaging_Tool>`__) rather
than `Anaconda <https://www.continuum.io/downloads>`__, BETSE may also
be manually installed in a platform-specific manner. This approach has
the obvious advantage of cleanly integrating with existing packaging
regimes but the non-obvious disadvantage of typically installing a
single-core version of BETSE with *no* multicore-aware hardware
optimizations. \ *That's bad.*\ 

Due to the difficulty of manually installing BETSE in a multicore-aware
manner, the `simple installation instructions </README.rst>`__ are
*strongly* recommended. For completeness, this document nonetheless
details the manual approach for several popular package managers.

POSIX
-----

BETSE is manually installable with *most* Linux-centric package
managers.

Debian
~~~~~~

Under `Debian <https://www.debian.org>`__-based Linux distributions
(e.g., `Linux Mint <https://www.linuxmint.com>`__,
`Ubuntu <https://www.ubuntu.com>`__), all mandatory dependencies are
installable in a system-wide manner as follows:

::

    $ sudo apt-get install python3-dev python3-dill python3-matplotlib \
      python3-numpy python3-pil python3-pip python3-scipy python3-setuptools \
      python3-six python3-yaml tcl tk

Under some (especially older) `Debian <https://www.debian.org>`__-based
Linux distributions, the above instructions may not entirely suffice to
satisfy all installation-time or runtime requirements. Under these
distributions, dependencies may require some form of recompilation,
relinking, or reinstallation.

Updated Matplotlib
''''''''''''''''''

BETSE requires a fairly recent version of matplotlib. If the newest
version of matplotlib installed by your distribution is insufficient,
the newest version of matplotlib is installable in a system-wide manner
as follows:

::

    $ sudo apt-get uninstall python3-matplotlib &&
      sudo apt-get install gcc gfortran libfreetype6-dev libpng-dev \
        libpython3-all-dev tcl-dev tk-dev &&
      sudo pip3 install matplotlib[all]

Optimized BLAS and LAPACK
'''''''''''''''''''''''''

BETSE strongly recommends that optimized (rather than the unoptimized
default) implementations of the BLAS and LAPACK APIs for linear algebra
be used. While there exist numerous alternatives both open-source (e.g.,
CBLAS) and proprietary (e.g., MKL), the following instructions assume
use of either ATLAS or OpenBLAS.

ATLAS
     

Automatically Tuned Linear Algebra Software (ATLAS) is the standard
baseline for all optimized BLAS and LAPACK implementations. ATLAS is
installable in a system-wide manner as follows:

::

    $ sudo apt-get install build-essential libatlas-base-dev

Note that OpenBLAS and ATLAS *cannot* be installed at the same time.

OpenBLAS
        

OpenBLAS is a more performant (*but arguably less stable*) optimized
BLAS and LAPACK implementation. While ATLAS is recommended for new
users, experienced users requiring improved performance may benefit from
installing OpenBLAS instead. OpenBLAS is installable in a system-wide
manner as follows:

::

    $ sudo apt-get install build-essential libopenblas-dev

Note that OpenBLAS and ATLAS *cannot* be installed at the same time.

Gentoo
~~~~~~

Under `Gentoo <https://www.gentoo.org>`__-based Linux distributions
(e.g., `Chrome OS <https://en.wikipedia.org/wiki/Chrome_OS>`__,
`Sabayon <https://www.sabayon.org>`__), all mandatory and optional
dependencies are installable in a system-wide manner as follows:

1. **Install ```layman`` <https://wiki.gentoo.org/wiki/Layman>`__,**
   Gentoo's official overlay manager.

   ::

       $ sudo emerge layman
       $ sudo echo 'source /var/lib/layman/make.conf' >> /etc/portage/make.conf

-  **Add the ```raiagent`` <https://github.com/leycec/raiagent>`__
   overlay,** religiously maintained by a `BETSE
   co-maintainer <https://github.com/leycec>`__.

   ::

       $ sudo layman -a raiagent

-  **Synchronize overlays.**

   ::

       $ sudo layman -S

-  Either:

-  ***(Recommended)*** Install the `optimized BLAS and LAPACK
   stack <https://wiki.gentoo.org/wiki/User_talk:Houseofsuns>`__
   published by the
   ```science`` <https://github.com/gentoo-science/sci>`__ overlay.
   While technically optional, failing to do so *will* reduce BETSE to
   unoptimized single-core behavior. To properly install this stack, see
   these `authoritative
   instructions <https://wiki.gentoo.org/wiki/User_talk:Houseofsuns>`__.
-  Disable the ``smp`` USE flag enabled by default for BETSE. In this
   case, the default unoptimized BLAS and LAPACK stack will be linked
   against instead.

   ::

           $ sudo echo 'sci-biology/betse -smp' >> /etc/portage/package.use

-  **Unmask BETSE.** Either:

-  ***(Recommended)*** Unmask the most recent stable release of BETSE.

   ::

           $ sudo echo '>=sci-biology/betse-0.4.1' >> /etc/portage/package.accept_keywords

-  Unmask the most recent unstable commit to the BETSE ``git``
   repository.

   ::

           $ sudo echo '>=sci-biology/betse-0.4.1 **' >> /etc/portage/package.accept_keywords

-  **Install BETSE.**

   ::

       $ sudo emerge betse

macOS (old)
-----------

Under Apple macOS, all mandatory dependencies are installable in a
system-wide manner with either:

-  ***(Recommended)*** **`Homebrew <http://brew.sh>`__,** an unofficial
   package manager for macOS. Homebrew provides robust support for
   features commonly required by BETSE developers, including the
   capacity to install older rather than merely the newest versions of
   packages.
-  **`MacPorts <https://www.macports.org>`__,** another unofficial
   package manager for macOS. MacPorts lacks robust support for features
   commonly required by BETSE developers, as described above. Since
   Homebrew and MacPorts install packages into different system
   directories (i.e., ``/usr/local`` for Homebrew and ``/opt`` for
   MacPorts), the two *can* technically be used on the same system.
   However, this is generally discouraged. If you currently use and
   prefer MacPorts, we recommend adopting the following instructions to
   use MacPorts rather than Homebrew.

For simplicity, the following instructions assume use of Homebrew:

1.  **Register as an `Apple
    Developer <https://developer.apple.com>`__.** While free,
    registration requires an existing Apple ID and hence ownership of an
    existing Apple product. \ *We don't make the awful rules. We only
    complain about them.*\ 
2.  **Upgrade your system** to the most recently released minor version
    for your currently installed major version of macOS. For example, if
    your system is macOS **10.8.3** (*Mountain Lion*), upgrade to
    **10.8.5** (*Mountain Lion*). Homebrew requires recent command-line
    tools (e.g., ``clang``, ``gcc``), requiring requires recent XCode
    Command Line Tools (CLT), requiring a recent version of XCode,
    requiring a recent version of macOS. Provided your system meets the
    minimum requirements noted above, it should *not* be necessary to
    upgrade your system to a newer major version of macOS (e.g., from
    10.8.5 to 10.9.5).
3.  **Open a terminal window** (e.g., by running the pre-bundled
    ``Applications/Utilities/Terminal.app`` application). All commands
    prefixed by ``$`` below *must* be run from within a terminal window.
    Note that, by Unix convention, the ``$`` prefix only denotes the
    default Bash shell prompt and should *not* actually be typed (e.g.,
    type ``xcode-select –install`` rather than
    ``$ xcode-select –install`` when asked to do so below). Likewise,
    the ``<return>`` key should be typed after each such command.
4.  If an older version of the XCode Command Line Tools (CLT) has
    already been installed, **`manually
    uninstall <https://stackoverflow.com/questions/27438457/xcode-6-1-how-to-uninstall-command-line-tools>`__
    the CLT.** While XCode itself is safely upgradable merely by
    installing a new version, the CLT generally is *not*. \ *You can
    thank Apple for that.*\ 
5.  **Download and install the most recent version of
    `XCode <https://developer.apple.com/downloads>`__** available for
    your version of macOS. While free, this download requires an Apple
    Developer login.
6.  ***(Optional)*** After installing Xcode, perform the following
    *before* running Xcode:
7.  **Instruct Gatekeeper to implicitly trust Xcode.** Gatekeeper is the
    macOS application security manager. By default, Gatekeeper uselessly
    verifies Xcode via a labouriously time-consuming and safely
    skippable bureaucratic process requiring in upwards of twenty
    minutes on lower-end laptops. Note that verification is *not* safely
    skippable for "dubious" applications downloaded from third-party
    sources.

    ::

            $ sudo xattr -d com.apple.quarantine /Applications/Xcode.app

8.  **Run Xcode** (e.g., by double-clicking ``Applications/Xcode`` from
    the Finder). If you did *not* instruct Gatekeeper to implicitly
    trust this application as described above, grab a bag of greasy
    popcorn and `*Blade Runner (The Final
    Cut)* <https://en.wikipedia.org/wiki/Versions_of_Blade_Runner>`__.
    You'll need both.
9.  **Agree to the Xcode license.** This *must* be done before
    attempting to run any Xcode-bundled commands from the terminal
    (e.g., ``clang``, ``gcc``, ``git``).
10. ***(Optional)* Close Xcode.**
11. **Download and install the exact same version of the `XCode Command
    Line Tools <https://developer.apple.com/downloads>`__ (CLT)** as the
    installed version of XCode. Attempting to install an older or newer
    version of the CLT may superficially succeed but *will* result in
    obscure and difficult-to-debug issues on attempting to install
    dependencies with Homebrew or MacPorts. There are various approaches
    to installing the correct version of the CLT – some inherently safer
    than others. Either:
12. ***(Recommended)*** **Manually download and install the CLT:**

    1. Browse to the `Apple Developer
       Downloads <https://developer.apple.com/downloads>`__ site.
    2. Enter ``xcode`` into the search bar.
    3. Manually search the resulting hits for the installed version of
       XCode.
    4. Note the official date of this version's release (e.g., June 12,
       2013 for XCode 4.6.3).
    5. Manually search the resulting hits for the most recent version of
       the CLT *preceding* this date (e.g., April 11, 2013 for the CLT
       corresponding to XCode 4.6.3).
    6. Download and install this version.

13. ***(Not recommended)*** **Automatically download and install the
    CLT.** While error-prone and hence discouraged, automatically
    downloading and installing the CLT with Apple-based automation *is*
    technically feasible in common edge cases. If your system has been
    upgraded to both the most recently released minor version of the
    currently installed major version of macOS *and* the most recently
    released version of XCode for that version of macOS, the following
    command *should* suffice. If in doubt, prefer the manual approach
    listed above instead.

    ::

            $ xcode-select –install

14. **Download and install `Homebrew <http://brew.sh>`__.** While these
    dependencies are also technically installable via
    `MacPorts <https://www.macports.org>`__, Homebrew provides
    significantly more robust support for features of interest to BETSE
    users. Critically, this includes the capacity to install alternative
    versions of dependencies rather than merely the newest.

    ::

        $ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

15. **Manually prepend the current ``${PATH}``** by the absoute paths of
    all directories to which Homebrew installs packages. To do so
    permanently, append the following line to the appropriate startup
    dotfile in your home directory for your preferred shell (e.g.,
    ``.bashrc`` for Bash, the default macOS shell).

    ::

        export PATH="/usr/local/bin:/usr/local/sbin:${PATH}"

16. **Activate this ``${PATH}`` change.** Either:

-  ***(Recommended)*** Close the current terminal window and open a new
   terminal window.
-  ***(Not recommended)*** Manually source the modified dotfile: e.g.,

   ::

           $ source ~/.bashrc

1. ***(Optional)*** Inspect your Homebrew installation for potential
   issues. The following command should report that
   ``"Your system is ready to brew."`` If it does *not*, consider
   resolving all reported issues before continuing.

   ::

       $ brew doctor

2. **Install all dependencies.**

   ::

       $ brew tap homebrew/python &&
         brew install python3 &&
         pip3 install --upgrade pip setuptools wheel &&
         brew install matplotlib --with-python3 --without-python &&
         brew install numpy --with-python3 --without-python &&
         brew install pillow --with-python3 --without-python &&
         brew install scipy --with-python3 --without-python &&
         brew install libyaml &&
         pip3 install dill pyyaml

Note that Homebrew is a source-based package manager and hence
relatively slow. Expect the installation process to require anywhere
from several hours to several days, depending on hardware performance.
We wish we were kidding.

Note also that these instructions link ``numpy`` against the most
optimized multicore implementation of the BLAS and LAPACK APIs available
under macOS as of this writing: Apple's **`Accelerate
Framework <https://developer.apple.com/reference/accelerate/1668466-blas>`__.**
No further BLAS or LAPACK configuration is required or recommended.

Windows (old)
-------------

    **Note:** these instructions are *woefully* inadequate at present,
    encouraging installation of the
    `Cygwin <https://www.cygwin.com>`__-based
    **`Babun <https://babun.github.io>`__** wrapper rather than usage of
    the existing **`Bash on ubuntu on
    Windows <https://msdn.microsoft.com/en-us/commandline/wsl/about>`__**
    environment bundled with the Windows 10's Anniversary Update. Until
    these instructions are updated accordingly, Windows users are
    *strongly* encouraged to follow the `simple installation
    instructions </README.rst>`__ instead.

Under Microsoft Windows, all mandatory dependencies are installable in a
system-wide manner via any number of POSIX compatibility layers. For
simplicity, the following instructions assume use of the
`Miniconda <http://conda.pydata.org/miniconda.html>`__ Python
distribution *and* `Babun <http://babun.github.io>`__ POSIX
compatibility layer under 64-bit Windows:

1. Download and install **`Babun <https://babun.github.io>`__**, an
   open-source `Cygwin <https://www.cygwin.com>`__ convenience wrapper
   complete with ``pact``, a CLI-based package manager for Windows. Due
   to unreconcilable flaws in Windows' non-POSIX-compatible process
   model, Cygwin and hence Babun is incompatible with all Windows
   applications on the `Big List of Dodgy Apps
   (BLODA) <https://cygwin.com/faq/faq.html#faq.using.bloda>`__.
   Unfortunately, this includes most antivirus software. If Babun begins
   behaving erratically, consider temporarily disabling such software
   for the duration of Babun usage. (This is the fault of neither Babun
   nor Cygwin!)
2. Download and install the 64-bit Python 3 Windows version of
   **Miniconda**. (See the "Wine" subsection below for further details.)
3. Double click the desktop shortcut ``babun`` to open a new terminal
   window.
4. Prioritize Miniconda- over Babun-installed Python packages. By
   default, Babun prioritizes Babun- over Miniconda-installed Python
   packages. Since Babun packages only a subset of the dependencies
   required by BETSE, Miniconda's ``conda`` rather than Babun's ``pact``
   package manager must be used to install such dependencies. To permit
   this, modify the ``${PATH}`` global exported at Babun startup by
   editing the ``.zshrc`` file in your home directory as follows:

   ::

       # Alter this...
       export PATH=$HOME/bin:/usr/local/bin:$PATH

       # ...to this.
       export MINICONDA_HOME="/cygdrive/c/Miniconda3"
       export PATH="${MINICONDA_HOME}:${MINICONDA_HOME}/Scripts:${HOME}/bin:${PATH}"

5. Apply such changes to the current shell session.

   ::

       $ source ~/.zshrc

6. Install Python dependencies via ``conda``, Miniconda's package
   manager:

   ::

       $ conda install dill numpy matplotlib pyside pyyaml pywin32 scipy

7. Symbolically link the Python 3 executable ``python.exe`` installed by
   Miniconda to ``python3``. For disambiguity, numerous core scripts
   including BETSE's ``setup.py`` installer run Python 3 as ``python3``
   rather than ``python``. For unknown reasons, the Python 3-specific
   version of Miniconda under Windows does *not* appear to provide a
   ``python3`` executable. To fix this:

   ::

       $ ln -s /c/Miniconda3/python.exe /c/Miniconda3/python3

Optional
~~~~~~~~

BETSE optionally leverages (but does *not* strictly require) the
following dependencies where available at runtime:

-  `NetworkX <https://networkx.github.io>`__ >= 1.11, for optionally
   analyzing BETSE networks.
-  `pprofile <https://github.com/vpelletier/pprofile>`__ >= 1.8, for
   optionally profiling BETSE in a line-granular manner.
-  `ptpython <https://github.com/jonathanslenders/ptpython>`__ >= 0.29,
   for optionally wrapping the BETSE REPL with an improved interface. By
   default, the BETSE REPL leverages the stock Python REPL.
-  `py.test <http://pytest.org>`__ >= 2.8.0, for optionally running unit
   tests.
-  `PyDot <https://github.com/erocarrera/pydot>`__ >= 1.0.28 and
   `GraphViz <http://www.graphviz.org>`__ >= 2.38, for optionally
   visualizing BETSE networks.
-  `PyInstaller <http://www.pyinstaller.org>`__ >= 3.0, for optionally
   freezing BETSE.
-  `UPX <http://upx.sourceforge.net>`__ (any version), for optionally
   compressing frozen BETSE executables.

These dependencies are installable as follows.

NetworkX
^^^^^^^^

To optionally analyze networks (e.g., gene regulatory, biochemical
reaction), BETSE requires NetworkX, a pure-Python graph theoretic
framework. This dependency is installable in a system-wide manner as
follows:

-  Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

   ::

       $ sudo apt-get install python3-networkx

-  Under all other supported platforms: Under Linux, additionally prefix
   this command by ``sudo``.

   ::

       $ pip3 install networkx

``pprofile``
^^^^^^^^^^^^

To optionally profile the BETSE codebase with line-granularity into
`callgrind <http://kcachegrind.sourceforge.net/>`__-compatible profile
files, BETSE requires ``pprofile``, an advanced pure-Python line
profiler. This dependency is installable in a system-wide manner as
follows:

-  Under all supported platforms: Under Linux, additionally prefix this
   command by ``sudo``.

   ::

       $ pip3 install pprofile

``ptpython``
^^^^^^^^^^^^

To optionally wrap the BETSE REPL with an improved interface providing
syntax highlighting, multiline editing, autocompletion, and presumably
more, BETSE requires ``ptpython``, an advanced pure-Python REPL. This
dependency is installable in a system-wide manner as follows:

-  Under all supported platforms: Under Linux, additionally prefix this
   command by ``sudo``.

   ::

       $ pip3 install ptpython

``py.test``
^^^^^^^^^^^

To optionally `run tests <#testing>`__, BETSE requires ``py.test``, a
pure-Python test harness. This dependency is installable in a
system-wide manner as follows:

-  Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

   ::

       $ sudo apt-get install python3-pytest

-  Under all other supported platforms: Under Linux, additionally prefix
   this command by ``sudo``.

   ::

       $ pip3 install pytest

``py.test`` Plugins
'''''''''''''''''''

While optional, BETSE provides out-of-the-box support for the following
third-party ``py.test`` plugins:

-  ``pytest-xdist``, parallelizing test runs across all available
   processors. ``py.test`` itself provides *no* built-in support for
   parallelization! Since BETSE's test suite is computationally
   expensive (if not prohibitive), this plugin is a hard prerequisite
   for sanity preservation.

Contributors are strongly encouraged to install these optional
dependencies, which BETSE's test suite will then implicitly detect and
set accordingly. These dependencies are installable in a system-wide
manner as follows:

-  Under all other supported platforms: Under Linux, additionally prefix
   this command by ``sudo``.

   ::

       $ pip3 install pytest-xdist

PyDot + GraphViz
^^^^^^^^^^^^^^^^

To optionally visualize networks (e.g., gene regulatory, biochemical
reaction), BETSE requires:

-  PyDot, a high-level pure-Python GraphViz wrapper.
-  GraphViz, a low-level C-based graph theoretic visualizer.

These dependencies are installable in a system-wide manner as follows:

-  For PyDot:
-  Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

   ::

           $ sudo apt-get install python3-pydot

-  Under all other supported platforms: Under Linux, additionally prefix
   this command by ``sudo``.

   ::

           $ pip3 install pydot

-  For GraphViz:
-  Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

   ::

           $ sudo apt-get install graphviz

-  Under Apple OS X:

   ::

           $ brew install graphviz

PyInstaller
^^^^^^^^^^^

To optionally `freeze BETSE <#freezing>`__, BETSE requires PyInstaller,
a non-pure-Python cross-platform command-line utility for freezing
Python applications. This dependency is installable in a system-wide
manner as follows:

-  Under all supported platforms: Under Linux, additionally prefix this
   command by ``sudo``.

   ::

       $ sudo pip3 install pyinstaller

UPX (Ultimate Packer for eXecutables)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To optionally compress executables while `freezing BETSE <#freezing>`__,
BETSE requires UPX, a non-Python cross-platform command-line utility for
compressing arbitrary executables. This dependency is installable in a
system-wide manner as follows:

-  Under Debian-based Linux distributions (e.g., Linux Mint, Ubuntu):

   ::

       $ sudo apt-get install upx-ucl

-  Under Apple OS X:

   ::

       $ brew install upx

Installation
------------

BETSE itself is installable into either:

-  A system-wide directory accessible to all users of the current
   system.
-  A venv (i.e., virtual environment) isolated to the current user.

The latter has the advantage of avoiding conflicts with already
installed system-wide Python and non-Python packages (e.g., in the event
that BETSE requires different versions of such packages), but the
corresponding disadvantage of requiring reinstallation of such packages
and all transitive dependencies of such packages. Since several
dependencies are heavy-weight (e.g., Qt4) and hence costly to reinstall,
this is a notable disadvantage.

Note that the string ``${BETSE\_DIR}`` should be replaced everywhere
below by the absolute path of the top-level directory containing the
source for BETSE.

System-wide
~~~~~~~~~~~

BETSE is installable into a system-wide directory as follows:

-  Compile BETSE.

   ::

       $ cd "${BETSE_DIR}"
       $ python3 setup.py build

-  Install BETSE.

   ::

       $ sudo python3 setup.py easy_install --no-deps .

Curiously, although the ``develop`` command for ``setuptools`` provides
a ``--no-deps`` option, the ``install`` command does not. Hence, the
``easy_install`` command is called above instead.

BETSE is subsequently uninstallable via ``pip`` as follows:

::

    $ sudo pip uninstall betse

User-specific
~~~~~~~~~~~~~

BETSE is installable into a user-specific
`venv <https://docs.python.org/3/library/venv.html>`__ by running the
following commands **from within that venv**:

::

    $ cd "${BETSE_DIR}"
    $ ./setup.py install

This command should *not* be run outside of a venv. Doing so will
reinstall all dependencies of BETSE already installed by the system-wide
package manager (e.g., ``apt-get``). This may superficially appear to
work but invites obscure and difficult to debug conflicts at BETSE
runtime between dependencies reinstalled by ``setuptools`` and
dependencies already installed by such package maneger.

BETSE is subsequently uninstallable via ``pip`` as follows:

::

    $ pip uninstall betse

Docker
~~~~~~

BETSE is also installable into a
`Docker <https://www.docker.com>`__-hosted Linux distribution contained
within an existing Linux distribution, circumventing the need to install
BETSE directly into an existing system. For simplicity, the following
instructions assume usage of Docker's official Ubuntu image:

1. Install `**Docker** <https://docs.docker.com/engine/installation>`__.
2. Instruct the Xauthority security mechanism to ignore hostnames,
   permitting Docker containers with different hostnames than that of
   the local host to access the current X11 socket. \_Thanks to `Jürgen
   Weigert <https://stackoverflow.com/users/3936284/j%c3%bcrgen-weigert>`__
   for his observant `Stackoverflow
   answer <https://stackoverflow.com/a/25280523>`__ inspiring this
   snippet.

   ::

       $ touch /tmp/.docker.xauth && xauth nlist :0 |
             sed -e 's/^..../ffff/' |
             xauth -f /tmp/.docker.xauth nmerge -

3. Download the latest version of `Continuum
   Analytics <https://www.continuum.io/downloads>`__' official `Anaconda
   3 Docker image <https://hub.docker.com/r/continuumio/anaconda3>`__
   and instantiate this image as a new Docker container named ``betse``
   running an interactive Bash session mounting the X11 socket of the
   host's current X11 session.

   ::

       $ docker run -it\
             --name betse\
             -v /tmp/.X11-unix:/tmp/.X11-unix\
             -v /tmp/.docker.xauth:/tmp/.docker.xauth\
             -e DISPLAY=$DISPLAY\
             -e XAUTHORITY=/tmp/.docker.xauth\
             continuumio/anaconda3 bash

4. Run the following commands from within this container:
5. ***(Optional)*** Test the X11 connection by running ``xeyes``.

   ::

           $ apt-get update && apt-get install -y x11-apps && xeyes

6. Download the live version of BETSE into the ``${HOME}`` directory of
   the current user (i.e., ``root``).

   ::

           $ cd ~ && git clone https://gitlab.com/betse/betse.git

7. Install BETSE.

   ::

           $ cd betse && python3 setup.py install

8. ***(Optional)*** Test BETSE by running a sample simulation.

   ::

           $ cd /tmp && betse try && rm -rf sample_sim

9. Exit this session.

   ::

           $ exit

To resume the previously instantiated container:

1. Restart this container.

   ::

       $ docker start betse

2. Reenter this container by running another interactive Bash session.

   ::

       $ docker attach betse

.. # ------------------( LINKS ~ os : linux                 )------------------
.. _APT:
   https://en.wikipedia.org/wiki/Advanced_Packaging_Tool
.. _POSIX:
   https://en.wikipedia.org/wiki/POSIX
.. _Ubuntu:
.. _Ubuntu Linux:
   https://www.ubuntu.com
.. _Ubuntu Linux 16.04 (Xenial Xerus):
   http://releases.ubuntu.com/16.04

.. # ------------------( LINKS ~ os : macos                 )------------------
.. _Homebrew:
   http://brew.sh
.. _MacPorts:
   https://www.macports.org

.. # ------------------( LINKS ~ os : windows               )------------------
.. _WSL:
   https://msdn.microsoft.com/en-us/commandline/wsl/install-win10

.. # ------------------( LINKS ~ soft                       )------------------
.. _Atom:
   https://atom.io
.. _dill:
   https://pypi.python.org/pypi/dill
.. _FFmpeg:
   https://ffmpeg.org
.. _Git:
   https://git-scm.com/downloads
.. _GitLab-CI:
   https://about.gitlab.com/gitlab-ci
.. _Graphviz:
   http://www.graphviz.org
.. _Libav:
   https://libav.org
.. _MEncoder:
   https://en.wikipedia.org/wiki/MEncoder
.. _VirtualBox:
   https://www.virtualbox.org
.. _YAML:
   http://yaml.org

.. # ------------------( LINKS ~ soft : py                  )------------------
.. _imageio:
   https://imageio.github.io
.. _Matplotlib:
   http://matplotlib.org
.. _NumPy:
   http://www.numpy.org
.. _Python 3:
   https://www.python.org
.. _pip:
   https://pip.pypa.io
.. _py.test:
   http://pytest.org
.. _SciPy:
   http://www.scipy.org

.. # ------------------( LINKS ~ soft : py : conda          )------------------
.. _Anaconda:
   https://www.anaconda.com/download
.. _Anaconda packages:
   https://anaconda.org
.. _conda-forge:
   https://conda-forge.org
