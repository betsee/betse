#!/usr/bin/env python3
# -*- mode: python -*-
# ====================[ betse.core.spec                    ]====================
#
# --------------------( SYNOPSIS                           )--------------------
# Low-level dependency of all high-level PyInstaller spec files freezing the CLI
# script "betse" into an executable file or directory.
#
# --------------------( CAVEATS                            )--------------------
# This file expects "betse" to have been editably installed via the custom
# setuptools command "symlink" rather than non-editably installed via the core
# setuptools command "install". For details on how to adjust this file to
# support the latter, see:
#     https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Setuptools-Entry-Point

#FIXME: The produced executables are rather large, at the moment. The reason why
#appears to be indiscriminate inclusion of *ALL* possible matplotlib backends.
#We should probably correct this explicitly excluding the largest such backends
#-- particularly "PyQt4" and subsidiaries. For details on how to exclude modules
#and libraries, see this great stackoverflow answer:
#    https://stackoverflow.com/questions/4890159/python-excluding-modules-pyinstaller/17595149#17595149

#FIXME: When we begin work on freezing "betse-qt", note that the correct
#solution is almost certainly to leverage a so-called "multipackage bundle."
#This prevents us from duplicating embedded libraries in "betse-qt" that are
#already bundled with "betse", dramatically reducing executable filesize. See:
#    http://pythonhosted.org/PyInstaller/#id67

# ....................{ IMPORTS                            }....................
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: PyInstaller spec files are *NOT* valid Python code. Specifically,
# attempting to import modules in such files does nothing at best and silently
# corrupts the current global namespace at worst, preventing further processing.
# The following modules are "imported" only to avoid lint errors during editing.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from os import path
import os, platform, sys

# ....................{ SANITY                             }....................
# If a requisite environment variable exported by a "freeze_*" setuptools
# command does *NOT* exist, this file was manually passed to a PyInstaller
# command without the help of setuptools. Since we are in uncharted territory,
# raise an exception.
if '__FREEZE_MODULE_FILENAME' not in os.environ:
    raise RuntimeError(
        'Setuptools-assisted PyInstaller environment not found.\n'
        'Please freeze under setuptools by running either:\n'
        '    ./setup.py freeze_dir\n'
        '    ./setup.py freeze_file\n'
    )

# ....................{ CONSTANTS ~ os                     }....................
# Cryptographic cipher with which to encrypt the frozen executable.
BLOCK_CIPHER = None

# ....................{ CONSTANTS ~ os                     }....................
# True if the current operating system is Linux.
IS_OS_LINUX = platform.system() == 'Linux'

# True if the current operating system is OS X.
IS_OS_OS_X = platform.system() == 'Darwin'

# True if the current operating system is Windows.
IS_OS_WINDOWS = sys.platform == 'win32'

# ....................{ CONSTANTS ~ path                   }....................
# Absolute path of the current directory.
CURRENT_DIRNAME = os.getcwd()

# Absolute path of the root Python module running such script. For portability,
# such path references a module in the Python package tree for this application
# rather than a script wrapper previously installed by setuptools (e.g.,
# "C:\Miniconda3\Scripts\betse.exe'). While such wrappers are valid Python
# scripts under POSIX-compatible platforms, such wrappers are binary blobs
# *NOT* analyzable by PyInstaller under other platforms (e.g., Windows).
MODULE_ROOT_FILENAME = os.environ['__FREEZE_MODULE_FILENAME']

# Absolute path of the top-level directory containing all non-Python data files
# to be bundled with such executable.
DATA_ROOT_DIRNAME = path.join(CURRENT_DIRNAME, 'betse', 'data')

# ....................{ CONSTANTS ~ exe                    }....................
# Basename of such executable.
EXE_basename = 'betse'
if IS_OS_WINDOWS:
    EXE_basename += '.exe'

# Keyword arguments to be subsequently passed to EXE() by parent specifications.
EXE_options = {
    'name': EXE_basename,
    'debug': False,
    'strip': None,
    'upx': True,
    'console': True,

    #FIXME: Uncomment *AFTER* adding such file. Note that this is probably
    #Linux- and Windows-specific, however; OS X appears to require a completely
    #different icon file format. (Naturally. Why not?)
    # 'icon': path.join(APPPATH, "..", "res", "Icon.ico"),
}

# ....................{ MAIN ~ analysis                    }....................
# Since PyInstaller curiously fails to do so, inform the current user of:
#
# * The version of the active Python interpreter.
# * The absolute path of the entry module.
print('Freezing Python {} entry point:\n\t{}'.format(
    platform.python_version(), MODULE_ROOT_FILENAME), file=sys.stderr)

# Analyze such top-level Python module.
a = Analysis(
    [MODULE_ROOT_FILENAME],
    pathex = [CURRENT_DIRNAME],
    # cipher = BLOCK_CIPHER,

    # List of the names of all imported and hence required modules *NOT*
    # automatically detectable by PyInstaller.
    hiddenimports = [],

    # List of the relative paths of all directories containing non-runtime
    # PyInstaller hooks.
    hookspath = [path.join('freeze', 'hooks')],

    # List of the relative paths of all directories containing runtime
    # PyInstaller hooks.
    runtime_hooks = None,

    # List of the fully-qualified names of all importables (e.g., modules,
    # packages, extensions) to *NOT* be frozen into such executable.
    excludes = None,
)

# a = Analysis(['/usr/local/bin/betse'],
#              pathex=['freeze'],
#              hiddenimports=[],
#              hookspath=['freeze/hooks'],
#              runtime_hooks=None)

# ....................{ MAIN ~ datas                       }....................
# Record all non-Python data files to be frozen into such executable. For
# details, see http://pythonhosted.org/PyInstaller/#id38.
a.datas += Tree(
    DATA_ROOT_DIRNAME,

    # Relative path of the output top-level directory containing such files.
    prefix = 'data',
)

# ....................{ MAIN ~ pyz                         }....................
# Create an intermediate archive containing only Python modules and scripts.
pyz = PYZ(
    a.pure,
    # cipher = BLOCK_CIPHER,
)

# ....................{ MAIN ~ exe, collect                }....................
# If freezing in "one-file" mode, do so.
if os.environ['__FREEZE_MODE'] == 'file':
    exe = EXE(
        pyz,
        a.scripts,
        a.binaries,
        a.zipfiles,
        a.datas,
        **EXE_options)
# Else, freeze in "one-directory" mode.
else:
    exe = EXE(
        pyz,
        a.scripts,
        exclude_binaries=True,
        **EXE_options
    )
    coll = COLLECT(
        exe,
        a.binaries,
        a.zipfiles,
        a.datas,
        strip=None,
        upx=True,

        #FIXME: Incorrect. This should be "betse-qt" and hence conditionally
        #depend on the current entry point.
        name='betse'
    )

# --------------------( WASTELANDS                         )--------------------
#FUXME: This hard-codes the CLI-based BETSE entry point and hence obviously
#fails to generalize to the GUI-based BETSE entry point (e.g., "betse-qt".) To
#support both, we'll need to generalize this. Unfortunately, it doesn't appear
#to be feasible to pass arguments to specification files. Possibly? Hmm.
#
#Actually, we could probably do it the old-fashioned way: set an environment
#variable in the "setup/freeze.py" setuptools command -- named, say,
#${_FREEZE_MODULE_ROOT_FILENAME} -- whose value is the absolute path of the
#desired module here. Assuming that works, such command should conditionally set
#such value in a platform-specific way. Under vanilla Windows, such path will
#need to be munged from the passed executable path to the desired Python wrapper
#(e.g., from
#"C:\Miniconda3\Scripts\betse.exe' to
#"C:\Miniconda3\Scripts\betse-script.py'); on all other platforms, such path
#should simply be set to the current entry point path unmodified.
#
#If passing environment variables does *NOT* work, see if there's a way that we
#can force a global module attribute to be set by PyInstaller (e.g., via a
#"pyinstaller" CLI argument). If even that doesn't work, things get considerably
#more complex. We'll need to perform such munging here (trivial) *AND* to split
#both "betse.file.spec" and "betse.dir.spec" into two new specs each: the first
#specific to the CLI-based BETSE entry point and the second specific to the
#GUI-based BETSE entry point (non-trivial).
#
#Actually, that's sufficiently horrible that we should do the following as a
#fallback: use a temporary file whose filename is guaranteed *NOT* to change but
#be sufficiently unique as to conflict with nothing (e.g.,
#"${system_temp_dir}/__MEI_PASS_betse_entry_point_path") and whose contents are
#simply the absolute path of the desired entry point. Simple, if annoying.
#
#Actually, I'd rather not munge around with external files. Too fragile, really.
#So let's just go the route of splitting up our specifications if we have to.
# MODULE_ROOT_FILENAME = path.join(
#     CURRENT_DIRNAME, 'betse', 'cli', '__main__.py')

# IS_OS_WINDOWS = platform.system() == 'Windows'
        # 'scipy.sparse.sparsetools._csr',

# ....................{ MAIN ~ binaries : scipy            }....................
# Include all `scipy`-specific shared libraries.

# import PyInstaller.log as logging
# from PyInstaller.utils.hooks import hookutils
# from PyInstaller.utils import misc
#
# #FUXME: The following code is unused. PyInstaller appears to be unable to
# #tolerate functions in spec files. (We're not kidding.) That said, the following
# #
#
# #FUXME: The following should probably be pushed to the PyInstaller codebase
# #itself as official "scipy" hooks. For a roughly analogous hook, see
# #"PyInstaller/hooks/hook-PyQt4.QtGui.py". Perhaps add such hooks by adding a
# #similar hook() method to a new "PyInstaller/hooks/hook-scipy.py" script.
# # hook_scipy(a)
#
# def hook_scipy(a):
#     '''
#     Include all `scipy`-specific shared libraries.
#     '''
#     a.binaries.extend(collect_package_binaries('scipy'))
#
# #FUXME: This is sufficiently useful that it should probably be generalized
# #into a new "hookutils" utility function.
# def collect_package_binaries(package_name):
#     '''
#     Get a `mod.binaries`-formatted list of 3-tuples installing all dynamic
#     libraries recursively found in the source directory containing the package
#     with the passed name to corresponding subdirectories of a target directory
#     with the same name.
#     '''
#     # Absolute path of such package's "__init__.py" script.
#     package_init_file = hookutils.eval_statement('''
#         import {package_name}
#         print({package_name}.__file__)
#     '''.format(package_name = package_name))
#
#     return collect_binaries(
#         os.path.dirname(package_init_file), package_name)
#
# #FUXME: This is sufficiently useful that it should probably be generalized
# #into a new "hookutils" utility function.
# #FUXME: qt5_qml_plugins_binaries() should be refactored to call this instead.
#
# def collect_binaries(src_root_dir, trg_root_dir):
#     '''
#     Get a `mod.binaries`-formatted list of 3-tuples installing all dynamic
#     libraries recursively found in the passed source directory to corresponding
#     subdirectories of the passed target directory.
#
#     Such source directory should be absolute; such target directory should be
#     relative. For example, to list all shared libraries for `scipy`:
#
#         >>> from PyInstaller.utils.hooks.hookutils import collect_binaries
#         >>> collect_binaries(
#         ...     '/usr/lib64/python3.3/site-packages/scipy', 'scipy')
#     '''
#     logger = logging.getLogger(__name__)
#     binaries = []
#     src_dlls = misc.dlls_in_subdirs(src_root_dir)
#
#     for src_dll in src_dlls:
#         src_dll_relative = os.path.relpath(src_dll, src_root_dir)
#         src_dll_relative_dir, src_dll_base = os.path.split(src_dll_relative)
#         trg_dll = os.path.join(trg_root_dir, src_dll_relative_dir, src_dll_base)
#
#         binaries.append((trg_dll, src_dll, 'BINARY'))
#         logger.info(
#             'Dynamic library "%s" installing to "%s".', src_dll, trg_dll)
#
#     return binaries

    # import scipy
    # a.binaries.extend(
    #     collect_binaries(os.path.dirname(scipy.__file__), 'scipy'))
    # Absolute paths of shared libraries in such directory.
        # Relative path of such library relative to such top-level directory.
            #FUXME: Just use "src_dll_relative", no?
    # # Relative path of the top-level "site-packages/scipy" directory.
    # trg_root_dir = 'scipy'
    #
    # # Absolute path of the top-level "site-packages/scipy" directory.
    # src_root_dir = os.path.dirname(scipy.__file__)

    # assert isinstance(module_object, module), \
    #     '"{}" not a module.'.format(module_object)
    # List of TOC-style tuples describing such libraries.
    # scipy_binaries = []

    # excludes = ['PyQt4', 'PyQt4.QtCore', 'PyQt4.QtGui'],
# PyInstaller-specific collection of such files, which may be passed as is to
# the EXE(), COLLECT(), and BUNDLE() functions. For further details, see:
    # Absolute path of the input top-level directory containing such files.
    # root = path.join(CURRENT_DIRNAME, 'betse', 'data'),
# a = Analysis(['/usr/bin/betse'],
# Absolute path of the root Python package for this application.
# PACKAGE_ROOT_DIRNAME = path.join(os.getcwd(), ')
