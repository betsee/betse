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

# ....................{ IMPORTS                            }....................
from os import path
import os, platform

# ....................{ CONSTANTS ~ os                     }....................
# True if the current operating system is Linux.
IS_OS_LINUX = platform.system() == 'Linux'

# True if the current operating system is OS X.
IS_OS_OS_X = platform.system() == 'Darwin'

# True if the current operating system is Windows.
IS_OS_WINDOWS = platform.system() == 'Windows'

# ....................{ CONSTANTS ~ path                   }....................
# Absolute path of the current directory.
CURRENT_DIRNAME = os.getcwd()

# Absolute path of the root Python module running such script. For portability,
# such path references a module in the Python package tree for this application
# rather than a script wrapper previously installed by setuptools (e.g.,
# "/usr/bin/betse").
MODULE_ROOT_FILENAME = path.join(
    CURRENT_DIRNAME, 'betse', 'cli', '__main__.py')

# Absolute path of the top-level directory containing all non-Python data files
# to be bundled with such executable.
DATA_ROOT_DIRNAME = path.join(CURRENT_DIRNAME, 'betse', 'data')

# ....................{ CONSTANTS ~ exe                    }....................
# Basename of the executable to be output.
EXE_basename = 'betse'
if IS_OS_WINDOWS:
    EXE_basename += '.exe'

# Keyword arguments to be subsequently passed to EXE() by parent files.
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

# ....................{ MAIN                               }....................
# Analyze such top-level Python module.
a = Analysis(
    [MODULE_ROOT_FILENAME],
    pathex = [CURRENT_DIRNAME],
    hiddenimports = [],
    hookspath = None,
    runtime_hooks = None,
)

# Record all non-Python data files to be bundled with such executable. For
# further details, see http://pythonhosted.org/PyInstaller/#id38.
a.datas += Tree(
    DATA_ROOT_DIRNAME,

    # Relative path of the output top-level directory containing such files.
    prefix = 'data',
)

# Create an intermediate archive containing only Python modules and scripts.
pyz = PYZ(a.pure)

# --------------------( WASTELANDS                         )--------------------
# PyInstaller-specific collection of such files, which may be passed as is to
# the EXE(), COLLECT(), and BUNDLE() functions. For further details, see:
    # Absolute path of the input top-level directory containing such files.
    # root = path.join(CURRENT_DIRNAME, 'betse', 'data'),
# a = Analysis(['/usr/bin/betse'],
# Absolute path of the root Python package for this application.
# PACKAGE_ROOT_DIRNAME = path.join(os.getcwd(), ')
