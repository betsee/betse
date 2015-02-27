# -*- mode: python -*-
from os import path
import os

a = Analysis(['/usr/bin/betse'],
             pathex=['/home/leycec/py/betse'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)

# Absolute path of the current directory.
CURRENT_DIRNAME = os.getcwd()

# Absolute path of the top-level directory containing all non-Python data files
# to be bundled with such executable.
DATA_ROOT_DIRNAME = path.join(CURRENT_DIRNAME, 'betse', 'data')

# Record all non-Python data files to be bundled with such executable. For
# further details, see http://pythonhosted.org/PyInstaller/#id38.
a.datas += Tree(
    DATA_ROOT_DIRNAME,

    # Relative path of the output top-level directory containing such files.
    prefix = 'data',
)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='betse',
          debug=False,
          strip=None,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='betse')
