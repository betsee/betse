# -*- mode: python -*-
# ====================[ betse.dir.spec                     ]====================
#
# --------------------( SYNOPSIS                           )--------------------
# PyInstaller spec file freezing the CLI script "betse" into a one-directory
# executable. PyInstaller evaluates this file as Python code when the
# "freeze_dir" setuptools command is run.

# ....................{ IMPORTS                            }....................
# Manually read and evaluate 'betse.core.spec', a mandatory dependency
# performing logic common to both this file and "betse.file.spec".
with open('betse.core.spec') as core_spec:
    exec(core_spec.read())

# ....................{ MAIN                               }....................
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          **EXE_options)
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='betse')

# --------------------( WASTELANDS                         )--------------------
