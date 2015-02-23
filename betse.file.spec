# -*- mode: python -*-
# ====================[ betse.file.spec                    ]====================
#
# --------------------( SYNOPSIS                           )--------------------
# PyInstaller spec file freezing the CLI script "betse" into a one-file
# executable. PyInstaller evaluates this file as Python code when the
# "freeze_file" setuptools command is run.

# ....................{ IMPORTS                            }....................
# Manually read and evaluate 'betse.core.spec', a mandatory dependency
# performing logic common to both this file and "betse.dir.spec".
with open('betse.core.spec') as core_spec:
    exec(core_spec.read())

# ....................{ MAIN                               }....................
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          **EXE_options)

# --------------------( WASTELANDS                         )--------------------
