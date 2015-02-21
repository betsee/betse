# -*- mode: python -*-
a = Analysis(['/usr/bin/betse'],
             pathex=['/home/leycec/py/betse'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='betse',
          debug=False,
          strip=None,
          upx=True,
          console=True )
