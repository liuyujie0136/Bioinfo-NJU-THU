@echo off
color f9
latex MCM
latex MCM
latex MCM
Call clean.bat
dvipdfm MCM
del *.dvi
exit