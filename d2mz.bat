@echo off
rem d2mz.bat — double-click launcher for the d2mz GUI (Windows).
rem Requires Python 3.8+ installed (the "py" launcher or "python" on PATH).
setlocal
where py >nul 2>nul
if %errorlevel%==0 (
    start "" py -3 "%~dp0d2mz_gui.pyw"
    goto :eof
)
where pythonw >nul 2>nul
if %errorlevel%==0 (
    start "" pythonw "%~dp0d2mz_gui.pyw"
    goto :eof
)
where python >nul 2>nul
if %errorlevel%==0 (
    start "" python "%~dp0d2mz_gui.pyw"
    goto :eof
)
echo Python 3.8+ is required (install it from python.org or the Microsoft Store,
echo then double-click this file again).
pause
