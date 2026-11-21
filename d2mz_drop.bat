@echo off
rem d2mz drop-and-run: drag ONE or SEVERAL folders onto this
rem file -- or a single .d folder. They are converted to mzML
rem 32-bit with the same options as the batch run. Output goes
rem to <dropped folder>\d2mz_output.
if "%~1"=="" (
  echo Nothing to convert.
  echo Drag one or several folders -- or a single .d folder -- onto d2mz_drop.bat.
  echo.
  pause
  goto :eof
)
py -3 "%~dp0d2mz.py" %* --format mzml --precision 32
echo.
echo === done. Press a key to close ===
pause >nul
