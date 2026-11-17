@echo off
REM ============================================================
REM  Double-click this file to start the FACS app.
REM  It finds Python, checks the needed packages are installed,
REM  starts the app and opens it in your web browser.
REM ============================================================

REM ---------------- SETTINGS (you can change these) ----------------
REM PORT       = the "door number" the app uses on this computer.
REM              Change it if 5050 is already in use (e.g. to 5051).
REM CONDA_ENV  = name of a conda environment to try first, if you have one.
set PORT=5050
set CONDA_ENV=scrna
REM -----------------------------------------------------------------

REM Always work from the folder this .bat file is in (works on X: too)
cd /d "%~dp0"

echo.
echo  FACS app - starting up
echo  ----------------------

REM ---- 1. Find Python ----
REM Try the usual Miniforge / Anaconda / Miniconda / Python install places.
set "PY="
for %%P in (
    "%USERPROFILE%\miniforge3\envs\%CONDA_ENV%\python.exe"
    "%USERPROFILE%\anaconda3\envs\%CONDA_ENV%\python.exe"
    "%USERPROFILE%\miniconda3\envs\%CONDA_ENV%\python.exe"
    "%USERPROFILE%\miniforge3\python.exe"
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
    "%LOCALAPPDATA%\miniforge3\python.exe"
    "%LOCALAPPDATA%\anaconda3\python.exe"
    "C:\ProgramData\miniforge3\python.exe"
    "C:\ProgramData\anaconda3\python.exe"
    "C:\ProgramData\miniconda3\python.exe"
) do (
    if not defined PY if exist %%P set "PY=%%~P"
)

REM Plain Python installs: the "py" launcher, then "python" on the PATH
if not defined PY (
    py -3 --version >nul 2>nul && set "PY=py -3"
)
if not defined PY (
    python --version >nul 2>nul && set "PY=python"
)

if not defined PY goto no_python

echo  Using Python: %PY%

REM ---- 2. Check the needed packages ----
%PY% -c "import flask, pandas, numpy, scipy, matplotlib, openpyxl" >nul 2>nul
if not errorlevel 1 goto run_app

echo.
echo  Some packages the app needs are not installed:
echo  flask, pandas, numpy, scipy, matplotlib, openpyxl
set "ANSWER="
set /p ANSWER= Install them now? Type Y and press Enter, or just Enter to cancel:
if /i not "%ANSWER%"=="Y" goto cancelled

%PY% -m pip install -r requirements.txt
%PY% -c "import flask, pandas, numpy, scipy, matplotlib, openpyxl" >nul 2>nul
if errorlevel 1 goto install_failed

REM ---- 3. Start the app and open the browser ----
:run_app
echo.
echo  The app is opening at  http://localhost:%PORT%
echo  Keep this window open while you use the app.
echo  To stop the app: close this window, or press Ctrl+C.
echo.
REM Open the browser after 3 seconds, so the app has time to start
start "" cmd /c "timeout /t 3 >nul & start http://localhost:%PORT%"
%PY% -m flask --app app run --port %PORT%
echo.
echo  The app has stopped.
pause
exit /b 0

:no_python
echo.
echo  Could not find Python on this computer.
echo  Install Python from https://www.python.org/downloads/
echo  and tick "Add python.exe to PATH" on the first installer screen.
echo  Then double-click this file again.
pause
exit /b 1

:install_failed
echo.
echo  Installing the packages did not work - see the messages above.
echo  You can try by hand:  %PY% -m pip install -r requirements.txt
pause
exit /b 1

:cancelled
echo.
echo  Cancelled. The app needs these packages to run.
pause
exit /b 1
