@echo off
REM Start the CCC Platform (Windows). Run from this folder after installation.
cd /d "%~dp0"
if exist .venv\Scripts\activate.bat call .venv\Scripts\activate.bat
python web\run.py
pause
