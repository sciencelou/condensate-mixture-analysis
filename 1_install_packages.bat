@echo off
chcp 65001 >nul
cd /d "%~dp0"
where python >nul 2>nul && (set PY=python) || (set PY=py)
%PY% -m pip install -r requirements.txt
echo.
echo Python packages installed. R packages are installed automatically on the first run.
pause
