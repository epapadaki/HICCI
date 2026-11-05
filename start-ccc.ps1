# Start the HICCI web interface (Windows PowerShell), from the project folder.
# Settings such as GEMINI_API_KEY and NCBI_EMAIL are read from the .env file.
Set-Location $PSScriptRoot
if (Test-Path .\.venv\Scripts\Activate.ps1) { & .\.venv\Scripts\Activate.ps1 }
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"
python web\run.py
