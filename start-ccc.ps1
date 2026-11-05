Set-Location C:\Users\eugpa\webccc
& .\.venv\Scripts\Activate.ps1
$env:NCBI_EMAIL = "eugpap02@gmail.com"
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"
python web\run.py
