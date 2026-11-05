"""
Start the CCC Platform web app.

    python web/run.py            # from the project folder
    python run.py                # from inside web/

Then open http://127.0.0.1:8000
"""

import os
import sys
from pathlib import Path

import uvicorn

WEB_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(WEB_ROOT))

if __name__ == "__main__":
    port = int(os.environ.get("PORT", "8000"))
    reload = os.environ.get("CCC_RELOAD", "0") == "1"
    print(f"\n  CCC Platform → http://127.0.0.1:{port}\n")
    uvicorn.run(
        "backend.main:app",
        host="127.0.0.1",
        port=port,
        reload=reload,
        reload_dirs=[str(WEB_ROOT)] if reload else None,
        app_dir=str(WEB_ROOT),
    )
