"""Run the unmodified upstream NATMI with current pandas.

NATMI (github.com/forrest-lab/NATMI, last change 2021) finishes its Excel
summary with ``ExcelWriter.save()``. pandas 2.0 removed ``save()``; ``close()``
is its documented replacement and writes the same file. This folder is added to
PYTHONPATH only for the NATMI process (see run_methods.py), so Python imports
this module at start-up and NATMI itself stays untouched.
"""

try:
    import pandas as _pd

    if not hasattr(_pd.ExcelWriter, "save"):
        _pd.ExcelWriter.save = _pd.ExcelWriter.close
except Exception:  # pandas missing or changed: let NATMI report its own error
    pass
