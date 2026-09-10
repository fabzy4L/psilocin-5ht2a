"""
Shared test fixtures for md/scripts/. The stage scripts are numbered
(01_import_charmm_gui.py, etc.) so they aren't importable with a plain
`import` statement -- load_script() below loads them by file path via
importlib, the same approach used to regenerate the corrected redock
report during the RMSD-bug fix (docs/REDOCK_BUG_HANDOFF.md).
"""
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MD_SCRIPTS = REPO_ROOT / "md" / "scripts"

# gate_utils.py has a normal (non-numeric) module name, so the scripts'
# own `sys.path.insert(0, ...)` + `from gate_utils import ...` works once
# md/scripts/ is importable -- make sure it is for the test process too.
sys.path.insert(0, str(MD_SCRIPTS))


def load_script(filename: str):
    """Load a md/scripts/NN_name.py file as a module, bypassing the
    numeric-leading-digit import restriction."""
    path = MD_SCRIPTS / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def load_md_script():
    return load_script
