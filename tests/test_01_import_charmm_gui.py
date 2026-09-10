"""
Covers the actual bug class this pipeline is trying to prevent: a
CHARMM-GUI export that's missing something (a topology, an index, a
coordinate file, or all .mdp files) should FAIL loudly here, not pass
through and break 02_minimize.py three steps later with a confusing gmx
error.
"""
from conftest import load_script


def _make_complete_export(root):
    (root / "topol.top").write_text("; dummy topology\n")
    (root / "index.ndx").write_text("[ System ]\n1 2 3\n")
    (root / "step5_input.gro").write_text("dummy\n1\n1.0 1.0 1.0\n")
    (root / "step6.1_equilibration.mdp").write_text("integrator = md\n")
    (root / "step7_production.mdp").write_text("integrator = md\n")


def test_validate_complete_export_passes(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["missing"] == []
    assert found["topol"] is not None
    assert found["index"] is not None
    assert found["coordinates"] is not None
    assert len(found["mdp_files"]) == 2


def test_validate_missing_topology_fails(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)
    (tmp_path / "topol.top").unlink()

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["topol"] is None
    assert any("topology" in m for m in found["missing"])


def test_validate_missing_mdp_files_fails(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)
    for mdp in tmp_path.glob("*.mdp"):
        mdp.unlink()

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["mdp_files"] == []
    assert any("mdp" in m for m in found["missing"])


def test_validate_accepts_alternate_coordinate_name(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    (tmp_path / "topol.top").write_text("; dummy\n")
    (tmp_path / "index.ndx").write_text("[ System ]\n")
    (tmp_path / "system.gro").write_text("dummy\n")
    (tmp_path / "step7_production.mdp").write_text("integrator = md\n")

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["coordinates"] is not None
    assert found["missing"] == []


def test_stage_into_copies_found_files(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    src = tmp_path / "src"
    src.mkdir()
    _make_complete_export(src)
    found = mod.validate_charmm_gui_dir(src)

    out_dir = tmp_path / "staged"
    mod.stage_into(out_dir, found)

    assert (out_dir / "topol.top").exists()
    assert (out_dir / "index.ndx").exists()
    assert (out_dir / "step5_input.gro").exists()
    assert (out_dir / "step6.1_equilibration.mdp").exists()
    assert (out_dir / "step7_production.mdp").exists()
