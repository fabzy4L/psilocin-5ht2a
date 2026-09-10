"""
05_analysis.py's MDAnalysis-dependent run_analysis() isn't covered here --
it's isolated behind a local import specifically so the pure gate/geometry
functions stay testable without MDAnalysis installed (it isn't in this
test environment). These tests cover is_plateaued() and
salt_bridge_distance(), which are the parts with actual logic.
"""
import numpy as np

from conftest import load_script


def test_is_plateaued_true_for_flat_tail():
    mod = load_script("05_analysis.py")
    # Drifts early, flat for the back half.
    series = np.array([5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.05, 0.95, 1.0, 1.02])
    plateaued, std = mod.is_plateaued(series, tail_frac=0.5, tolerance=0.5)
    assert plateaued is True


def test_is_plateaued_false_for_still_drifting():
    mod = load_script("05_analysis.py")
    series = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    plateaued, std = mod.is_plateaued(series, tail_frac=0.5, tolerance=0.5)
    assert plateaued is False


def test_is_plateaued_empty_series():
    mod = load_script("05_analysis.py")
    plateaued, std = mod.is_plateaued(np.array([]))
    assert plateaued is False
    assert std != std  # NaN


def test_salt_bridge_distance_picks_nearer_oxygen():
    mod = load_script("05_analysis.py")
    amine = np.array([0.0, 0.0, 0.0])
    od1 = np.array([3.0, 0.0, 0.0])   # distance 3.0
    od2 = np.array([0.0, 4.0, 0.0])   # distance 4.0

    d = mod.salt_bridge_distance(amine, od1, od2)

    assert d == 3.0


def test_salt_bridge_distance_matches_manual_investigation_scale():
    """Sanity check against the real crystal-structure geometry found
    during the redock-bug investigation (docs/REDOCK_BUG_HANDOFF.md):
    Asp155 OD1/OD2 vs. 7LD's N2, ~2.8-3.3 A in the crystal pose."""
    mod = load_script("05_analysis.py")
    amine = np.array([24.891, 39.201, 52.668])
    od1 = np.array([26.873, 38.215, 50.894])
    od2 = np.array([26.549, 36.421, 52.124])

    d = mod.salt_bridge_distance(amine, od1, od2)

    assert 2.5 < d < 3.5
