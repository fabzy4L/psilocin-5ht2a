from conftest import load_script

SAMPLE_XVG = """\
# This file was created by gmx energy
# gmx energy is part of G R O M A C S
@    title "Temperature"
@    xaxis  label "Time (ps)"
@    yaxis  label "Temperature (K)"
@TYPE xy
0.0 298.5
10.0 299.1
20.0 300.2
80.0 300.0
90.0 299.9
100.0 300.1
"""


def test_parse_xvg_from_file(tmp_path):
    mod = load_script("03_equilibrate.py")
    xvg_path = tmp_path / "temp.xvg"
    xvg_path.write_text(SAMPLE_XVG)

    points = mod.parse_xvg(xvg_path)

    assert points == [
        (0.0, 298.5), (10.0, 299.1), (20.0, 300.2),
        (80.0, 300.0), (90.0, 299.9), (100.0, 300.1),
    ]


def test_tail_stats_windows_by_time(tmp_path):
    mod = load_script("03_equilibrate.py")
    xvg_path = tmp_path / "temp.xvg"
    xvg_path.write_text(SAMPLE_XVG)
    points = mod.parse_xvg(xvg_path)

    mean, std, n = mod.tail_stats(points, tail_ps=30.0)

    # last 30 ps of a 100 ps trace = t >= 70 -> the 80/90/100 points only
    assert n == 3
    assert abs(mean - ((300.0 + 299.9 + 300.1) / 3)) < 1e-9


def test_tail_stats_empty_series_returns_nan():
    mod = load_script("03_equilibrate.py")
    mean, std, n = mod.tail_stats([], tail_ps=100.0)
    assert n == 0
    assert mean != mean  # NaN != NaN
