"""B01's selector flags and complete-matrix CLI use the maintained declaration."""
from types import SimpleNamespace

import numpy as np
import pytest

from tools import robot_run


def test_matrix_cli_uses_the_saved_complete_preset(monkeypatch, tmp_path, capsys):
    from met import robot_budget
    calls = []
    monkeypatch.setattr(robot_budget, "run", lambda path, out, progress: (
        calls.append((path, out, progress)) or (tmp_path / "saved", {"rows": ["all"]})))
    monkeypatch.setattr(robot_budget, "summary", lambda record: "complete table\n")
    assert robot_run.main(["--matrix", "--out", str(tmp_path)]) == 0
    assert calls[0][0].name == "robot_budget.yaml"
    assert calls[0][1] == tmp_path and calls[0][2] is None
    assert "complete table" in capsys.readouterr().out


@pytest.mark.parametrize("extra", (["--no-floor"], ["--equation", "ours29"],
                                   ["--runs", "1-baseline"], ["--white", "autocorrelated"],
                                   ["--white", "independent"], ["--rotation", "nominal"]))
def test_matrix_cli_refuses_silently_narrowed_comparisons(extra):
    with pytest.raises(SystemExit) as exc:
        robot_run.main(["--matrix", *extra])
    assert exc.value.code == 2


def test_selector_flags_reach_the_windows_without_changing_static_controls(monkeypatch, capsys):
    calls = []
    cal = SimpleNamespace(mass=.9323)
    declaration = SimpleNamespace(metadata={"seed": 20261030})
    result = dict(median=.93, interval=(.90,.96), raw_interval=(.91,.95),
                  readout=.93, eta=1., g=np.array([.5]))
    monkeypatch.setattr(robot_run.robot, "calibrate", lambda: cal)
    monkeypatch.setattr(robot_run.robot, "calibration_declaration", lambda given: declaration)
    monkeypatch.setattr(robot_run.robot, "windows", lambda *args, **kwargs: (
        calls.append(kwargs) or ([object()], [])))
    monkeypatch.setattr(robot_run.robot, "static_heldout", lambda *args, **kwargs: ([object()], .9323))
    monkeypatch.setattr(robot_run.robot, "series_ratio", lambda readings: 1.)
    monkeypatch.setitem(robot_run.EQUATIONS, "ours30", lambda readings, c: result)
    assert robot_run.main(["--runs", "1-baseline", "--widths", ".2",
                           "--white", "autocorrelated", "--rotation", "static"]) == 0
    assert calls == [dict(floor=True, declaration=declaration, white="autocorrelated", rotation="static")]
    assert "static controls retain their original nominal" in capsys.readouterr().out
