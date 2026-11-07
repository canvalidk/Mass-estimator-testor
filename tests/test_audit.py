"""The audit keeps comparison-theorem findings visible with their declared scope."""
from tools import audit

import pytest


@pytest.fixture
def living_layer(tmp_path, monkeypatch):
    definition = tmp_path / "definition"
    theorems = definition / "theorems"
    problems = tmp_path / "problems"
    theorems.mkdir(parents=True)
    problems.mkdir()
    (definition / "premises.md").write_text(
        "| id | name | statement | status |\n"
        "| L1 | one object | one mass | bedrock |\n"
        "| L2 | noisy type | positive noise | bedrock |\n",
        encoding="utf-8",
    )
    (definition / "equation.md").write_text("L1\n", encoding="utf-8")
    (problems / "README.md").write_text(
        "| [O2](O2.md) | proof | missing file |\n", encoding="utf-8"
    )
    (tmp_path / "rulings.md").write_text(
        "### R-2026-10-28-1\n", encoding="utf-8"
    )
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    monkeypatch.setattr(audit, "DEF", definition)
    monkeypatch.setattr(audit, "PROB", problems)
    monkeypatch.setattr(audit, "RUL", tmp_path / "rulings.md")

    def theorem(number, status, *, source="L1", statement="a statement"):
        path = theorems / f"T{number:02}.md"
        path.write_text(
            f"# T{number}: Example\n\n"
            f"**From:** {source}\n"
            "**Test:** none yet\n"
            f"**Status (d4):** {status}\n\n"
            f"## Statement\n\n{statement}\n\n"
            "## Proof or support\n\nUnwritten.\n",
            encoding="utf-8",
        )
        return path

    return theorem


@pytest.mark.parametrize("marker", ["Not load-bearing", "**not load-bearing**"])
def test_comparison_scope_keeps_untested_and_other_findings(living_layer, capsys, marker):
    living_layer(1, "derived")
    living_layer(2, f"{marker} for ours30; comparison readout only", source="L99")

    assert audit.main([]) == 0
    output = capsys.readouterr().out

    assert "untested theorem (2)" in output
    assert "T1 (T01.md): none yet" in output
    assert "T2 (T02.md) [not load-bearing for ours30]: none yet" in output
    assert "theorems/T02.md refers to L99" in output
    assert "premise no theorem uses (1)\n  - L2" in output
    assert "open item with no file (1)\n  - O2" in output


def test_table_scope_comes_from_status_only(living_layer, capsys):
    living_layer(1, "derived", statement="not load-bearing is mentioned in the statement")
    living_layer(2, "Not load-bearing for ours30; comparison readout only")

    assert audit.main(["--table"]) == 0
    rows = capsys.readouterr().out.splitlines()
    active = next(row for row in rows if row.startswith("| T1 |"))
    comparison = next(row for row in rows if row.startswith("| T2 |"))

    assert active.endswith("| from L1; test: none yet |")
    assert comparison.endswith(
        "| from L1; test: none yet; not load-bearing for ours30 (comparison readout only) |"
    )


def test_removing_status_marker_restores_ordinary_finding(living_layer, capsys):
    living_layer(2, "Not load-bearing for ours30; comparison readout only")
    audit.main([])
    assert "T2 (T02.md) [not load-bearing for ours30]" in capsys.readouterr().out

    living_layer(2, "derived")
    audit.main([])
    output = capsys.readouterr().out
    assert "T2 (T02.md): none yet" in output
    assert "not load-bearing" not in output


def test_ids_left_on_exploration_branches_are_known_not_dangling(tmp_path, monkeypatch):
    """explore/README.md's 'Left on exploration branches' table names ids whose files are not on main (7 Nov)."""
    (tmp_path / "explore").mkdir()
    (tmp_path / "explore" / "README.md").write_text(
        "## The board\n\n| E1 | T99 is not in the table below |\n\n"
        "## Left on exploration branches\n\n| Ids | Branch | What |\n| --- | --- | --- |\n"
        "| T54–T56 | `codex/x` | a line |\n| T66 | `codex/y` | one |\n\n## Next\n| T70 | no |\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    assert audit.left_on_branches() == {"T54", "T55", "T56", "T66"}


def test_the_repo_table_covers_the_line_left_on_codex_ours32_consistency():
    left = audit.left_on_branches()
    assert {f"T{n}" for n in range(54, 69)} <= left
    assert not left & {f"T{n}" for n in range(39, 54)}  # harvested, on main
