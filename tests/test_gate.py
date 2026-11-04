"""tools/gate.py: the intake gate for agent branches (findings/README.md)."""
import re
import subprocess
from pathlib import Path

import pytest

from tools import gate

ROOT = Path(__file__).resolve().parent.parent


def _example_card():
    text = (ROOT / "findings" / "README.md").read_text()
    block = text.split("## Example")[1].split("```")[1]
    return "# Findings: x\n\n**Run:** B99\n\n## Headline\n- F1\n\n## Cards\n" + block


def test_the_readme_example_card_passes():
    assert gate.check_findings(_example_card(), exploration=False) == []


def test_a_card_missing_a_field_or_with_an_unknown_kind_fails():
    card = _example_card()
    assert any("missing Does not show" in p for p in gate.check_findings(re.sub(r"- \*\*Does not show:\*\*.*\n", "", card), False))
    assert any("unknown kind" in p for p in gate.check_findings(card.replace("**Kind:** no-go", "**Kind:** vibes"), False))
    assert any("Stop" in p for p in gate.check_findings(card, exploration=True))


def _repo(tmp_path):
    def run(*a):
        subprocess.run(["git", *a], cwd=tmp_path, check=True, capture_output=True)
    run("init", "-q", "-b", "main")
    run("config", "user.email", "t@t"); run("config", "user.name", "t")
    (tmp_path / "docs").mkdir(); (tmp_path / "definition").mkdir()
    (tmp_path / "docs/STATE.md").write_text("# State\n\nmain page\n\n## Log\n\n- old line\n")
    (tmp_path / "definition/proof.md").write_text("proof\n")
    run("add", "-A"); run("commit", "-qm", "base")
    run("checkout", "-qb", "codex/b99")
    return run


def test_gate_on_a_clean_and_a_dirty_branch(tmp_path, monkeypatch):
    run = _repo(tmp_path)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "findings").mkdir()
    (tmp_path / "findings/codex-b99.md").write_text(_example_card())
    s = (tmp_path / "docs/STATE.md").read_text().replace("## Log\n\n", "## Log\n\n- new line\n")
    (tmp_path / "docs/STATE.md").write_text(s)
    (tmp_path / "met").mkdir(); (tmp_path / "met/x.py").write_text("x = 1\n")
    run("add", "-A"); run("commit", "-qm", "clean")
    assert gate.main(["codex/b99", "--base", "main"]) == 0
    (tmp_path / "definition/proof.md").write_text("proof, rewritten by the agent\n")
    (tmp_path / "docs/STATE.md").write_text(s.replace("main page", "main page, narrated"))
    run("add", "-A"); run("commit", "-qm", "dirty")
    assert gate.main(["codex/b99", "--base", "main"]) == 1
