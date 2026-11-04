"""The intake gate for agent branches (findings/README.md).

    python tools/gate.py <branch> [--base origin/main]

Lists (1) files the branch changes outside the agent's slots and (2) whether findings/<branch>.md exists
with a Headline, a Stop line for explorations, and every card's fields. Exit 0 if clean, 1 otherwise.
"""
import argparse
import re
import subprocess
import sys

AGENT_NEW_ONLY = (re.compile(r"^definition/theorems/T\d+\.md$"),)
AGENT_SLOTS = (re.compile(r"^findings/"), re.compile(r"^met/"), re.compile(r"^tests/"), re.compile(r"^tools/"),
               re.compile(r"^presets/"), re.compile(r"^problems/(O\d+|README)\.md$"), re.compile(r"^briefs/B\d+.*\.md$"))
STATE = "docs/STATE.md"
FIELDS = ("Kind", "Bears on", "Status", "Check", "Standard name", "Does not show", "Where", "Absorb as", "Retracts")
KINDS = {"break", "no-go", "identity", "construction", "measurement", "correction"}


def git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def changed(base, branch):
    out = git("diff", "--name-status", f"{base}...{branch}")
    return [line.split("\t") for line in out.splitlines() if line]


def state_ok(base, branch):
    """Only lines added to STATE.md's log (after '## Log'), nothing removed or changed elsewhere."""
    diff = git("diff", "-U0", f"{base}...{branch}", "--", STATE)
    old = git("show", f"{base}:{STATE}").splitlines()
    log_line = next(i for i, l in enumerate(old) if l.strip() == "## Log") + 1   # 1-indexed
    for hunk in re.finditer(r"^@@ -(\d+)(?:,(\d+))? \+\d+(?:,\d+)? @@", diff, re.M):
        start, count = int(hunk.group(1)), int(hunk.group(2) if hunk.group(2) is not None else 1)
        if count > 0:              # an old line removed or changed
            return False
        if start < log_line:       # an insertion above the log
            return False
    return True


def check_findings(text, exploration):
    problems = []
    if "## Headline" not in text:
        problems.append("no Headline section")
    if exploration and "## Stop" not in text:
        problems.append("exploration run with no Stop line")
    cards = re.split(r"^### ", text, flags=re.M)[1:]
    if not cards:
        problems.append("no cards")
    for card in cards:
        title = card.splitlines()[0][:60]
        for f in FIELDS:
            if not re.search(rf"^- \*\*{re.escape(f)}[^*]*:\*\*\s*\S", card, re.M):
                problems.append(f"card '{title}': missing {f}")
        m = re.search(r"^- \*\*Kind:\*\*\s*([a-z-]+)", card, re.M)
        if m and m.group(1) not in KINDS:
            problems.append(f"card '{title}': unknown kind {m.group(1)!r}")
        if len(card.strip().splitlines()) > 14:
            problems.append(f"card '{title}': longer than twelve lines")
    return problems


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("branch")
    p.add_argument("--base", default="origin/main")
    a = p.parse_args(argv)
    problems = []
    for status, *paths in changed(a.base, a.branch):
        path = paths[-1]
        if path == STATE:
            if not state_ok(a.base, a.branch):
                problems.append(f"{STATE}: changed outside its log (agents add one log line only)")
            continue
        if any(r.match(path) for r in AGENT_NEW_ONLY):
            if status != "A":
                problems.append(f"{path}: existing theorem file changed (agents add new theorem files only)")
            continue
        if any(r.match(path) for r in AGENT_SLOTS):
            continue
        problems.append(f"{path}: outside the agent's slots ({status})")
    name = "findings/" + a.branch.replace("origin/", "").replace("/", "-") + ".md"
    try:
        text = git("show", f"{a.branch}:{name}")
    except subprocess.CalledProcessError:
        problems.append(f"{name}: missing")
    else:
        exploration = bool(re.search(r"\*\*Run:\*\*\s*exploration", text))
        problems += check_findings(text, exploration)
    for line in problems:
        print(line)
    print(f"gate: {'clean' if not problems else f'{len(problems)} problem(s)'}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
