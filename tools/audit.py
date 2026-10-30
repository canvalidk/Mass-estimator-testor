"""Audit the living layer: definition/, problems/, rulings.md.

Reports, without changing anything:

  1. theorems whose Test line is "none yet" or "partial", retaining the finding
     but marking those whose Status says "not load-bearing";
  2. premises whose status neither cites a dated ruling that exists in rulings.md
     nor carries a grade (bedrock, strong, medium, mild), nor is struck;
  3. premises that no theorem's From line uses;
  4. open items (problems/README.md table) with no file, and problem files not in the table;
  5. ids referenced in From lines, problem files or equation.md that do not exist;
  6. dates written as "ruled <date>" in premises or theorems with no entry in rulings.md.

Exit status 0 always: this is a report, not a gate. `--table` prints the premises and
theorems as one markdown table, the human view of the folder.

Run from the repo root: python tools/audit.py
"""
import re
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEF = ROOT / "definition"
PROB = ROOT / "problems"
RUL = ROOT / "rulings.md"

ID = re.compile(r"(?<![A-Za-z])([FLTO]\d{1,2}′?)(?!\d)")
GRADES = ("bedrock", "strong", "medium", "mild", "struck", "withdrawn")
MONTHS = "Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec"
RULED = re.compile(r"ruled(?:\s+\w+){0,3}?\s+(\d{1,2})\s+(" + MONTHS + r")\b", re.I)
RULING_DATE = re.compile(r"^### R-(\d{4})-(\d{2})-(\d{2})", re.M)


def table_rows(text, prefix):
    rows = {}
    for ln in text.splitlines():
        m = re.match(r"^\|\s*(~~)?(" + prefix + r"\d{1,2}′?)(~~)?\s*\|", ln)
        if m:
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            struck = bool(m.group(1)) or (len(cells) > 1 and cells[1].startswith("~~"))
            rows[m.group(2)] = (struck, cells)
    return rows


def main(argv):
    problems = []
    prem_text = (DEF / "premises.md").read_text(encoding="utf-8")
    prem = table_rows(prem_text, "[FL]")
    theorems = {}
    for f in sorted((DEF / "theorems").glob("T*.md")):
        t = f.read_text(encoding="utf-8")
        head = re.match(r"# (T\d+)", t)
        if not head:
            continue
        tid = head.group(1)
        frm = re.search(r"\*\*From:\*\*\s*(.*)", t)
        test = re.search(r"\*\*Test:\*\*\s*(.*)", t)
        status = re.search(r"\*\*Status[^*]*\*\*\s*(.*)", t)
        status_text = status.group(1).strip() if status else ""
        theorems[tid] = {
            "file": f.name,
            "from": frm.group(1).strip() if frm else "",
            "test": test.group(1).strip() if test else "",
            "status": status_text,
            "not_load_bearing": "not load-bearing" in status_text.casefold(),
            "text": t,
        }

    ruling_days = set()
    rulings_text = ""
    if RUL.exists():
        rulings_text = RUL.read_text(encoding="utf-8")
        for y, mo, d in RULING_DATE.findall(rulings_text):
            ruling_days.add((int(mo), int(d)))
    else:
        problems.append(("rulings", "rulings.md is missing"))

    def month_num(name):
        return ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"].index(name[:3].lower()) + 1

    # 1. untested theorems
    for tid, th in theorems.items():
        if th["test"].startswith("none") or th["test"].startswith("partial"):
            scope = " [not load-bearing for ours30]" if th["not_load_bearing"] else ""
            problems.append(("untested theorem", f"{tid} ({th['file']}){scope}: {th['test'][:80]}"))

    # 2. uncited premises; 6. ruled dates with no ruling
    all_text_for_dates = prem_text + "".join(th["text"] for th in theorems.values())
    for pid, (struck, cells) in prem.items():
        if struck or pid.startswith("F"):
            status = cells[2] if pid.startswith("F") else cells[3]
        else:
            status = cells[3] if len(cells) > 3 else ""
        low = status.lower()
        has_grade = any(g in low for g in GRADES) or "inherited" in low or "proposed" in low or "declared" in low
        has_ruling = bool(RULED.search(status)) or "ruled" in low
        if not struck and not has_grade and not has_ruling:
            problems.append(("uncited premise", f"{pid}: status = {status[:90]}"))
    if ruling_days:
        for m in RULED.finditer(all_text_for_dates):
            day, mon = int(m.group(1)), month_num(m.group(2))
            if (mon, day) not in ruling_days:
                problems.append(("ruling not in rulings.md", f"'{m.group(0)}' cited but no R-2026-{mon:02d}-{day:02d}-n entry"))

    # 3. unused premises
    used = set()
    for th in theorems.values():
        used |= set(ID.findall(th["from"]))
    for pid, (struck, cells) in prem.items():
        if pid.startswith("L") and not struck and pid not in used:
            problems.append(("premise no theorem uses", pid))

    # 4. open items vs problem files
    readme = (PROB / "README.md").read_text(encoding="utf-8") if (PROB / "README.md").exists() else ""
    listed = set(re.findall(r"\[(O\d+)\]\(O\d+\.md\)", readme))
    files = {f.stem for f in PROB.glob("O*.md")}
    for o in sorted(listed - files):
        problems.append(("open item with no file", o))
    for o in sorted(files - listed):
        problems.append(("problem file not in problems/README.md", o))

    # 5. dangling ids
    known = set(prem) | set(theorems) | files | {"O1", "O3", "O8", "O10", "O12"}
    for name, text in [("equation.md", (DEF / "equation.md").read_text(encoding="utf-8"))] + \
                      [(f"problems/{f.name}", f.read_text(encoding="utf-8")) for f in PROB.glob("O*.md")] + \
                      [(f"theorems/{th['file']}", th["from"]) for th in theorems.values()]:
        for ref in sorted(set(ID.findall(text))):
            if ref not in known:
                problems.append(("dangling id", f"{name} refers to {ref}"))

    if "--table" in argv:
        print("| id | name / statement | status |")
        print("| --- | --- | --- |")
        for pid, (struck, cells) in prem.items():
            name = cells[1] if len(cells) > 3 else ""
            stmt = cells[2] if len(cells) > 3 else cells[1]
            status = cells[3] if len(cells) > 3 else cells[2]
            print(f"| {'~~' if struck else ''}{pid}{'~~' if struck else ''} | {name} {stmt} | {status} |")
        for tid, th in theorems.items():
            stmt = re.search(r"## Statement\s+(.*?)\s+## ", th["text"], re.S)
            scope = "; not load-bearing for ours30 (comparison readout only)" if th["not_load_bearing"] else ""
            print(f"| {tid} | {stmt.group(1).strip() if stmt else ''} | from {th['from']}; test: {th['test']}{scope} |")
        return 0

    print(f"definition/: {len(prem)} premise rows, {len(theorems)} theorems; problems/: {len(files)} files; "
          f"rulings.md: {len(RULING_DATE.findall(rulings_text))} rulings")
    if not problems:
        print("audit: nothing out of order")
        return 0
    by_kind = {}
    for kind, msg in problems:
        by_kind.setdefault(kind, []).append(msg)
    for kind, msgs in by_kind.items():
        print(f"\n{kind} ({len(msgs)})")
        for m in msgs:
            print(f"  - {m}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
