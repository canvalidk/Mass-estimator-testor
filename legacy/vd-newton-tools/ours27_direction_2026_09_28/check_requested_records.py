"""Structural checks on the four user-requested Markdown records."""
from pathlib import Path
import json
import re
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "06_newton_analysis" / "mass_estimator_records"
PIN = "b4bfc2040031ef010416334d0c5f4aa2ad05cbd8"
FILES = sorted(DEST.glob("*.md"))
assert len(FILES) == 4
errors = []
rows = []
commands = set()
for path in FILES:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    assert "\ufffd" not in text
    for n, line in enumerate(text.splitlines(), 1):
        if line.rstrip() != line:
            errors.append(f"{path.name}:{n}: trailing whitespace")
    # Dollar delimiters must all belong to complete inline or display formulas.
    formula_re = re.compile(r"\$\$([\s\S]*?)\$\$|\$([^\n$]*?)\$")
    formulas = list(formula_re.finditer(text))
    leftover = formula_re.sub("", text)
    if "$" in leftover:
        errors.append(f"{path.name}: unmatched math delimiter")
    for match in formulas:
        math = match.group(1) if match.group(1) is not None else match.group(2)
        depth = 0
        for char in re.sub(r"\\[{}]", "", math):
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth < 0:
                    errors.append(f"{path.name}: premature closing math brace")
        if depth:
            errors.append(f"{path.name}: unbalanced math braces")
        begins = re.findall(r"\\begin\{([^}]+)\}", math)
        ends = re.findall(r"\\end\{([^}]+)\}", math)
        if begins != ends:
            errors.append(f"{path.name}: environment mismatch")
        if len(re.findall(r"\\left\b", math)) != len(re.findall(r"\\right\b", math)):
            errors.append(f"{path.name}: left/right mismatch")
        commands.update(re.findall(r"\\([A-Za-z]+)", math))
    links = re.findall(r"\[[^\]]*\]\(([^)]+)\)", text)
    for link in links:
        if link.startswith("https://github.com/canvalidk/VD-docs/"):
            if f"/blob/{PIN}/" not in link:
                errors.append(f"{path.name}: unpinned managed link {link}")
        elif not link.startswith("https://"):
            if not (path.parent / link).is_file():
                errors.append(f"{path.name}: missing sibling link {link}")
    if re.search(r"\b(?:TODO|TBD|FIXME)\b", text):
        errors.append(f"{path.name}: placeholder")
    rows.append({
        "file": path.name,
        "words_approx": len(text.split()),
        "formulas": len(formulas),
        "links": len(links),
    })

# Check the exact arithmetic actually reproduced in the records.
J = Fraction(1, 10) * Fraction(1, 2)**2 + Fraction(99, 10) * Fraction(10, 101)**2
G = Fraction(1, 2)**2 + Fraction(10, 101)**2
assert J == Fraction(49801, 408040)
assert G == Fraction(10601, 40804)
assert J / (J + G) == Fraction(49801, 155811)
assert Fraction(5, 6) * (J + G) / J == Fraction(259685, 99602)
assert Fraction(1, 2) * (1 + 16) - Fraction(5, 2)**2 == Fraction(9, 4)

summary = {
    "records": rows,
    "math_commands": sorted(commands),
    "errors": errors,
    "scope": "Markdown/link/math-delimiter checks and quoted exact arithmetic; not a rendered-math check.",
}
print(json.dumps(summary, indent=2))
assert not errors
