"""Build final Markdown from canonical source and perform structural QA."""
from pathlib import Path
from urllib.parse import urlparse
import hashlib
import json
import re

base = Path(__file__).resolve().parent
workspace = base.parents[2]
source = base / 'report-source.md'
destination = workspace / 'pass to VD-docs' / 'mass_estimator_originality_research_2026-09-05.md'
content = source.read_text(encoding='utf-8')
assert destination.parent.resolve().is_relative_to(workspace)
if destination.exists():
    assert destination.read_text(encoding='utf-8') == content, 'Existing deliverable has different content.'
lines = content.splitlines()
assert sum(line.strip() == '$$' for line in lines) % 2 == 0
assert sum(line.startswith('```') for line in lines) % 2 == 0
assert content.count('$') % 2 == 0
assert not re.search(r'\b(?:TODO|TBD|FIXME)\b|turn\d+(?:view|search)|\ue200', content)
headings = [line for line in lines if line.startswith('#')]
assert len(headings) == len(set(headings))
links = re.findall(r'\[([^\]]+)\]\((https?://[^\s)]+)\)', content)
assert len(links) >= 12
for label, url in links:
    assert label.strip() and urlparse(url).netloc
for i, line in enumerate(lines):
    if line.startswith('#') and i+1 < len(lines):
        assert not lines[i+1].strip(), ('Missing blank line after heading',i)
    if line.startswith('|') and i and lines[i-1].startswith('|'):
        assert line.count('|') == lines[i-1].count('|'), ('Uneven table columns',i)
destination.write_bytes(source.read_bytes())
assert destination.read_bytes() == source.read_bytes()
result = {
    'status':'passed', 'artifact':str(destination),
    'words':len(content.split()),'headings':len(headings),'web_citations':len(links),
    'display_math_blocks':sum(line.strip() == '$$' for line in lines)//2,
    'sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),
    'checks':['source/artifact equality','math and code delimiters','headings','table columns','citation structure','no placeholders'],
    'visual_rendering':'not performed; Markdown structural QA only',
    'mathematics':'check_math.py passed; analytic divergence proof reviewed separately'
}
(base/'qa-result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
