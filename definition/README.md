# definition/: the living deliverable

The goal of the project is the function that is the definition of inertial mass, with its
proof. This folder is that deliverable as it stands, and nothing else:

| File | Holds |
| --- | --- |
| `type.md` | what the function takes and returns |
| `premises.md` | the frame (F), the properties of mass (L, graded), the test-constraining rulings (R) |
| `theorems/` | one file per theorem (T): statement, what it rests on, the test that checks it, status |
| `equation.md` | the equation itself, stated once |
| `ladder.md` | the version's commitments in order, what each removes from the body, and the status of every open freedom by that order (TASK.md §7) |

Open items live in `../problems/`, one file each. Rulings live in `../rulings.md`.

## Rules

1. **Edit in place.** There are no d1/d2 drafts here; `git log` is the history, and each file
   ends with a short History list for the changes that matter.
2. **Every status cites something.** A premise is *ruled* (then it names a date that is in
   `../rulings.md`) or *graded* (bedrock, strong, medium, mild). A theorem is *proved* only if
   its file names the test that checks it, or says why no test can.
3. **Nothing is edited to agree with the code, and no code is edited to agree with this
   folder** (tester README, rule 5). A disagreement is a finding: a strict xfail in `tests/`
   naming both numbers, and a line here.
4. **No records.** A conversation's results go into the file they concern: a theorem, a
   problem file's attempts list, a ruling appended to `../rulings.md`. If a result fits none of
   those, it is not kept here.
5. `python tools/audit.py` lists what is out of order: theorems with no test, premises with no
   citation, premises nothing uses, open items with no problem file, references to ids that do
   not exist. Run it before committing.

## Provenance

Created 2026-10-24 by splitting lemma list d4
(`../docs/vd-docs-context/mass_estimator_lemmas_2026-10-10.md`, provisional, 10 Oct) with every
cell carried verbatim. The d4 banner still applies: everything marked 10 Oct there is written as
accepted on Can's instruction and is provisional as a whole; in particular L14's status is
unchanged and the L14/swap trade is unruled.
