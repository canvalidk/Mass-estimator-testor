# For agents (Codex and the like)

Read `TASK.md`, then `findings/README.md`. Everything else follows from those two.

- **Every run ends with `findings/<branch>.md`**, in the card format of `findings/README.md`. That file is
  what gets read; the rest of your work is the evidence it points to. A result not on a card is not found.
- **Write only in your slots:** new theorem files, new dated sections in problem files, `met/`, `tests/`,
  `tools/`, `presets/`, your findings file, one line in `docs/STATE.md`'s log. Never `definition/proof.md`,
  `equation.md`, `premises.md`, `ladder.md`, `GOAL.md`, `TASK.md`, STATE's main page or `rulings.md`; put what
  you would change there in a card's "Absorb as" line. `python tools/gate.py <branch>` checks this; a branch
  that fails it is returned unread.
- **Briefed runs** follow `briefs/README.md` and the brief's "Do not" list.
- **Exploration runs** have a target (a probe, a premise, a hole). Pressure it until it breaks, card each
  break with the premise it names, and **stop when a break names part of your own construction instead of a
  premise of the definition**; say so in the findings file's Stop line.
- Before working a fork, name its standard mathematics and cite it (`GOAL.md`, standing rule 6).
- Full suite and `python tools/audit.py` before you finish; give both results in the findings header.
