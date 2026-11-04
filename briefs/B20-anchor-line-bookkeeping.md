# B20: the score-law and anchor constructions marked as benchmarks, not probes

**Status:** open.
**Item:** theorems T56–T65 (`codex/ours32-consistency`); `definition/theorems/README.md`; `definition/proof.md` G24, G26, G27; `tools/audit.py`.
**Branch:** `codex/b20`, from the head of `codex/ours32-consistency`.
**Written:** 2026-11-04, Claude.

## The sorted fork (Claude's, provisional; do not re-sort it)

From T56 on, the branch builds a law of the form $Q(dm)\propto\exp[-T(m)^2/2+\lambda L(m)]\,\Pi(dm)$: the exact
conditional score statistic (T45/T50) exponentiated, a vanishing profile anchor $\lambda=(1+R)^{-\alpha}$ to
suppress the false roots T56 found, and carried state (T60–T62). **Standard name:** a Laplace-type /
quasi-Bayesian estimator (Chernozhukov & Hong 2003) with a Godambe sandwich; T57 cites it. That class is
calibrated as a *class*: any criterion with the right first-order behaviour gives a member with the same
guarantees, so the construction can never select one function. Its unforced pieces are its own status
lines: the anchor class, the exponent $\alpha$, the permitted state. Under TASK §1 the anchor is a bolted-on
clause (the score alone cannot tell a profile maximum from a minimum, so a bit of profile is glued on), and
its breaks (T58, T63) name the anchor, not a premise. So the line is **good-estimator work**: valuable as a
benchmark of attainable width and for its no-go theorems (T56, T58, T63), not a probe of the definition.
This is the bookkeeping B09 did for the cloud readout.

## Deliverable

1. In `definition/theorems/README.md`, mark T54–T65 with scope "benchmark: quasi-Bayesian score law;
   comparison only, not on the path" (T54/T55's joint-protocol constructions included; T45 and T50, the
   exact conditional score *tests*, stay as they are: B15 uses T45). Keep every theorem; remove nothing.
2. In `definition/proof.md`, retype G24 and G27 so their "where" column says the constructions are
   benchmarks and the open question is the one the probes carry (G3/G26: η's estimate and form); leave
   G25/O19 to B19.
3. `tools/audit.py`: tag the benchmark theorems as not load-bearing (as B09 did for T21), checked by
   `tests/test_audit.py`.
4. One line in `docs/STATE.md`'s log. Then stop extending this line: no T66.

## Do not

- Delete or rewrite any theorem's content; the no-go results (T56, T58, T63) are findings and stay.
- Build a further anchor, exponent or state class.
- Touch the production equations or any premise row.
