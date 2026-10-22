# Why mass can't be undefined

**2026-09-29.** Record, written by Claude at the user's request ("take this
argument. Think about it thoroughly and make a record on it"). The argument
is his, quoted verbatim in §7. The reconstruction (§1), the sharpening
(§§2–3) and the consequences (§4) are Claude's.

**What it records.** The argument behind his ruling of 29 Sept on O1 of
`mass_estimator_lemmas_2026-09-29.md`: on measured inputs the definition never
returns "undefined", including when every reading is zero–zero. From the
ruling: "the information we gather from 0/0 doesn't rule any mass value out,
meanwhile 0/a or F/0 could".

**Relations.** It scopes the records named here and binds nothing. It answers Case I of
`mass_estimation_zero_zero_two_cases.md` and reaches Case II's conclusion by a
different route. It gives the argument for the O1 ruling and a forced value
for O1/O2 of the lemma record (§4). On measured inputs only, it displaces the
"mass undefined from this trial" clause of
`disagreement_rule_mass_estimator_trace.md` §5 and the spec d1. T1, the exact
branch, is untouched.

---

## 0. The claim

On every measured input, meaning every declared $\sigma>0$ (L2), the
definition of mass returns its full output: a law of mass on $(0,\infty)$ and
a point $\widehat m\in(0,\infty)$. It never returns "undefined". That includes
the **all-zero experiment**, a series in which every reading is
$\tilde{\mathbf F}=\mathbf 0,\ \tilde{\mathbf a}=\mathbf 0$.

## 1. His argument, in steps

- **A1.** The estimator is the definition (F1).
- **A2.** Other routes give the same object a mass at the same moment: a
  spring; a density route (photograph the object, get its volume and
  material, compute 5 ± 0.5 kg); weighing (named in his ruling).
- **A3.** "What did we just calculate?" If it was mass, then the definition's
  "undefined" and the density route's 5 ± 0.5 kg disagree. A definition
  cannot disagree with a determination of the thing it defines.
- **A4.** Disagreement is not tested with `==`. Two good solar-mass
  measurements, 0.95 ± 0.5 and 1 ± 0.1, fail `==` too, and nobody calls
  Newtonian mechanics contradictory. Contradiction is judged by a
  compatibility relation, `==^`.
- **A5.** We are free to define `==^` so that `undefined ==^ 5 ± 0.5` passes.
  The price is that every equation using $m$ must accept $m=$ undefined (what
  is kinetic energy $=$ undefined?), plus rules for when the
  undefined-accepting versions apply. That work needs a purpose.
- **A6.** No purpose is visible. So mass is not undefined.

He left one door open: "maybe this is the answer. I didn't try doing all
that." §2.3 shows what trying finds.

## 2. Why it holds

### 2.1 There is no per-moment mass to be undefined

The standing reply is Case I §1 of the 5 Sept record: "We need not claim that
the object has an undefined physical property. We claim that this route has
not supplied a value for it." This is the retreat A3 anticipates ("it's
undefined for exactly this moment"). Two things close it.

- Under F1 there is no mass of the object sitting behind the routes, for a
  route to fail to supply ("there is no super meta mass"). What the
  definition outputs *is* the mass.
- Mass is not indexed by moment. By L1 one $m$ is shared by every reading of a
  series, and by T24 a null reading's factor is constant. A stationary moment
  inside a series is pooled and contributes nothing. The definition never
  asks it a mass question of its own.

So "undefined" could only live in a series whose readings are all null: the
all-zero experiment, which is O1's case.

### 2.2 What the density route calculated

Either answer to A3's question gives a number.

- **It calculated mass.** The route is lumping (L5) over volume elements,
  plus a homogeneity premise (equal volumes of one material have equal mass)
  and the persistence of the calibrating sample's mass (O10). Every step
  belongs to the definition. So an "undefined" would be the definition
  contradicting itself through another of its own routes, not the definition
  against an outside method.
- **It calculated something else, $Q$.** Then Newtonian practice needs a
  bridge $m=g(Q)$, because density measurements are used to get masses. The
  bridge gives a number at that moment, and "undefined" contradicts the
  bridge instead.

The three examples are not equally strong. The spring needs its own push, so
it answers a different trace and needs persistence to reach the stationary
moment. Weighing reads the object as it stands (a stationary object on a
scale is itself a $0/0$ net-force trace, read off one member of the force
sum, $W=mg$). But it needs inertial mass to equal gravitational mass, which
is a separate law in NM, so a defender can push it to the second answer.
The density route reads the object as it stands and uses only lumping. It is
the example that carries the argument.

### 2.3 The cost step, made exact

`==^` *compares*. The definition's own relation to other information is
**combination**. Readings join by multiplying factors (T12), and other
knowledge of the object enters the prior slot
(`old_information_in_the_prior_slot_2026-09-24.md`; L3′, proposed). A null
series' factor is constant (T24), so it is the identity of combination. The
definition therefore already represents "this rules nothing out", and it
represents it as a value, not a gap.

An "undefined" output would have to behave in one of three ways when combined
with other information:

1. **Absorbing**, as NaN is: combine anything with it and get undefined.
   Watching the body sit still would erase the mass the density route gave
   it. This contradicts 2.2.
2. **Identity**: combining with it changes nothing. That is the constant
   factor the null series already contributes, under a new name. The special
   `==^` and the $m=$ undefined versions of every equation are then the price
   of the renaming, and they buy nothing. Read this way, "kinetic energy $=$
   undefined" only means "kinetic energy under the unchanged law", which
   ordinary propagation already gives.
3. **Neither**: it favours or excludes some masses with no data doing so.
   That is arbitrary.

So an undefined output is contradictory, redundant or arbitrary. That is his
"ok buddy" detector turned into a proof, and it is what trying A5's extension
finds. Case I takes the identity reading when another route exists ("any
later operation requiring that value must obtain it through another licensed
route") and halts when none does. The halt is answered next.

### 2.4 With nothing else known, the law is already there

Case I's remaining position is to halt if no other route has been run. But at
the all-zero experiment the definition's law is the declared starting law
$\pi$ (T24), which is proper on $(0,\infty)$. So the definition has already
said what the mass is, given this information: it is distributed as $\pi$.
An undefined point beside that law withholds a summary of a law already
stated, and every proper law has one. It would also make the output type
depend on the data: a law and a point for every input but one, and a law and
nothing at that one. That is A5's cost in type form, since every consumer
would have to branch.

### 2.5 The best case for undefined, and why it fails

At a standstill $\mathbf p=m\mathbf v=\mathbf 0$, $K=0$ and
$m\mathbf a=\mathbf 0$, whatever $m$ is. Mass looks idle at that moment, so
"undefined" looks harmless. But $0/0$ says $\mathbf F=\mathbf 0$ and
$\mathbf a=\mathbf 0$, not $\mathbf v=\mathbf 0$. A body coasting at
constant velocity is a $0/0$ trace, and it carries momentum $m\mathbf v$ that
it will deliver to the next wall. A body at rest on a table presses on it
with $mg$. Mass is not idle at a null moment. That trace simply does not
exercise it.

## 3. The trap on the other side

Returning a number is not enough on its own. The all-zero law $\pi$ has
finite central intervals: for the half-Cauchy, 95% lies in
$[0.039c,\ 25.45c]$ (Case II §9). If `==^` compared those intervals, the
all-zero result would "contradict" a density reading of 5000 kg when
$c=1$ kg. That is A3's contradiction again, at a different number.

Those intervals belong to the declared starting law, not to the trial. The
trial excludes nothing: its profile compatibility set is $(0,\infty)$ (Case I
§8). So the definition's relation to another route has to be combination
(2.3), or a comparison of what the data exclude. It must never be a
comparison against the starting law's width. Case II §9 already requires the
output to keep this distinction, and his argument depends on it.

## 4. What it forces now

**O1 is decided.** O1 names three things that cannot all stand. His ruling
keeps the readout at the all-zero experiment. The best-guess repair of T29
is out as it stands, because it returns $0/0$ there.

**The value is forced to be $c$.** Z at the point (L20) lets any null reading
be dropped or added without changing the answer. The combination is
order-free (T12), so the all-zero answer depends neither on the number of
readings nor on their $\Sigma_i$. It is a function of the declared starting
law (and the declared shapes) alone. Its only mass scale is $c$, so by units
(L3) it is $kc$. The swap maps $\pi$ to itself in inverse-mass units (L4). If
the point keeps the swap property the readout has elsewhere (T5), then
$kc=c^2/(kc)$, so $k=1$:

$$
\widehat m(\text{all-zero})=c,
$$

the centre of the declared starting law, whatever $c$ is declared to be.
T14's value $s$ is the case $c=s$ (L3′, unruled). The old cloud readout
(T10) gave this value at $c=s$, and fails Z elsewhere.

**The best-guess repair cannot be patched to reach $c$.** Its weight is
homogeneous of degree one in the data, and the law tends to $\pi$ as the data
shrink. So as the readings shrink to zero along a fixed direction, the
readout tends to a limit that depends on the direction.

| 1D, unit noise, $c=1$, readings $\varepsilon(\tilde F,\tilde a)$, $\varepsilon\to0$ | limit of $\widehat m$ |
|---|---|
| force only | $\pi/2\approx1.571$ |
| acceleration only | $2/\pi\approx0.637$ |
| equal sizes, aligned or anti-aligned | $1$ |
| over all directions | every value from $0.546$ to $1.830$ (extremes at anti-aligned readings with sizes about 3.4 : 1; their product is $1$, by swap) |

The table is the same under ours27 and ours28, because the size law's factor
tends to a constant. So "return $c$ at the all-zero experiment" would be a
bolted-on clause. Just off zero, a reading carrying almost no information
would move the point by up to a factor of 1.83, while moving the law only by
$O(\varepsilon^2)$. That does not breach Z as ruled, since Z speaks of exact
nulls, but it is Z's intent failing at its limit.

**O2, sharpened.** What is wanted is a readout that composes (L5, L14),
passes Z, and returns $c$ at the all-zero experiment. Whether it must also be
continuous there is a question for him.

## 5. What it does not do

- **Exact inputs.** T1 ($\sigma=0$) is untouched. L2 rules $\sigma>0$, and his
  ruling adds that uncertainties cannot be zero, so no measured input reaches
  T1. The argument would run the same way there, but that is not ruled.
- **Rejection.** There is no gate (L1, ruled 7 Sept), so "never undefined"
  does not compete with a rejection output.
- **The choice of $c$** (L3′, O11) and **the readout** (O2) stay open.
- **Persistence.** The density route leans on O10 across series. Within a
  series, L1 is enough.

## 6. Relations, in detail

- `mass_estimation_zero_zero_two_cases.md` (5 Sept). Answers Case I (§§1, 7,
  8). Reaches Case II's conclusion by a different route: the plurality of
  routes and the cost of an undefined value, not continuity. Case II §6's
  continuity holds for the cloud readout, not for the best-guess repair (§4).
  The ours26 sweep had already marked Case I broken, against L2.
- `mass_estimator_lemmas_2026-09-29.md`. The argument for his O1 ruling. §4's
  value and limits are offered for O2, not placed. L2's "as d2 reads it" is
  now his ruling.
- `disagreement_rule_mass_estimator_trace.md` §5, and the spec d1's "mass
  undefined from this trial": displaced on measured inputs; exact inputs
  untouched.

Nothing is bound.

## 7. His words, 29 Sept (verbatim)

> You might be missing my point. Pay attention. Let me set up the
> contradiction. If I say that 'the estimator is the definition!' and the
> estimator gives out undefined, here is what could happen. First you argue,
> no its undefined for exactly this moment. I say, it is possible to
> determine the mass for that moment using the spring. But other mass
> determining processes exist also for that material at that time. For
> example, we use a density method. We take pictures of the object, assess
> its volume and material, and determine the 'mass' value. Now the question
> is, what did we just calculate? Did we calculate mass? If we did, we have a
> seeming contradiction. One method of determining mass gives us:
> 'undefined'. The other gives us '5 +- 0.5 kg'. Now if we do 'undefined ==
> 5 +- 0.5 kg' we will get False if I'm correct. And that's where we have a
> problem. If both methods really determine 'mass' and the estimator is the
> ultimate definition, we can't have them disagree. Because that would be a
> contradiction. Now there is an alternative case that clears this up. Say we
> make a measurement of the sun. We make 2 measurements: "0.95 +-0.5 solar
> mass" and "1 +-0.1 solar mass". Now if we run the same python equality
> test, these results will also contradict. It'll give us 'False'. So is
> Newtonian mechanics contradictory? No of course we don't accept this
> contradiction. We say there is some different process "==^" that
> determines whether there is a contradiction or not. Now we need the
> undefined mass to pass this particular type of equality. And we have
> freedom to define it. But if we do that especially and concoct some
> sensile equality condition for specifically the undefined mass case, that
> gives us SOOOO MUCH MORE work to do. Because now we need all our equations
> to accept m = undefined as well as understanding the rules for when to pick
> that equivalence over the non-undefined accepting versions. And we need to
> do all that for a fucking reason. It needs to have some point. Some purpose
> that it would be that way. What would it gain us to hypothetically figure
> out the interpretation of kinetic energy = undefined... I mean maybe...
> maybe this is the answer. I didn't try doing all that. But doesn't your "ok
> buddy" detectors go off here?

## 8. Checks

§4's limits were computed with scipy quadrature at
$\varepsilon=10^{-1},10^{-2},10^{-3},10^{-4}$ (converged by $10^{-3}$), plus a
1441-angle direction scan. The models are 1D with unit noise and $c=1$, so
$\pi(m)=1/(1+m^2)$ up to a constant, and $H(m)=(m\tilde F+\tilde a)^2/(m^2+1)$.

- **ours27:** law $\pi\,e^{H/2}$; best-guess candidate
  $E[a^\ast\mid m]=(m\tilde F+\tilde a)/(m^2+1)$.
- **ours28** (arcsine size law, $k=1$): law
  $\pi\,{}_1F_1(\tfrac12;\tfrac32;\tfrac H2)$; candidate multiplied by
  $E[\beta\mid m]={}_1F_1(\tfrac32;\tfrac52;\tfrac H2)\big/3\,{}_1F_1(\tfrac12;\tfrac32;\tfrac H2)$.
- **Readout:** $\int m\,A\,p\big/\!\int A\,p$ with $A=\lvert E[a^\ast\mid m]\rvert$.

The swap check $\widehat m(\tilde F,\tilde a)\,\widehat m(\tilde a,\tilde F)=1$
held to $10^{-15}$. The scripts were not kept.

## Sources

- `mass_estimator_lemmas_2026-09-29.md`: F1, L1–L4, L20, T1, T5, T10, T12,
  T14, T24, T29, O1, O2, O10, O11.
- `mass_estimation_zero_zero_two_cases.md`: Case I §§1, 7, 8; Case II §§6, 9.
- `old_information_in_the_prior_slot_2026-09-24.md`.
- `disagreement_rule_mass_estimator_trace.md` §5.
- `ours26/ours26_sweep_of_the_flat_records_2026-09-27.md` (Case I marked
  broken).
- His ruling of 29 Sept on O1, and his argument of 29 Sept (§7).
