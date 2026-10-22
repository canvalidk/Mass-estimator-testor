# Butler (2004): bounded mass-candidate audit

Date: 2026-09-08. Status: local evidence note, awaiting transfer; NOT submitted to VD-docs. Scope: one primary publication and the supplied comparison contract. No managed VD-docs source is interpreted here; the combined review records the managed target snapshot.

Source: Hans Butler, *Controlling a position of a mass, especially in a lithographic apparatus*, EP1455231A2, published 2004-09-08. [Primary patent text](https://patents.google.com/patent/EP1455231A2/en).

## Published prescription

Sections 2.1-2.2 estimate mass from controller force and acceleration derived from position, using recursive least squares:

$$
\widehat m_k=\widehat m_{k-1}+\Gamma_k a_k
(F_k-\widehat m_{k-1}a_k).
$$

Section 3.1 gives the batch normal equation

$$
A^T A\widehat m=A^TF.
$$

The estimator retains signed acceleration-force products. Section 2.3 adds high-pass filtering for force offsets. Section 2.4 switches adaptation off at zero setpoint acceleration, or alternatively limits the adaptation-gain trace. These are stated engineering provisions. The audited equations specify no full joint measurement covariance input or mass uncertainty distribution. The recursive estimator depends on previous state; retaining that state can retain a finite estimate during zero readings.

## Our input mapping and counterexamples

The batch normal equation can be tested without confusing time samples with spatial directions: set the scalar sample entries to the corresponding components of the two given vectors, using equal weights. This is an algebraic specialization, not a claim that the patent originally performs vector measurement inference. Where $a^Ta>0$, the solution is

$$
\widehat m=\frac{a^TF}{a^Ta}.
$$

For measured vectors in two coordinates, give every example joint error covariance $I_4$, which is positive definite:

| Measured inputs | Consequence of the normal equation |
|---|---|
| $a=(1,0), F=(1,0)$ | $\widehat m=1$ |
| $a=(1,0), F=(-1,0)$ | $\widehat m=-1$ |
| $a=(1,0), F=(0,1)$ | $\widehat m=0$ |
| $a=(1,0), F=(0,0)$ | $\widehat m=0$ |
| $a=(0,0)$, any finite $F$ | The normal equation is $0=0$, giving no unique mass; the quotient is undefined. |

No assumption of an exact, error-free contradiction is being made: these are finite noisy measurements within the test domain. For the continuous path $a=(t,0), F=(1,0)$, the nonzero-$t$ estimate is $1/t$, with a singularity and sign change at zero. Positive weights or an ordinary forgetting factor cannot make the estimate positive when every nonzero sample satisfies $F_i=-a_i$.

Keeping an initialized recursive value at zero is not itself a failure of finiteness. It does mean that the result requires historical information in addition to the two current measured vectors; the mass uncertainty requested by the contract remains unspecified. The batch counterexamples must not be misreported as every recursive implementation literally dividing by zero.

## Verdict and contribution

**Demonstrated limitation of the published mass-estimation prescription under the comparison contract.** Direction/sign information is retained, but an everywhere finite, strictly positive scale and its uncertainty are not supplied. The native controller uses additional time history and a stated adaptation gate. Adding positive priors, bounded updates, a joint error model, or uncertainty propagation would create a separately auditable adaptation; no impossibility claim is made about those adaptations. Exact identity with the workspace equation is not established.

What this adds: a direct force-acceleration precedent and a clear distinction between operationally holding an estimate through a zero phase and performing uncertainty-aware inference from a zero measurement. The signed least-squares counterexample itself repeats an already established failure mechanism, so more angle/error sweeps on this prescription would have low value.

Access: primary HTML sections 2.1-2.4, 2.10-2.11, 3.1-3.2 and 4 read. Text-extracted normal and recursive mass equations support the result. Image-only gain/forgetting-factor formulas were not visually recovered and are not transcribed here. No legal-status or patent-scope conclusion is intended.
