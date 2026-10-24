# Read-only audit of area_conditioning_checks.py

Independent mathematical and numerical-method review, 2026-09-15. Reviewed
the production research script and its imported gaussian_conditioning.summary
normalization. No reviewed code was changed or main run repeated.

## Verdict

The finite-band reduction, normalization constant, and channel-swap recovery
of E||Y|| are correct for the stated scope: three dimensions, independent
unit-covariance Gaussian channels, and collinear (possibly opposed) centers.
The finite-band norm ratio is correctly labelled a convergence diagnostic,
not a mass policy for non-collinear latent pairs. No blocking issue found.

## 1. Conditional reduction

Write X=x u, x>0, t=u_3, with X~N(F e_3,I_3) and Y~N(A e_3,I_3).
Conditional on X,

    Y_parallel ~ N(A t,1),
    Y_perpendicular ~ N(mu_perp,I_2),
    ||mu_perp||=|A| sqrt(1-t^2),

and these components are independent. Since

    X dot Y = x Y_parallel,
    ||X cross Y|| = x ||Y_perpendicular||,

the event is exactly Y_parallel>0 together with a perpendicular disk of
radius b=epsilon/x. Its conditional probability is

    Phi(A t) D(b,h),   h=|A| sqrt(1-t^2).

Integrating the azimuth of X gives its radial/angular volume-density factor

    (1/sqrt(2pi)) x^2 exp[-(x^2+F^2-2xFt)/2] dx dt.

After division by epsilon^2, x^2 D/epsilon^2 is D/b^2. This is exactly the
script's weighted integrand; neither a radial factor nor an angular factor
is missing.

## 2. Disk probability

For a two-dimensional unit-covariance Gaussian displaced by h,

    D(b,h)=integral_0^b r exp[-(r^2+h^2)/2] I_0(hr) dr.

Substitution r=bq gives

    D(b,h)/b^2=integral_0^1 q exp[-(b^2q^2+h^2)/2] I_0(hbq) dq,

which matches disk_probability_over_radius_squared. The 2pi from the angular
integral cancels the bivariate normal normalization.

If b>=h+12, the omitted disk complement is at most exp(-72), since
||Z+mu||>b implies ||Z||>b-h>=12 for Z~N(0,I_2). The script reports this
conditional absolute probability bound; it should not be described as a
certified relative error of the final ratio for arbitrary inputs.

Within the actual centers used here, the direct exponential/Bessel product
is safely within its declared numerical range. The argument guard does not
claim a universally stable noncentral-disk routine.

## 3. Limiting normalization

The imported summary() returns

    Z_flat=integral_0^infty integral_0^infty
             average_u exp[-Q(f,alpha,u)/2] df d(alpha),

using normalized uniform sphere measure. The original ambient 6D Gaussian
normalization is (2pi)^(-3); a small aligned angular cap has area
pi epsilon^2/(f alpha)^2; the remaining shared sphere area is 4pi.
Consequently

    lim P(event)/epsilon^2
      = (2pi)^(-3) pi (4pi) Z_flat
      = Z_flat/(2pi).

The factor used in run() is correct.

At measured zero-zero, Z_flat=pi/2, giving the limit 1/4. The exact finite
event probability is

    P(event)=[1-(1+epsilon)exp(-epsilon)]/2.

This independently agrees with the script's null oracle and its small-
epsilon limit.

## 4. Obtaining the other mean by swapping channels

The event is symmetric under X<->Y: the dot product and cross-product norm
are unchanged. Swapping the two Gaussian centers therefore leaves event
probability unchanged and turns E||X|| into E||Y|| under the original
conditioning. Thus the reverse call's mean_x is exactly the desired mean_y.

The separate radial cutoffs differ after swapping, so equality holds exactly
for the untruncated integrals. The recorded normalization discrepancy and
padding extension check appropriately diagnose the numerical approximation.

## 5. Radial truncation and numerical verification scope

With upper=|F|+padding, the excluded X probability is bounded by

    P(||Z_3||>padding)
      = erfc(padding/sqrt(2))
        + sqrt(2/pi) padding exp(-padding^2/2).

At padding=12 this is about 5.2e-31. Even after dividing by the smallest
tested epsilon^2=9e-6, the absolute probability error bound remains tiny.
A corresponding first-moment bound follows from

    E[||Z_3|| 1_{||Z_3||>p}]
      = sqrt(2/pi)(p^2+2)exp(-p^2/2)

and ||X||<=|F|+||Z_3||. The empirical 12-to-16 padding comparison is
therefore consistent with an independently controlled tail.

Quadrature-order refinement and channel-swap equality are numerical checks,
not uniform error guarantees. The exact null probability provides a separate
analytic oracle. The final 1% proximity assertion for the norm ratio is a
finite-tolerance convergence check, not an assertion that epsilon=0.003 is
the exact limiting law. The tabulated probabilities and means should likewise
remain labelled finite-tolerance values.

## 6. Noise-sampling identities

For the script's collinear true centers (.6 strength, .8 strength), an
orthogonal channel rotation puts all signal in one column and preserves
wedge area. Thus area^2 is a product of independent variables

    U~noncentral-chi-square_d(strength^2),
    V~chi-square_(d-1).

Its mean and variance are precisely those used in noise_checks(). At null,
the area itself has Gamma(d-1,scale=1) law, and changing Y to -Y proves the
positive-dot half-event has half its area-CDF probability. The d=2 and d=3
Gamma CDF implementations are correct.

The noise checks test properties of the ambient measurement experiment.
They do not turn the finite wedge tolerance into a selected estimator rule.
