"""Independent numerical checks of finite identities; not a tail-divergence test."""
import json
import math


def integrate(fn, lo, hi, n=40000):
    h = (hi-lo)/n
    return h/3 * (fn(lo)+fn(hi)+sum((4 if i % 2 else 2)*fn(lo+i*h) for i in range(1,n)))


def close(actual, expected, tol=2e-9):
    assert abs(actual-expected) < tol*max(1, abs(expected)), (actual,expected)


c = math.sqrt(2/math.pi)
close(integrate(lambda x: x*c*math.exp(-x*x/2), 0, 12), c)
s = 4.0
quantiles = {}
for p in (.025, .16, .5, .84, .975):
    q = s*math.tan(math.pi*p/2)
    close(integrate(lambda m: 2*s/(math.pi*(s*s+m*m)), 0, q), p)
    quantiles[str(p)] = q/s

# Alternate null measure: chi_3 density c*x²*exp(-x²/2).
e_x = integrate(lambda x: c*x**3*math.exp(-x*x/2), 0, 12)
e_inv = integrate(lambda x: c*x*math.exp(-x*x/2), 0, 12)
e_x2 = integrate(lambda x: c*x**4*math.exp(-x*x/2), 0, 12)
e_inv2 = integrate(lambda x: c*math.exp(-x*x/2), 0, 12)
close(e_x*e_inv, 4/math.pi)
close(e_x2*e_inv2, 3)

# Resolved one-zero profile thresholds, with unequal channel scales.
k, sf, sa, A, F = 3.84, 2.0, .5, 5.0, 10.0
upper = math.sqrt(k)*sf/math.sqrt(A*A-k*sa*sa)
lower = math.sqrt(F*F-k*sf*sf)/(math.sqrt(k)*sa)
q_0a = lambda m: m*m*A*A/(sf*sf+m*m*sa*sa)
q_f0 = lambda m: F*F/(sf*sf+m*m*sa*sa)
close(q_0a(upper), k)
close(q_f0(lower), k)
assert q_0a(.5*upper) < k < q_0a(2*upper)
assert q_f0(2*lower) < k < q_f0(.5*lower)

print(json.dumps({"status":"passed", "half_normal_mean":c,
 "half_cauchy_quantiles_in_scale_units":quantiles,
 "chi3_ratio_mean_in_scale_units":e_x*e_inv,
 "chi3_ratio_second_moment_in_scale_squared_units":e_x2*e_inv2,
 "one_zero_thresholds":{"upper":upper,"lower":lower}}, indent=2))
