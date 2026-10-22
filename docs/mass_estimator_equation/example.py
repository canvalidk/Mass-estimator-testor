"""Run with: python example.py (from this folder). Inputs below are SI."""
import numpy as np
from mass_estimator import estimate_mass

force = [10.0, 0.0, 0.0]         # N: measured net force vector
acceleration = [2.0, 0.1, 0.0]   # m/s^2: measured acceleration vector
sigma_force = 0.5               # N: per-coordinate standard uncertainty
sigma_acceleration = 0.1        # m/s^2: per-coordinate standard uncertainty
covariance = np.diag([sigma_force**2] * 3 + [sigma_acceleration**2] * 3)

result = estimate_mass(force, acceleration, covariance)
print(f"Mass point: {result.mass:.6f} kg")
print(f"95% conditional interval: [{result.interval[0]:.6f}, {result.interval[1]:.6f}] kg")
print(f"SD of log mass: {result.log_standard_deviation:.6f}")
print(f"Inverse mass point: {result.inverse_mass:.6f} 1/kg")
print(f"Last numerical refinement change: {result.numerical_change:.3g}")
