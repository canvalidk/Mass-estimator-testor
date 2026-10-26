"""Composable case inspection: data, independent geometry, adapters, checks, reports.

The small functions in these modules can be used separately. No registration,
base class or fixed report layout is required to add an adapter or a check.
"""
from .case import Case, Inputs, Observation, load_cases

__all__ = ["Case", "Inputs", "Observation", "load_cases"]
