"""Shared pytest settings."""


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: takes minutes; runs only when the environment variable MET_SLOW is set")
