import pytest
from typing import Union
import numpy as np


def get_cpu_backends():
    return ["numpy", "dace:cpu"]

@pytest.fixture(name="dtypes", scope="module")
def dtypes_fixtures(precision: Union["single", "double"] = 'single'):
    if precision == "single":
        return {
            float: np.float32,
            int: np.int32,
        }
    elif precision == "double":
        return {
            float: np.float64,
            int: np.int64,
        }
    else:
        raise KeyError("precision not in 'single' or 'double'")

@pytest.fixture(name="origin_with_halo", scope="module")
def origin_with_halo_fixture(h: int = 5):
    return h, h, 0

@pytest.fixture(name="origin", scope="module")
def origin_fixture():
    return 0, 0, 0

@pytest.fixture(name="inner_domain_origin", scope="module")
def inner_domain_origin_fixture():
    return 1, 0, 0

@pytest.fixture(name="origin_1_padding", scope="module")
def o1_padding_origin_fixture():
    return 1, 1, 0

@pytest.fixture(name="origin_2_padding", scope="module")
def o2_padding_origin_fixture():
    return 2, 2, 0

@pytest.fixture(name="origin_5_padding", scope="module")
def o5_padding_origin_fixture():
    return 5, 5, 0

@pytest.fixture(name="domain", scope="module")
def domain(nx: int = 50, ny: int = 50, nz: int = 10):
    return nx, ny, nz

@pytest.fixture(name="inner_domain", scope="module")
def inner_domain_fixture(domain):
    return tuple(dim - 1 for dim in domain)

@pytest.fixture(name="large_domain", scope="module")
def large_domain_fixture(nx: int = 255, ny: int = 255, nz: int = 120):
    return nx, ny, nz
