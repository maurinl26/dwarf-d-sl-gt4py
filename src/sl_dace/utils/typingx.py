# -*- coding: utf-8 -*-
import logging
from typing import TypeVar, Union

import numpy as np

try:
    import cupy as cp
except ImportError:
    cp = None


T = TypeVar("T")

Pair = tuple[T, T]
Triple = tuple[T, T, T]

if cp is not None:
    ArrayLike = Union[np.ndarray, cp.ndarray]
else:
    ArrayLike = np.ndarray

import os

import dace

try:
    precision = os.environ["PRECISION"]
except KeyError as ke:
    logging.error(f"{ke}")
    precision = "double"

match precision:
    case "simple":
        dtype_float = dace.float32
        dtype_int32 = dace.int32
    case "double":
        dtype_float = dace.float64
        dtype_int = dace.int64
    case _:
        dtype_float = np.float64
        dtype_int = np.int64

