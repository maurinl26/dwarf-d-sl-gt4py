import dace
from typing import Callable, Literal


def build_sdfg(function: Callable, mode: Literal["aot", "jit"], device: Literal["gpu", "cpu"]):

    # dace
    sdfg = (
        dace.program(function)
        .to_sdfg()
    )
    match device:
        case "gpu":
            sdfg.apply_gpu_transformations()

    match mode:
        case "aot":
            compiled_sdfg = sdfg.compile()
            return compiled_sdfg
        case "jit":
            return sdfg