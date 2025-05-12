from contextlib import contextmanager
import numpy as np
from ifs_physics_common.framework.grid import ComputationalGrid

# note : similar to ifs_physics_common
@contextmanager
def temporaries(
        temporary_fields: dict,
        grid: ComputationalGrid
):

    yield {
        name: np.zeros(grid.grids[field_desc["grid"]].shape, dtype=field_desc[type])
        for name, field_desc in temporary_fields.items()
    }