# sl_jax - JAX Implementation of Semi-Lagrangian Advection

This module is a JAX translation of the `sl_python` module, providing GPU-accelerated semi-lagrangian advection schemes.

## Overview

The `sl_jax` module implements 2D semi-lagrangian advection using JAX for automatic differentiation and GPU acceleration. The code is translated from the NumPy-based `sl_python` implementation.

## Key Features

- **JAX Arrays**: All operations use `jax.numpy` (jnp) instead of NumPy for GPU compatibility
- **Automatic Differentiation**: JAX enables gradient computation through the entire advection scheme
- **JIT Compilation**: Functions can be JIT-compiled for performance improvements
- **Vectorized Operations**: Leverages JAX's vectorization capabilities

## Module Structure

```
sl_jax/
├── __init__.py
├── sl_2D.py                    # Main 2D semi-lagrangian advection
├── boundaries.py               # Boundary condition handling
├── diagnostics.py              # Diagnostic functions (Lipschitz, overshoots)
├── filter.py                   # Over/undershoot filters (non-periodic)
├── periodic_filters.py         # Over/undershoot filters (periodic boundaries)
└── interpolation/
    ├── __init__.py
    └── interpolation_2d.py     # 2D interpolation (linear, cubic, min/max)
```

## Main Components

### sl_2D.py

Main module containing:
- `dep_search_1d()`: 1D departure point search
- `sl_init()`: Initialize velocity fields (LSETTLS or LNESC methods)
- `backup()`: Copy fields for next iteration
- `lagrangian_search()`: Iterative departure point search
- `sl_xy()`: Complete 2D semi-lagrangian advection

### interpolation/interpolation_2d.py

Interpolation functions:
- `interpolate_lin_2d()`: Bilinear interpolation
- `interpolate_cub_2d()`: Bicubic interpolation
- `max_interpolator_2d()`: Maximum value in stencil
- `min_interpolator_2d()`: Minimum value in stencil

### boundaries.py

- `boundaries()`: Apply boundary conditions (periodic or fixed)

### diagnostics.py

- `diagnostic_lipschitz()`: Compute Lipschitz stability condition
- `diagnostic_overshoot()`: Compute overshoots
- `diagnostic_undershoot()`: Compute undershoots

### filter.py & periodic_filters.py

Over/undershoot filtering for monotonicity preservation. Note: The filter implementations are simplified placeholders as the original sequential update logic is challenging to vectorize efficiently in JAX.

## Key Differences from sl_python

1. **Array Types**: `numpy.ndarray` → `jax.numpy.ndarray`
2. **Import Changes**: `import numpy as np` → `import jax.numpy as jnp`
3. **Copy Operations**: JAX arrays are immutable by default, `.copy()` creates new arrays
4. **Indexing**: JAX uses similar indexing to NumPy but with stricter rules for JIT compilation
5. **Control Flow**: Some sequential operations (filters) may need `jax.lax.fori_loop` or `jax.lax.scan` for full implementation

## Usage Example

```python
import jax.numpy as jnp
from config import Config
from sl_jax.sl_2D import sl_xy

# Initialize configuration
config = Config(...)

# Initialize fields (as JAX arrays)
vx = jnp.array(...)
vy = jnp.array(...)
tracer = jnp.array(...)

# Perform advection
tracer_new = sl_xy(
    config=config,
    vx=vx,
    vy=vy,
    vx_e=vx_e,
    vy_e=vy_e,
    tracer=tracer,
    tracer_e=tracer_e,
    nitmp=4
)
```

## JIT Compilation

For performance, functions can be JIT-compiled:

```python
from jax import jit

# JIT compile the advection function
sl_xy_jit = jit(sl_xy, static_argnums=(0, 6))  # config and nitmp are static

# Use JIT-compiled version
tracer_new = sl_xy_jit(config, vx, vy, vx_e, vy_e, tracer, tracer_e, nitmp)
```

## Limitations and TODOs

1. **Filter Functions**: The over/undershoot filter implementations in `filter.py` and `periodic_filters.py` are simplified placeholders. The original implementations use sequential updates that affect neighboring cells, which is challenging to vectorize in JAX while maintaining exact equivalence. Full implementation would require:
   - Using `jax.lax.fori_loop` or `jax.lax.scan` for sequential state updates
   - Careful handling of boundary conditions and neighbor dependencies
   - Potentially restructuring the algorithm for better parallelization

2. **Performance Optimization**: Additional optimizations may be needed:
   - Strategic use of `@jit` decorators
   - Handling of static vs traced arguments
   - Memory layout optimization for GPU execution

3. **Testing**: Comprehensive testing against the NumPy implementation is needed to verify correctness.

## Dependencies

- JAX
- jaxlib (for GPU support)
- The project's Config class

## Notes

- JAX arrays are immutable, so operations create new arrays rather than modifying in place
- For best performance, ensure JIT compilation is used for hot loops
- GPU execution requires jaxlib with CUDA support
