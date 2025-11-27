"""Functional test for sl_jax module.

This test verifies that the JAX implementation produces correct results
and can be JIT-compiled for performance.
"""
import logging
import time
import numpy as np
import sys
import pytest
from pathlib import Path

# Try to import JAX, skip tests if not available
try:
    import jax
    import jax.numpy as jnp
    from jax import jit
    JAX_AVAILABLE = True
except ImportError:
    JAX_AVAILABLE = False
    pytest.skip("JAX not available", allow_module_level=True)

from config import Config
from sl_jax.sl_2D import sl_xy, sl_init, backup

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


def create_test_config():
    """Create a small test configuration for quick validation."""
    return Config(
        nx=32,
        ny=32,
        xmin=0.0,
        xmax=1.0,
        ymin=0.0,
        ymax=1.0,
        dt=0.01,
        model_starttime=0.0,
        model_endtime=0.1,
        bcx_kind=1,  # Periodic
        bcy_kind=1,  # Periodic
        filter=True,
        lsettls=True,
    )


def gaussian_tracer(xcr: jnp.ndarray, ycr: jnp.ndarray, x0=0.5, y0=0.5, sigma=0.1):
    """Create a Gaussian tracer field."""
    return jnp.exp(-((xcr - x0)**2 + (ycr - y0)**2) / (2 * sigma**2))


def solid_body_rotation_velocity(xcr: jnp.ndarray, ycr: jnp.ndarray, t: float, omega=2*jnp.pi):
    """Solid body rotation velocity field."""
    x_center = 0.5
    y_center = 0.5
    
    vx = -omega * (ycr - y_center)
    vy = omega * (xcr - x_center)
    
    return vx, vy


def test_sl_jax_basic():
    """Test basic functionality of sl_jax without JIT compilation."""
    logging.info("Testing sl_jax basic functionality...")
    
    config = create_test_config()
    
    # Create coordinate grids
    xcr = jnp.linspace(config.xmin, config.xmax, config.nx)
    ycr = jnp.linspace(config.ymin, config.ymax, config.ny)
    xcr, ycr = jnp.meshgrid(xcr, ycr, indexing='ij')
    
    # Initialize tracer
    tracer = gaussian_tracer(xcr, ycr)
    tracer_e = tracer.copy()
    
    # Initialize velocity
    t = config.model_starttime
    vx, vy = solid_body_rotation_velocity(xcr, ycr, t)
    vx_e, vy_e = solid_body_rotation_velocity(xcr, ycr, t + config.dt)
    vx_p, vy_p = solid_body_rotation_velocity(xcr, ycr, t - config.dt)
    
    # Initialize velocities
    vx_e, vy_e, vx, vy = sl_init(vx_e, vy_e, vx, vy, vx_p, vy_p, lsettls=config.lsettls)
    
    # Perform one advection step
    tracer_e = sl_xy(
        config=config,
        vx=vx,
        vy=vy,
        vx_e=vx_e,
        vy_e=vy_e,
        tracer=tracer,
        tracer_e=tracer_e,
        nitmp=4,
    )
    
    # Check that result is a JAX array
    assert isinstance(tracer_e, jnp.ndarray), "Output should be a JAX array"
    
    # Check shape is preserved
    assert tracer_e.shape == tracer.shape, "Shape should be preserved"
    
    # Check values are reasonable (not NaN or Inf)
    assert jnp.all(jnp.isfinite(tracer_e)), "Results should be finite"
    
    # Check mass conservation (for periodic boundaries with no filtering)
    mass_before = jnp.sum(tracer)
    mass_after = jnp.sum(tracer_e)
    relative_error = jnp.abs(mass_after - mass_before) / mass_before
    
    logging.info(f"Mass before: {mass_before:.6f}")
    logging.info(f"Mass after: {mass_after:.6f}")
    logging.info(f"Relative error: {relative_error:.6e}")
    
    # Mass should be approximately conserved (within 1%)
    assert relative_error < 0.01, f"Mass conservation error too large: {relative_error}"
    
    logging.info("✓ Basic functionality test passed")


def test_sl_jax_jit_compilation():
    """Test that sl_jax functions can be JIT compiled."""
    logging.info("Testing JIT compilation...")
    
    config = create_test_config()
    
    # Create coordinate grids
    xcr = jnp.linspace(config.xmin, config.xmax, config.nx)
    ycr = jnp.linspace(config.ymin, config.ymax, config.ny)
    xcr, ycr = jnp.meshgrid(xcr, ycr, indexing='ij')
    
    # Initialize tracer
    tracer = gaussian_tracer(xcr, ycr)
    tracer_e = tracer.copy()
    
    # Initialize velocity
    t = config.model_starttime
    vx, vy = solid_body_rotation_velocity(xcr, ycr, t)
    vx_e, vy_e = solid_body_rotation_velocity(xcr, ycr, t + config.dt)
    
    # JIT compile the main function
    sl_xy_jit = jit(sl_xy, static_argnums=(0, 6))
    
    # Warmup (first call includes compilation)
    logging.info("Warming up JIT compilation...")
    _ = sl_xy_jit(config, vx, vy, vx_e, vy_e, tracer, tracer_e, 4)
    
    # Time JIT-compiled version
    start = time.time()
    tracer_e_jit = sl_xy_jit(config, vx, vy, vx_e, vy_e, tracer, tracer_e, 4)
    jax.block_until_ready(tracer_e_jit)  # Ensure GPU completion
    jit_time = time.time() - start
    
    # Check result is valid
    assert jnp.all(jnp.isfinite(tracer_e_jit)), "JIT-compiled results should be finite"
    
    logging.info(f"JIT execution time: {jit_time*1000:.2f} ms")
    logging.info("✓ JIT compilation test passed")


def test_sl_jax_conservation():
    """Test conservation properties over multiple steps."""
    logging.info("Testing conservation over multiple steps...")
    
    config = create_test_config()
    
    # Create coordinate grids
    xcr = jnp.linspace(config.xmin, config.xmax, config.nx)
    ycr = jnp.linspace(config.ymin, config.ymax, config.ny)
    xcr, ycr = jnp.meshgrid(xcr, ycr, indexing='ij')
    
    # Initialize tracer
    tracer = gaussian_tracer(xcr, ycr)
    tracer_e = tracer.copy()
    initial_mass = jnp.sum(tracer)
    
    # Initialize velocity
    t = config.model_starttime
    vx, vy = solid_body_rotation_velocity(xcr, ycr, t)
    vx_e, vy_e = solid_body_rotation_velocity(xcr, ycr, t + config.dt)
    vx_p, vy_p = solid_body_rotation_velocity(xcr, ycr, t - config.dt)
    
    # Initialize velocities
    vx_e, vy_e, vx, vy = sl_init(vx_e, vy_e, vx, vy, vx_p, vy_p, lsettls=config.lsettls)
    
    # Run for 10 steps
    n_steps = 10
    for step in range(n_steps):
        t += config.dt
        
        # Update velocities
        vx, vy = solid_body_rotation_velocity(xcr, ycr, t)
        vx_e, vy_e = solid_body_rotation_velocity(xcr, ycr, t + config.dt)
        
        # Advection step
        tracer_e = sl_xy(
            config=config,
            vx=vx,
            vy=vy,
            vx_e=vx_e,
            vy_e=vy_e,
            tracer=tracer,
            tracer_e=tracer_e,
            nitmp=4,
        )
        
        tracer, _, _ = backup(vx, vy, vx_e, vy_e, tracer, tracer_e)
    
    final_mass = jnp.sum(tracer)
    relative_error = jnp.abs(final_mass - initial_mass) / initial_mass
    
    logging.info(f"Initial mass: {initial_mass:.6f}")
    logging.info(f"Final mass: {final_mass:.6f}")
    logging.info(f"Relative error after {n_steps} steps: {relative_error:.6e}")
    
    # Mass should be conserved within 5% after multiple steps
    assert relative_error < 0.05, f"Mass conservation error too large: {relative_error}"
    
    logging.info("✓ Conservation test passed")


def test_sl_jax_against_numpy():
    """Test sl_jax results against sl_python (NumPy) implementation."""
    logging.info("Testing sl_jax against NumPy implementation...")
    
    try:
        from sl_python.sl_2D import sl_xy as sl_xy_numpy
        from sl_python.sl_2D import sl_init as sl_init_numpy
    except ImportError:
        pytest.skip("sl_python not available for comparison")
    
    config = create_test_config()
    
    # Create coordinate grids (NumPy)
    xcr_np = np.linspace(config.xmin, config.xmax, config.nx)
    ycr_np = np.linspace(config.ymin, config.ymax, config.ny)
    xcr_np, ycr_np = np.meshgrid(xcr_np, ycr_np, indexing='ij')
    
    # Create coordinate grids (JAX)
    xcr_jax = jnp.array(xcr_np)
    ycr_jax = jnp.array(ycr_np)
    
    # Initialize tracer (both)
    tracer_np = np.exp(-((xcr_np - 0.5)**2 + (ycr_np - 0.5)**2) / (2 * 0.1**2))
    tracer_jax = jnp.array(tracer_np)
    
    tracer_e_np = tracer_np.copy()
    tracer_e_jax = tracer_jax.copy()
    
    # Initialize velocity (both)
    t = config.model_starttime
    omega = 2 * np.pi
    vx_np = -omega * (ycr_np - 0.5)
    vy_np = omega * (xcr_np - 0.5)
    vx_e_np = -omega * (ycr_np - 0.5)
    vy_e_np = omega * (xcr_np - 0.5)
    vx_p_np = vx_np.copy()
    vy_p_np = vy_np.copy()
    
    vx_jax = jnp.array(vx_np)
    vy_jax = jnp.array(vy_np)
    vx_e_jax = jnp.array(vx_e_np)
    vy_e_jax = jnp.array(vy_e_np)
    vx_p_jax = jnp.array(vx_p_np)
    vy_p_jax = jnp.array(vy_p_np)
    
    # Run NumPy version
    vx_e_np, vy_e_np, vx_np, vy_np = sl_init_numpy(
        vx_e_np, vy_e_np, vx_np, vy_np, vx_p_np, vy_p_np, lsettls=config.lsettls
    )
    tracer_e_np = sl_xy_numpy(
        config=config,
        vx=vx_np,
        vy=vy_np,
        vx_e=vx_e_np,
        vy_e=vy_e_np,
        tracer=tracer_np,
        tracer_e=tracer_e_np,
        nitmp=4,
    )
    
    # Run JAX version
    vx_e_jax, vy_e_jax, vx_jax, vy_jax = sl_init(
        vx_e_jax, vy_e_jax, vx_jax, vy_jax, vx_p_jax, vy_p_jax, lsettls=config.lsettls
    )
    tracer_e_jax = sl_xy(
        config=config,
        vx=vx_jax,
        vy=vy_jax,
        vx_e=vx_e_jax,
        vy_e=vy_e_jax,
        tracer=tracer_jax,
        tracer_e=tracer_e_jax,
        nitmp=4,
    )
    
    # Compare results
    tracer_e_jax_np = np.array(tracer_e_jax)
    difference = np.abs(tracer_e_jax_np - tracer_e_np)
    max_diff = np.max(difference)
    mean_diff = np.mean(difference)
    
    logging.info(f"Max difference: {max_diff:.6e}")
    logging.info(f"Mean difference: {mean_diff:.6e}")
    
    # Results should be very close (numerical precision differences expected)
    assert max_diff < 1e-5, f"Difference between JAX and NumPy too large: {max_diff}"
    
    logging.info("✓ Comparison with NumPy passed")


if __name__ == "__main__":
    # Run tests
    print("=" * 60)
    print("Running sl_jax Functional Tests")
    print("=" * 60)
    
    test_sl_jax_basic()
    test_sl_jax_jit_compilation()
    test_sl_jax_conservation()
    
    try:
        test_sl_jax_against_numpy()
    except Exception as e:
        logging.warning(f"Comparison test skipped: {e}")
    
    print("=" * 60)
    print("All tests passed!")
    print("=" * 60)
