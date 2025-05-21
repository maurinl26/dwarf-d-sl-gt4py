# sl-gt4py

Semi-Lagrangian advection schemes with python and DaCe

## sl_python

sl_python implements the classical sl scheme in 2d.

## sl_dace

- interpolation :
  - interpolation_2d.py : 2d linear interpolation for classical SL
  - flux_integral.py : flux integrals along the trajectory for FFSL
  
- stencils :
  - ppm.py : ppm reconstruction and limiter 
  - ffsl.py : stencils for 1d FFSL
  - dep_search_1d.py : depature search for classical SL

- ffsl_x.py / ffsl_y.py : ffsl 1d orchestration
- (WIP) ffsl_xy.py : ffsl 2d with swift splitting 
- elarche.py : departure search for classical SL
- (WIP) sl_init.py : SETTLS or NESC init orchestration
- sl_xy.py : classical SL orchestration
- (WIP) sl_driver.py : driver for full sl or ffsl schemes in dace

## tests

Tests are build using pytests.
All required test objects (grids, dimensions, backends, etc) are implemented as fixtures in conftest.py

To run the tests :

```bash
    pytest
```

Tests structure :
- unit : for unit tests / stencils compilation / dace generation
- functional : 
  - test_blossey.py : implements blossey 2d advection test on a plate
  - test_one_step.py : implements 1 step run for an advection scheme
  - test_uniform.py : implements uniform velocity advection

## Setup

The project dependencies are managed with uv

```bash
    uv init                     # init project
    uv venv                     # create virtual environment
    source .venv/bin/activate   # activate virtual environment
    uv sync                     # load and synchronize project dependencies
```

If you don't have uv, it can be downloaded from :

```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Build the doc 

```bash
   uv run sphinx-autobuild -M html docs/ docbuild/
```