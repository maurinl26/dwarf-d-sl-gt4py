import numpy as np

p_1 = lambda l: (1 / 6) * l * (l - 1) * (2 - l)
p0 = lambda l: (1 - l**2) * (1 - l / 2)
p1 = lambda l: (1 / 2) * l * (l + 1) * (2 - l)
p2 = lambda l: (1 / 6) * l * (l**2 - 1)


def interpolate_linear_2d(
    lx: np.float64,
    ly: np.float64,
    ii: np.int64,
    jj: np.int64,
    field: np.ndarray,
    bc_kind: int,
):
    """Interpolate sequentially on x axis and on y axis

    Args:
        l (np.float64): _description_
        ii (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    p0 = lambda l: 1 - l
    p1 = lambda l: l

    px = np.array([p0(lx), p1(lx)])
    py = np.array([p0(ly), p1(ly)])

    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    return np.dot(np.matmul(padded_field[ii + 1 : ii + 3, jj + 1 : jj + 3], px), py)


def interpolate_cubic_2d(
    lx: np.float64,
    ly: np.float64,
    ii: np.int64,
    jj: np.int64,
    field: np.ndarray,
    bc_kind: int,
) -> np.float64:
    """Interpolate sequentially on x axis and on y axis.

    Args:
        lx (np.float64): _description_
        ly (np.float64): _description_
        ii (int): _description_
        jj (int): _description_
        field (np.ndarray): _description_

    Returns:
        _type_: _description_
    """

    # Padding on interpolation field
    # Fixed
    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    # Polynomes de lagrange d'ordre 3
    p_1 = lambda l: (1 / 6) * l * (l - 1) * (2 - l)
    p0 = lambda l: (1 - l**2) * (1 - l / 2)
    p1 = lambda l: (1 / 2) * l * (l + 1) * (2 - l)
    p2 = lambda l: (1 / 6) * l * (l**2 - 1)

    px = np.array([p_1(lx), p0(lx), p1(lx), p2(lx)])
    py = np.array([p_1(ly), p0(ly), p1(ly), p2(ly)])

    assert round(sum(px)) == 1
    assert round(sum(py)) == 1

    return np.dot(np.matmul(padded_field[ii : ii + 4, jj : jj + 4].T, px), py)


def cubic_interpolation(
    lx: np.ndarray,
    ly: np.ndarray,
    ii: np.int64,
    jj: np.int64,
    field: np.int64,
    bc_kind: int,
):
    if bc_kind == 0:
        padded_field = np.pad(field, (1, 2), "edge")

    # Periodic
    else:
        padded_field = np.pad(field, (1, 2), "wrap")

    py_1 = (
        p_1(lx) * padded_field[ii, jj]
        + p0(lx) * padded_field[ii + 1, jj]
        + p1(lx) * padded_field[ii + 2, jj]
        + p2(lx) * padded_field[ii + 3, jj]
    )

    py0 = (
        p_1(lx) * padded_field[ii, jj + 1]
        + p0(lx) * padded_field[ii + 1, jj + 1]
        + p1(lx) * padded_field[ii + 2, jj + 1]
        + p2(lx) * padded_field[ii + 3, jj + 1]
    )

    py1 = (
        p_1(lx) * padded_field[ii, jj + 2]
        + p0(lx) * padded_field[ii + 1, jj + 2]
        + p1(lx) * padded_field[ii + 2, jj + 2]
        + p2(lx) * padded_field[ii + 3, jj + 2]
    )

    py2 = (
        p_1(lx) * padded_field[ii, jj + 3]
        + p0(lx) * padded_field[ii + 1, jj + 3]
        + p1(lx) * padded_field[ii + 2, jj + 3]
        + p2(lx) * padded_field[ii + 3, jj + 3]
    )

    p = py_1 * p_1(ly) + py0 * p0(ly) + py1 * p1(ly) + py2 * p2(ly)

    return p