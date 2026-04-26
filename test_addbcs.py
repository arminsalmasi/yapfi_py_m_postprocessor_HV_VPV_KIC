import pytest
from addbcs import add_scalars_limits, add_coords_limits

np = pytest.importorskip("numpy")

def test_add_scalars_limits_1d():
    in_array = np.array([1, 2, 3])
    nxyz = [3]
    result = add_scalars_limits(in_array, nxyz)
    expected = np.array([1, 1, 2, 3, 3])
    np.testing.assert_array_equal(result, expected)

def test_add_scalars_limits_2d():
    in_array = np.array([1, 2, 3, 4])
    nxyz = [2, 2]
    # It reshapes to 2x2: [[1, 2], [3, 4]]
    # Then inserts row 1 at idx 1: [[1, 2], [1, 2], [3, 4]]
    # Then inserts last row at idx -1 (before last): [[1, 2], [1, 2], [3, 4], [3, 4]]
    # Then col 1: [[1, 1, 2], [1, 1, 2], [3, 3, 4], [3, 3, 4]]
    # Then col -1: [[1, 1, 2, 2], [1, 1, 2, 2], [3, 3, 4, 4], [3, 3, 4, 4]]
    # Reshaped to (1, 16): [[1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4]]

    result = add_scalars_limits(in_array, nxyz)
    expected = np.array([[1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4]])
    np.testing.assert_array_equal(result, expected)

def test_add_scalars_limits_3d():
    in_array = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    nxyz = [2, 2, 2]
    with pytest.raises(ValueError, match="Unsupported dimension: 3. Only 1D and 2D arrays are supported."):
        add_scalars_limits(in_array, nxyz)

def test_add_scalars_limits_0d():
    in_array = np.array([1])
    nxyz = []
    with pytest.raises(ValueError, match="Unsupported dimension: 0. Only 1D and 2D arrays are supported."):
        add_scalars_limits(in_array, nxyz)

def test_add_coords_limits_1d():
    cc_old = np.array([10, 20, 30])
    nxyz = [3]
    lxyz = [40]
    result = add_coords_limits(cc_old, nxyz, lxyz, 1)
    expected = np.array([0, 10, 20, 30, 40])
    np.testing.assert_array_equal(result, expected)

def test_add_coords_limits_2d():
    cc_old = np.array([10, 5, 10, 15, 20, 5, 20, 15])
    nxyz = [2, 2]
    lxyz = [30, 40]
    result = add_coords_limits(cc_old, nxyz, lxyz, 2)
    expected = np.array([
         0.,  0.,  0.,  5.,  0., 15.,  0., 40.,
        10.,  0., 10.,  5., 10., 15., 10., 40.,
        20.,  0., 20.,  5., 20., 15., 20., 40.,
        30.,  0., 30.,  5., 30., 15., 30., 40.
    ])
    np.testing.assert_array_equal(result, expected)
