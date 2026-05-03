import unittest
import numpy as np
import pytest
from benchmark import get_itemwise_scalars_optimized

np = pytest.importorskip('numpy')

class TestBenchmark(unittest.TestCase):
    def test_get_itemwise_scalars_optimized_standard(self):
        """Test with a standard interleaved array."""
        in_array = np.array([1, 2, 3, 4, 5, 6])
        # 2 items, 3 points each
        out_array = get_itemwise_scalars_optimized(in_array, 2, 3)
        expected_array = np.array([[1, 3, 5], [2, 4, 6]])
        np.testing.assert_array_equal(out_array, expected_array)

    def test_get_itemwise_scalars_optimized_single_item(self):
        """Test with a single item (num_items=1)."""
        in_array = np.array([10, 20, 30])
        out_array = get_itemwise_scalars_optimized(in_array, 1, 3)
        expected_array = np.array([[10, 20, 30]])
        np.testing.assert_array_equal(out_array, expected_array)

    def test_get_itemwise_scalars_optimized_empty(self):
        """Test with an empty input array."""
        in_array = np.array([])
        # 2 items, 0 points each
        out_array = get_itemwise_scalars_optimized(in_array, 2, 0)
        expected_array = np.zeros((2, 0))
        np.testing.assert_array_equal(out_array, expected_array)

if __name__ == '__main__':
    unittest.main()
