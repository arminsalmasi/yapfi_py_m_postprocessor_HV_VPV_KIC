import unittest
import numpy as np
import loaddata

class TestLoadData(unittest.TestCase):
    def test_get_itemwise_scalars_standard(self):
        """Test with a standard interleaved array."""
        in_array = np.array([1, 2, 3, 4, 5, 6])
        # 2 items, 3 points each
        out_array = loaddata.get_itemwise_scalars(in_array, 2, 3)
        expected_array = np.array([[1, 3, 5], [2, 4, 6]])
        np.testing.assert_array_equal(out_array, expected_array)

    def test_get_itemwise_scalars_single_item(self):
        """Test with a single item (num_items=1)."""
        in_array = np.array([10, 20, 30])
        out_array = loaddata.get_itemwise_scalars(in_array, 1, 3)
        expected_array = np.array([[10, 20, 30]])
        np.testing.assert_array_equal(out_array, expected_array)

    def test_get_itemwise_scalars_empty(self):
        """Test with an empty input array."""
        in_array = np.array([])
        # 2 items, 0 points each
        out_array = loaddata.get_itemwise_scalars(in_array, 2, 0)
        expected_array = np.zeros((2, 0))
        np.testing.assert_array_equal(out_array, expected_array)

if __name__ == '__main__':
    unittest.main()
