import unittest
from unittest.mock import MagicMock, patch, call
import sys
import os

class TestReadTxtToList(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Mock numpy to handle missing dependency before anything is imported
        cls.mock_np = MagicMock()
        cls.patcher = patch.dict('sys.modules', {'numpy': cls.mock_np})
        cls.patcher.start()

    @classmethod
    def tearDownClass(cls):
        cls.patcher.stop()

    @patch('os.makedirs')
    @patch('os.path.exists')
    @patch('builtins.open', new_callable=MagicMock)
    @patch('read_txt_to_list.write_vtk')
    @patch('read_txt_to_list.add_scalars_limits')
    @patch('read_txt_to_list.add_coords_limits')
    @patch('read_txt_to_list.get_itemwise_scalars')
    @patch('read_txt_to_list.read_files')
    @patch('os.getcwd')
    def test_main_default_path(self, mock_getcwd, mock_read_files, mock_get_itemwise, mock_add_coords,
                               mock_add_scalars, mock_write_vtk, mock_open, mock_exists, mock_makedirs):
        # Import after mocking numpy
        import read_txt_to_list
        from read_txt_to_list import main

        # Setup mocks
        mock_getcwd.return_value = '/mock/cwd'
        expected_default_path = os.path.join('/mock/cwd', '10-2D-F2275-TCFE-AIMD-FittedWithYapfi')

        # Mock returns for read_files
        # Return signature: fin_volcentr_coord, chem_potentials, domain_size, grad_energy_contr, mole_fractions, n_elements, n_gridpoints, n_phases, permeabilities, ph_field, ph_fractions, time, el_names, ph_names, n_dimensions, hcc, k1c
        mock_read_files.return_value = [
            MagicMock(), # fin_volcentr_coord
            MagicMock(), # chem_potentials
            [1.0, 1.0, 1.0], # domain_size
            MagicMock(), # grad_energy_contr
            MagicMock(), # mole_fractions
            [2], # n_elements
            [5, 5, 5], # n_gridpoints
            [2], # n_phases
            MagicMock(), # permeabilities
            MagicMock(), # ph_field
            MagicMock(), # ph_fractions
            [0.1, 0.2], # time (2 timesteps)
            ['El1', 'El2'], # el_names
            ['Ph1', 'Ph2'], # ph_names
            [3], # n_dimensions
            MagicMock(), # hcc
            MagicMock() # k1c
        ]

        # For simplicity, make lengths non-zero to test conditionals
        mock_read_files.return_value[15].__len__.return_value = 1
        mock_read_files.return_value[16].__len__.return_value = 1

        mock_exists.return_value = False
        mock_write_vtk.return_value = "Mock Header "

        # Mock numpy behaviors
        read_txt_to_list.np.prod.return_value = 125 # 5*5*5
        read_txt_to_list.np.array_str.return_value = "[1 2 3]"

        # Make mocked lists look like indexable/sliceable structures
        mock_get_itemwise.return_value = MagicMock()
        mock_add_scalars.return_value = MagicMock()
        mock_add_coords.return_value = MagicMock()

        # Run
        main()

        # Assertions
        mock_read_files.assert_called_once_with(expected_default_path)
        mock_exists.assert_called_once_with(os.path.join(expected_default_path, 'BCvtk'))
        mock_makedirs.assert_called_once_with(os.path.join(expected_default_path, 'BCvtk'))

        # Check files were written
        self.assertEqual(mock_open.call_count, 2)
        mock_open.assert_has_calls([
            call(os.path.join(expected_default_path, 'BCvtk', 'tstp_bc_0.vtk'), 'w'),
            call(os.path.join(expected_default_path, 'BCvtk', 'tstp_bc_1.vtk'), 'w')
        ], any_order=True)

    @patch('os.makedirs')
    @patch('os.path.exists')
    @patch('builtins.open', new_callable=MagicMock)
    @patch('read_txt_to_list.write_vtk')
    @patch('read_txt_to_list.add_scalars_limits')
    @patch('read_txt_to_list.add_coords_limits')
    @patch('read_txt_to_list.get_itemwise_scalars')
    @patch('read_txt_to_list.read_files')
    def test_main_custom_path_and_empty_arrays(self, mock_read_files, mock_get_itemwise, mock_add_coords,
                                              mock_add_scalars, mock_write_vtk, mock_open, mock_exists, mock_makedirs):
        import read_txt_to_list
        from read_txt_to_list import main

        custom_path = '/custom/path/to/data'

        # Mock read_files returning empty hcc and k1c
        hcc_mock = MagicMock()
        hcc_mock.__len__.return_value = 0
        k1c_mock = MagicMock()
        k1c_mock.__len__.return_value = 0

        mock_read_files.return_value = [
            MagicMock(), MagicMock(), [1.0, 1.0, 1.0], MagicMock(), MagicMock(), [1], [2, 2, 2], [1],
            MagicMock(), MagicMock(), MagicMock(), [0.5], ['El1'], ['Ph1'], [3], hcc_mock, k1c_mock
        ]

        mock_exists.return_value = True # Directory exists
        mock_write_vtk.return_value = "Mock Header "
        read_txt_to_list.np.prod.return_value = 8

        # Run
        main(custom_path)

        # Assertions
        mock_read_files.assert_called_once_with(custom_path)
        mock_makedirs.assert_not_called() # Should not be called since directory exists

        # Verify open was called
        mock_open.assert_called_once_with(os.path.join(custom_path, 'BCvtk', 'tstp_bc_0.vtk'), 'w')

if __name__ == '__main__':
    unittest.main()
