import unittest
import numpy as np
from writevtk import write_vtk

class TestWriteVTK(unittest.TestCase):
    def test_write_vtk_1d(self):
        coord = np.array([[0.0], [1.0], [2.0]])
        num_ngd = [3]
        tstp = 10
        elnames = ['Fe']
        phnames = ['Alpha']
        mf = [np.array([[0.1, 0.2, 0.3]])]
        mur = [np.array([[-1.0, -2.0, -3.0]])]
        phf = [np.array([[0.5, 0.6, 0.7]])]

        result = write_vtk(coord, num_ngd, tstp, elnames, phnames, mf, mur, phf)

        self.assertIn('# vtk DataFile Version 2.0', result)
        self.assertIn('All information at time step: 10', result)
        self.assertIn('DATASET UNSTRUCTURED_GRID', result)
        self.assertIn('POINTS 3 Double', result)

        # Check cell representation for 1D
        self.assertIn('CELLS 2 6', result) # nx-1 = 2, (nx-1)*(2**1+1) = 2*(3) = 6
        self.assertIn('2 0 1\n2 1 2\n', result)
        self.assertIn('CELL_TYPES 2', result)
        self.assertIn('4 4\n', result) # 1D cells are line type (VTK type 4)

        self.assertIn('POINT_DATA 3', result)
        self.assertIn('SCALARS mole-fraction(Fe) Double 1', result)
        self.assertIn('SCALARS chemical-potential(Fe) Double 1', result)
        self.assertIn('SCALARS phase-fraction(Alpha) Double 1', result)

        # Check coordinates format: n_dim to 3 so padded with 0.0s
        self.assertIn('0.000000000000000E+00 0.000000000000000E+00 0.000000000000000E+00', result)
        self.assertIn('1.000000000000000E+00 0.000000000000000E+00 0.000000000000000E+00', result)
        self.assertIn('2.000000000000000E+00 0.000000000000000E+00 0.000000000000000E+00', result)

    def test_write_vtk_2d(self):
        # 2x2 grid = 4 points
        coord = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        num_ngd = [2, 2]
        tstp = 20
        elnames = ['C']
        phnames = ['Gamma']
        mf = [np.array([[0.01, 0.02, 0.03, 0.04]])]
        mur = [np.array([[-10.0, -20.0, -30.0, -40.0]])]
        phf = [np.array([[1.0, 1.0, 0.0, 0.0]])]

        result = write_vtk(coord, num_ngd, tstp, elnames, phnames, mf, mur, phf)

        self.assertIn('# vtk DataFile Version 2.0', result)
        self.assertIn('All information at time step: 20', result)
        self.assertIn('DATASET UNSTRUCTURED_GRID', result)
        self.assertIn('POINTS 4 Double', result)

        # Check cell representation for 2D
        # (nx-1)*(ny-1) = 1*1 = 1
        # (nx-1)*(ny-1)*(2**2+1) = 1*5 = 5
        self.assertIn('CELLS 1 5', result)
        self.assertIn('4 0 1 2 3\n', result)
        self.assertIn('CELL_TYPES 1', result)
        self.assertIn('8\n', result) # 2D cells are pixel/quad type (VTK type 8)

        self.assertIn('POINT_DATA 4', result)
        self.assertIn('SCALARS mole-fraction(C) Double 1', result)
        self.assertIn('SCALARS phase-fraction(Gamma) Double 1', result)

if __name__ == '__main__':
    unittest.main()
