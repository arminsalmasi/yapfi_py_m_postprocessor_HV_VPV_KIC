import numpy as np
import re
import sys


def write_vtk(coord, num_ngd, tstp, elnames, phnames, mf, mur, phf):
    n_dim = len(num_ngd)
    h_str = '# vtk DataFile Version 2.0' + '\n'\
            + 'All information at time step: '+str(tstp) + ' ' + "\\n" + '\n'\
            + 'ASCII' + '\n' \
            + 'DATASET UNSTRUCTURED_GRID' + '\n' \
            + 'POINTS ' + str(np.prod(num_ngd)) + ' Double' + '\n'
    coord = coord.reshape(np.prod(num_ngd), n_dim)
    app_array = np.zeros((np.prod(num_ngd), 1))
    for i in range(n_dim, 3):
        coord = np.append(coord, app_array, axis=1)
    h_str = h_str + ('\n '.join(' '.join(str('{:1.15E}'.format(cell)) for cell in row) for row in coord)) + '\n'


    if n_dim == 1:
        nx, ny, nz = num_ngd[0], 1, 1
        h_str = h_str + 'CELLS' + ' ' + str(nx-1) + ' ' + str((nx-1)*(2**n_dim+1)) + '\n'
        x = np.arange(nx)
        n_dim_pow = 2**n_dim
        p0 = x[:-1]
        p1 = x[1:]
        cells_str = "".join(f"{n_dim_pow} {a} {b}\n" for a, b in zip(p0, p1))
        h_str = h_str + cells_str
        h_str = h_str+'CELL_TYPES ' + str(nx-1) + '\n'
        type_array = np.ones((1, (nx-1)), dtype=int)*4
        h_str = h_str + ' '.join(map(str, type_array.ravel())) + '\n'
    if n_dim == 2:
        nx, ny, nz = num_ngd[0], num_ngd[1], 1
        h_str = h_str + 'CELLS' + ' ' + str((nx-1)*(ny-1)) + ' ' + str((nx-1)*(ny-1)*(2**n_dim+1)) + '\n'
        x = np.arange(nx*ny*nz).reshape(nz, ny, nx)
        n_dim_pow = 2**n_dim
        p0 = x[0, :-1, :-1].ravel()
        p1 = x[0, :-1, 1:].ravel()
        p2 = x[0, 1:, :-1].ravel()
        p3 = x[0, 1:, 1:].ravel()
        cells_str = "".join(f"{n_dim_pow} {a} {b} {c} {d}\n" for a, b, c, d in zip(p0, p1, p2, p3))
        h_str = h_str + cells_str
        h_str = h_str+'CELL_TYPES ' + str((nx-1)*(ny-1)) + '\n'
        type_array = np.ones((1, (nx-1)*(ny-1)), dtype=int)*8
        h_str = h_str + ' '.join(map(str, type_array.ravel())) + '\n'
    parts = [h_str, 'POINT_DATA ' + str(np.prod(num_ngd)) + '\n']
    for k, i in enumerate(elnames):
        parts.append('SCALARS ' + 'mole-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n'  \
                     + ' '.join(map(str, mf[k].ravel()))+'\n' \
                     + 'SCALARS ' + 'chemical-potential(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default'+'\n' \
                     + ' '.join(map(str, mur[k].ravel()))+'\n')
    for k, i in enumerate(phnames):
        parts.append('SCALARS ' + 'phase-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n' \
                     + ' '.join(map(str, phf[k].ravel()))+'\n')
    return "".join(parts)
