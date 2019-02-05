import numpy as np
import re
import sys


def write_vtk(coord, num_ngd, tstp, elnames, phnames, mf, mur, phf):
    np.set_printoptions(threshold=np.inf, linewidth=sys.maxint)
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
        for ix in range(nx-1):
            h_str = h_str + str(2**n_dim) + ' ' + str(x[ix]) + ' ' + str(x[ix+1])+ '\n'
        h_str = h_str+'CELL_TYPES ' + str(nx-1) + '\n'
        type_array = np.ones((1, (nx-1)), dtype=np.int)*4
        h_str = h_str+re.sub('[\[\]]', '', np.array_str(type_array)) + '\n'
    if n_dim == 2:
        nx, ny, nz = num_ngd[0], num_ngd[1], 1
        h_str = h_str + 'CELLS' + ' ' + str((nx-1)*(ny-1)) + ' ' + str((nx-1)*(ny-1)*(2**n_dim+1)) + '\n'
        x = np.arange(nx*ny*nz).reshape(nz, ny, nx)
        for iz in range(nz):
            for iy in range(ny-1):
                for ix in range(nx-1):
                    h_str = h_str + str(2**n_dim) + ' ' \
                            + str(x[0][iy][ix]) + ' ' \
                            + str(x[0][iy][ix+1]) + ' ' \
                            + str(x[0][iy+1][ix]) + ' ' \
                            + str(x[0][iy+1][ix+1]) + '\n'
        h_str = h_str+'CELL_TYPES ' + str((nx-1)*(ny-1)) + '\n'
        type_array = np.ones((1, (nx-1)*(ny-1)), dtype=np.int)*8
        h_str = h_str+re.sub('[\[\]]', '', np.array_str(type_array)) + '\n'
    h_str = h_str+'POINT_DATA ' + str(np.prod(num_ngd)) + '\n'
    k = 0
    for i in elnames:
        h_str = h_str + 'SCALARS ' + 'mole-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n'  \
                        + re.sub('[\[\]]', '', np.array_str(mf[k][:]))+'\n' \
                        + 'SCALARS ' + 'chemical-potential(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default'+'\n' \
                        + re.sub('[\[\]]', '', np.array_str(mur[k][:]))+'\n'
        k += 1
    k = 0
    for i in phnames:
        h_str = h_str + 'SCALARS ' + 'phase-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n' \
                        + re.sub('[\[\]]', '', np.array_str(phf[k][:]))+'\n'
        k += 1
    return h_str
