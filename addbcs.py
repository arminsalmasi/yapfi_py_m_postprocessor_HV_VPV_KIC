import numpy as np

def add_scalars_limits(in_array, nxyz):
    # type: (object, object) -> object
    n_dim=len(nxyz)
    if n_dim == 1:
        out_array = np.insert(in_array, 0, in_array[0])
        out_array = np.insert(out_array, nxyz[0]+1, in_array[-1])
    if n_dim == 2:
        temp_array = np.reshape(in_array, (nxyz[0], nxyz[1]))
        x_app_first = temp_array[0, :]
        x_app_last = temp_array[-1, :]
        temp_array = np.insert(temp_array, 1, x_app_first, axis=0)
        temp_array = np.insert(temp_array, -1, x_app_last, axis=0)
        y_app_first = temp_array[:, 0]
        y_app_last = temp_array[:, -1]
        temp_array = np.insert(temp_array, 1, y_app_first, axis=1)
        temp_array = np.insert(temp_array, -1, y_app_last, axis=1)
        out_array = np.reshape(temp_array, (1, (nxyz[0]+2)*(nxyz[1]+2)))
    return out_array

def add_coords_limits(cc_old, nxyz , lxyz, n_dim):
    if n_dim == 1:
        nx, lx, = nxyz[0], lxyz[0]
        cc_new = np.insert(cc_old, 0, [0], axis=0)
        cc_new = np.insert(cc_new, nx+1, lx, axis=0)
    if n_dim == 2:
        nx, ny, lx, ly = nxyz[0], nxyz[1], lxyz[0], lxyz[1]
        cc_old_y = cc_old[1:ny*2:2]
        cc_old_x = cc_old[0:-1:nx*2]
        #print(" ".join(str(y) for y in cc_old_y))
        #print(" ".join(str(x) for x in cc_old_x))
        cc_old_x = np.insert(cc_old_x, 0, [0], axis=0)
        cc_old_y = np.insert(cc_old_y, 0, [0], axis=0)
        cc_old_x = np.insert(cc_old_x, nx+1, lx, axis=0)
        cc_old_y = np.insert(cc_old_y, ny+1, ly, axis=0)
        cc_new = np.zeros((nx+2)*(ny+2)*2)
        k = 0
        for i in range(nx+2):
            for j in range(ny+2):
                cc_new[k] = cc_old_x[i]
                cc_new[k+1] = cc_old_y[j]
                k += 2
        #print(" ".join(str(y) for y in cc_new))
    if n_dim == 3:
        nx, ny, nz, lx, ly, lz = nxyz[0], nxyz[1], nxyz[2], lxyz[0], lxyz[1], lxyz[2]
    # Return Values
    return cc_new


