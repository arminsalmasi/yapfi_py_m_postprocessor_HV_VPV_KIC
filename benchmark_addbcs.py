import time
import numpy as np

def original_add_coords_limits(cc_old, nxyz, lxyz, n_dim):
    if n_dim == 2:
        nx, ny, lx, ly = nxyz[0], nxyz[1], lxyz[0], lxyz[1]
        cc_old_y = cc_old[1:ny*2:2]
        cc_old_x = cc_old[0:-1:nx*2]
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
        return cc_new

def optimized_add_coords_limits(cc_old, nxyz, lxyz, n_dim):
    if n_dim == 2:
        nx, ny, lx, ly = nxyz[0], nxyz[1], lxyz[0], lxyz[1]
        cc_old_y = cc_old[1:ny*2:2]
        cc_old_x = cc_old[0:-1:nx*2]
        cc_old_x = np.insert(cc_old_x, 0, [0], axis=0)
        cc_old_y = np.insert(cc_old_y, 0, [0], axis=0)
        cc_old_x = np.insert(cc_old_x, nx+1, lx, axis=0)
        cc_old_y = np.insert(cc_old_y, ny+1, ly, axis=0)

        cc_new = np.empty((nx+2) * (ny+2) * 2)
        cc_new[0::2] = np.repeat(cc_old_x, ny+2)
        cc_new[1::2] = np.tile(cc_old_y, nx+2)
        return cc_new

nx, ny = 1000, 1000
lx, ly = 10.0, 10.0
# Dummy cc_old: size needs to be large enough.
# cc_old_y will be cc_old[1:ny*2:2] -> needs ny elements
# cc_old_x will be cc_old[0:-1:nx*2] -> wait, it skips nx*2 elements?
# Ah, cc_old is for coordinates, so it's probably size 2*nx*ny or something?
# Let's just create a large array of size 2*nx*ny
cc_old = np.random.rand(2 * nx * ny)

nxyz = [nx, ny]
lxyz = [lx, ly]

# Check correctness
res1 = original_add_coords_limits(cc_old, nxyz, lxyz, 2)
res2 = optimized_add_coords_limits(cc_old, nxyz, lxyz, 2)
print("Correctness check:", np.allclose(res1, res2))

# Benchmark
import timeit
t1 = timeit.timeit(lambda: original_add_coords_limits(cc_old, nxyz, lxyz, 2), number=10)
t2 = timeit.timeit(lambda: optimized_add_coords_limits(cc_old, nxyz, lxyz, 2), number=10)

print(f"Original: {t1:.4f} s")
print(f"Optimized: {t2:.4f} s")
print(f"Speedup: {t1/t2:.2f}x")
