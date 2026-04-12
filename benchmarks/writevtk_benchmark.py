import time

nx, ny, nz = 100, 100, 1
n_dim = 2

# Create a fake x 3D list
x = [[[0 for _ in range(nx)] for _ in range(ny)] for _ in range(nz)]
val = 0
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            x[iz][iy][ix] = val
            val += 1

def old_way():
    h_str = ''
    for iz in range(nz):
        for iy in range(ny-1):
            for ix in range(nx-1):
                h_str = h_str + str(2**n_dim) + ' ' \
                        + str(x[0][iy][ix]) + ' ' \
                        + str(x[0][iy][ix+1]) + ' ' \
                        + str(x[0][iy+1][ix]) + ' ' \
                        + str(x[0][iy+1][ix+1]) + '\n'
    return h_str

start_time = time.time()
for _ in range(50):
    res1 = old_way()
end_time = time.time()
print(f"Old way (50 runs): {end_time - start_time:.4f} seconds")

def new_way():
    parts = []
    prefix = str(2**n_dim) + ' '
    for iz in range(nz):
        for iy in range(ny-1):
            for ix in range(nx-1):
                parts.append(prefix \
                        + str(x[0][iy][ix]) + ' ' \
                        + str(x[0][iy][ix+1]) + ' ' \
                        + str(x[0][iy+1][ix]) + ' ' \
                        + str(x[0][iy+1][ix+1]) + '\n')
    return ''.join(parts)

start_time = time.time()
for _ in range(50):
    res2 = new_way()
end_time = time.time()
print(f"New way (50 runs): {end_time - start_time:.4f} seconds")

assert res1 == res2
print("Results match!")
