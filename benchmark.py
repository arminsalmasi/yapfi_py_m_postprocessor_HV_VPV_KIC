import time
import numpy as np
from loaddata import get_itemwise_scalars

def get_itemwise_scalars_optimized(in_array, num_items, num_point):
    return in_array.reshape(num_point, num_items).T

num_items = 100
num_point = 100000
in_array = np.random.rand(num_items * num_point)

# Baseline
start = time.time()
for _ in range(100):
    res1 = get_itemwise_scalars(in_array, num_items, num_point)
base_time = time.time() - start

# Optimized
start = time.time()
for _ in range(100):
    res2 = get_itemwise_scalars_optimized(in_array, num_items, num_point)
opt_time = time.time() - start

print("Baseline:", base_time)
print("Optimized:", opt_time)
print("Speedup:", base_time / opt_time if opt_time > 0 else "inf")
print("Correctness:", np.allclose(res1, res2))
