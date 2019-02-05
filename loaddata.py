import numpy as np


def get_itemwise_scalars(in_array, num_items, num_point):
    out_array = np.zeros((num_items, num_point))
    for i in range(num_items):
        out_array[i][:] = in_array[i::num_items]
    return out_array
