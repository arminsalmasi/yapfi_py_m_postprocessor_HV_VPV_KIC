import numpy as np


def get_itemwise_scalars(in_array, num_items, num_point):
    return in_array.reshape(num_point, num_items).T
