import functools 
import numpy as np 


@functools.lru_cache(maxsize=None)
def compute_turn (oldnode, currentnode, newnode):
    """
    Method to calculate the turn / angle in the path moving from 
    currentnode to newnode.

    NOTE: We use caching to speed up the computation.

    :param oldnode, currentnode, newnode: The coordinates (x, y) of last visited node
                                        current node, and newnode respectively.
    :return: The turn in radians
    """
    vector1 = np.array([currentnode[0] - oldnode[0], currentnode[1] - oldnode[1]])
    vector2 = np.array([newnode[0] - currentnode[0], newnode[1] - currentnode[1]])
    
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)

    dot_product = np.dot(unit_vector1, unit_vector2)

    return np.arccos(dot_product)