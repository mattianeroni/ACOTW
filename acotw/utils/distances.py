import math 


def euclidean_from_coordinates (c1, c2):
    """ 
    Method to calculate the Euclidean distance between two positions / nodes.
    
    :param c1: (x, y) coordinates of node 1.
    :param c2: (x, y) coordinates of node 2.
    :return: The euclidean distance.
    """
    x1, y1 = c1 
    x2, y2 = c2 
    return math.sqrt(math.pow(x1 - x2, 2) + math.pow(y1 - y2, 2))


def euclidean_from_ids (grid, n1, n2):
    """ 
    Method to calculate the Euclidean distance between two nodes knowing 
    their ids and the grid they belong to.
    
    :param grid: The grid instance.
    :param n1: First node id.
    :param n2: Second node id.
    :return: The euclidean distance.
    """
    pos = grid.pos 
    return euclidean_from_coordinates(pos[n1], pos[n2])
