import math 
import collections 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 


def _distance (c1, c2):
    """ Method to calculate the distance between two cells """
    return math.sqrt(math.pow(c1.x - c2.x, 2) + math.pow(c1.y - c2.y, 2))



class Cell:
    """ A cell of the Grid """
    def __init__(self, x, y):
        self.x = x 
        self.y = y  

        # Reference to cells on top, bottom, left, right
        self.top = None 
        self.bottom = None 
        self.left = None 
        self.right = None 

        # The time windows in which the cell is busy
        self.windows = collections.deque()


    def __repr__(self):
        return f"Cell({self.x},{self.y})"


class Grid:

    """ The grid where the AGVs moves """

    def __init__(self, size, cell_sizex, cell_sizey):
        """
        Initialise. 

        :param size: A tuple representing the size of the grid --i.e., the number of cells 
        :param cell_sizex: The size of cells along x-axis 
        :param cell_sizey: The size of cells along y-axis
        :param times: A matrix of times to move from a cell to the others (if not passed by the user 
                    a simple matrix of distances is built)
        """
        G = nx.Graph()
        node_id = 0
        pos = dict()

        Y, X = size

        for x in range(X):
            for y in range(Y):

                # Update the graph
                G.add_node(node_id)
                pos[node_id] = (x * cell_sizex + cell_sizex / 2, y * cell_sizey + cell_sizey / 2)

                if y > 0:
                    G.add_edge(node_id, node_id - 1, weight=cell_sizey)
                
                if x > 0:
                    G.add_edge(node_id, node_id - X, weight=cell_sizex)

                node_id += 1

        self.G = G 
        self.pos = pos 

        # The set of cells that make the Grid
        # This can go into the previous loops (it is here just for code beauty)
        self._cells = [[Cell(x * cell_sizex + cell_sizex / 2, y * cell_sizey + cell_sizey / 2) for y in range(Y) ] for x in range(X)]


        # Pass to each cell the instances of the next ones
        for i, _row in enumerate(self._cells):
            for j, cell in enumerate(_row):
                
                if j > 0:
                    cell.top = self._cells[i][j - 1]
                
                if i > 0:
                    cell.left = self._cells[i - 1][j]

                if j < size[1] - 1:
                    cell.bottom = self._cells[i][j + 1]

                if i < size[0] - 1:
                    cell.right = self._cells[i + 1][j]

        # Save the moving times 
        self.times = np.asarray([ [ _distance(c1, c2) for c2 in self.itercells()] for c1 in self.itercells()]).astype(np.float32).round(3)


    def itercells (self):
        """ To iterate over all the cells """
        return (i for row in self._cells for i in row)

    
    def plot (self):
        """ Plot the grid """
        nx.draw(self.G, pos=self.pos, with_labels=True, font_size=8, font_weight="bold")
        plt.show()

    
    def __getitem__(self, index):
        """ Get a Cell """
        row, col = index
        return self._cells[row][col]






if __name__ == '__main__':
    grid = Grid( (10, 10), 1, 1)
    cell = grid[1, 2]
    print( cell )
    print(cell.top, cell.bottom, cell.left, cell.right )
    print(grid.times[0])
    grid.plot()