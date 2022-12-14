import math 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 

from acotw.utils.distances import euclidean_from_coordinates






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
        pos = dict()
        node_id = 0
        nodes = np.zeros(size).astype(np.int32)

        X, Y = size

        for x in range(X):
            for y in range(Y):

                # Update the graph
                G.add_node(node_id, windows=set(), conflicting=None)
                nodes[x,y] = node_id
                pos[node_id] = (x * cell_sizex + cell_sizex / 2, y * cell_sizey + cell_sizey / 2)

                if y > 0:
                    G.add_edge(node_id, node_id - 1, weight=cell_sizey)
                
                if x > 0:
                    G.add_edge(node_id, node_id - Y, weight=cell_sizex)
                    if y > 0:
                        G.add_edge(node_id, node_id - Y - 1, weight=math.sqrt(math.pow(cell_sizex,2) + math.pow(cell_sizey,2)))
                    if y < Y - 1:
                        G.add_edge(node_id, node_id - Y + 1, weight=math.sqrt(math.pow(cell_sizex,2) + math.pow(cell_sizey,2)))
            
                node_id += 1

        self.G = G 
        self.pos = pos 
        self.nodes = nodes 

        # Save the moving times 
        self.times = np.asarray([[ 
            G[i][j]["weight"] if G.get_edge_data(i, j) else 0
            for i in G.nodes] 
                for j in G.nodes]
        ).astype(np.float32).round(3)

    
    def plot (self, path=None):
        """ Plot the grid """
        nx.draw(self.G, pos=self.pos, with_labels=True, font_size=7, font_weight='bold')
        if path:
            path_edges = list(zip(path[:-1],path[1:]))
            nx.draw_networkx_nodes(self.G, self.pos,nodelist=path,node_color='r')
            nx.draw_networkx_edges(self.G, self.pos,edgelist=path_edges,edge_color='r')
        plt.show()

    
    def __getitem__(self, index):
        """ Get a node by position on the chessboard or by id """
        if isinstance(index, tuple):
            row, col = index
            return self.G.nodes[self.nodes[row, col]]
        return self.G.nodes[index]






if __name__ == '__main__':
    grid = Grid( (10, 20), 1, 1)
    node = grid[1, 2]
    print( node )
    #print(grid.G.edges((0,0), (1,1)))
    print(grid.times)
    print(grid.G[0][21])
    print(grid.G[21][0])
    print(grid.times[0, 21], grid.times[21, 0])

    grid.plot()