from acotw.aco import AntColony
from acotw.grid import Grid 

import matplotlib.pyplot as plt 
from networkx.algorithms.shortest_paths.astar import astar_path, astar_path_length



if __name__ == '__main__':
    grid = Grid( (10, 20), 1, 1)
    #grid.plot()

    aco = AntColony(
        grid=grid, 
        source=0, 
        target=139,
        pher_init = 0.1, ro = 0.5, Q = 5.0, alpha = 1.0, beta = 5.0, phi=8.0,
        evaporate = False, max_iter = 1000, max_noimp = 500
    
    )
    path, t = aco.run(verbose=False)

    print(f"Computational time: {aco.computational_time}s")
    print("Best solution cost: ", t)
    print("Best solution path: ", path)
    print("Iterations required: ", aco.computations)
    grid.plot(path=path)
    plt.plot(aco.history)
    plt.show()

    print("A* best time is: ", astar_path_length(grid.G, 0, 139))
    #grid.plot(astar_path(grid.G, 0, 139))