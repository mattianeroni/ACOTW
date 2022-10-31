from acotw.aco import AntColony
from acotw.grid import Grid 


from networkx.algorithms.shortest_paths.astar import astar_path, astar_path_length



if __name__ == '__main__':
    grid = Grid( (10, 20), 1, 1)
    #grid.plot()

    aco = AntColony(
        grid=grid, 
        source=0, 
        target=139,
        pher_init = 0.1, ro = 0.5, Q = 5.0, alpha = 1.0, beta = 5.0, phi=8.0,
        evaporate = False, max_iter = 3000, max_noimp = 3000
    
    )
    path, t = aco.run(verbose=False)


    print("Best time found: ", t)
    grid.plot(path=path)

    print("A* best time is: ", astar_path_length(grid.G, 0, 139))
    #grid.plot(astar_path(grid.G, 0, 139))