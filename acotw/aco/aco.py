import random
import numpy as np
import time
import operator
import collections 
import itertools
import matplotlib.pyplot as plt 

from acotw.utils.distances import euclidean_from_ids
from acotw.utils.angle import compute_turn 

from .ants import Ant 



class AntColony:
    """
    This is the Ant Colony Optimization algorithm with Warm-Up.
    """
    def __init__ (self, grid, source, target,
                pher_init = 0.1, ro = 0.5, Q = 5.0, alpha = 1.0, beta = 5.0, phi=0.0,
                evaporate = False, max_iter = 3000, max_noimp = 1000):
        """
        Initialize.

        :param grid: The grid.
        :param source: The source node the AGV comes from.
        :param target: The target node where the AGV is going.

        :param ro: A parameter that defines the evaporation of the pheromone.

        :param Q: A parameter that defines the increment of the pheromone on
                the new best path.

        :param alpha, beta: Parameters of the empirical distribution used to
                            select the next node at each step.

        :param phi: Penalty given to turns during the selection of next node.

        :param evaporate: If TRUE the pheromone evaporates at every iteration,
                        otherwise only when a better best solution is found.

        :param max_iter: The number of iterations.

        :param max_noimp: Maximum number of iterations without improvement.

        :param pher_init: The initial pheromone (it is always 0 on arcs (i,j) where i == j).

        """
        self.grid = grid 
        self.times = grid.times
        self.source = source 
        self.target = target

        self.ro = ro
        self.Q = Q
        self.phi = phi
        self.alpha = alpha
        self.beta = beta
        self.evaporate = evaporate
        self.max_iter = max_iter
        self.max_noimp = max_noimp
        self.pher_init = pher_init

        # Initialize the pheromone
        self.pheromone = np.full(grid.times.shape, pher_init)
        np.fill_diagonal(self.pheromone, 0)

        # Initialise a NULL best solution and NULL picking list for the moment
        self.best_ant = None 
        self.best_time = float("inf")       

        # Initialize the history and the number of iterations needed to find the best
        # and other statistics.
        self.history = collections.deque()
        self.computations = 0
        self.computational_time = 0.0


    def reset (self):
        """ This method resets the algorithm """
        # Initialize the best
        self.best_ant = None 
        self.best_time = float('inf')

        # Initialize the pheromone
        self.pheromone = np.full(self.grid.times.shape, self.pher_init)
        np.fill_diagonal(self.pheromone, 0)

        # Initialize the history and the number of iterations needed to find the best
        # and other statistics.
        self.history = collections.deque()
        self.computations = 0
        self.computational_time = 0.0


    def plot_history (self):
        """ Plot evolution of best solution found """
        plt.plot(self.history)
        plt.xlabel("Iterations")
        plt.ylabel("Best solution")
        plt.show()


    def _evap (self):
        """ This method evaporates the pheromone """
        self.pheromone *= self.ro



    def _update (self, ant):
        """ 
        Method to update the pheromone on the best path 
        
        As in the Max-Min Ant System, we update the pheromone only on the best 
        path, when a new best path is found.
        """
        self.pheromone += ant.pathmap * ( self.Q / ant.time )



    def _probability (self, ant, nextnode):
        """
        Method to compute the desirability of nextnode.

        Given j a possible next node, and i the current node, the desirability p(i,j) is 
        computed as:

        p(i,j) = ph(i, j)^alpha * n(i,j)^beta / sum[k](ph(i, k)^alpha * n(i,k)^beta)

        where n(i, j) is:

        n(i, j) = Q / ( d(j, s) + d(j, t) + w(i, j) + phi * turn(i, j) )

        where d(j, s) and d(j, t) are the distances as the crow flies between j and source
        and target nodes relatively, while w(i, j) is the idle time the AGV would have 
        to wait before entering node j (because of a busy time window), turn (i, j)
        is the turn the AGV has to do to reach node j, and phi is the penalty assigned 
        tu turns.
        """
        # Move useful variables to the stack
        pheromone, grid, alpha, beta, Q, phi = self.pheromone, self.grid, self.alpha, self.beta, self.Q, self.phi
        source, target, ctime, cnode, oldnode = ant.source, ant.target, ant.time, ant.cnode, ant.oldnode

        # Compute euclidean distances between next node and source / target
        ds = euclidean_from_ids(grid, nextnode, source)
        dt = euclidean_from_ids(grid, nextnode, target)

        # Estimate arrival time at next node
        arrival_at_nextnode = ctime + grid.times[cnode, nextnode]

        # Get eventually interested busy time window on next node
        time_windows = grid.G.nodes[nextnode]["windows"]
        window = next((i for i in time_windows if i[0] <= arrival_at_nextnode and i[1] > arrival_at_nextnode), None)

        # Time the ant has to wait to see next node available
        wait = 0.0 if not window else window[1] - arrival_at_nextnode

        # Compute the turn 
        turn = 0.0 if not oldnode else compute_turn(grid.pos[oldnode], grid.pos[cnode], grid.pos[nextnode])

        # Compute n(i, j)
        ni = Q / ( ds + dt + wait + (phi * turn))

        # Return the desirability
        return pheromone[cnode, nextnode]**alpha * ni**beta


    def _next_node (self, ant):
        """
        This method returns the next node during the constructing process that
        brings to a new solution.

        The node is chosen from the set of neighbour nodes according to the 
        current position of the ant. Nodes already covered by the ant and part 
        of the path are excluded.

        Given the set of possible next nodes/cells (i.e., i, j), the next node is 
        chosen through a roulette wheel according to their desirability. 

        :param ant: The Ant instance currently used to explore a new path.
        :return: The selected node.

        """
        # Move useful variables to the stack
        grid, G, tabu = self.grid, self.grid.G, set(ant.path) 

        # Define possible next nodes
        options = {i : self._probability(ant, i) for _, i in G.edges(ant.cnode) if i not in tabu}

        # If the ant stacks, we interrupt the process and go to the next ant
        if len(options) == 0:
            return None
        # Choose the next node
        cumulate, total_probability = 0.0, sum(options.values()) 
        r = random.random() 
        for node, probability in options.items():
            cumulate += probability / total_probability
            if r < cumulate:
                return node
        


    def _new_solution (self, ant):
        """ Method to explore a new pathconstructed node-by-node """
        grid, times = self.grid, self.grid.times 
        cnode, target, _next_node = ant.source, ant.target, self._next_node

        # Until the considered ant do not reach the target or an error is raised...
        while cnode != target:
            
            # Pick next node
            nextnode = _next_node(ant)
            if not nextnode:
                return True

            # Update ant's path, position, and travel times
            ant.path.append(nextnode)
            ant.pathmap[cnode, nextnode] = 1

            arrival_at_nextnode = ant.time + times[cnode, nextnode]
            window = next((i for i in grid.G.nodes[nextnode]["windows"] if i[0] <= arrival_at_nextnode and i[1] > arrival_at_nextnode), (0, 0))
            ant.time = max(arrival_at_nextnode, window[1])

            # Create the time window for the cnode
            # NOTE: In this way we are not going to have a time window for the target node
            _tw = ant.windows[-1] if len(ant.windows) > 0 else (0, 0)  # Time Windows of oldnode (i.e., the one visited before cnode)
            ant.windows.append( (_tw[1], ant.time) )
            
            cnode = nextnode

        return False



    def run (self, verbose = False):
        """ Main execution method for the ACO """
        # Move useful variables to the stack for a little more performant execution
        max_iter, max_noimp, _new_solution = self.max_iter, self.max_noimp, self._new_solution
        source, target = self.source, self.target 
        grid, times = self.grid, self.times
        history = self.history

        # Initialize variables and ants
        ants = (Ant(times.shape, source, target) for _ in range(max_iter))
        best_ant = None 
        best_time = float('inf')
        start = time.time()
        noimp = 0

        # Main exploration loop
        for i, ant in enumerate(ants):

            error = _new_solution(ant)
            if error:
                continue

            # Evaporate if required at each iteration
            if self.evaporate is True:
                self._evap ()

            if ant.time < best_time:
                # Update the best
                best_ant, best_time = ant, ant.time

                # Evaporate if only required when a new best is found 
                if self.evaporate is False:
                    self._evap ()
                
                # Update pheromone
                self._update(best_ant)

                # Update metrics
                noimp = 0
                self.computations = i
            
            else:
                # Check if stopping for no improvement since many iterations
                noimp += 1
                if noimp > self.max_noimp:
                    break
            
            history.append(best_time)

            # Print the best path every 100 iterations if verbose
            if verbose and i % 100 == 0:
                grid.plot(path=tuple(best_ant.path))
        
        # Save best solution found and computational time
        self.computational_time = round(time.time() - start, 3)
        self.best_ant = best_ant 
        self.best_time = round(best_time, 3) 
        # Return best path and solution cost
        return tuple(best_ant.path), round(best_time, 3)
