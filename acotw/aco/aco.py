import random
import numpy as np
import time
import operator
import itertools

from acotw.utils.distances import euclidean_from_ids
from acotw.utils.angle import turn 

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

        self.ro = ro
        self.Q = Q
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
        self.best = None
        self.vbest = float("inf")       

        # Initialize the history and the number of iterations needed to find the best
        # and other statistics.
        self.history = collections.deque()
        self.computations = 0
        self.computational_time = 0.0


    def reset (self):
        """ This method resets the algorithm """
        # Initialize the best
        self.best = None
        self.vbest = float('inf')

        # Initialize the pheromone
        self.pheromone = np.full(self.grid.times.shape, self.pher_init)
        np.fill_diagonal(self.pheromone, 0)

        # Initialize the history and the number of iterations needed to find the best
        # and other statistics.
        self.history = collections.deque()
        self.computations = 0
        self.computational_time = 0.0



    def _evap (self):
        """ This method evaporates the pheromone """
        self.pheromone *= self.ro



    def _update (self):
        """ This method updates the pheromone on the best path """
        for i in range (len(self.picking_list) - 1):
            self.pheromone[self.best[i], self.best[i + 1]] += (self.Q / self.distances[self.best[i], self.best[i + 1]])
        self.pheromone[0, self.best[0]] += (self.Q / self.distances[0, self.best[0]])
        self.pheromone[self.best[-1], 0] += (self.Q / self.distances[self.best[-1], 0])


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
        turn = 0.0 if not oldnode else turn(grid.pos[oldnode], grid.pos[cnode], grid.pos[nextnode])

        # Compute n(i, j)
        ni = Q / ( ds + dt + wait + phi * turn)

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
        # Choose next node
        options = {i: self._probability(ant, i) for _, i in G.edges(ant.cnode) if i not in tabu}
        return random.choices(options.keys(), weights=options.values(), k=1)[0]


    def _new_solution (self):
        """ Method to explore a new path """
        pass


    def run (self, picking_list, verbose = False):
        """ Main execution method for the ACO """
        
        # Initialise the picking list
        self.picking_list = list(picking_list)
        # Initialize the best
        self.best = list(self.picking_list)
        random.shuffle(self.best)
        self.vbest = _compute_distance (self.best, self.distances)

        # Start the effective execution
        start = time.time()
        noimp = 0
        for i in range (self.max_iter):
            # Build a new solution
            new_sol, vnew_sol = self._new_solution ()
            # Eventually evaporate pheromone
            if self.evaporate is True:
                self._evap ()
            # Eventually update best, iterations with no improvement
            # and computations needed to find the best.
            if vnew_sol < self.vbest:
                self.best, self.vbest = new_sol, vnew_sol
                if self.evaporate is False:
                    self._evap ()
                self._update ()
                noimp = 0
                self.computations = i
            else:
                noimp += 1
                if noimp > self.max_noimp:
                    break

            # Update history
            self.history.append(self.vbest)
            # Logs
            if verbose is True and i % self.print_every == 0:
                print('Epoch: ', i, ', Best: ', self.vbest)

        # Set computational time
        self.computational_time = time.time() - start
        # Return the best solution found
        return self.best, self.vbest
