import random
import numpy as np
import time
import operator
import collections
import itertools





class AntColony (object):
    """
    This is the Ant Colony Optimization algorithm with Warm-Up.
    """
    def __init__ (self, grid, source, target,
                pher_init = 0.1, ro = 0.5, Q = 5.0, alpha = 1.0, beta = 5.0,
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
        """
        This method resets the algorithm.
        """
        # Initialize the best
        self.best = list(self.picking_list)
        random.shuffle(self.best)
        self.vbest = _compute_distance (self.best, self.distances)

        # Initialize the pheromone
        self.pheromone = self.saved_pheromone.copy()

        # Initialize the history and the number of iterations needed to find the best
        # and other statistics.
        self.history = collections.deque((self.vbest,))
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


    def _next_node (self, path):
        """
        This method returns the next node during the constructing process that
        brings to a new solution.

        The node is chosen from the set of neighbour nodes according to the 
        current position of the ant. Nodes already covered by the ant and part 
        of the path are excluded.

        Given the set of possible next nodes/cells (i.e., i, j), the next node is 
        chosen through a roulette wheel according to their desirability. Given j a 
        possible next node, and i the current node, the desirability p(i,j) is 
        computed as:

        p(i,j) = ph(i, j)^alpha * n(i,j)^beta / sum[k](ph(i, k)^alpha * n(i,k)^beta)

        where n(i, j) is:

        n(i, j) = Q / ( d(s) + d(t) + )

        where d(s) and d(t) are the distances as the crow flies between j and source
        and target nodes relatively, while 

        :param path: Set of nodes already covered by the current ant.
        :return: The selected node.

        """

        prob = 0.0
        r = random.random()
        options.sort(key=operator.itemgetter(1), reverse=True)
        total = sum(desirability for _, desirability in options)

        for op, desirability in options:
            prob += desirability / total
            if r < prob:
                return op
        return -1


    def _new_solution (self):
        """
        This method construct node by node a new solution.
        """
        c_node = 0
        new_sol, vnew_sol = [], 0
        options = list(self.picking_list)

        for i in range (len(self.picking_list)):
            options_params = [(op, self.pheromone[c_node, op]**self.alpha / self.distances[c_node, op]**self.beta) for op in options]
            n_node = self._next_node (options_params)
            new_sol.append(n_node)
            options.remove (n_node)
            vnew_sol += self.distances[c_node, n_node]
            c_node = n_node
        vnew_sol += self.distances[c_node, 0]

        return new_sol, vnew_sol


    def run (self, picking_list, verbose = False):
        """
        This method represents the execution of the algorithm.

        :param picking_list: The tour of nodes for which the problem must be solved.
        :param verbose: If TRUE a log takes place every <print_every> iterations.
        :return: The best solution and its cost.

        """
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
