import collections
import numpy as np 



class Ant:
    """ An ant used by the ACO to explore a possible path """
    def __init__(self, mapsize, source, target):
        """ 
        :param mapsize: the size of the times matrix.
        :param source: Id of source node.
        :param target: Id of target node.

        :attr time: The current time seen by the ant.
        :attr path: The path under construction.
        """
        self.path = collections.deque((source,))
        self.pathmap = np.zeros(mapsize).astype(np.int32)
        self.windows = collections.deque()
        self.source = source 
        self.target = target 
        self.time = 0.0 
    
    @property 
    def cnode (self):
        """ The node where the ant currently is """
        return self.path[-1]

    @property 
    def oldnode (self):
        """ Return the node visited before the current one (if any) """
        if len(self.path) > 1:
            return self.path[-2]
        return None