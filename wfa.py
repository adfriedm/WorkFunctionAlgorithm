import numpy as np
import itertools
import bisect
from munkres import Munkres


class WFA:
    def __init__(self, space_pts, space_dist, init_conf, exp_len=20):
        """ Initializer for Work Function Algorithm instance

        Keyword arguments:
        space_pts -- points from a finite metric space
        space_dist -- pseduometric defined on the points of the space
        init_conf -- the initial configuration of k points
        exp_len -- the expected sequence length (default 20)
        """
        self.space_pts = space_pts
        self.space_dist = space_dist
        self.init_conf = tuple( self.__pts_to_inds(init_conf) )
        self.n = 0
        self.n_max = exp_len
        self.k = len(init_conf)
        self.s = len(space_pts)

        self.hungarian = Munkres()
        self.confs_dist_dict = {}
        self.pts_dist_dict = {}

        # All configurations are represented by tuples
        self.confs = list( itertools.combinations(xrange(self.s), self.k) )
        # Lookup table to find configurations by request
        self.confs_dict = {}
        for conf in self.confs:
            for ind in conf:
                if ind not in self.confs_dict:
                    self.confs_dict[ind] = []
                self.confs_dict[ind].append(conf)

        self.reqs = []
        self.prev_conf = self.init_conf

        # Predesignate space for work costs
        self.n_conf_res = len(self.confs_dict[ind])
        self.work_cost = np.zeros( (self.n_conf_res, self.n_max) )

    def __dist(self, pt1, pt2):
        """ Finds the distance between points by index, checking for stored
        values first

        Keyword arguments:
        pt1 -- index of first point
        pt2 -- index of second point

        Returns:
        pt_dist -- distance between the points in the space
        """
        if pt1 > pt2:
            key_pair = (pt2,pt1)
        else:
            key_pair = (pt1,pt2)

        if key_pair in self.pts_dist_dict:
            pt_dist = self.pts_dist_dict[ key_pair ]
        else:
            pt_dist = self.space_dist( self.space_pts[pt1], self.space_pts[pt2] )
            self.pts_dist_dict[ key_pair ] = pt_dist
        return pt_dist

    def __conf_dist(self, conf1, conf2):
        """ Finds the minimum bipartite matching distance between
        configurations by index, checking for stored values first

        Keyword arguments:
        conf1 -- index of first configuration
        conf2 -- index of second configuration

        Returns:
        conf_dist -- distance between the configurations in the space
        """
        if (conf1, conf2) in self.confs_dist_dict:
            conf_dist = self.confs_dist_dict[ (conf1, conf2) ]
        elif (conf2, conf1) in self.confs_dist_dict:
            conf_dist = self.confs_dist_dict[ (conf2, conf1) ]
        else:
            # Calculate minimum bipartite matching
            adj_mat = [ [self.__dist(a,b) for b in conf2] for a in conf1 ]
            min_indices = self.hungarian.compute( adj_mat )
            conf_dist = sum( adj_mat[i][j] for (i,j) in min_indices )
            self.confs_dist_dict[ (conf1, conf2) ] = conf_dist

        return conf_dist

    def serve(self, pt):
        # Lookup index of requested point. TO DO ERROR CHECK
        pt_idx = self.space_pts.index(pt)

        # Expand work_cost array if full
        if self.n == self.n_max:
            new_mat = np.zeros( (self.n_conf_res, self.n_max*2) )
            new_mat[0:self.n_conf_res, 0:self.n_max] = self.work_cost
            self.work_cost = new_mat
            self.n_max *= 2

        # Get lists of all current possible configs
        cur_confs = self.confs_dict[pt_idx]

        # For the first request calculate the minimum cost matching
        if self.n == 0:
            for i,conf in enumerate(cur_confs):
                #min cost matching init_conf to conf
                #adj_mat = [ [self.__dist(a,b) for b in conf] for a in self.init_conf ]
                #min_indx = self.hungarian.compute( adj_mat )
                #calc_dist = sum( adj_mat[ix][jx] for (ix,jx) in min_indx )
                #self.work_cost[i,0] = calc_dist
                self.work_cost[i,0] = self.__conf_dist(self.init_conf, conf)
            min_idx = self.work_cost[:,0].argmin()
            best_conf = cur_confs[ min_idx ]
        else:
            # Get the list of the last iterations configurations
            last_confs = self.confs_dict[ self.reqs[-1] ]

            # Generate the minimum cost path from previous to next configs
            path_dict = {c_conf:float('inf') for c_conf in cur_confs}
            for i,conf in enumerate(last_confs):
                # Find all configurations we can end up in by moving just one
                # server, as well as the cost of such a movement
                pos_perm = self.__perm_swap(conf, pt_idx)
                # For every possible next config, see if the current path beats
                # the current best, and if so then store the cost
                for cost,perm in pos_perm:
                    path_cost = self.work_cost[i,self.n-1] + cost
                    if path_dict[perm] > path_cost:
                        path_dict[perm] = path_cost

            best_conf = None
            best_score = float('inf')
            for i,conf in enumerate(cur_confs):
                # Store the minimum path costs
                self.work_cost[i,self.n] = path_dict[conf]
                # Calculate work function objective function
                obj_score = path_dict[conf] + self.__conf_dist(self.prev_conf, conf)
                if obj_score < best_score:
                    best_conf = conf
                    best_score = obj_score


        # Store request
        self.reqs.append(pt_idx)
        self.n += 1
        # Store chosen config according to the work function algorithm
        self.prev_conf = best_conf

        return tuple(self.space_pts[ix] for ix in best_conf)

    def __pts_to_inds(self, pts):
        """ Take list of pts and return generator of indices

        Keyword arguments:
        pts -- list of points

        Returns:
        pt_gen -- generator of indices
        """
        pt_gen = (self.space_pts.index(pt) for pt in pts)
        return pt_gen

    def __perm_swap(self, perm, idx):
        """ Generates a list of ordered permutations by inserting val and
        removing one of each element

        Keyword argument:
        perm -- permutation as tuple
        val  -- value to insert

        Note: val must be different to all elements in perm, the check has been
        omitted for speed

        Returns:
        list of tuples with the first entry containing the cost of transition,
        the second is the subsequent permutation of indices
        """
        if idx in perm:
            return [(0.,perm)]
        perml = list(perm)
        bisect.insort_left(perml, idx)
        return [(self.__dist(idx,y), tuple(x for x in perml if x!=y)) for y in reversed(perm)]
