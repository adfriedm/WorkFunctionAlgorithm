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



def d_euc(pt1, pt2):
    """ Calculates the euclidean distance between Rd points 

    Keyword arguments:
    pt1 -- coordinate of the first point
    pt2 -- coordinate of the second point

    Returns:
    euc_dist -- 2-norm distance between pt1 and pt2
    """
    euc_dist = np.sqrt( np.sum( (c1-c2)**2 for (c1,c2) in zip(pt1,pt2) ) )
    return euc_dist


def gen_grid(n=3, nx=None, ny=None, xs=0., xf=1., ys=0., yf=1.):
    """Generate a list of grid points in R2

    Keyword arguments:
    n -- defines nx and ny if they arent
    nx -- number of points along x axis
    ny -- number of points along y axis
    xs,xf -- start and end x grid component
    ys,yf -- start and end y grid component

    Return:
    grid_pts -- list of grid coordinates
    """
    if not nx: nx = n
    if not ny: ny = n

    grid_pts = [ ((xf-xs)*xi/(nx-1.), (yf-ys)*yi/(ny-1.)) \
                 for xi in xrange(nx) \
                 for yi in xrange(ny) ]
    return grid_pts



def wfa_test(verbose=False, draw_test=False, wait_for_user=False, save_figs=False, save_folder=None):
    ### Parameter setup ###
    k = 4   # number of servers
    s = 100  # number of requests
    n = 5   # use an nxn grid

    # Generate a grid
    space_pts = gen_grid(n)
    # Generate sequence of requests
    seq = [space_pts[i] for i in np.random.randint(len(space_pts), size=s)]
    # Initial configuration (for now this requires all points to be different)
    init_conf = tuple( space_pts[i] for i in np.random.choice(len(space_pts), k, replace=False) )

    if verbose:
        print "Initial configuration: %s" % (init_conf,)

    stored_reqs = []
    stored_configs = [init_conf]

    # Setup work function object
    wfa = WFA(space_pts, d_euc, init_conf, 1)
    # Service requests
    for i,req in enumerate(seq):
        # Serve request
        next_conf = wfa.serve( req )
        stored_reqs.append( req )
        stored_configs.append( next_conf )

        if verbose:
            print "Request (%d): %s" % (i,req)
            print "Next Configuration: %s" % (next_conf,)


    if draw_test:
        import matplotlib.pyplot as plt
        if wait_for_user:
            plt.ion()
        fig,ax = plt.subplots(1,1)

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        for i in xrange( len(stored_reqs) ):
            ### DRAW CONFIGURATION ###
            plt.cla()
            ax.axis([-0.1,1.1,-0.1,1.1])
            ax.scatter( *zip(*stored_configs[i]), s=150, color='b')
            ax.text(-0.05, -0.08, "Request #: %d"%(i,))

            # Save figures if necessary
            if save_figs:
                if save_folder:
                    plt.savefig("%sconf%03d.png"%(save_folder,i), bbox_inches='tight')
                else:
                    # If no folder specified then save to the current
                    # directory. This is not recommended
                    plt.savefig("conf%03dr.png"%(i), bbox_inches='tight')

            # Delay for mouse click
            if wait_for_user:
                plt.waitforbuttonpress()

            ### DRAW CONFIGURATION + REQUEST ###
            plt.cla()
            ax.axis([-0.1,1.1,-0.1,1.1])
            ax.scatter( *zip(*stored_configs[i]), s=150, color='b')
            ax.scatter([stored_reqs[i][0]], [stored_reqs[i][1]], s=50, color='r')
            ax.text(-0.05, -0.08, "Request #: %d"%(i,))

            # Save figures if necessary
            if save_figs:
                if save_folder:
                    plt.savefig("%sconf%03dr.png"%(save_folder,i), bbox_inches='tight')
                else:
                    # If no folder specified then save to the current
                    # directory. This is not recommended
                    plt.savefig("conf%03dr.png"%(i), bbox_inches='tight')

            # Delay for mouse click
            if wait_for_user:
                plt.waitforbuttonpress()


        ### DRAW FINAL CONFIGURATION ###
        plt.cla()
        ax.axis([-0.1,1.1,-0.1,1.1])
        ax.scatter( *zip(*stored_configs[i+1]), s=150, color='b')

        # Save figures if necessary
        if save_figs:
            if save_folder:
                plt.savefig("%sconf%03dr.png"%(save_folder,i), bbox_inches='tight')
            else:
                # If no folder specified then save to the current
                # directory. This is not recommended
                plt.savefig("conf%03dr.png"%(i), bbox_inches='tight')

        # Delay for mouse click
        if wait_for_user:
            plt.waitforbuttonpress()




if __name__ == '__main__':
    wfa_test(verbose=True, draw_test=True, wait_for_user=False, save_figs=True, save_folder="figs/")
