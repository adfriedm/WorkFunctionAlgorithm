""" Example of using WFA on a grid in Euclidean space """
from wfa import WFA
import numpy as np

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
