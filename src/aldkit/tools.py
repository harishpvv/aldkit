def get_mp_data(source, **kwargs):

    '''
    source : path for ase readable file | ase.Atoms object | mp-id as string
    returns: mp_api.Summary object
    '''
    
    from ase.io import read
    import pymatgen as pg
    from mp_api import MPRester
    import ase
    
    mpr = MPRester("R5mXfMsamJjSWfVMcFh0SrP8IrUS4Xy8")
    
    if type(source)==str:
        try:
            at = read(source)
            a  = pg.io.ase.AseAtomsAdaptor().get_structure(at)
            mpid = mpr.find_structure(a, **kwargs)
        except:
            mpid = source
            
    elif type(source)==ase.atoms.Atoms:
        a  = pg.io.ase.AseAtomsAdaptor().get_structure(source)
        mpid = mpr.find_structure(a, **kwargs)
        
    else:
        print("Please enter structure file or Atoms object or mp id\n")   
        
    
    return mpr.summary.get_data_by_id(mpid)


def log_it(func):
    import time
    def wrapper(*args, **kwargs):
        print(f"[{time.asctime()}]:  {func.__name__} called.")
        r = func(*args, **kwargs)
        print(f"[{time.asctime()}]:  {func.__name__} ended.")
        return r
    return wrapper


def is_layered(at, cutoff=2.5):
    
    """
    Returns 0 if not layered. 1 for layered in X-dir, 2 for layered in Y-dir etc.
    Input: Atoms object.
    """
    import numpy as np

    cell = at.cell.cellpar()
    px = at.get_scaled_positions()[:,0]; py = at.positions[:,1]; pz = at.positions[:,2]
    
    if 0 in np.histogram(px, bins=int(cell[0]/cutoff))[0]:
        return 1
    elif 0 in np.histogram(py, bins=int(cell[1]/cutoff))[0]:
        return 2
    elif 0 in np.histogram(pz, bins=int(cell[2]/cutoff))[0]:
        return 3
    
    else:
        return 0

def read_flfrc(flfrc_path):
    
    '''
    Returns; fc, positions, masses_by_positions, cell2, cart_cell
    Harmonic force constants as array with indices atom1, atom2, dir1, dir2, cell1, cell2, cell3
    '''
    
    import numpy as np

    with open(flfrc_path, 'r') as f:
        lines = f.readlines()

    natoms = int(lines[0].split()[1])
    nspecies = int(lines[0].split()[0])
    alat = float(lines[0].split()[3])
    direct_cell = np.array([[float(lines[i].split()[j]) for j in range(3)] for i in range(1,4)])
    cart_cell = direct_cell*alat
    masses_by_species = np.array([float(lines[i].split()[3]) for i in range(4,4+nspecies)])
    positions = np.array([[float(lines[i].split()[j]) for j in range(2, 5)] for i in range(4+nspecies, 4+nspecies+natoms)])
    masses_by_list = np.array([masses_by_species[int(lines[i].split()[1])-1] for i in range(4+nspecies, 4+nspecies+natoms)])
    
    start = natoms+nspecies+5

    cell2 = np.array([int(lines[start].split()[i]) for i in range(3)])
    tcells = cell2[0]*cell2[1]*cell2[2]

    fcs = np.zeros([natoms, natoms, 3, 3, cell2[0], cell2[1], cell2[2]])
    
    dividend = 0
    for n, l in enumerate(lines[start+1:]):
        
        if (n==0) or (n//(tcells+1)==dividend):
            a1 = int(l.split()[3]); a2 = int(l.split()[2])
            d1 = int(l.split()[1]); d2 = int(l.split()[0])
            dividend += 1

        else:
            c1 = int(l.split()[0]); c2 = int(l.split()[1]); c3 = int(l.split()[2])
            fcs[a1-1][a2-1][d1-1][d2-1][c1-1][c2-1][c3-1] = float(l.split()[-1])
            
    return fcs, positions, masses_by_list, cell2, cart_cell

def plot_contour_histogram(x, y, scaling = 1.0, nx=500, ny=500,
                           point_size = 0.1, norm=None, vmin=None, vmax=None,
                           least_value = 1, cmap='plasma'):
    
    '''
    This function plots contour histogram. x, y are 1d arrays, z-axis will count of point within dx*dy square.
    If scaling is specified plotted z will be equal to z/scaling. norm='log' is supported.
    least_value = 1e-2 set all values of z < 1e-2 to 1e-2. 
    After calling this function, we can edit the plot with the commands like plt.xlabel() etc., 
    '''
    
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    z, xedges, yedges = np.histogram2d(x, y, bins=[nx, ny])

    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    
    z = z.T.ravel()/scaling    
    
    z[z<least_value] = least_value
    
    if norm == 'log':
        norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        
    plt.scatter(X.ravel(), Y.ravel(), c=z, s=point_size, norm=norm, cmap=cmap)
    plt.xlim(np.min(x), np.max(x))
    plt.ylim(np.min(y), np.max(y))
    
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 10
    cbar.ax.set_ylabel('Count')

def load_data(filename,head=0):
        import numpy as np
        f = open(filename, 'r')
        lines = [i.split() for i in f.readlines()[head:]]
        f.close()

        for n, i in enumerate(lines):
            for o, j in enumerate(i):
                try:
                    lines[n][o] = float(j)
                except:
                    lines[n][o] = j

        length = max(map(len, lines))
        y=np.array([xi+[None]*(length-len(xi)) for xi in lines])

        return y

def get_pdf(x, w, bandwidth=0.05, xmin=0, xmax=None, npts=1000):
        
        from scipy.stats import gaussian_kde
        import numpy as np

        gk = gaussian_kde(x, weights=w)
        gk.covariance_factor = lambda : bandwidth
        gk._compute_covariance()

        if xmax == None:
            xmax = np.max(x)

        f = np.linspace(xmin, xmax, npts)
        g = gk(f)
        g[-1] = 0

        return f, gk(f)

def move_axes(ax, fig, pos=[0.1,0.1,0.5,0.5]):
    
        '''
        To move an axes from one figure to another from: ax, to: fig
        '''
        old_fig = ax.figure
        ax.remove()
        ax.figure = fig
        fig.axes.append(ax)
        fig.add_axes(ax)
        dummy_ax = fig.add_subplot()
        ax.set_position(pos)
        dummy_ax.remove()
        plt.close(old_fig)

def pair_sort(x,y):
        '''
        Returns sorted_x, y_values_wrto_sorted_x
        '''
        import numpy as np
        a = np.argsort(x)
        sorted_y = []

        for i in range(len(a)):
            sorted_y.append(y[a[i]])

        return np.sort(x), np.array(sorted_y)

def get_close_arrays(a, norm_limit=1e-2):
        '''
        returns Args of close arrays, args of similar array, lists of close arrays.
        '''
        import numpy as np
        from  itertools import chain
        args = []
        similar_args = []
        similar_arrays = []
        
        for i in range(len(a)):
            
            if i not in chain(*similar_args):
                args.append(i)
                sargs = []
            else:
                continue        
            
            for j in range(len(a))[i:]:
                if np.linalg.norm(a[i]-a[j]) < norm_limit:
                    sargs.append(j)
            similar_args.append(sargs)
            similar_arrays.append(a[i])
        return np.array(args), np.array(similar_args), np.array(similar_arrays)

def print_status(m,n, clear_screen=True):
        from IPython.display import clear_output
        
        if clear_screen:
            clear_output(wait=True)
        l = round(m/n*100)
        print(f"[{'*'*l}{'.'*(100-l)}]({m}/{n})")

def repeated(atoms, vdr=None):
    '''
    Return repeated atoms object, with original atoms
    in the center

    '''
    from ase import Atoms    
    import numpy as np

    d_pos = atoms.get_scaled_positions()
    c_pos = atoms.get_positions()
    cell = np.matrix(atoms.get_cell().array)
    sym = atoms.get_chemical_symbols()
    sym2 = sym.copy()
    n_pos = c_pos.copy()
    
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                if [i,j,k] != [0,0,0]:
                    p = c_pos.copy()
                    inc = np.array([i,j,k])
                    if (vdr != None and inc[vdr] != 0):
                        continue
                    disp = np.matmul(np.matrix(inc),cell)
                    p += disp
                    n_pos = np.append(n_pos, p, axis=0)
                    sym2 = np.append(sym2, sym)
    
    times = np.array([3,3,3])
    if vdr: times[vdr] = 1
    cell = atoms.get_cell().array
    cell *= times    
    at = Atoms(symbols=sym2, positions=n_pos, cell=cell)
    at.center()
    
    return at

def get_voronoi_volumes(atoms, vdr=None):
    
    from freud.box import Box
    from freud.locality import Voronoi

    n = len(atoms)
    c = atoms.get_cell().array
    a2 = repeated(atoms, vdr=vdr)

    pos = a2.get_positions()
    com = a2.get_center_of_mass()
    pos -= com

    c = a2.get_cell()
    b = Box.from_matrix(c.array)
    v = Voronoi()
    vc = v.compute((b, pos))

    return vc.volumes[:n]

def get_mp_atoms(idd):
    
    from aldkit.tools import get_mp_data
    from pymatgen.io.ase import AseAtomsAdaptor
    
    a = get_mp_data(idd)
    return AseAtomsAdaptor.get_atoms(a.structure) 

def get_neighbors( *args, atoms, cutoff=6, pbc=True):
    
    from ase.neighborlist import neighbor_list as nl
    from ase.build.tools import sort as asesort
    
    a = nl('ijd', atoms, cutoff=8)
    k = args
    
    atoms.pbc=pbc
    
    if len(k)==0:
        nn = np.array(a).T
        return asesort(nn, nn[:,-1])
    elif len(k)==1:
        nn = np.array(a).T[a[0]==args[0]]
        return asesort(nn, nn[:,-1])
    elif len(k)==2:
        nn = np.array(a).T[(a[0]==args[0]) & (a[1]==args[1])]
        return asesort(nn, nn[:,-1])
    else:
        return "Maximum two positional arguments are allowed, first atom id, second atom id."
