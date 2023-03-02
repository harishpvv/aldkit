def main():
    from ase import Atoms
    from ase.io import read
    from ase.dft.kpoints import bandpath
    import argparse
    import os
    import pickle
    import ase.data

    ################

    path_isotope_data = os.path.dirname(__file__)+"/aldcodes/"
    ald_codes = os.path.dirname(__file__)+"/aldcodes/"
    pseudo_path = os.path.dirname(__file__)+"/pseudos/"

    ################
    
    os.system("mkdir -p scratch")
    parse = argparse.ArgumentParser()
    parse.add_argument('-f', '--file', help='Poscar or cif or pw file')
    parse.add_argument('-ppd', '--ppd', help='PseudoDirectory, you can use pbe, lda, pbesol or path where Cu.UPF like files are present.', default='pbe')
    parse.add_argument('-bp', '--band_path', help='Bandpath', default=None, type=str)

    args = vars(parse.parse_args())
    file = args['file']
    ppd = args['ppd']
    band_path = args['band_path']


    with open("%s/isotopesData.pickle" % (path_isotope_data), "rb") as iso:
        isotopes_data = pickle.load(iso)


    atoms = ase.io.read(file, index=-1)
    pos = atoms.get_scaled_positions()
    symbols = atoms.get_chemical_symbols()
    masses = atoms.get_masses()
    cell = atoms.get_cell().array

    try:
        spoints = bandpath(band_path, atoms.get_cell(), len(band_path)).kpts
    except:
        bp = atoms.cell.bandpath()
        band_path = bp.path.split(',')[0]
        sps = bp.special_points
        spoints = [sps[i] for i in band_path]
        npoints = len(spoints)
        natoms = atoms.get_global_number_of_atoms()
        print('Using ASE default bandpath: '+band_path)
    
    bandpathstring = ''
    for i, j in enumerate(spoints[:-1]):
        k = spoints[i+1]
        bandpathstring += f'{j[0]} {j[1]} {j[2]} {k[0]}  {k[1]}  {k[2]}  100 \n'

    lattice_vecs = ''
    positions = ''
    species = ''
    cutfrc = ''

    for i in cell:
        lattice_vecs += f'{i[0]} {i[1]} {i[2]}\n'

    positions += f'natoms = {len(symbols)}\n'
    for i, j in enumerate(symbols):
        positions += f'{j}  {round(pos[i][0],6)} {round(pos[i][1],6)} {round(pos[i][2],6)} {masses[i]}\n' 

    sp = []
    for i, j in enumerate(symbols):
        if j not in sp:
            sp.append(j)
            species += f'{j}  '
            a = ase.data.atomic_numbers[j]
            entry = isotopes_data[a]
            for key in entry.keys():
              if entry[key]['composition'] > 1e-6:
                species += f"{round(entry[key]['mass'], 3)} {entry[key]['composition']} "
            species += ' \n'

    els = list(dict.fromkeys(symbols))
    for i, j in enumerate(els):
        for k in els[i:]:
            cutfrc += f'{j}  {k}  6 6 0\n'

    pwd = os.getcwd()
    save = 'save'

    f = open('input.dat', 'w')
    f.write(f'''# SYSTEM
path_restart = restart
path_save = {save}
path_output = output
calculation = phph
symprecision = 6
precision = 10
temperature = 300
ntime = 1.0e-12
nlength = 1.0e-10
nmass = 1.66054e-27
nproc_calculator = 10

# STRUCTURE
alat
{lattice_vecs}
coordinates = direct
{positions}
species
{species}

# FORCE_CONSTANTS
ftype = finite-difference
iprim = yes
ddisp = 0.025
ndisp = 200
cutfrc
{cutfrc}
potential = empirical
cell4 = 4 4 4
cell3 = 4 4 4
cell2 = 6 6 6

# HARMONIC_PROPERTIES
fermi_invariance = yes
renormalize = no
helmholtz_quartic = no
helmholtz_cubic = yes
igruneissen = yes

# BAND STRUCTURE ({band_path})
idirdisp = yes
dirtype = direct
ndirdisp = {len(spoints)-1}
{bandpathstring}

# SCATTERINGS
normal_umklapp = both
ifourphonon = no
isoscat = no
ibdry = no
lbdry = 1e+9
ielph = no
fkcell = 1 1 1
ckcell = 1 1 1
cqcell = 1 1 1

# SOLVER
fsolve = taylor-series
fulldisp = pre
iterate_bte = no
coherent_bte = no
cutfermi
{cutfrc}
udisp = 200
fqcell = 20 20 20
occupation = quantum
isensitivity3 = no
iscattering_phase_space_only = no
start_q = 0
end_q = -1
fourphonon_maxfreq = 50
itime_reversal_symmetry = yes
broadening_function = gaussian
split_fit = no
equilibrium_force_fit = yes
maxeqn = 10000
itrans = yes

# NUMERICAL_PARAMETERS
scale_force_fit = 100
deltaq = 0.00001
epsil_anharmonic_renormalize = 0.01
epsil_renormalize = 0.01
betamix_renormalization = 0.5
energy_delta_cutoff3 = 3
energy_delta_cutoff4 = 2
energy_delta_scale3 = 0.1
energy_delta_scale4 = 0.1
fc_threshold3 = 0.0001
fc_threshold4 = 0.01
energy_delta_cutoff_isotope = 4
energy_delta_scale_isotope = 0.1
epsil_iter_bte = 1.0
betamix_iter_bte = 0.5
    ''')
    f.close()

    # Copy python codes

    os.system(f'cp {ald_codes}/*.py ./')

    # Change Cluster settings
    pwd = os.getcwd()
    f = open('cluster_settings.py', 'r')
    lines = f.readlines()
    f.close()
    if ppd in ['pbe', 'lda', 'pbesol']:
      lines[2] = f'pseudo_dir = "{pwd}/"\n'
    else:
      lines[2] = f'pseudo_dir = "{os.path.abspath(ppd)}/"\n'
    f = open('cluster_settings.py', 'w')
    f.writelines(lines)
    f.close()

    # Copy Pseudopotentials
    for i in els:
        if ppd=='pbe':
            string = f'cp {pseudo_path}/PBE_ONCV/{i}.UPF ./{i}.UPF'
            os.system(string)
        elif ppd=='lda':
            string = f'cp {pseudo_path}/LDA_ONCV/{i}.UPF ./{i}.UPF'
            os.system(string)
        elif ppd=='pbesol':
            string = f'cp {pseudo_path}/PBESOL_ONCV/{i}.UPF ./{i}.UPF'
            os.system(string)
