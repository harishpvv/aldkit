def main():

    from ase import Atoms
    from ase.io import read
    from ase.dft.kpoints import bandpath
    import argparse
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    parse = argparse.ArgumentParser()
    parse.add_argument('-s', '--structure', help='Structure file', default='./vc-relax/vc-relax.out')
    parse.add_argument('-n', '--ncores', help = 'Number of cpus', type=int, default=2)
    parse.add_argument('-np', '--nqpoints', help = 'Number q points per path', type=int, default=100)
    parse.add_argument('-in', '--input', help='input file for ald', default='./input.dat')
    parse.add_argument('-c', '--calc', help='Does ALD calculation needed?', default=False, type=bool)
    parse.add_argument('-bp', '--band_path', help='Bandpath', default=None, type=str)
    parse.add_argument('-p', '--parallel', help='parallelizer mpirun or aprun', default='mpirun', type=str, dest='parallel')
    parse.add_argument('-ald', '--ald_cmd', help='ald executable command', default='ald', type=str, dest='ald_cmd')
    args = vars(parse.parse_args())
    ncores = args['ncores']
    structure = args['structure']
    nqpoints = args['nqpoints']
    input_file = args['input']
    calc = args['calc']
    band_path = args['band_path']
    parallel = args['parallel']
    ald_cmd = args['ald_cmd']

    def read_file(file):
        f = open(file)
        lines = [ i.split() for i in f.readlines() ]
        try:
            lines = np.array([[float(i) for i in j] for j in lines]) 
        except:
            lines = np.array([[float(i) for i in j] for j in lines[1:]])
        f.close() 
        return lines


    atoms = read(structure)
    atoms.set_pbc([1, 1, 1])
    if band_path is not None:
        spoints = bandpath(band_path, atoms.get_cell(), len(band_path)).kpts
    else:
        bp = atoms.cell.bandpath()
        band_path = bp.path.split(',')[0]
        sps = bp.special_points
        spoints = [sps[i] for i in band_path]
    npoints = len(spoints)
    natoms = atoms.get_global_number_of_atoms()

    if calc == True:
        os.system('if [ -f log.dat ]; then mv log.dat log1.dat; fi')
        input = open(input_file, 'r')
        lines = input.readlines()
        input.close()

        for i, j in enumerate(lines):
            if 'calculation' in j:
                lines[i] = 'calculation = harmonic\n'
            if 'path_output' in j:
                lines[i] = 'path_output = bandsoutput\n'

            if 'ndirdisp' in j:
                n = int(j.split()[-1])
                for k in range(n+1):
                    lines[i+k] = ''
            if 'idirdisp' in j:
              lines[i] = ''

        lines.append('idirdisp = yes\n')
        lines.append(f'ndirdisp = {npoints-1}\n')

        path = ''
        for i, j in enumerate(spoints[:-1]):
            k = spoints[i+1]
            path += f'{j[0]} {j[1]} {j[2]} {k[0]}  {k[1]}  {k[2]}  {nqpoints}\n'

        lines.append(path)

        bands = open('bandsin.dat', 'w')
        bands.writelines(lines)
        bands.close()

        string = f'{parallel} -n {ncores} {ald_cmd} bandsin.dat'
        os.system(string)
        os.system('if [ -f log.dat ]; then mv log.dat log_d.dat; fi')
        os.system('if [ -f log1.dat ]; then mv log1.dat log.dat; fi')
        dlines = read_file('./bandsoutput/directional_harmonic_properties.dat')
        hlines = read_file('./bandsoutput/harmonic_properties.dat')

    else:
        dlines = read_file('./output/directional_harmonic_properties.dat')
        hlines = read_file('./output/harmonic_properties.dat')


    dist = []
    for i in range(len(dlines[:, 8])):
        if i%6 == 0:
            dist.append(dlines[:, 8][i])


    os.system('mkdir -p plots')
    # Band Structure
    nmodes = natoms*3

    freqs = [[] for i in range(nmodes)]
    dgps = [[] for i in range(nmodes)]

    for i in dlines:
        mode = int(i[1])
        freqs[mode].append(i[9])
        dgps[mode].append(i[14] + i[15] + i[16])

    plt.figure(1)
    for i in freqs:
        plt.plot(dist, i, 'k-', lw=2)
    plt.xlim([0, dist[-1]])
    plt.savefig('./plots/bands.svg', dpi=600)
    plt.savefig('./plots/bands.png', dpi=300)


