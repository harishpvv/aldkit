import os, sys
import numpy
import ase
import ase.io
import argparse
import copy

import cluster_settings
outdir = cluster_settings.outdir
pseudo_dir = cluster_settings.pseudo_dir
home_bashrc = cluster_settings.home_bashrc
submission_command = cluster_settings.submission_command
OMP_NUM_THREADS = cluster_settings.OMP_NUM_THREADS
modules = cluster_settings.modules
parallel = cluster_settings.parallel

nproc = os.cpu_count()

# ------- read user inputs -----------------------------------------------------
def get_args():

  parser = argparse.ArgumentParser(description='dft: wrapper for dft"\
                                                          " supercell-calc')
  
  parser.add_argument("-sfd", "--structure_file_directory", \
                          default = None, type = str, \
                          help = "directory from where structure"
                                    " files are read. ",
                          dest = "structure_file_directory")
  parser.add_argument("-sft", "--structure_file_type", \
                          default = "ald", type = str, \
                          choices = ["ald"], \
                          help = "file format from which structure is read.",
                          dest = "structure_file_type")
  parser.add_argument("-ppn", default = nproc, type = int, \
                          help = "processes-per-node", 
                          dest = "ppn")
  parser.add_argument("-N", "--nodes", default = 1, type = int, \
                          help = "number-of-nodes", 
                          dest = "nodes")

  parser.add_argument("-e", "--encut", default = None, \
                          type = float, 
                          help = "Energy-cutoff for pw.x", \
                          dest = "encut")
  parser.add_argument("-k", "--kgrid", default = None, \
                          type = int, \
                          nargs = 3, \
                          help = "scf kgrid for pw.x", \
                          dest = "kgrid")
  parser.add_argument("-v", "--vacuum", default = None, \
                          type = float, 
                          help = "vacuum to-be adjusted", \
                          dest = "vacuum")
  parser.add_argument("-t", "--title", \
                          default = "MyTitle", type = str, \
                          help = "title for job",
                          dest = "title")
  parser.add_argument("-c", "--calc", \
                          default = None, type = str, \
                          required = True, \
                          help = "type-of-calc to be performed", \
                          choices = ["pre", "post"],
                          dest = "calc")

  args = parser.parse_args(sys.argv[1:])
  
  if args.structure_file_directory is None:
    args.structure_file_directory = "save/"
    print("Using %s as the structure_file directory." % \
                                            (args.structure_file_directory))
  
  if args.structure_file_type is None:
    args.structure_file_type = "ald"
    print("Using %s as the structure_file_type." % (args.structure_file_type))

  if args.calc == "pre":
    if args.encut == None:
      parser.error("Encut needed!")
  
  if args.calc == "pre":
    if args.kgrid == None:
      parser.error("kgrid needed!")
  

  return args



#--read structure ------------------------------------------
def read_structure(filename, args):
  'returns atoms, alat'  
 
  filedir = args.structure_file_directory 
  filetype = args.structure_file_type
  if filetype == "ald":
   
    scale = 1.0 
    # read atoms now
    atoms_filename = filedir + "/" + filename   
    with open(atoms_filename, "r") as f:  
      lines = f.readlines()
      natom = int(lines[0].split()[0])
      a1 = [float(_) for _ in lines[1].split()[0:3]] 
      a2 = [float(_) for _ in lines[1].split()[3:6]] 
      a3 = [float(_) for _ in lines[1].split()[6:9]] 
      alat = [a1,a2,a3]
      atoms = []
      inv_alat = numpy.linalg.inv(alat)
      for i in range(natom):
        vals = lines[2+i].split() 
        symbol = vals[0]
        cart = [float(_) for _ in vals[2:5]]
        direct = numpy.matmul(cart, inv_alat)
        atoms.append({'symbol' : symbol,\
                      'direct' : direct, \
                      'cart' : cart})

  else:
    a = ase.io.read(filename, format=filetype, index=-1)
    alat = a.get_cell()
    symbols = a.get_chemical_symbols()
    scaled_positions = a.get_scaled_positions()
    positions = a.get_positions()
    atoms = []
    scale = 1.0
    for symbol, direct, cart in zip(symbols, scaled_positions, positions):
      atoms.append({'symbol' : symbol, \
                    'direct' : direct, \
                    'cart' : cart})


  if args.vacuum is not None:
    
    # get the vacuum dircetion
    vac_dim = -1
    for d in range(3):
      if args.kgrid[d] == 1:
        vac_dim = d
    
    # find max distance from center in this direction
    max_dist = -1e+6
    min_dist = 1e+6
    for atom in atoms:
      max_dist = max(max_dist, abs(atom['direct'][vac_dim]- 0.5))
      min_dist = min(min_dist, abs(atom['direct'][vac_dim]- 0.5))
  
    shift = min_dist + 0.5*(max_dist - min_dist)
    for atom in atoms:
      atom['direct'][vac_dim] -= shift
  
    # now bring close to center as per nearest image convention
    for atom in atoms:
      atom['direct'][vac_dim] = atom['direct'][vac_dim] % 1
   
    # fix the cart coordinates now
    for atom in atoms:
      atom['cart'] = numpy.matmul(atom['direct'], alat)

  return scale, alat, atoms


# create run script ----------------------------------------------------------
def generate_run_file(prefix, path, args, jobid, nstructure):

  with open("%s/run_%d.sh" % (path, jobid), "w") as f:
    f.write("#!/bin/bash\n")
    for module in modules:
      f.write("module load %s\n" % (module))
    f.write("export OMP_NUM_THREADS=%d\n" % (OMP_NUM_THREADS))
    f.write("\n")

    f.write("\nfor j in {0..0}\n")
    f.write("do\n")
    f.write("  for i in {0..%d}\n" % (nstructure-1))
    f.write("  do\n")
    f.write("    cd $i\n")
    f.write('    if [ ! -f "start_$j" ]\n')
    f.write('    then\n')  
    f.write('      touch "start_$j"\n')  
    f.write("      /usr/bin/rm -rf " + outdir + "/%s_$i\n" % (prefix))
    f.write("      mkdir -p " + outdir + "/%s_$i\n" % (prefix))
    f.write(f"     {parallel} -n {args.ppn} pw.x -in pw.in | tee pw.out\n")
    f.write("      /usr/bin/rm -rf " + outdir + "/%s_$i\n" % (prefix))
    f.write('    fi\n')  
    f.write("    cd ..\n")
    f.write("  done\n")
    f.write("done\n\n")
    
    f.write("\nfor i in {0..%d}\n" % (nstructure-1))
    f.write("do\n")
    f.write("  cd $i\n")
    f.write('  if ! grep -q "JOB DONE" pw.out\n')
    f.write('  then\n')  
    f.write('    if test `find "pw.out" -mmin +100`\n')
    f.write('    then\n')  
    f.write(f"      {parallel} -n {args.ppn} pw.x -in pw.in | tee pw.out\n")
    f.write('    fi\n')  
    f.write('  fi\n')  
    f.write("  cd ..\n")
    f.write("done\n\n")


# create id file  ------------------------------------------------------
def generate_id_file(path, id):
  
  with open("%s/id" % (path), "w") as f:
    f.write(args.structure_file_directory + "/" + id)


# create pw.x input file  ------------------------------------------------------
def generate_pw_file(prefix, path, args):
  
  species = set([_['symbol'] for _ in args.atoms])
  with open("%s/pw.in" % (path), "w") as f:
    f.write(" &control\n")  
    f.write("  title='%s'\n" % (args.title))  
    f.write("  calculation='scf'\n")  
    f.write("  restart_mode='from_scratch'\n")  
    f.write("  outdir='" + outdir + "/%s'\n" % (str(prefix)))  
    f.write("  pseudo_dir='" + pseudo_dir + "'\n")  
    f.write("  prefix='%s'\n" % (str(prefix)))  
    f.write("  tprnfor=.true.\n")  
    f.write(" /\n\n")
    
    f.write(" &system\n")  
    f.write("  ibrav=0\n")  
    f.write("  nat=%d\n" % (len(args.atoms)))  
    f.write("  ntyp=%d\n" % (len(species)))  
    f.write("  ecutwfc=%5.2f\n" % (args.encut))  
    #f.write("  ecutrho=%5.2f\n" % (8*args.encut))  
    f.write("  occupations='smearing'\n")  
    f.write("  smearing='gaussian'\n")  
    f.write("  degauss=0.05\n")  
    f.write(" /\n\n")  
  
    f.write(" &electrons\n")  
    #f.write("  conv_thr=%e\n" % (1e-9 * len(args.atoms)))  
    f.write("  conv_thr=%e\n" % (1e-9))  
    f.write("  mixing_beta=0.7\n")  
    f.write(" /\n\n")  
  
    f.write(" &IONS\n")  
    f.write("  \n")  
    f.write(" /\n\n")  
  
    f.write(" &CELL\n")  
    f.write(" /\n\n")  

    f.write(" ATOMIC_SPECIES\n")  
    for specie in species:
      mass = ase.atom.Atom(specie).mass 
      f.write("  %s  %5.2f %s.UPF\n" % (specie, mass, specie))  
      os.system("cp " + pseudo_dir + "%s.UPF %s/" % (specie, path))
    f.write(" \n\n")  

    f.write(" ATOMIC_POSITIONS Angstrom\n")  
    for atom in args.atoms:
      f.write("  %s  %12.10f %12.10f %12.10f\n" % (atom['symbol'], \
                                                 args.scale*atom['cart'][0], \
                                                 args.scale*atom['cart'][1], \
                                                 args.scale*atom['cart'][2]))  
    f.write(" \n\n")  

    if(args.kgrid[0] == 1 and args.kgrid[1] == 1 and args.kgrid[2] == 1):
      f.write(" K_POINTS gamma\n")  
      f.write("  %d %d %d 0 0 0\n" % 
                          (args.kgrid[0], args.kgrid[1], args.kgrid[2]))  
      f.write(" \n\n")  
    
    else:
      f.write(" K_POINTS automatic\n")  
      f.write("  %d %d %d 0 0 0\n" % 
                          (args.kgrid[0], args.kgrid[1], args.kgrid[2]))  
      f.write(" \n\n")  

    f.write(" CELL_PARAMETERS Angstrom\n")  
    for i in range(3):
      f.write("  %12.10f %12.10f %12.10f\n" % \
                  (args.scale*args.alat[i][0], \
                   args.scale*args.alat[i][1], \
                   args.scale*args.alat[i][2]))  
    f.write(" \n\n")  


#------- read id from fiven path
def read_id(path):
  
  with open("%s/id" % path, "r") as f:
    line = f.readline()
    id = line.split()[0]
  return id

#------- read froces from fiven path
def read_forces(path, natom):

  forces = []
  with open(path + "/pw.out", "r") as f:
    lines = f.readlines()
    
    # check in max_iter
    converged = False
    for line in lines:
      if "convergence has been achieved in" in line:
        converged = True
    if not converged:
      print("Job not done %s" % (path))
      return None

    # check if terminated
    terminated  = False
    for line in lines:
      if "This run was terminated on" in line:
        terminated = True
    if not terminated:
      print("Job not done %s" % (path))
      return None

    # read and append forces now
    for l in range(len(lines)):
      line = lines[l]
      if "Forces acting on atoms (cartesian axes, Ry/au)" in line:
        for j in range(natom):
          line = lines[l + j + 2]
          force = [float(_) for _ in line.split()[-3:]]
          forces.append(force)
  
  return forces


#-------------------------------------------------------------------------------
######### ------ MAIN ----- ####################################################

# read args 
args = get_args()
pwd = os.getcwd()
  
script = open("script.sh", "w")
script.write("#!/bin/bash\n")
#script.write("source " + home_bashrc + "/.bashrc\n")

# count num files
if args.structure_file_type == "ald":
  files_ = os.listdir(args.structure_file_directory)
  files = []
  for file in files_:
    if ".xyz" == file[-4:]:
      files.append(file)

  scale, alat, atoms = read_structure(files[0], args)
  natom = len(atoms)

# set the path 
path = "dft_%d/" % (natom)
  

if args.calc == "pre":

  # generate prefix
  prefix = str(numpy.random.randint(0,1e+9))
  prefix = path[:-1] + prefix
  os.system("rm -rf %s" % (path))

  for i, file in enumerate(files):
    job_path = "%s/%d/" % (path, i)
    job_prefix = "%s_%d" % (prefix, i)
 
    # create job dir
    os.system("mkdir -p %s" % (job_path))
  
    # generate pw file
    args.scale, args.alat, args.atoms = read_structure(file, args)
    generate_pw_file(job_prefix, job_path, args)
    generate_id_file(job_path, file)
  

  generate_run_file(prefix, path, args, 0, len(files))
   
      
  script.write("cd %s/\n" % (path))
  script.write(submission_command + " ./run_0.sh -f %s\n" % (path))
  script.write("cd %s\n" % (pwd))


if args.calc == "post":
  
  # read all forces
  factor = (13.605698066/0.52917720859)
  for s in range(len(files)):
    job_path = "%s/%d/" % (path, s)
    
    with open(job_path + "/pw.out", "r") as f:                                   
      lines = f.readlines()                                                      
      for ll in lines:                                                           
        if 'number of atoms/cell' in ll:                                         
          natom = int(ll.split()[-1])
          break
    
    forces = read_forces(job_path, natom)
    id = read_id(job_path)
    if forces == None:
      continue
    for j in range(natom):
      for d in range(3):
        forces[j][d] = forces[j][d] * factor

    id = id[:-4] + ".force"
    with open(id, "w") as f:
      for j in range(natom):
        for d in range(3):
          f.write("%12.10f " % (forces[j][d]))
        f.write("\n")  

