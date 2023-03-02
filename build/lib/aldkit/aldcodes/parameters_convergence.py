#!/opt/python/3.6.1.1/bin/python3
import os, sys
import numpy
import ase
import ase.io
import argparse
import copy


#------ these parameters need to be changed depending on lcuster
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

  parser = argparse.ArgumentParser(description='QE: pw.x simple jobs wrappers')
  
  parser.add_argument("-sf", "--structure_file", \
                          default = None, type = str, \
                          help = "file from which structure is read. ",
                          dest = "structure_file")
  parser.add_argument("-sft", "--structure_file_type", \
                          default = "ald", type = str, \
                          choices = ["ald", "espresso-in", "espresso-out", \
                                     "vasp"], \
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
                          required = True, \
                          help = "Energy-cutoff for pw.x", \
                          dest = "encut")
  parser.add_argument("-v", "--vacuum", default = None, \
                          type = float, 
                          help = "vacuum to-be adjusted", \
                          dest = "vacuum")
  parser.add_argument("-k", "--kgrid", default = None, \
                          type = int, \
                          nargs = 3, \
                          required = True, \
                          help = "scf kgrid for pw.x", \
                          dest = "kgrid")
  parser.add_argument("-t", "--title", \
                          default = "MyTitle", type = str, \
                          help = "title for job",
                          dest = "title")
  parser.add_argument("-df", "--dofree", \
                          default = 'all', type = str, \
                          help = "degrees allowed to change in vc-relax", \
                          choices = ["all", 'ibrav', 'x', 'y', 'z', \
                                     'xy', 'xz', 'yz', 'xyz', \
                                     "shape", 'volume', '2Dxy', "2Dshape"], \
                          dest = "dofree")
  parser.add_argument("-c", "--calc", \
                          default = None, type = str, \
                          required = True, \
                          help = "type-of-calc to be performed", \
                          choices = ["scf", "relax", "vc-relax", \
                                     "kgrid", "encut", "vacuum"],
                          dest = "calc")

  args = parser.parse_args(sys.argv[1:])
  
  if args.structure_file is None:
    args.structure_file = "save/dfpt"
    print("Using %s as the structure_file." % (args.structure_file))
  
  if args.structure_file_type is None:
    args.structure_file_type = "ald"
    print("Using %s as the structure_file_type." % (args.structure_file_type))

  if args.calc == "vc-relax" and args.dofree is None:
    parser.error("dofree required for vc-relax")
     
  return args



#--read structure ------------------------------------------
def read_structure(filename, filetype, args):
  'returns atoms, alat'  
  
  if filetype == "ald":
    with open(filename, "r") as f:  
      lines = f.readlines()
      scale = float(lines[0].split()[0])
      scale *= 1e+10
      a1 = [float(_) for _ in lines[1].split()] 
      a2 = [float(_) for _ in lines[2].split()] 
      a3 = [float(_) for _ in lines[3].split()] 
      alat = [a1,a2,a3]
      natom = int(lines[4].split()[-1])
      atoms = []
      for i in range(natom):
        vals = lines[5+i].split() 
        symbol = vals[0]
        direct = [float(_) for _ in vals[1:]]
        cart = numpy.matmul(direct, alat)
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


  if args.calc == "vacuum" or args.vacuum is not None:
    
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


# create batch script ----------------------------------------------------------
def generate_batch_file(prefix, path, args):

  with open(path + "/run.sh", "w") as f:
    f.write("#!/bin/bash\n")
    for module in modules:
      f.write("module load %s\n" % (module))
    f.write("export OMP_NUM_THREADS=%d\n" % (OMP_NUM_THREADS))   
    f.write("\n")
    f.write("\n/usr/bin/rm -rf " + outdir + "/%s\n" % (prefix))
    f.write("mkdir -p " + outdir + "/%s\n" % (prefix))
    f.write(f"\n{parallel} -n {args.ppn} pw.x -in {args.calc}.in | tee {args.calc}.out\n")
    f.write("\n/usr/bin/rm -rf " + outdir + "/%s\n" % (prefix))



# create pw.x input file  ------------------------------------------------------
def generate_qe_file(prefix, path, args):
  
  species = set([_['symbol'] for _ in args.atoms])
  with open("%s/%s.in" % (path, args.calc), "w") as f:
    f.write(" &control\n")  
    f.write("  title='%s'\n" % (args.title))  
    if args.calc not in ["kgrid", "encut", "vacuum"]:
      f.write("  calculation='%s'\n" % (args.calc))  
    else:
      f.write("  calculation='%s'\n" % ("scf"))  
    f.write("  restart_mode='from_scratch'\n")  
    f.write("  outdir='" + outdir + "/%s'\n" % (str(prefix)))  
    f.write("  pseudo_dir='" + pseudo_dir + "'\n")  
    f.write("  prefix='%s'\n" % (str(prefix)))  
    f.write("  tprnfor=.true.\n")  
    f.write("  tstress=.true.\n")  
    f.write("  etot_conv_thr=1.0d-8\n")  
    f.write("  forc_conv_thr=1.0d-5\n")  
    f.write(" /\n\n")
    
    f.write(" &system\n")  
    f.write("  ibrav=0\n")  
    #f.write("  celldm(1)=%12.10f\n" % (args.scale/0.52917720859))  
    f.write("  nat=%d\n" % (len(args.atoms)))  
    f.write("  ntyp=%d\n" % (len(species)))  
    f.write("  ecutwfc=%5.2f\n" % (args.encut))  
    #f.write("  ecutrho=%5.2f\n" % (8*args.encut))  
    f.write("  occupations='smearing'\n")  
    f.write("  smearing='gaussian'\n")  
    f.write("  degauss=0.05\n")  
    f.write(" /\n\n")  
  
    f.write(" &electrons\n")  
    f.write("  conv_thr=1.0d-12\n")  
    f.write("  mixing_beta=0.7\n")  
    f.write(" /\n\n")  
  
    f.write(" &IONS\n")  
    f.write("  \n")  
    f.write(" /\n\n")  
  
    f.write(" &CELL\n")  
    f.write("  cell_dofree='%s'\n" % (args.dofree))  
    f.write(" /\n\n")  

    f.write(" ATOMIC_SPECIES\n")  
    for specie in species:
      mass = ase.atom.Atom(specie).mass 
      f.write("  %s  %5.2f %s.UPF\n" % (specie, mass, specie))  
      os.system("cp " + pseudo_dir + "%s.UPF %s/" % (specie, path))
    f.write(" \n\n")  

    f.write(" ATOMIC_POSITIONS angstrom\n")  
    for atom in args.atoms:
      f.write("  %s  %12.10f %12.10f %12.10f\n" % (atom['symbol'], \
                                                 args.scale*atom['cart'][0], \
                                                 args.scale*atom['cart'][1], \
                                                 args.scale*atom['cart'][2]))  
    f.write(" \n\n")  

    f.write(" K_POINTS automatic\n")  
    f.write("  %d %d %d 0 0 0\n" % (args.kgrid[0], args.kgrid[1], args.kgrid[2]))  
    f.write(" \n\n")  

    f.write(" CELL_PARAMETERS Angstrom\n")  
    for i in range(3):
      f.write("  %12.10f %12.10f %12.10f\n" % \
                  (args.scale*args.alat[i][0], \
                   args.scale*args.alat[i][1], \
                   args.scale*args.alat[i][2]))  
    f.write(" \n\n")  



#-------------------------------------------------------------------------------
######### ------ MAIN ----- ####################################################

# read args 
args = get_args()
  
# read structure
args.scale, args.alat, args.atoms = read_structure(\
                          args.structure_file, args.structure_file_type, args)


f = open("script.sh", "w")
f.write("#!/bin/bash\n")
f.write("source " + home_bashrc + "\n")
prefix = str(numpy.random.randint(0,1e+9))

pwd = os.getcwd()
if args.calc not  in ["kgrid", "encut", "vacuum"]:

    
  path = args.calc + "/"
  full_prefix = args.title + "_" + prefix + "_" + path
  full_prefix = full_prefix[:-1] 
  os.system("rm -rf %s && mkdir -p %s" % (path, path))
  
  f.write("cd %s\n" % (path))
  f.write(submission_command + " ./run.sh\n")
  f.write("cd %s\n" % (pwd))

  vac_dim = -1
  for d in range(3):
    if args.kgrid[d] == 1:
      vac_dim = d
  if args.vacuum is not None:
    args.alat[vac_dim][vac_dim] += args.vacuum

  generate_qe_file(full_prefix, path, args)
  generate_batch_file(full_prefix, path, args)


else:

  original_kgrid = None
  original_encut = None
  original_vacuum = None

  if args.calc == "encut":
    for delta in range(0,60,10):

      if original_encut is None:
        original_encut = args.encut
      
      path = args.calc  +  "/" + str(original_encut + delta)
      full_prefix = args.title + "_" + prefix + "_" + path
      full_prefix = "-".join(full_prefix.split("/"))
      
      os.system("rm -rf %s && mkdir -p %s" % (path, path))
  
      f.write("cd %s\n" % (path))
      f.write(submission_command + " ./run.sh\n")
      f.write("cd %s\n" % (pwd))
     
      if original_encut + delta > 20: 
        args.encut = original_encut + delta 
      else:
        args.encut = 20

      generate_qe_file(full_prefix, path, args)
      generate_batch_file(full_prefix, path, args)

  if args.calc == "kgrid":
    for delta in range(0,7,2):

      if original_kgrid is None:
        original_kgrid = copy.deepcopy(args.kgrid)
        name_kgrid = max(original_kgrid)
      
      path = args.calc  +  "/" + str(name_kgrid + delta)
      full_prefix = args.title + "_" + prefix + "_" + path
      full_prefix = "-".join(full_prefix.split("/"))
      
      os.system("rm -rf %s && mkdir -p %s" % (path, path))
  
      f.write("cd %s\n" % (path))
      f.write(submission_command + " ./run.sh\n")
      f.write("cd %s\n" % (pwd))
     
      for d in range(3):
        if original_kgrid[d] != 1:
          if original_kgrid[d] + delta > 1:
            args.kgrid[d] = original_kgrid[d] + delta
          else:
            args.kgrid[d] = 1      
   
      generate_qe_file(full_prefix, path, args)
      generate_batch_file(full_prefix, path, args)


  if args.calc == "vacuum":
    for delta in range(0,7,2):

      if original_vacuum is None:
        vac_dim = -1
        for d in range(3):
          if args.kgrid[d] == 1:
            vac_dim = d
        original_vacuum = args.alat[vac_dim][vac_dim]
      
      path = args.calc  +  "/" + str(int(original_vacuum + delta))
      full_prefix = args.title + "_" + prefix + "_" + path
      full_prefix = "-".join(full_prefix.split("/"))
      
      os.system("rm -rf %s && mkdir -p %s" % (path, path))
  
      f.write("cd %s\n" % (path))
      f.write(submission_command + " ./run.sh\n")
      f.write("cd %s\n" % (pwd))
      
      args.alat[vac_dim][vac_dim] = original_vacuum + delta 
    
      generate_qe_file(full_prefix, path, args)
      generate_batch_file(full_prefix, path, args)
        
f.close()

