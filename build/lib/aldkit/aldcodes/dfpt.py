
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

  parser = argparse.ArgumentParser(description='dfpt: wrapper for dfpt calc')
  
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
                          help = "Energy-cutoff for pw.x", \
                          dest = "encut")
  parser.add_argument("-irr", "--irr_per_job", default = None, \
                          type = int, 
                          help = "num irr per job", \
                          dest = "irr_per_job")
  parser.add_argument("-v", "--vacuum", default = None, \
                          type = float, 
                          help = "vacuum to-be adjusted", \
                          dest = "vacuum")
  parser.add_argument("-k", "--kgrid", default = None, \
                          type = int, \
                          nargs = 3, \
                          help = "scf kgrid for pw.x", \
                          dest = "kgrid")
  parser.add_argument("-qpt", "--qgrid", default = None, \
                          type = int, \
                          nargs = 3, \
                          help = "phonon dfpt qgrid for ph.x", \
                          dest = "qgrid")
  parser.add_argument("-t", "--title", \
                          default = "MyTitle", type = str, \
                          help = "title for job",
                          dest = "title")
  parser.add_argument("-c", "--calc", \
                          default = None, type = str, \
                          required = True, \
                          help = "type-of-calc to be performed", \
                          choices = ["pre", "mid", "post", "irr"],
                          dest = "calc")

  args = parser.parse_args(sys.argv[1:])
  
  if args.structure_file is None:
    args.structure_file = "save/dfpt"
    print("Using %s as the structure_file." % (args.structure_file))
  
  if args.structure_file_type is None:
    args.structure_file_type = "ald"
    print("Using %s as the structure_file_type." % (args.structure_file_type))

  if args.structure_file_type != "ald" and args.qgrid is None:
    parser.error("missing qgrid! Please specify -qpt/--qgrid argument!")

  if args.calc == "pre":
    if args.encut == None:
      parser.error("Encut needed!")
  
  if args.calc == "pre":
    if args.kgrid == None:
      parser.error("kgrid needed!")
  
  if args.calc == "mid":
    if args.irr_per_job == None:
      parser.error("irr needed!")
  
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

      # read qgrid now
      args.qgrid = [int(_) for _ in lines[5+natom].split()]

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
def generate_run_file(prefix, path, args, sq, lq, sr, lr, jobstring):

  if args.calc == "pre":
    filename = "%s/run_pre.sh" % (path)
  if args.calc == "mid":
    filename = "%s/run_%d_%d_%d_%d.sh" % (path, sq, lq, sr, lr)
  if args.calc == "post":
    filename = "%s/run_post.sh" % (path)
  
  with open(filename, "w") as f:
    f.write("#!/bin/bash\n")
    for module in modules:
      f.write("module load %s\n" % (module))
    f.write("export OMP_NUM_THREADS=%d\n" % (OMP_NUM_THREADS))
    f.write("\n")

    if args.calc == "pre":
      f.write("\n/usr/bin/rm -rf " + outdir + "/%s\n" % (prefix))
      f.write("mkdir -p " + outdir + "/%s\n" % (prefix))
      f.write("\n")
    
      f.write(f"\n{parallel} -n {args.ppn} pw.x -in pw.in | tee pw.out\n")
      f.write(f"\n{parallel} -n {args.ppn} ph.x -in ph_pre.in | tee ph_pre.out\n")
      f.write("\n")
      f.write("mkdir -p data\n") 
      f.write("mkdir -p data/dvscf\n") 
      f.write("cp -r %s/%s/* data/\n" % (outdir, prefix))
      f.write("\n")
     
      f.write("\n/usr/bin/rm -rf %s/%s\n" % (outdir,prefix))


    if args.calc == "mid":
      f.write(jobstring)
      f.write("\n\npython /home/MatSim/a_jain/softwares/"
                                        "SpaceTime/mydel.py -f %s\n" % (path))  
          
    if args.calc == "post":
      f.write("mkdir -p %s/%s/\n" % (outdir, prefix))
      f.write("cp -r data/* %s/%s/\n\n" % (outdir, prefix))
      
      f.write(f"\n{parallel} -n {args.ppn} ph.x -in ph_post.in | tee ph_post.out\n\n")
      f.write(f"\n{parallel} -n 1 q2r.x -in q2r.in | tee q2r.out\n\n")
      f.write("/usr/bin/rm -rf %s/%s\n" % (outdir, prefix))
      


#----------- generate job string combining all jobs --------------------------
def create_job_string(prev_string, sq, lq, sr, lr):
  string = prev_string

  string += '#' * 80
  string += '\ntodo=false\n'
  string += 'if [ ! -f "start_%d_%d_%d_%d" ]\n' % (sq,lq,sr,lr)
  string += "then\n"
  string += "  todo=true\n"
  string += "else\n"
  string += '  if ! grep -q "JOB DONE" ph_%d_%d_%d_%d.out\n' % (sq,lq,sr,lr)
  string += "  then\n"
  string += '    if test `find "ph_%d_%d_%d_%d.out" -mmin +100`\n' \
                                                          % (sq,lq,sr,lr)
  string += "    then\n"
  string += "      todo=true\n"
  string += "    fi\n"
  string += "  fi\n"
  string += "fi\n"

  string += 'if [ "$todo" = true ] ; then\n'
  string += '  touch start_%d_%d_%d_%d\n' % (sq, lq, sr, lr)
  string += "  mkdir -p %s/%s/%d_%d_%d_%d\n" % (outdir, prefix, sq, lq, sr, lr)
  string += "  cp -r data/* %s/%s/%d_%d_%d_%d/\n" % \
                                        (outdir, prefix, sq, lq, sr, lr)
  string += f"  {parallel} -n {args.ppn} ph.x -in ph_{sq}_{lq}_{sr}_{lr}.in | tee ph_{sq}_{lq}_{sr}_{lr}.out\n\n")
  
  for iq in range(sq, lq+1):
    for ir in range(sr, lr+1): 
      string += "  cp -r %s/%s/%d_%d_%d_%d/_ph0/%s.phsave/dynmat.%d.%d.xml data/_ph0/%s.phsave/\n" % \
               (outdir, prefix, sq, lq, sr, lr, prefix, iq, ir, prefix)
          
      # 0th rep induced charge contribution
      if ir == 1:
        string += "  cp -r %s/%s/%d_%d_%d_%d/_ph0/%s.phsave/dynmat.%d.%d.xml data/_ph0/%s.phsave/\n" % \
               (outdir, prefix, sq, lq, sr, lr, prefix, iq, 0, prefix)
            
        #  get electric filed part
        if iq == 1:
          string += "  cp -r %s/%s/%d_%d_%d_%d/_ph0/%s.phsave/tensors.xml data/_ph0/%s.phsave/\n" % \
               (outdir, prefix, sq, lq, sr, lr, prefix, prefix)
                
    if iq == 1:
      string += "  cp %s/%s/%d_%d_%d_%d/_ph0/%s.dfpt.dvscf1 data/dvscf/dvscf_%d_%d_%d_%d\n" % \
                                          (outdir, prefix, sq, lq, sr, lr, \
                                                        prefix, sq, lq, sr, lr)
    else:
      string += "  cp %s/%s/%d_%d_%d_%d/_ph0/%s.q_%d/%s.dfpt.dvscf1 data/dvscf/dvscf_%d_%d_%d_%d\n" % \
                                   (outdir, prefix, sq, lq, sr, lr, prefix, iq,\
                                                        prefix, sq, lq, sr, lr)
                     
  string += "  /usr/bin/rm -rf %s/%s/%d_%d_%d_%d\n" % \
                                               (outdir, prefix, sq, lq, sr, lr)
  
  string += "fi\n\n"
  
  return string            


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
    f.write("  tstress=.true.\n")  
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
    f.write("  conv_thr=1.0d-12\n")  
    f.write("  mixing_beta=0.3\n")  
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


# create q2rx input file  ------------------------------------------------------
def generate_q2r_file(prefix, path, args):

  with open("%s/q2r.in" % (path), "w") as f:
    f.write(" &input\n")  
    f.write("  zasr='simple'\n")
    f.write("  flfrc='harmonic_flfrc.dat'\n")
    f.write("  fildyn='dfpt.dyn'\n")
    f.write(" /\n\n")
    


# create ph.x input file  ------------------------------------------------------
def generate_ph_file(prefix, path, args, sq, lq, sr, lr):
  
  if args.calc == "pre":
    filename = "%s/ph_pre.in" % (path)
  if args.calc == "mid":
    filename = "%s/ph_%d_%d_%d_%d.in" % (path, sq, lq, sr, lr)
  if args.calc == "post":
    filename = "%s/ph_post.in" % (path)
 
  with open(filename, "w") as f:
    f.write("title='%s'\n" % (args.title))  
    f.write(" &inputph\n")  
    f.write("  tr2_ph=1.0d-22\n")
    f.write("  alpha_mix(1)=0.1\n")
    f.write("  ldisp=.true.\n")
    f.write("  search_sym=.false.\n")
    
    f.write("  trans=.true.\n")
    #f.write("  epsil=.true.\n")
    
    if args.calc != "pre":
      f.write("  recover=.true.\n")
    
    if args.calc != "post": 
      f.write("  start_irr=%d\n" % (sr))
      f.write("  last_irr=%d\n" % (lr))
      f.write("  start_q=%d\n" % (sq))
      f.write("  last_q=%d\n" % (lq))
    
    f.write("  nq1=%d\n" % (args.qgrid[0]))
    f.write("  nq2=%d\n" % (args.qgrid[1]))
    f.write("  nq3=%d\n" % (args.qgrid[2]))
    f.write("  prefix='" + prefix + "'\n")
    f.write("  fildyn='dfpt.dyn'\n")
    f.write("  fildvscf='dfpt.dvscf'\n")

    if args.calc == "pre":
      f.write("  outdir='%s/%s'\n" % (outdir, prefix))
    if args.calc == "mid":
      f.write("  outdir='%s/%s/%d_%d_%d_%d'\n" % \
                                          (outdir, prefix, sq, lq, sr, lr))
    if args.calc == "post":
      f.write("  outdir='%s/%s'\n" % (outdir, prefix))

    f.write(" /\n\n")
 

# -------------- read num irreducible reps
def read_num_irr(prefix, path):

  # read number of representations for different cases
  import xml
  import xml.etree
  import xml.etree.ElementTree
  doc = xml.etree.ElementTree.parse("%s/data/_ph0/%s.phsave/control_ph.xml" % \
                                                                (path, prefix))
  text = doc.find('Q_POINTS').find('Q-POINT_COORDINATES').text.split() 
  text = [float(_) for _ in text]
  nqpts = int(len(text)/3)
  qpts = []
  for i in range(nqpts):
    qpts.append([text[3*i + 0], text[3*i+1], text[3*i+2]])
  num_irr = []
  for i in range(len(qpts)):  
    doc = xml.etree.ElementTree.parse(\
                                    "%s/data/_ph0/%s.phsave/patterns.%d.xml" % \
                                                           (path, prefix, i+1))
    num_irr.append(int(doc.find('IRREPS_INFO').\
                                find('NUMBER_IRR_REP').text.split()[0]))

  return num_irr
     

#-------------------------------------------------------------------------------
######### ------ MAIN ----- ####################################################

# read args 
args = get_args()
  
# read structure
args.scale, args.alat, args.atoms = read_structure(\
                          args.structure_file, args.structure_file_type, args)


script = open("script.sh", "w")
script.write("#!/bin/bash\n")
#script.write("source " + home_bashrc + "/.bashrc\n")

pwd = os.getcwd()
path = "dfpt_%d_%d_%d/" % (args.qgrid[0], args.qgrid[1], args.qgrid[2])
if args.calc == "pre":
  os.system(f"rm -rf {pwd}/{path}")
  os.mkdir(f'{pwd}/{path}')

  # generate prefix
  prefix = str(numpy.random.randint(0,1e+9))
  prefix = path[:-1] + prefix
  with open(path + "/prefix", "w") as f:
    f.write(prefix)

  # generate pw file
  generate_pw_file(prefix, path, args)
  generate_ph_file(prefix, path, args, 1, 1, 0, 0)
  generate_run_file(prefix, path, args, 1, 1, 0, 0, "")
      
  script.write("cd %s\n" % (path))
  script.write(submission_command + " ./run_pre.sh\n")
  script.write("cd %s\n" % (pwd))


if args.calc == "irr":
  
  # read prefix
  with open(path + "/prefix", "r") as f:
    prefix = f.readline().split()[0]

  num_irr = read_num_irr(prefix, path)
  for irr in num_irr:
    print(irr)


if args.calc == "mid":

  # read prefix
  with open(path + "/prefix", "r") as f:
    prefix = f.readline().split()[0]
  
  num_irr = read_num_irr(prefix, path)

  job_string = ""
  # now create job files
  for iq in range(len(num_irr)):
    for ir in range(0, num_irr[iq], args.irr_per_job):
      
      sr = ir + 1
      qpt = iq + 1
      if (ir + args.irr_per_job) <= num_irr[iq]:
        lr = ir + args.irr_per_job
      else:
        lr = num_irr[iq]  

      generate_ph_file(prefix, path, args, qpt, qpt, sr, lr)
      job_string = create_job_string(job_string, qpt, qpt, sr, lr)
  
  # generate run files now
  for iq in range(len(num_irr)):
    for ir in range(0, num_irr[iq], args.irr_per_job):
      
      sr = ir + 1
      qpt = iq + 1
      if (ir + args.irr_per_job) <= num_irr[iq]:
        lr = ir + args.irr_per_job
      else:
        lr = num_irr[iq]  

      generate_run_file(prefix, path, args, qpt, qpt, sr, lr, job_string)

      script.write("cd %s\n" % (path))
      script.write(submission_command + " ./run_%d_%d_%d_%d.sh -f %s\n" % \
                                        (qpt, qpt, sr, lr, path))
      script.write("cd %s\n" % (pwd))
      break
    break
  


if args.calc == "post":

  # read prefix
  with open(path + "/prefix", "r") as f:
    prefix = f.readline().split()[0]

  
  generate_ph_file(prefix, path, args, 1, 1, 0, 0)
  generate_q2r_file(prefix, path, args)
  generate_run_file(prefix, path, args, 1, 1, 0, 0, "")
  
  script.write("cd %s\n" % (path))
  script.write(submission_command + " ./run_post.sh\n")
  script.write("cd %s\n" % (pwd))
  
  # create dvscf files now
  num_irr = read_num_irr(prefix, path)
  num_modes = 3*len(args.atoms)
  
  # --- list of all dvscf files
  dvscf_files = os.listdir("%s/data/dvscf/" % (path))
  dvscf_files.sort(key = lambda entry: [int(_) for _ in entry.split("_")[1:]]) 
 

  #----- filesize of each mode
  statinfo = os.stat("%s/data/dvscf/%s" % (path, dvscf_files[-1]))
  per_mode_size = statinfo.st_size / num_modes
  
  # --- combine now
  prev_q = 0
  modes_done = 0
  for dvscf_file in dvscf_files:
    q = int(dvscf_file.split('_')[1])
    
    if prev_q != q:
      prev_q = q
      modes_done = 0
      
    statinfo = os.stat("%s/data/dvscf/%s" % (path, dvscf_file))
    modes_todo = int(statinfo.st_size / per_mode_size)
    modes_todo -=  modes_done
    
    print(dvscf_file)
    os.system("dd bs=%d count=%d skip=%d seek=%d if=%s/data/dvscf/%s "
                                                "of=%s/dfpt.dvscf_q_%d" % \
                   (per_mode_size, modes_todo, modes_done, modes_done, \
                                                    path, dvscf_file, path, q))
    modes_done += modes_todo

script.close()

