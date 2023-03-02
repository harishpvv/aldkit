
def main():
  print('''--------  ALD Parameters Description  -----------------------
path_restart
  Required: Optional
  Type: open_string
  Default: restart
  Description: path where restart files will be read/write on disk

path_save
  Required: Optional
  Type: open_string
  Default: save
  Description: path where save files will be read/write on disk

path_output
  Required: Optional
  Type: open_string
  Default: output
  Description: path where output files will be written on disk

calculation
  Required: Optional
  Type: string
  Options: [forceconstant, harmonic, phph, helmholtz]
  Default: phph
  Description: Type of calculation to be done

symprecision
  Required: Optional
  Type: int
  Default: 6
  Description: number of precision digits in symmetry calculations

precision
  Required: Optional
  Type: int
  Default: 10
  Description: number of precision digits in calculations

temperature
  Required: Yes
  Type: double
  Default: None
  Description: Temperature in Kelvin

ntime
  Required: Optional
  Type: double
  Default: 1.0e-12
  Description: non-dimensional time. All time in input file should be normalized using this value.

nlength
  Required: Optional
  Type: double
  Default: 1.0e-10
  Description: non-dimensional length. All length in input file should be normalized using this value.

nmass
  Required: Optional
  Type: double
  Default: 1.66054e-27
  Description: non-dimensional mass. All mass in input file should be normalized using this value.

alat
  Required: Yes
  Type: composite
  Default: None
  Description: Lattice vectors. Single keyword 'alat' on first line followed by three lines listing a1, a2, and a3 [each line listing x,y,z component of lattice vectors] resp.

coordinates
  Required: Optional
  Type: string
  Options: [direct, reduced, cartesian, absolute]
  Default: direct
  Description: atomic basis is in reduced/direct coordinates or cartesian/absoulte coordinates

natom
  Required: Yes
  Type: composite
  Default: None
  Description: basis for structure. The units are specified by coordinates tag. Format is: on first line 'natom' = number_of_atoms, followed by natom atom entries on subsequent lines. Each atom entry is [ Symbol_i {string}, px_i {double}, py_i {double}, pz_i {double}, mass_i {double}, tag_i {int, optional} ], where symbol_i, px_i, py_i, pz_i, mass_i, and tag_i are symbol, coordinates, mass and tag of concerned atom. 
PLEASE NOTE that only symbol/pos are involved in the symmetry search of lattice/FCs, etc

iprim
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if use primitive cell or not for all calculations that does not involve supercell [prim is not used for DFT force calculation]

cutfrc
  Required: Yes
  Type: composite
  Default: None
  Description: Cutoff used in force contant extraction. Format is [specie-1 specie-2 2nd-order, 3rd-order, 4th-order].

cutfermi
  Required: Optional
  Type: composite
  Default: None
  Description: Cutoff used in actual dispersion, scattering, etc calcualtions. Format is [specie-1 specie-2 2nd-order, 3rd-order, 4th-order].

potential
  Required: Yes
  Type: string
  Options: [espresso, vasp, empirical, internal]
  Default: None
  Description: Type of calculator used in getting forces.cell2, cell3, & cell4 are needed for 'espresso', 'vasp', and 'empirical'.if ftype='md', then only empirical allowed.

cell4
  Required: Optional
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: Suercell-size for 4th order FC extraction. Specified numbers should are integral copies of 'alat' in each lattice vector direction

cell3
  Required: Optional
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: Suercell-size for 3rd order FC extraction. Specified numbers should are integral copies of 'alat' in each lattice vector direction

cell2
  Required: Optional
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: supercell size for harmonic force calculations. Specified numbers should are integral copies of 'alat' in each lattice vector direction

ftype
  Required: Optional
  Type: string
  Options: [finite-difference, distribution, md]
  Default: finite-difference
  Description: type of forces; from finite-difference,  distribution, or MD-traj

fsolve
  Required: Optional
  Type: string
  Options: [simple, taylor-series]
  Default: simple
  Description: solution method for force-constants. one eqn at a time (simple) or taylor-series

fulldisp
  Required: Yes
  Type: string
  Options: [pre, full, disp]
  Default: None
  Description: Use: 'pre', 'disp', 'full' steps with ftype='distribituon' and 'disp', 'full' for other ftypes. With 'pre', only harmonic-files are generated so that can calculate harmonic forces to get  frequenceis as need in ftype='distribution'

normal_umklapp
  Required: Optional
  Type: string
  Options: [normal, umklapp, both]
  Default: both
  Description: Use: 'normal', 'umklapp', 'both' types of scattering in the computation of scattering rates.

nproc_calculator
  Required: Optional
  Type: int
  Default: nPE [all processes]
  Description: number of procs to be used in the calculator class [memory heavy operations]

species
  Required: Optional
  Type: composite
  Default: * 1.0 1.0
  Description: isotope information for species. Format is nspecie lines of: [Specie_symbol_i mass_i fraction_i], where i loops over all specie types.

ddisp
  Required: Optional
  Type: double
  Default: 0.025
  Description: delta-disp for displacement of atoms in Force-Constant calculations [suggested value: 0.025 A]

maxeqn
  Required: Optional
  Type: int
  Default: 10000
  Description: Maxium num-eqns for fsolve='taylor-series'

itrans
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if enforce Translational Invariance on FCs

harmonic_invariance
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if apply translation invariance on harmonic FCs?

fermi_invariance
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if apply translation invariance on FCs within Fermi-cutoff? Should have itrans='yes' in order to work.

deltaq
  Required: Optional
  Type: double
  Default: 0.00001
  Description: dq to be used in the harmonic velocity calc by  a numerical derivative of dynamical matrix

ndisp
  Required: Optional
  Type: int
  Default: 100
  Description: number of displaced-cells for getting FCs using ftype='distribution'

udisp
  Required: Optional
  Type: int
  Default: 100
  Description: number of displaced-cells for actually used for extracting 3/4-order FC using ftype='distribution/md'. Should be <= ndisp for ftype='distribution'

split_fit
  Required: Optional
  Type: yes/no
  Default: no
  Description: if fit [2nd], 3rd and 4th order FC separately in force-disp data?

equilibrium_force_fit
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if fit equilibrium forces while extracting [2nd], 3rd, and 4th order FCs from force-disp data?

scale_force_fit
  Required: Optional
  Type: double
  Default: 100
  Description: 4th-order will be scaled by this scale^3.

fourshift
  Required: Optional
  Type: yes/no
  Default: no
  Description: if include quartic contribution in the calculation of freqshift from Dyson/Green-function approach. If using renomalize, make sure tht freqshift4 is no. Else would count quartic correction twice!

freqshift
  Required: Optional
  Type: yes/no
  Default: no
  Description: if calculate freqshift from Dyson/Green-function approach. Would include cubic and quartic shifts

renormalize
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if renormalize FCs

epsil_anharmonic_renormalize
  Required: Optional
  Type: double
  Default: 0.01
  Description: convergence threshold for renormalize-anharmonic-FC. (maximum % change in total norm of any order-FC).

epsil_renormalize
  Required: Optional
  Type: double
  Default: 0.01
  Description: convergence threshold for renormalize-FC. (maximum % change in any harmonic-FC).

betamix_renormalization
  Required: Optional
  Type: double
  Default: 0.5
  Description: mix-fraction between prev and current step for renormalization. given value is frac of current step.

fqcell
  Required: Yes
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: phonon-grid size (fine-grid) for full-phonon calculutions.

occupation
  Required: Optional
  Type: string
  Options: [quantum, classical]
  Default: quantum
  Description: Statistic to be used for phonon occupation

igruneissen
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if compute Gruneissen parameters

ifourphonon
  Required: Optional
  Type: yes/no
  Default: no
  Description: if compute 4-phonon scattering

iscattering_phase_space_only
  Required: Optional
  Type: yes/no
  Default: no
  Description: if calc. only scattering phase-space?

start_q
  Required: Optional
  Type: int
  Default: 0
  Description: phonon-id to which start the phph calc from. Useful insplitting trivially parallel computations of symmetry  reduced qpoints over multiplt jobs.

end_q
  Required: Optional
  Type: int
  Default: -1
  Description: phonon-id at which end the phph calc from. Useful insplitting trivially parallel computations of symmetry  reduced qpoints over multiplt jobs.

parallel3ph
  Required: Optional
  Type: string
  Options: [inner, outer]
  Default: outer
  Description: in 3-ph scattering, whether to do paralleization over q1 [outer] or q2 [inner] loop for  q1 + q2 -> q3. Outer use less memory and is fast but save data only after finish for q1. So inner is beneficial for systems with large atoms.

parallel4ph
  Required: Optional
  Type: string
  Options: [inner, outer]
  Default: inner
  Description: in 4-ph scattering, whether to do paralleization over q1 [outer] or q2 [inner] loop for  q1 + q2 + q3 -> q4. Outer use less memory and is fast but save data only after finish for q1. So inner is beneficial for systems with large atoms.

fourphonon_maxfreq
  Required: Optional
  Type: double
  Default: 1000.0 [no cutoff]
  Description: Maximum freq. mode calculated in 4ph-scattering (THz units)

coherent_bte
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if calculate coherent (off-diagonal) contribution to thermal conductivity

iterate_bte
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if do iterative (full, i.e., beyong RTA) solution of BTE

epsil_iter_bte
  Required: Optional
  Type: double
  Default: 1.0
  Description: convergence threshold for iterative-BTE. (in percent change of thermal-K.

betamix_iter_bte
  Required: Optional
  Type: double
  Default: 0.5
  Description: mix-fraction between prev and current step for iterative-BTE. given value is frac of current step.

itime_reversal_symmetry
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if use time-reversal symmetry in reducing BZ 

broadening_function
  Required: Optional
  Type: string
  Options: [gaussian, lorentzian]
  Default: gaussian
  Description: Type of broadening to be used for delta functions in ph-ph scatterings

energy_delta_cutoff3
  Required: Optional
  Type: double
  Default: 3
  Description: Range of energies (around) mean/target upto which energy delta function in 3rd order ph-ph processes is considered non-zero. Specified in terms of standard deviation.

energy_delta_cutoff4
  Required: Optional
  Type: double
  Default: 2
  Description: Range of energies (around) mean/target upto which energy delta function in 4th order ph-ph processes is considered non-zero. Specified in terms of standard deviation.

energy_delta_scale3
  Required: Optional
  Type: double
  Default: 0.1
  Description:  for 3-ph processes

energy_delta_scale4
  Required: Optional
  Type: double
  Default: 0.1
  Description:  for 4-ph processes

fc_threshold2
  Required: Optional
  Type: double
  Default: 0.00000001
  Description: For 2nd order FCs, if the valus of FC is less than  fc_threshold2 times maximum(abs(FC_2)) value, then value of given FC is set to zero.

fc_threshold3
  Required: Optional
  Type: double
  Default: 0.0001
  Description: For 3rd order FCs, if the valus of FC is less than  fc_threshold3 times maximum(abs(FC_3)) value, then value of given FC is set to zero.

fc_threshold4
  Required: Optional
  Type: double
  Default: 0.01
  Description: For 4th order FCs, if the valus of FC is less than  fc_threshold4 times maximum(abs(FC_4)) value, then value of given FC is set to zero.

isoscat
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if consider isotope scattering

energy_delta_cutoff_isotope
  Required: Optional
  Type: double
  Default: 4
  Description: Range of energies (around) mean/target upto which energy delta function in isotope-scattering processes is considered non-zero. Specified in terms of standard deviation.

energy_delta_scale_isotope
  Required: Optional
  Type: double
  Default: 0.1
  Description:  for ph-isotope scattering processes

ibdry
  Required: Optional
  Type: yes/no
  Default: no
  Description: if consider boundary scattering

lbdry
  Required: Optional
  Type: double
  Default: 1e+9
  Description: characteristic length for phonon-bdry scattering

ielph
  Required: Optional
  Type: yes/no
  Default: no
  Description: if computer el-ph scattering

fkcell
  Required: Optional
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: electron-grid size (fine-grid) for full-electron calculutions [with ielph=True].

ckcell
  Required: Optional
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: electron-grid size (coarse-grid) for elph calculutions [with ielph=True].

cqcell
  Required: Optional
  Type: int_array
  Format: [int, int, int]
  Default: None
  Description: phonon-grid size (coarse-grid) for elph calculutions [with ielph=True].

helmholtz_quartic
  Required: Optional
  Type: yes/no
  Default: no
  Description: if calculate quartic contribution to Helmholtz Energy

helmholtz_cubic
  Required: Optional
  Type: yes/no
  Default: yes
  Description: if calculate cubic contribution to Helmholtz Energy

idirdisp
  Required: Optional
  Type: yes/no
  Default: no
  Description: if compute directional dispersion

dirtype
  Required: Optional
  Type: string
  Options: [direct, cartesian]
  Default: direct
  Description: directions k-points are in reduced or cartesian coordinates

ndirdisp
  Required: Optional
  Type: composite
  Default: 0
  Description: Number of directions to be computed! The format is:
   ndirdisp number of lines with each line having format:
     start_x start_y start_z end_x end_y end_x num
   where start_x, ..., end_z are qpts coordinates and num is the number of points in that direction.

---------------------------------------------------------
---------------------------------------------------------
''')
