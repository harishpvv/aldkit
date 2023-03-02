# ALD

ALD package is a tool to parse the output of ALD code which is written by Prof. Ankit Jain, IIT Bombay.
ALD code is to solve the linearized Boltzman transport equation for three, four phonon and electron-phonon scattering rates from abinitio methods.
For more details contact a_jain@iitb.ac.in, harishpvv@gmail.com 

import package as `from aldkit.ald import ALD`
ALD class instance has to be created as `ald = ALD(outdir='path')`
Then methods can be called as `ald.get_frequencies()`

Cite: Phonon properties and thermal conductivity from first principles, lattice dynamics, and the Boltzmann transport equation, A Jain et al, https://doi.org/10.1063/1.5064602 
