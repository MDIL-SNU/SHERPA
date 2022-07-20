# LAMMPS_dimer
Dimer method combined with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator).  
Neural Network Potentials (NNPs) made by [SIMPLE-NN](https://github.com/MDIL-SNU/SIMPLE-NN_v2) are also possible.  

## Requirement
- CMake >= 2.8.12
- LAMMPS == 29Oct2020

## Installation
1. Build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd /path/to/lammps/src
make mpi
make mode=shlib mpi
```
2. Set the environment variable. [[link](https://docs.lammps.org/Build_link.html)]
```bash
export LD_LIBRARY_PATH=/path/to/lammps/src:$LD_LIBRARY_PATH
```
3. Modify `LMP_PATH` in `CMakeLists.txt`.
```text
SET ( LMP_PATH /path/to/lammps )
```
4. Make `Makefile` and build
``` bash
cmake CMakeLists.txt
make
```

## INPUT
```bash
# potential parameter #
NELEMENT    = 2
ATOM_TYPE   = O Pt
PAIR_STYLE  = nn
PAIR_COEFF  = * * potential_saved O Pt
CUTOFF      = 6.0

# dimer parameter #
INIT_CONFIG = ./POSCAR
TARGET      = ./TARGET
DIMER_DIST  = 0.001
FTOL        = 0.01
F_ROT_MIN   = 0.01
F_ROT_MAX   = 1.0
MAX_NUM_ROT = 1
TRIAL_ANGLE = 45
DISP_CUTOFF = 4.4
STDDEV      = 0.1
MAX_STEP    = 0.1
TRIAL_STEP  = 0.001
INIT_MODE   = 0

# random parameter #
RANDOM_SEED = 0
```

|Tag|Description|Units|
|:---|:---|:---|
|NELEMENT|||
|ATOM_TYPE|||
|PAIR_STYLE|||
|PAIR_COEFF|||
|CUTOFF|||
|INIT_CONFIG|||
|TARGET_LIST|||
|DIMER_DIST|||
|F_TOL|||
|F_ROT_MIN|||
|F_ROT_MAX|||
|MAX_NUM_ROT|||
|TRIAL_ANGLE|||
|DISP_CUTOFF|||
|STDDEV|||
|MAX_STEP|||
|TRIAL_STEP|||
|INIT_MODE|||
|RANDOM_SEED|||

## TARGET
```bash
I 0 1 2 3
```

TBA

## Command
```bash
mpirun -np $numproc ./LAMMPS_dimer
```
$numproc stands for the number of CPU cores in parallel computation.

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. 
2. `numproc` in command should be the multiple of nimages. 
