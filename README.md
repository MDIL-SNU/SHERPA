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
DIMER_DIST  = 0.0001
F_TOL       = 0.01
F_ROT_MIN   = 0.1
F_ROT_MAX   = 1.0
MAX_NUM_ROT = 100
TRIAL_ANGLE = 45
DISP_CUTOFF = 5.1
STDDEV      = 0.1
MAX_STEP    = 0.1
TRIAL_STEP  = 0.001
INIT_RELAX  = 1
INIT_MODE   = 0

# random parameter #
RANDOM_SEED = 0

# output parameter #
OUTPUT_DIR  = ./OUTPUT
```

|Tag|Description|Units|
|:---|:---|:---|
|NELEMENT|The number of elements||
|ATOM_TYPE|Atomic symbols of elements||
|PAIR_STYLE|Pair style for LAMMPS input||
|PAIR_COEFF|Pair coeff for LAMMPS input||
|CUTOFF|Cutoff radius of potential file|Angstrom|
|INIT_CONFIG|Initial configuration file||
|TARGET_LIST|File containing target information||
|DIMER_DIST|Dimer distance from the center|Angstrom|
|F_TOL|Force tolerance for dimer method|eV/Angstrom|
|F_ROT_MIN|Minimum force criteria for rotation|eV/Angstrom|
|F_ROT_MAX|Maximum force criteria for rotation|eV/Angstrom|
|MAX_NUM_ROT|The maximum number of rotation||
|TRIAL_ANGLE|Trial rotation angle|radian|
|DISP_CUTOFF|Radius of displacement sphere||
|STDDEV|Standard deviation of gaussian displacement||
|MAX_STEP|Maximum step size of translation|Angstrom|
|TRIAL_STEP|Trial step size of translation|Angstrom|
|INIT_MODE|User defined initial eigenmode||
|INIT_RELAX|Initial structure optimization||
|RANDOM_SEED|Seed for random number||
|OUTPUT_DIR|Directory for output files||

## TARGET
```bash
I 0 1 2 3
```

* I: Index
* A: All (not supported yet)
* R: Random (not supported yet)

## Command
```bash
mpirun -np $numproc ./LAMMPS_dimer
```
$numproc stands for the number of CPU cores in parallel computation.

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. Selective dynamics are also supported.
2. `numproc` in command should be the multiple of nimages. 
