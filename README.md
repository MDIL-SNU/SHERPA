# Saddle_point_search
Dimer method combined with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator).  
Neural Network Potentials (NNPs) made by [SIMPLE-NN](https://github.com/MDIL-SNU/SIMPLE-NN_v2) are also possible.  

## Requirement
- CMake >= 3.10
- LAMMPS >= 23Jun2022

## Installation
1. Build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd lammps
mkdir build; cd build
cmake ../cmake -D BUILD_SHARED_LIBS=yes
cmake --build . --target install
```
2. Check compiler type in `CMakeLists.txt`
```bash
cd Saddle_point_search/src
```
```text
SET (CMAKE_C_COMPILER "mpiicc" CACHE PATH "")
```
If you don't have intel compiler, change "mpiicc" to "mpicc".


3. Build and install
``` bash
cmake CMakeLists.txt
make
```

## INPUT
```text
# potential parameter #
NELEMENT    = 2
ATOM_TYPE   = O Pt
PAIR_STYLE  = nn
PAIR_COEFF  = * * potential_saved O Pt
CUTOFF      = 6.0

# kMC parameter #
TEMPERATURE = 353.0
ATT_FREQ    = 1e12
END_TIME    = 1e4
END_STEP    = 1000000

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
CONFIDENCE  = 0.9

# random parameter #
RANDOM_SEED = -1

# output parameter #
OUTPUT_DIR  = ./OUTPUT

# parallelism parameter #
NCORE       = 8
```

|Tag|Description|Units|
|:---|:---|:---|
|NELEMENT|The number of elements||
|ATOM_TYPE|Atomic symbols of elements||
|PAIR_STYLE|Pair style for LAMMPS input||
|PAIR_COEFF|Pair coeff for LAMMPS input||
|CUTOFF|Cutoff radius of potential file|Angstrom|
|TEMPERATURE|System temperature|Kelvin|
|ATT_FREQ|Attempt frequency of reaction|1/s|
|END_TIME|Termination condition for kMC time|s|
|END_STEP|Termination condition for kMC step||
|INIT_CONFIG|Initial configuration file||
|TARGET_LIST|File containing target information||
|DIMER_DIST|Dimer distance from the center|Angstrom|
|F_TOL|Force tolerance for dimer method|eV/Angstrom|
|F_ROT_MIN|Minimum force criteria for rotation|eV/Angstrom|
|F_ROT_MAX|Maximum force criteria for rotation|eV/Angstrom|
|MAX_NUM_ROT|The maximum number of rotation||
|TRIAL_ANGLE|Trial rotation angle|radian|
|DISP_CUTOFF|Radius of displacement sphere|Angstrom|
|STDDEV|Standard deviation of gaussian displacement||
|MAX_STEP|Maximum step size of translation|Angstrom|
|TRIAL_STEP|Trial step size of translation|Angstrom|
|INIT_RELAX|Initial structure optimization||
|CONFIDENCE|Confidence level of event table||
|RANDOM_SEED|Seed for random number||
|OUTPUT_DIR|Directory for output files||
|NCORE|The number of cores for each dimer method||

## TARGET
```text
I 0 1 2 3
```

* I: Index
* T: Type (not supported yet)
* A: All (not supported yet)
* R: Random (not supported yet)

## Command
```bash
mpirun -np $numproc ./SPS
```
$numproc stands for the number of CPU cores in parallel computation.

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. Selective dynamics are also supported.
2. `NCORE` of 4-8 is recommended. 
3. Set DISP_CUTOFF shorter than PAIR_CUTOFF.
