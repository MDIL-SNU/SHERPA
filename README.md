# Saddle_point_search
Saddle point search combined with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator) and VASP (Vienna Ab initio Simulation Package).
Neural Network Potentials (NNPs) made by [SIMPLE-NN](https://github.com/MDIL-SNU/SIMPLE-NN_v2) are also possible.  

## Requirement
- Intel C compiler
- Intel MPI
- Intel oneAPI MKL
- CMake >= 3.13
- LAMMPS >= 23Jun2022

## Installation
1. Build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd lammps
mkdir build; cd build
cmake ../cmake -D BUILD_SHARED_LIBS=yes
cmake --build . --target install
```
2. Build Check compiler type in `CMakeLists.txt`
```bash
cd Saddle_point_search
mkdir build; cd build
cmake ../src
cmake --build .
```


## INPUT
```text
# general parameter #
NELEMENT    = 2
ATOM_TYPE   = O Pt
INIT_CONFIG = ./POSCAR
TARGET_LIST = ./TARGET
DISP_DIST   = 0.001
ACTI_CUTOFF = 5.1
F_TOL       = 0.01
STDDEV      = 0.1
MAX_STEP    = 0.1
TRIAL_STEP  = 0.001
INIT_RELAX  = 1
INIT_DISP   = 1
CONFIDENCE  = 0.9
MAX_SEARCH  = 1
WRITE_MODE  = 1

# LAMMPS parameter #
PAIR_STYLE  = nn
PAIR_COEFF  = * * potential_saved O Pt
PAIR_CUTOFF = 6.0
NCORE       = 8

# VASP parameter #
VASP_CMD    = mpirun -np $nprocs vasp_std
ISTART      = 1

# dimer parameter #
KAPPA_DIMER = 0
SNC_DIMER   = 0
F_ROT_MIN   = 0.1 
F_ROT_MAX   = 1.0
MAX_NUM_ROT = 4
TRIAL_ANGLE = 45

# art_nouveau parameter #
ART_NOUVEAU = 1
LAMBDA_CRIT = 0.0
LAMBDA_CONV = 0.01
MAX_NUM_RLX = 4
ART_DELAY   = 3
ART_MIXING  = 3

# system parameter #
FREQUENCY   = 1e12
TEMPERATURE = 353

# random parameter #
RANDOM_SEED = -1

# directory parameter #
OUTPUT_DIR  = ./gen_1

# restart parameter #
RESTART     = 0
RESTART_DIR = ./gen_0
```

|Tag|Description|Remark|
|:---|:---|:---|
|NELEMENT|The number of elements||
|ATOM_TYPE|Atomic symbols of elements||
|INIT_CONFIG|Initial configuration file||
|TARGET_LIST|File containing target information||
|DISP_DIST|Finite difference step (=dimer distance)|Angstrom|
|ACTI_CUTOFF|Cutoff radius of active volume|Angstrom|
|F_TOL|Force tolerance for dimer method|eV/Angstrom|
|STDDEV|Standard deviation of gaussian displacement||
|MAX_STEP|Maximum step size of optimization|Angstrom|
|TRIAL_STEP|Trial step size of optimization|Angstrom|
|INIT_RELAX|Initial structure optimization||
|INIT_DISP|Initial structure perturbation||
|CONFIDENCE|Confidence level of saddle point search||
|MAX_SEARCH|Maximum number of saddle point searches||
|WRITE_MODE|Write eigenmode for each step||
|PAIR_STYLE|Pair style for LAMMPS input||
|PAIR_COEFF|Pair coeff for LAMMPS input||
|PAIR_CUTOFF|Cutoff radius of potential file|Angstrom|
|NCORE|The number of cores for each LAMMPS instance||
|VASP_CMD|Command for running VASP||
|ISTART|ISTAT tag in INCAR||
|KAPPA_DIMER|Basin constrained dimer method|Ref.[1](https://doi.org/10.1063/1.4898664)|
|SNC_DIMER|Scaled normal coordinate dimer method|Ref.[2](https://doi.org/10.1016/j.commatsci.2021.110785)|
|F_ROT_MIN|Minimum force criteria for rotation|eV/Angstrom|
|F_ROT_MAX|Maximum force criteria for rotation|eV/Angstrom|
|MAX_NUM_ROT|Maximum number of rotation||
|TRIAL_ANGLE|Trial rotation angle|radian|
|ART_NOUVEAU|Activation and relaxation technique|Ref.[3](http://dx.doi.org/10.1103/PhysRevE.62.7723)|
|LAMBDA_CRIT|Criteria value of inflection region|eV/Angstrom^2|
|LAMBDA_CONV|Convergence criteria value for Lanczos method|eV/Angstrom^2|
|MAX_NUM_RLX|Maximum number of relaxation||
|ART_DELAY|Initial step without Lanczos method||
|ART_MIXING|Mixing step above inflection line||
|FREQUENCY|Attempt frequency of reaction|1/s|
|TEMPERATURE|System temperature|Kelvin|
|RANDOM_SEED|Seed for random number||
|OUTPUT_DIR|Directory for output files||
|RESTART|Restart from previous SPS||
|RESTART_DIR|Directory of previous SPS output||

## TARGET
It contains the target atom indices to be the center of active volume.
```text
I 0 1 2 3
T 1
```

* I: Index
* T: Type
* A: All (not supported yet)
* R: Random (not supported yet)

## Command
If you run SPS with LAMMPS, please use command following:
```bash
mpirun -np $numproc ./SPS_LMP
```
$numproc stands for the number of CPU cores in parallel computation.


Otherwise, execute code without parallel command following:
```bash
./SPS_VASP
```

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. Selective dynamics are also supported.
2. Selective dynamics have priority over `ACTI_CUTOFF`.
3. `NCORE` of 4-8 is recommended with LMP. 
4. The filename of final structure follows Final _ `count` _ `atomic index`.POSCAR
5. For VASP calculation, INCAR, KPOINTS, and POTCAR should be prepared in working directory.
