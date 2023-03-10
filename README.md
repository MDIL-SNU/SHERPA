# Saddle point search
Script for saddle point searches combined with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator) and VASP (Vienna Ab initio Simulation Package).

## Requirement
- Intel oneAPI >= 21.3
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
TARGET_LIST = ./TARGET
FINITE_DIFF = 0.01
ACTI_CUTOFF = 6.0
F_TOL       = 0.01
DIFF_TOL    = 0.02
MAX_MOVE    = 0.1
TRIAL_MOVE  = 0.01
CONFIDENCE  = 0.9
MAX_SEARCH  = 10
WRITE_MODE  = 0

# initial structure parameter #
NELEMENT    = 2
ATOM_TYPE   = O Pt
INIT_CONFIG = ./POSCAR
INIT_RELAX  = 1
INIT_DISP   = 1
DISP_CUTOFF = 3
DISP_STDDEV = 0.1

# LAMMPS parameter #
PAIR_STYLE  = nn
PAIR_COEFF  = * * potential_saved O Pt
NCORE       = 4

# VASP parameter #
VASP_CMD    = mpirun -np $nprocs $vasp_path

# dimer parameter #
KAPPA_DIMER = 0
F_ROT_MIN   = 0.1 
F_ROT_MAX   = 1.0
MAX_NUM_ROT = 4
MAX_NUM_TLS = 500

# art_nouveau parameter #
ART_NOUVEAU = 1
LAMBDA_CRIT = 0.0
LAMBDA_CONV = 0.01
MAX_NUM_RLX = 4
ART_DELAY   = 0
ART_MIXING  = 3

# random parameter #
RANDOM_SEED = -1

# directory parameter #
OUTPUT_DIR  = ./gen_1

# restart parameter #
RESTART     = 0
RESTART_DIR = ./gen_0
```

### General parameter
* **TARGET_LIST** [strings]
  - A path of file containing target atoms
* **FINITE_DIFF** [real]
  - A displacement for the finite difference method (Angs)
* **ACTI_CUTOFF** [real]
  - A cutoff radius of active volume (Angs)
* **F_TOL** [real]
  - A force tolerance for dimer method (eV/Angs)
* **DIFF_TOL** [real]
  - A distance tolerance for identical structures (Angs)
* **MAX_MOVE** [real]
  - A maximum step size of image movement (Angs)
* **TRIAL_MOVE** [real]
  - A trial step size of image for cg optimization (Angs)
* **CONFIDENCE** [real]
  - A confidence level of saddle point search (Ref. [1](https://doi.org/10.1063/1.2976010))
* **MAX_SEARCH** [integer]
  - A maximum number of saddle point searches
* **WRITE_MODE** [integer]
  - An eigenmode for each step
### Initial structure parameter
* **NELEMENT** [integer]
  - The number of elements
* **ATOM_TYPE** [strings]
  - Atomic symbols of elements
* **INIT_CONFIG** [strings]
  - A path of file containing initial atomic positions
* **INIT_RELAX** [integer]
  - An initial structure optimization
* **INIT_DISP** [integer]
  - An initial structure perturbation
* **DISP_CUTOFF** [real]
  - A cutoff radius of initial displacement (Angs)
* **DISP_STDDEV** [real]
  - A standard deviation of gaussian displacement
### LAMMPS parameter
* **PAIR_STYLE** [strings]
  - Pair style in LAMMPS input
* **PAIR_COEFF** [strings]
  - Pair coeff in LAMMPS input
* **NCORE** [integer]
  - The number of cores for each LAMMPS instance
### VASP parameter
* **VASP_CMD** [strings]
  - A command for running VASP
### Dimer parameter
* **KAPPA_DIMER** [integer]
  - Basin constrained dimer method (Ref.[2](https://doi.org/10.1063/1.4898664))
* **F_ROT_MIN** [real]
  - A minimum force criteria for rotation (eV/Angs)
* **F_ROT_MAX** [real]
  - A maximum force criteria for rotation (eV/Angs)
* **MAX_NUM_ROT** [integer]
  - A maximum number of rotation steps
* **MAX_NUM_TLS** [integer]
  - A maximum number of translation steps
* **TRIAL_ANGLE** [real]
  - A trial rotation angle (degree)
### ART nouveau parameter
* **ART_NOUVEAU** [integer]
  - Activation and relaxation technique (Ref.[3](http://dx.doi.org/10.1103/PhysRevE.62.7723))
* **LAMBDA_CRIT** [real]
  - A criteria value of inflection points (eV/Angs^2)
* **LAMBDA_CONV** [real]
  - A convergence criteria value for Lanczos method (eV/Angs^2)
* **MAX_NUM_RLX** [integer]
  - A maximum number of orthogonal relaxation steps below inflection points
* **ART_DELAY** [integer]
  - A number of initial steps without Lanczos method
* **ART_MIXING** [integer]
  - A number of mixing steps above inflection points
### Random parameter
* **RANDOM_SEED** [integer]
  - A seed for random number
### Directory parameter
* **OUTPUT_DIR** [strings]
  - A directory for output files
### Restart parameter
* **RESTART** [integer]
  - A restart from previous SPS
* **RESTART_DIR** [strings]
  - A directory of previous SPS output

## TARGET
The TARGET file lists either the target atom indices or types that serve as the center of the active volume, with each condition being added sequentially.
```text
I 0 1 2 3
T 1
A
```

* I: Index (starting from 0)
* T: Type (starting from 1)
* A: All

If `R` is appended after a single character, the indices of target atoms are randomly shuffled.
```text
IR 0 1 2 3
TR 1
AR
```

* +R: random shuffle

## Command
If you run SPS with LAMMPS, please use the command following:
```bash
mpirun -np $numproc ./SPS_LMP
```
$numproc stands for the number of CPU cores in parallel computation.


Otherwise, execute code with the command following:
```bash
./SPS_VASP
```

## Outputs
### SPS_`count`.log
It contains information such as steps, potential energy, eigenvalue, or curvature and so on.
Its format depends on the algorithm you used.
In addition, the type of the saddle point, barrier energy, reaction energy, and elapsed time are written.

### SPS_`count`.XDATCAR
The trajectory from the initial to the saddle point are appended into here.

### `count`.MODECAR
The eigenvector at the saddle point is written as N x 3 array, where N is the number of all atoms.

### Saddle_`count`_ `atomic_index`.POSCAR
The configuration of the saddle point is written as POSCAR format.

### Final_`count`_ `atomic_index`.POSCAR
If the saddle point is connected and unique, the final configuration is written as "Final_ `count`_ `atomic_index`.POSCAR".


If the saddle point is connected and not unique, the final configuration is written as "`previous_count`_ Final_ `count`_ `atomic_index`.POSCAR", where `previous_count` stands for the `count` of the same saddle point.


If the saddle point is not splited, the final configurations are written as "x0_Final_`count`_ `atomic_index`.POSCAR".


If the saddle point is disconnected, the final configurations are written as "x1_Final_`count`_ `atomic_index`.POSCAR" and "x2_Final_`count`_ `atomic_index`.POSCAR".

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format, which supports `Selective dynamics`.
2. For VASP calculation, INCAR, KPOINTS, and POTCAR should be located in the current directory.
3. To provide specific initial vector, set `RESTART` of 1 and provide "Final_ `count`_ `atomic index`.POSCAR" and "`count`.MODECAR" in `RESTART_DIR`.
4. To use hybrid potentials, write down all `pair_coeff`s within `PAIR_COEFF` with the delimiter (|) between them.
For example,
```text
PAIR_COEFF = 1 1 potential_saved_1 | 1 2 potential_saved_2
```
