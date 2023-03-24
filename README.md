# Saddle point search
The scripts for saddle point searches combined with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator) and VASP (Vienna Ab initio Simulation Package).

## Requirement
- Intel oneAPI >= 21.3
- CMake >= 3.13
- LAMMPS and VASP are tested by version of 23Jun2022 and 6.3.2, respectively.

## Build
Three executables can be built.
- SPS_LMP: saddle point searches script through LAMMPS C-API
- SPS_VASP: saddle point searches script by reading/writing VASP
- EXTRACTOR: script for postprocessing of outputs (see `Outputs`)

0. If `SPS_LMP` is needed, build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd lammps
mkdir build; cd build
cmake ../cmake -D BUILD_SHARED_LIBS=yes
cmake --build . --target install
```

1. Make `build` directory and configure with `CMakeLists.txt` in `src` directory.
```bash
cd Saddle_point_search
mkdir build; cd build
cmake ../src
```

2. Build targets through **one of following commands**. 
```bash
# all (SPS_LMP, SPS_VASP, and EXTRACTOR)
cmake --build .
# SPS_LMP
cmake --build . --target SPS_LMP
# SPS_VASP
cmake --build . --target SPS_VASP 
# EXTRACTOR
cmake --build . --target EXTRACTOR
```

## INPUT
### General parameter
* **TARGET_LIST** [strings]
  - *TARGET_LIST* indicates the file that contains the index of target atoms.
* **FINITE_DIFF** [real]
  - *FINITE_DIFF* sets the displacement in finite difference method (in Angst).
* **ACTI_CUTOFF** [real]
  - *ACTI_CUTOFF* sets the cutoff radius of active volume (in Angst).
* **ACTI_NEVERY** [int]
  - *ACTI_NEVERY* sets the interval step to check the active volume.
* **F_TOL** [real]
  - *F_TOL* sets the force tolerance for saddle point searches and relaxation (in eV/Angst).
* **DIFF_TOL** [real]
  - *DIFF_TOL* sets the distance tolerance for each atom in identical structures (in Angst).
* **MAX_MOVE** [real]
  - *MAX_MOVE* sets the maximum step size of image movement (in Angst).
* **TRIAL_MOVE** [real]
  - *TRIAL_MOVE* sets the trial step size of image for cg optimization (in Angst).
* **CONFIDENCE** [real]
  - *CONFIDENCE* sets termination condition through confidence level (Ref. [1](https://doi.org/10.1063/1.2976010)).
* **MAX_SEARCH** [integer]
  - *MAX_SEARCH* sets termination condition through the number of saddle point searches.
### Initial structure parameter
* **NELEMENT** [integer]
  - *NELEMENT* is the number of elements
* **ATOM_TYPE** [strings]
  - *ATOM_TYPE* is the symbols of elements
* **INIT_CONFIG** [strings]
  - *INIT_CONFIG* indicates the structure file written as POSCAR format.
* **INIT_RELAX** [0/1]
  - *INIT_RELAX* determines whether or not to relax *INIT_CONFIG* before saddle point searches.
* **INIT_DISP** [0/1]
  - *INIT_DISP* determines whether or not to displace structure before saddle point searches.
* **DISP_CUTOFF** [real]
  - *DISP_CUTOFF* sets the cutoff radius of displaced region (in Angst).
* **DISP_STDDEV** [real]
  - *DISP_STDDEV* sets the standard deviation of gaussian displacement vector (in Angst).
* **INIT_MODE** [0/1]
  - *INIT_MODE* determines whether or not to provide the initial eigenmode.
* **MODE_LIST** [strings]
  - *MODE_LIST* indicates the file that contains the eigenmodes.
### LAMMPS parameter
* **PAIR_STYLE** [strings]
  - *PAIR_STYLE* stands for the pair style in LAMMPS input.
* **PAIR_COEFF** [strings]
  - *PAIR_COEFF* stands for the pair coeff in LAMMPS input.
* **NCORE** [integer]
  - *NCORE* sets the number of cores for each LAMMPS instance.
### VASP parameter
* **VASP_CMD** [strings]
  - *VASP_CMD* is the command to run VASP.
### Dimer parameter
* **KAPPA_DIMER** [0/1]
  - *KAPPA_DIMER* activates the basin constrained dimer method (Ref.[2](https://doi.org/10.1063/1.4898664)).
* **F_ROT_MIN** [real]
  - *F_ROT_MIN* sets the minimum force criteria for rotation (in eV/Angst).
* **F_ROT_MAX** [real]
  - *F_ROT_MAX* sets the maximum force criteria for rotation (in eV/Angst).
* **MAX_NUM_ROT** [integer]
  - *MAX_NUM_ROT* sets the maximum number of rotation steps in dimer method.
* **MAX_NUM_TLS** [integer]
  - *MAX_NUM_TLS* sets the maximum number of translation steps in dimer method.
### ART nouveau parameter
* **ART_NOUVEAU** [0/1]
  - *ART_NOUVEAU* activates the activation and relaxation technique (Ref.[3](http://dx.doi.org/10.1103/PhysRevE.62.7723)).
* **LAMBDA_CONV** [real]
  - *LAMBDA_CONV* sets the convergence criteria value for Lanczos method (in eV/Angs^2).
* **MAX_NUM_RLX** [integer]
  - *MAX_NUM_RLX* sets the maximum number of perpendicular relaxation steps at eigenvalue less than *LAMBDA_CRIT*.
* **ART_DELAY** [integer]
  - *ART_DELAY* sets the number of initial steps without Lanczos method.
* **ART_MIXING** [integer]
  - *ART_MIXING* sets the number of mixing steps above inflection points.
* **HYPER_RLX** [integer]
  - *HYPER_RLX* sets the number of relaxation steps, where the configuration is on the hyperplane that is spanned by eigenmode and displacement vector.
### Random parameter
* **RANDOM_SEED** [-1/nonnegative integer]
  - *RANDOM_SEED* sets the seed for random number generation. *RANDOM_SEED* is generated automatically when provided -1.
### Directory parameter
* **OUTPUT_DIR** [strings]
  - *OUTPUT_DIR* indicates the output directory.

## TARGET
The TARGET file lists either the target atom indices or types that serve as the center of the active volume, with each condition being added sequentially.
If `R` is appended after a single character, the indices of target atoms are randomly shuffled.
```text
# indices of 0th, 1st, 2nd, and 3rd atoms
I 0 1 2 3
# randomly shuffled indices of type 1 atoms
TR 1
# indices of all atoms
A
```

* I: index (starting from 0)
* T: type (starting from 1)
* A: all
* +R: random shuffle


## Usage
To run scripts, use the following commands:
```bash
# run SPS_LMP
mpirun -np ${numproc} SPS_LMP
# run SPS_VASP
SPS_VASP
```
, where `${numproc}` stands for the number of processors.


## Outputs (will be updated soon)
||Unconverged|Not splited|Disconnected|Connected|
|:---|:---:|:---:|:---:|:---:|
|SPS.log|O|O|O|O|
|SPS.XDATCAR|O|O|O|O|
|SPS.MODECAR|X|O|O|O|
|Saddle.POSCAR|X|O|O|O|
|Final.POSCAR|X|X|X|O|


All outputs have header like `{count}_{index}`, meaning that `{count}` and `{index}` are the number of saddle point searches trials and the index of target atom, respectively. 
`EXTRACTOR` executable helps to specific `{count}` of configuration and trajectories from outputs.
```bash
EXTRACTOR ${output} ${count}
```
, where `${output}` indicates one of output files.


### SPS.log
*SPS.log* contains information of saddle point searches such as steps, potential energy, and curvature.
*SPS.log* also provides the type of the saddle point, barrier energy, reaction energy, and elapsed time.

### SPS.XDATCAR
*SPS.XDATCAR* is the trajectory of configuration during saddle point searches.

### SPS.MODECAR
*SPS.MODECAR* is the eigenvector at the saddle point.

### Saddle.POSCAR
*Saddle.POSCAR* is the configuration of the saddle point.

### Final.POSCAR
*Final.POSCAR* is the configuration of the final point (=next minimum).


## Miscellaneous
1. Fixed atoms in `INIT_CONFIG` will not move during saddle point searches.
2. INCAR, KPOINTS, and POTCAR should be prepared in the current directory for `SPS_VASP`.
3. To use hybrid potentials for LAMMPS, append `pair_coeff`s with the delimiter (|) among them.
For example,
```text
PAIR_COEFF = 1 1 potential_saved_1 | 1 2 potential_saved_2
```
