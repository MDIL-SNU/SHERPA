# SHERPA
**SHERPA** (**S**addle points **H**unting based on **E**nergy surface for **R**eaction **PA**thways) is the script for finding saddle points of atomic reaction with minimum-mode following methods such as the dimer method[1](https://doi.org/10.1063/1.480097) and activation-relaxation technique nouveau (ARTn)[2](http://dx.doi.org/10.1103/PhysRevE.62.7723). This script can interface with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator) and VASP (Vienna Ab initio Simulation Package).
<p align="center">
<img src="./src/logo.png" width="200"/>
</p>

## Requirement
- Intel oneAPI >= 21.3
- CMake >= 3.13
- LAMMPS and VASP are tested by version of 23Jun2022 and 6.3.2, respectively.

## Build
Three executables can be built.
- SHERPA_LMP: saddle point searches script through LAMMPS C-API
- SHERPA_VASP: saddle point searches script by reading/writing VASP
- EXTRACTOR: script for postprocessing of outputs (see `Outputs`)
- KMC: script for kMC simulation utilizing SHERPA

1. If `SHERPA_LMP` is needed, build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd lammps
mkdir build; cd build
cmake ../cmake -D BUILD_SHARED_LIBS=yes
cmake --build . --target install
```

2. Make `build` directory and configure with `CMakeLists.txt` in `src` directory.
```bash
cd SHERPA
mkdir build; cd build
cmake ../src
```
If cmake fails to find LAMMPS package automatically, please check environment variables of system or specify the path of LAMMPS in `src/LMP/CMakeLists.txt` manually.

3. Build targets through **one of following commands**.
```bash
# All (SHERPA_LMP, SHERPA_VASP, EXTRACTOR, and KMC)
cmake --build .
# Only SHERPA_LMP
cmake --build . --target SHERPA_LMP
# Only SHERPA_VASP
cmake --build . --target SHERPA_VASP
# Only EXTRACTOR
cmake --build . --target EXTRACTOR
# Only KMC
cmake --build . --target KMC
```

## Usage (SHERPA)
Before running scripts, *INPUT*, *POSCAR*, and *TARGET* should be prepared.
*POSCAR* is an initial structure file written in VASP5 format, which supports `Selective dynamics`.
To run scripts, use the following commands:
```bash
# SHERPA_LMP
mpirun -np ${numproc} SHERPA_LMP
# SHERPA_VASP
SHERPA_VASP
```
, where `${numproc}` stands for the number of processors.


### INPUT
#### General parameter
* **ACTI_CUTOFF** [real, 5.0 (default)]
  - *ACTI_CUTOFF* sets the cutoff radius of active volume (in Angst).
* **ACTI_NEVERY** [integer, 3 (default)]
  - *ACTI_NEVERY* sets the interval step to expand the active volume (in Angst).
* **FINITE_DIFF** [real, 0.01 (default)]
  - *FINITE_DIFF* sets the displacement in finite difference method (in Angst).
* **F_TOL** [real, 0.01 (default)]
  - *F_TOL* sets the force tolerance for saddle point searches and relaxation (in eV/Angst).
* **DIFF_TOL** [real, 0.4 (default)]
  - *DIFF_TOL* sets the distance tolerance for each atom in identical structures (in Angst).
* **MAX_MOVE** [real, 0.1 (default)]
  - *MAX_MOVE* sets the maximum step size of image movement (in Angst).
* **TRIAL_MOVE** [real, 0.01 (default)]
  - *TRIAL_MOVE* sets the trial step size of image for cg optimization (in Angst).
* **MAX_SEARCH** [integer, 100 (default)]
  - *MAX_SEARCH* sets termination condition through the number of saddle point searches.
* **WRITE_TRAJ** [0/1, 1 (default)]
  - *WRITE_TRAJ* determines whether or not to write the saddle point search trajectories.
* **CONTINUE** [0/1, 0 (default)]
  - *CONTINUE* determines whether or not to continue SHERPA from previous results. *Statistics.log* and *Event.log* should be prepared.
#### Initial structure parameter
* **NELEMENT** [integer]
  - *NELEMENT* is the number of elements
* **ATOM_TYPE** [strings]
  - *ATOM_TYPE* is the symbols of elements
* **INIT_RELAX** [0/1, 1 (default)]
  - *INIT_RELAX* determines whether or not to relax *INIT_CONFIG* before saddle point searches.
* **INIT_DISP** [0/1, 0 (default)]
  - *INIT_DISP* determines whether or not to displace structure before saddle point searches.
* **DISP_CUTOFF** [real, 5.0 (default)]
  - *DISP_CUTOFF* sets the cutoff radius, where initial displacement is defined (in Angst).
* **DISP_MOVE** [real, 0.0 (default)]
  - *DISP_MOVE* sets the magnitude of the initial displacement vector (in Angst).
* **INIT_MODE** [0/1, 0 (default)]
  - *INIT_MODE* determines whether or not to provide the initial eigenmode. The initial eigenmode can be provided by the file named *Initial.MODECAR*.
#### LAMMPS parameter (for SHERPA_LMP)
* **PAIR_STYLE** [strings]
  - *PAIR_STYLE* stands for the pair style in LAMMPS input.
* **PAIR_COEFF** [strings]
  - *PAIR_COEFF* stands for the pair coeff in LAMMPS input.
* **NCORE** [integer]
  - *NCORE* sets the number of cores for each LAMMPS instance.
#### VASP parameter (for SHERPA_VASP)
* **VASP_CMD** [strings]
  - *VASP_CMD* is the command to run VASP.
#### Dimer parameter
* **KAPPA_DIMER** [0/1, 0 (default)]
  - *KAPPA_DIMER* activates the basin-constrained dimer method ([3](https://doi.org/10.1063/1.4898664)).
* **F_ROT_MIN** [real, 0.1 (default)]
  - *F_ROT_MIN* sets the minimum force criteria for rotation (in eV/Angst).
* **F_ROT_MAX** [real, 1.0 (default)]
  - *F_ROT_MAX* sets the maximum force criteria for rotation (in eV/Angst).
* **MAX_NUM_ROT** [integer, 4 (default)]
  - *MAX_NUM_ROT* sets the maximum number of rotation steps in dimer method.
* **MAX_NUM_TLS** [integer, 500 (default)]
  - *MAX_NUM_TLS* sets the maximum number of translation steps in dimer method.
#### ART nouveau parameter
* **ART_NOUVEAU** [0/1, 1 (default)]
  - *ART_NOUVEAU* activates the activation and relaxation technique.
* **LAMBDA_CONV** [real, 0.01 (default)]
  - *LAMBDA_CONV* sets the convergence criteria value for Lanczos method (in eV/Angs^2).
* **MAX_NUM_ITR** [integer, 500 (default)]
  - *MAX_NUM_ITR* sets the maximum number of iteration in ARTn method.
* **MAX_NUM_RLX** [integer, 1 (default)]
  - *MAX_NUM_RLX* sets the maximum number of perpendicular relaxation steps at eigenvalue less than *LAMBDA_CRIT*.
* **DELAY_STEP** [integer, 0 (default)]
  - *DELAY_STEP* sets the number of initial steps without Lanczos method.
* **MIXING_STEP** [integer, 0 (default)]
  - *MIXING_STEP* sets the number of mixing steps above inflection points.
* **HYPER_STEP** [integer, 3 (default)]
  - *HYPER_STEP* sets the number of relaxation steps, where the configuration is on the hyperplane that is spanned by eigenmode and displacement vector.
#### Random parameter
* **RANDOM_SEED** [unsigned integer, $RANDOM (default)]
  - *RANDOM_SEED* sets the seed for random number generation. If *RANDOM_SEED* is not specified, random seed is generated randomly.

### TARGET
The TARGET file lists the conditions for the center atoms of the sphere. The conditions are appended sequentially and independently. Each one consists of numbers following characters. Three characters (I, T, A) and one additive (R) are supported currently. Three characters should be used once in one condition. Fixed atoms cannot be the center atoms.

* I: index (starting from 0)
* T: type (starting from 1)
* A: all
* +R: random shuffle

If the characters contain `R`, the indices of atoms are randomly shuffled.
For example,
```text
# indices of 0th, 1st, 2nd, and 3rd atoms
I 0 1 2 3
# randomly shuffled indices of type 1 atoms
TR 1
# indices of all atoms
A
```

### Outputs
#### Event.log
*Event.log* contains `{count}`, barrier energy and frequency of unique reactions.

#### Statistics.log
*Statistics.log* shows a current progress of reaction finding. The number of unique reaction, relevant (=connected) reaction, and trials are written.

#### Redundancy.log
*Redundancy.log* tells the which reactions are duplicated.

Not all saddle points are written in all output files. The rules of the outputs are summarized below table.
||Unconverged|Not splited|Disconnected|Connected|
|:---|:---:|:---:|:---:|:---:|
|SHERPA.log|O|O|O|O|
|SHERPA.XDATCAR|O|O|O|O|
|SHERPA.MODECAR|X|O|O|O|
|Saddle.POSCAR|X|O|O|O|
|Final.POSCAR|X|X|X|O|

#### SHERPA.log
*SHERPA.log* contains information of saddle point searches such as steps, potential energy, and curvature.
*SHERPA.log* also provides the type of the saddle point, barrier energy, reaction energy, and elapsed time.

#### Time.log
*Time.log* says how long the SHERPA takes to run.

#### SHERPA.XDATCAR
*SHERPA.XDATCAR* is the trajectory of configuration during saddle point searches.

#### SHERPA.MODECAR
*SHERPA.MODECAR* is the eigenvector at the saddle point.

#### Saddle.POSCAR
*Saddle.POSCAR* is the configuration of the saddle point.

#### Final.POSCAR
*Final.POSCAR* is the configuration of the final point (=next minimum).

## Usage (EXTRACTOR)
The below outputs have header like `{count}_{index}`, meaning that `{count}` and `{index}` are the number of reaction finding trials and the index of target atom, respectively.
`EXTRACTOR` executable helps to specific `{count}` of configuration and trajectories from outputs through the follow command.
```bash
EXTRACTOR ${output} ${count}
```
, where `${output}` indicates one of output files.

## Usage (KMC)
To run KMC script, several settings can be provided by in-command.
```bash
KMC --sherpa_cmd "{sherpa_command}" --att_freq 1e13 --temperature 298 --end_step 10000 --inputs_path "./INPUTS"
```
Note that "{sherpa_command}" is positional (and mandatory) parameter.

For example,
```bash
KMC "mpirun -np 16 SHERPA_LMP" --temperature 500 --end_step 1000
```

In *inputs_path*, input files for SHERPA should be located.
For example,
``` text
--- KMC
 |- INPUTS
      |- POSCAR
      |- TARGET
      |- potential files
      |- ...
```

### Options
* **sherpa_cmd** [string]
  - *sherpa_cmd* is the command to run SHERPA.
* **att_freq** [real, 1e13 (default)]
  - *att_freq* sets the attempt frequency of reaction rates (in Hertz).
* **temperature** [real, 298 (default)]
  - *temperature* sets the temperature of system (in Kelvin).
* **end_step** [int, 10000 (default)]
  - *end_step* sets the maximum step in kMC simulation.
* **inputs_path** [string, "./INPUTS" (default)]
  - *inputs_path* indicates the directory that contains input files for SHERPA.


## Miscellaneous
1. INCAR, KPOINTS, and POTCAR should be prepared in the current directory for `SHERPA_VASP`.
2. To use hybrid potentials for LAMMPS, append `pair_coeff`s with the delimiter (|) among them.
For example,
```text
PAIR_COEFF = 1 1 potential_saved_1 | 1 2 potential_saved_2
```
