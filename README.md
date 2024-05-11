# SHERPA
**SHERPA** (**S**addle points **H**unting based on **E**nergy surface for **R**eaction **PA**thways) is the script for finding saddle points of atomic reaction with minimum-mode following methods such as the dimer method[1](https://doi.org/10.1063/1.480097) and activation-relaxation technique nouveau (ARTn)[2](http://dx.doi.org/10.1103/PhysRevE.62.7723). This script can interface with LAMMPS (Large Atomic/Molecular Massively Parallel Simulator), VASP (Vienna Ab initio Simulation Package), and ASE[3](https://doi.org/10.1088/1361-648X/aa680e) (Atomic Simulation Environment) calculator.
* **DIMER** [0/1, 0 (default)]
  - *DIMER* activates the original dimer method ([4](https://doi.org/10.1063/1.480097)).
* **KAPPA_DIMER** [0/1, 0 (default)]
  - *KAPPA_DIMER* activates the basin-constrained dimer method ([5](https://doi.org/10.1063/1.4898664)).
* **ART_NOUVEAU** [0/1, 1 (default)]
  - *ART_NOUVEAU* activates the activation and relaxation technique. ([6](https://doi.org/10.1103/PhysRevE.62.7723))
<p align="center">
<img src="./assets/logo.png" width="200"/>
</p>

## Requirement
- CMake >= 3.13
- LAMMPS, VASP, ASE are tested by version of 23Jun2022, 6.3.2, and 3.22.1, respectively.

## Installation
Five executables are supported.
- sherpa_lmp: saddle point searches script through LAMMPS C-API
- sherpa_vasp: saddle point searches script via reading/writing VASP
- sherpa_ase: saddle point searches script with embedded ASE python calculator
- extractor: script for postprocessing of outputs (see `Outputs`)
- kmc: script for kMC simulation combined with SHERPA

1. Build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
Skip this step, if `sherpa_lmp` is not needed.
```bash
cd lammps
mkdir build; cd build
cmake ../cmake -D BUILD_SHARED_LIBS=yes
cmake --build . --target install
```

2. Make `build` directory and configure with `CMakeLists.txt`.
```bash
cd SHERPA
mkdir build; cd build
cmake ../
```
The path for Python package for `sherpa_ase` can be specified through `CMAKE_PREFIX_PATH`.
For example, if virtual environment is managed by Anaconda, python interpreter would be located in `/path/.conda/envs/{name}/bin/python`. Then,
```bash
cmake ../ -D CMAKE_PREFIX_PATH=/path/.conda/envs/{name}
```

3. Build and install.
```bash
# specific executable
cmake --build . --target {name of executable}
# all executables
cmake --build . --target install
```
The built executable will be copied to `SHERPA/bin`.

## Usage (sherpa)
*INPUT*, *POSCAR*, and *TARGET* are required to run the `sherpa` script.
*POSCAR* is an initial structure file written in VASP5 format, which supports `Selective dynamics`.

**Note**
`sherpa_vasp` requires additional files to run VASP such as INCAR, KPOINTS, and POTCAR. `sherpa_ase` also needs a python file that defines ase calculator. An example python file is provided as `ase_calc.py`.

Use the following commands:
```bash
# sherpa_lmp and sherpa_ase
mpirun -np ${numproc} sherpa_lmp
# sherpa_vasp
sherpa_vasp
```
, where `${numproc}` stands for the number of processors.
The number of processors for VASP can be defined in `VASP_CMD` in `INPUT`.

### INPUT
#### General
* **ALGORITHM** [Dimer, kappa-dimer, art_nouveau | art_nouveau (default)]
  - *ALGORITHM* sets which algorithm to be used for saddle point search. 
* **ACTI_CUTOFF** [real | 5.0 (default)]
  - *ACTI_CUTOFF* sets the cutoff radius of active volume (in Angst).
* **ACTI_NEVERY** [integer | 3 (default)]
  - *ACTI_NEVERY* sets the interval step to expand the active volume (in Angst).
* **FINITE_DIFF** [real | 0.01 (default)]
  - *FINITE_DIFF* sets the displacement in finite difference method (in Angst).
* **F_TOL** [real | 0.01 (default)]
  - *F_TOL* sets the force tolerance for saddle point searches and relaxation (in eV/Angst).
* **DIFF_TOL** [real | 0.4 (default)]
  - *DIFF_TOL* sets the distance tolerance for each atom in identical structures (in Angst).
* **MAX_MOVE** [real | 0.1 (default)]
  - *MAX_MOVE* sets the maximum step size of image movement (in Angst).
* **TRIAL_MOVE** [real | 0.01 (default)]
  - *TRIAL_MOVE* sets the trial step size of image for cg optimization (in Angst).
* **MAX_SEARCH** [integer | 100 (default)]
  - *MAX_SEARCH* sets termination condition through the number of saddle point searches.
* **WRITE_TRAJ** [True/False | True (default)]
  - *WRITE_TRAJ* determines whether or not to write the saddle point search trajectories.
* **CONTINUE** [True/False | False (default)]
  - *CONTINUE* determines whether or not to continue SHERPA from previous results. *Statistics.log* and *Event.log* should be prepared.
#### Initial structure
* **NELEMENT** [integer]
  - *NELEMENT* is the number of elements
* **ATOM_TYPE** [strings]
  - *ATOM_TYPE* is the symbols of elements
* **INIT_RELAX** [True/False | True (default)]
  - *INIT_RELAX* determines whether or not to relax *INIT_CONFIG* before saddle point searches.
* **INIT_DISP** [True/False | False (default)]
  - *INIT_DISP* determines whether or not to displace structure before saddle point searches.
* **DISP_CUTOFF** [real | 5.0 (default)]
  - *DISP_CUTOFF* sets the cutoff radius, where initial displacement is defined (in Angst).
* **DISP_MOVE** [real | 0.0 (default)]
  - *DISP_MOVE* sets the magnitude of the initial displacement vector (in Angst).
* **INIT_MODE** [True/False | False (default)]
  - *INIT_MODE* determines whether or not to provide the initial eigenmode. The initial eigenmode can be provided by the file named *Initial.MODECAR*.
#### LAMMPS (only work for sherpa_lmp)
* **PAIR_STYLE** [strings]
  - *PAIR_STYLE* stands for the pair style in LAMMPS input.
* **PAIR_COEFF** [strings]
  - *PAIR_COEFF* stands for the pair coeff in LAMMPS input.
* **NCORE** [integer]
  - *NCORE* sets the number of cores for each LAMMPS instance.
#### VASP (only work for sherpa_vasp)
* **VASP_CMD** [strings]
  - *VASP_CMD* is the command to run VASP.
#### ASE (only work for sherpa_ase)
* **ASE_CALC** [strings]
  - *ASE_CALC* indicates the python file containing ASE calculator.
* **MODEL_PATH** [strings]
  - *MODEL_PATH* indicates the potential file for ASE calculator.
#### Dimer (only work for dimer and kappa-dimer method)
* **F_ROT_MIN** [real | 0.1 (default)]
  - *F_ROT_MIN* sets the minimum force criteria for rotation (in eV/Angst).
* **F_ROT_MAX** [real | 1.0 (default)]
  - *F_ROT_MAX* sets the maximum force criteria for rotation (in eV/Angst).
* **MAX_NUM_ROT** [integer | 4 (default)]
  - *MAX_NUM_ROT* sets the maximum number of rotation steps in dimer method.
* **MAX_NUM_TLS** [integer | 500 (default)]
  - *MAX_NUM_TLS* sets the maximum number of translation steps in dimer method.
#### ART nouveau (only work for art_nouveau)
* **LAMBDA_CONV** [real | 0.01 (default)]
  - *LAMBDA_CONV* sets the convergence criteria value for Lanczos method (in eV/Angs^2).
* **MAX_NUM_ITR** [integer | 500 (default)]
  - *MAX_NUM_ITR* sets the maximum number of iteration in ARTn method.
* **MAX_NUM_RLX** [integer | 1 (default)]
  - *MAX_NUM_RLX* sets the maximum number of perpendicular relaxation steps at eigenvalue less than *LAMBDA_CRIT*.
* **DELAY_STEP** [integer | 0 (default)]
  - *DELAY_STEP* sets the number of initial steps without Lanczos method.
* **MIXING_STEP** [integer | 0 (default)]
  - *MIXING_STEP* sets the number of mixing steps above inflection points. ([7](https://doi.org/10.1021/acs.jctc.0c00541))
* **HYPER_STEP** [integer | 3 (default)]
  - *HYPER_STEP* sets the number of relaxation steps, where the relaxation trajectory is constrained on the hyperplane orthogonal to eigenmode and displacement vectors. (TBA)
#### Random
* **RANDOM_SEED** [unsigned integer | $RANDOM (default)]
  - *RANDOM_SEED* sets the seed for random number generation. If *RANDOM_SEED* is not specified, random seed is generated randomly.

**Note**
In the case of *ALGORITHM* and boolean tags, the first letter determines the applied algorithm.

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

## Usage (extractor)
The below outputs have header like `{count}_{index}`, meaning that `{count}` and `{index}` are the number of reaction finding trials and the index of target atom, respectively.
`EXTRACTOR` executable helps to specific `{count}` of configuration and trajectories from outputs through the follow command.
```bash
extractor ${output} ${count}
```
, where `${output}` indicates one of output files.

## Usage (kmc)
To run KMC script, several settings can be provided by in-command.
```bash
kmc --sherpa_cmd "{sherpa_command}" --att_freq 1e13 --temperature 298 --kmc_step 10000 --inputs_path "./INPUTS"
```
In *inputs_path*, input files for the `sherpa` should be located (INPUT, POSCAR, and TARGET).


### Options
* **sherpa_cmd** [string]
  - *sherpa_cmd* is the command to run SHERPA.
* **att_freq** [real | 1e13 (default)]
  - *att_freq* sets the attempt frequency of reaction rates (in Hertz).
* **temperature** [real | 298 (default)]
  - *temperature* sets the temperature of system (in Kelvin).
* **kmc_step** [int | 10000 (default)]
  - *kmc_step* sets the maximum step in kMC simulation.
* **low_cut** [real | 0.1 (default)]
  - *low_cut* sets the minimum activation energy to be in event table.
* **restart** [int | 0 (default)]
  - *restart* sets the initial step to continue from the previous kMC simulation.
* **inputs_path** [string | "./INPUTS" (default)]
  - *inputs_path* indicates the directory that contains input files for SHERPA.
* **random_seed** [unsigned int | $RANDOM (default)]
  - *random_seed* sets the seed for random number generation. If *random_seed* is not specified, random seed is generated randomly.
