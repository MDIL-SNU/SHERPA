To run SHERPA package with ASE Python module (`sherpa_ase`), five input files are required.

1. INPUT
2. TARGET
3. POSCAR
4. interatomic potential
5. ase_calc.py

Here, we provide the input files with checkpoint of pretrained [SevenNet](https://github.com/MDIL-SNU/SevenNet) model. 

> [!NOTE]
> ASE calculators module can be loaded by modifying `ase_calc.py` only if the module can be found in `CMAKE_PREFIX_PATH` defined in the build.

> [!TIP]
> Multiple GPU cards allows simultaneous multiple saddle point search with the same numbers of `ntasks-per-node` and `gres:gpu` (see `run.sh`).
