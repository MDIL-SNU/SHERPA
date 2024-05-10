from typing import Tuple, List, Optional

from ase import Atoms
from ase.optimize import LBFGS
from ase.constraints import FixAtoms


# modify here!
def ase_initialize(model_path: Optional[str] = None) -> None:
    global calculator

    from ase.calculators.abc import ABC

    calculator = ABC(model_path)


def oneshot(
    atomic_numbers: List[int], positions: List[List[float]], cell: List[List[float]]
) -> Tuple[float, List[List[float]]]:
    global calculator

    atoms = Atoms(
        numbers=atomic_numbers, positions=positions, cell=cell, pbc=[True] * 3
    )
    atoms.calc = calculator
    energy = atoms.get_potential_energy(force_consistent=True)
    forces = atoms.get_forces().tolist()

    return energy, forces


def atom_relax(
    atomic_numbers: List[int],
    positions: List[List[float]],
    cell: List[List[float]],
    fix: List[bool],
    ftol: float,
) -> Tuple[float, List[List[float]]]:
    global calculator

    atoms = Atoms(
        numbers=atomic_numbers, positions=positions, cell=cell, pbc=[True] * 3
    )
    atoms.constraints = FixAtoms(mask=fix)
    atoms.calc = calculator
    opt = LBFGS(atoms)
    opt.run(fmax=ftol, steps=10000)

    energy = atoms.get_potential_energy(force_consistent=True)
    pos = atoms.get_positions().tolist()

    return energy, pos
