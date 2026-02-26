# uv --directory /cluster/home/amsoni/g-xtb-ase run python /cluster/home/amsoni/gxTB_HOME/test_gxtb.py > test_gxtb.out 2>&1

import os
from pathlib import Path
from ase.build import molecule
from ase.io import read
from gxtb_ase import GxTB

# Resolve the directory where this script lives, so output files are written
# here regardless of the process working directory (e.g. when using uv --directory).
SCRIPT_DIR = str(Path(__file__).resolve().parent)

atoms = read("structures/Ag_max_203atoms.xyz")

# Attach g-xTB calculator
calc = GxTB(
    label="Ag_max_203atoms",
    charge=0,
    spin=0,
    write_log=True,
    directory=SCRIPT_DIR,
)  # Singlet state (closed shell)
atoms.calc = calc

energy = atoms.get_potential_energy()
print(f"Ag_max_203atoms energy: {energy:.6f} eV")

forces = atoms.get_forces()
print("Forces (eV/Å):")
for i, (symbol, force) in enumerate(zip(atoms.get_chemical_symbols(), forces)):
    print(
        f"  {i}: {symbol} " f"[{force[0]:8.5f}, {force[1]:8.5f}, {force[2]:8.5f}]"
    )

# Access gap directly from results
gap = calc.results.get("gap_eV")
print(f"\nHOMO-LUMO gap: {gap} eV")

# Access all electronic properties at once
elec_props = calc.get_electronic_properties()
print("\nElectronic properties:")
for key, value in sorted(elec_props.items()):
    print(f"  {key}: {value}")