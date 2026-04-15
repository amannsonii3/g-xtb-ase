"""
Example: accessing electronic properties (HOMO-LUMO gap, IP, EA, spin) via atoms.info.

After a g-xTB calculation, non-standard electronic properties are attached to
atoms.info under gxtb_* keys.  Standard ASE properties (energy, forces, dipole,
charges) remain in calc.results as usual.

Run from the repo root:
    python examples/electronic_properties.py
"""

from pathlib import Path
from ase.io import read
from gxtb_ase import GxTB

SCRIPT_DIR = Path(__file__).resolve().parent

atoms = read(SCRIPT_DIR / "structures" / "Ag_min_35atoms.xyz")

# write_log=True ensures the log file is parsed so that electronic
# properties (gap, IP, EA, ...) are populated into atoms.info.
calc = GxTB(
    label="Ag35",
    charge=0,
    spin=0,
    write_log=True,
    directory=str(SCRIPT_DIR),
)
atoms.calc = calc

energy = atoms.get_potential_energy()
print(f"Total energy : {energy:.6f} eV")

# Standard ASE properties live in calc.results
forces = atoms.get_forces()
print(f"Max force    : {abs(forces).max():.6f} eV/Å")

# Non-standard electronic properties are stored in atoms.info
print("\nElectronic properties (atoms.info):")
for key in sorted(k for k in atoms.info if k.startswith("gxtb_")):
    print(f"  {key}: {atoms.info[key]}")

# Access individual properties directly
gap = atoms.info.get("gxtb_gap_eV")
ip  = atoms.info.get("gxtb_IP_Janak_eV")
ea  = atoms.info.get("gxtb_EA_Janak_eV")
print(f"\nHOMO-LUMO gap : {gap} eV")
if ip is not None:
    print(f"IP (Janak)    : {ip} eV")
if ea is not None:
    print(f"EA (Janak)    : {ea} eV")
