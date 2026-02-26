"""
Parse g-xTB log files and store electronic properties into ASE Atoms.info.

Usage:
    python parse_gxtb.py                      # process all .log files in cwd
    python parse_gxtb.py path/to/file.log     # process a single log file
"""

import re
import sys
import glob
from pathlib import Path
from ase.io import read, write


def parse_gap_from_text(text: str) -> dict:
    """
    Extract gap values from g-xTB output text.

    Returns dict with 'is_uks' and 'gap_eV' (float or dict).
    """
    result = {}
    is_uks = bool(re.search(r'UKS \?\s+T', text))
    result['is_uks'] = is_uks

    if is_uks:
        gaps = {}
        # These appear in the summary section after population analysis
        m = re.search(r'gap \(eV\)\s+alpha\s*->\s*alpha\s*:\s+([\d.]+)', text)
        if m:
            gaps['alpha_alpha'] = float(m.group(1))
        m = re.search(r'beta\s*->\s*beta\s*:\s+([\d.]+)', text)
        if m:
            gaps['beta_beta'] = float(m.group(1))
        m = re.search(r'alpha\s*->\s*beta\s*:\s+([\d.]+)', text)
        if m:
            gaps['alpha_beta'] = float(m.group(1))
        result['gap_eV'] = gaps if gaps else None
    else:
        # RKS: "gap (eV)                :        3.28876431" in the summary
        # Find the LAST occurrence (summary, not SCF table)
        matches = re.findall(r'gap \(eV\)\s+:\s+([\d.]+)', text)
        if matches:
            result['gap_eV'] = float(matches[-1])
        else:
            # Fallback: last gap value from SCF iteration table
            scf_gaps = re.findall(r'[\d.]+\s+\d+\s+F\s*$', text, re.MULTILINE)
            if scf_gaps:
                parts = scf_gaps[-1].split()
                if len(parts) >= 3:
                    try:
                        result['gap_eV'] = float(parts[0])
                    except ValueError:
                        result['gap_eV'] = None
    return result


def parse_gxtb_log(logfile: str) -> dict:
    """
    Parse a g-xTB log file and return a dict of properties.

    Returns
    -------
    dict with keys:
        - 'xyz_file': str, the input coordinate file name
        - 'total_energy_Eh': float, total energy in Hartree
        - 'is_uks': bool, whether UKS (open-shell) calculation
        - 'gap_eV': float or dict or None
            RKS: single float
            UKS: dict with 'alpha_alpha', 'beta_beta', 'alpha_beta'
            None if not found
        - 'dipole_au': list[float], final SCF dipole [x, y, z]
        - 'dipole_total_au': float
        - 'dipole_total_Debye': float
        - 'nopen': int, number of unpaired electrons
        - 'charge': int
        - 'nel': int, number of electrons
        - 'IP_Janak_eV': float (RKS only)
        - 'EA_Janak_eV': float (RKS only)
        - 'spin_contamination': float (UKS only)
    """
    text = Path(logfile).read_text()
    result = {}

    # --- Input XYZ file ---
    m = re.search(r'coord file:\s*(\S+)', text)
    if m:
        result['xyz_file'] = m.group(1)

    # --- Charge, nel, nopen ---
    m = re.search(r'^\s*charge\s+(-?\d+)', text, re.MULTILINE)
    if m:
        result['charge'] = int(m.group(1))
    m = re.search(r'^\s*nel\s+(\d+)', text, re.MULTILINE)
    if m:
        result['nel'] = int(m.group(1))
    m = re.search(r'^\s*nopen\s+(\d+)', text, re.MULTILINE)
    if m:
        result['nopen'] = int(m.group(1))

    # --- Gap values (using shared helper) ---
    gap_info = parse_gap_from_text(text)
    result['is_uks'] = gap_info['is_uks']
    result['gap_eV'] = gap_info.get('gap_eV')

    # --- Total energy (Hartree) ---
    m = re.search(r'^\s*total\s+(-[\d.]+)\s*$', text, re.MULTILINE)
    if m:
        result['total_energy_Eh'] = float(m.group(1))

    # --- Final dipole (last occurrence, i.e. post-SCF) ---
    dipole_matches = re.findall(
        r'dipole moment\s+X\s+Y\s+Z\s*\n'
        r'\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+total\s+\(au/Debye\)\s*:\s+([-\d.]+)\s+([-\d.]+)',
        text
    )
    if dipole_matches:
        last = dipole_matches[-1]
        result['dipole_au'] = [float(last[0]), float(last[1]), float(last[2])]
        result['dipole_total_au'] = float(last[3])
        result['dipole_total_Debye'] = float(last[4])

    # --- Janak's theorem IP/EA (RKS only) ---
    m = re.search(r"Janaks theorem for IP\s*:\s+([\d.]+)", text)
    if m:
        result['IP_Janak_eV'] = float(m.group(1))
    m = re.search(r"Janaks theorem for EA\s*:\s+([\d.]+)", text)
    if m:
        result['EA_Janak_eV'] = float(m.group(1))

    # --- Spin contamination (UKS only) ---
    m = re.search(r'spin-contamination:\s+([\d.]+)', text)
    if m:
        result['spin_contamination'] = float(m.group(1))

    return result


def attach_to_xyz(logfile: str, outdir: str = None):
    """Parse log, read xyz, attach info, write new xyz."""
    props = parse_gxtb_log(logfile)
    logpath = Path(logfile)

    # Find the xyz file (look in same directory as log)
    xyz_name = props.get('xyz_file')
    if xyz_name is None:
        print(f"WARNING: Could not find xyz filename in {logfile}, skipping.")
        return

    xyz_path = logpath.parent / xyz_name
    if not xyz_path.exists():
        print(f"WARNING: {xyz_path} not found, skipping.")
        return

    atoms = read(str(xyz_path))

    # Attach all parsed properties to atoms.info
    # Nested dicts (like UKS gaps) are fine in atoms.info
    atoms.info['gxtb_total_energy_Eh'] = props.get('total_energy_Eh')
    atoms.info['gxtb_gap_eV'] = props.get('gap_eV')  # float or dict
    atoms.info['gxtb_is_uks'] = props.get('is_uks', False)
    atoms.info['gxtb_charge'] = props.get('charge', 0)
    atoms.info['gxtb_nel'] = props.get('nel')
    atoms.info['gxtb_nopen'] = props.get('nopen', 0)

    if 'dipole_au' in props:
        atoms.info['gxtb_dipole_au'] = props['dipole_au']
        atoms.info['gxtb_dipole_Debye'] = props['dipole_total_Debye']

    if 'IP_Janak_eV' in props:
        atoms.info['gxtb_IP_Janak_eV'] = props['IP_Janak_eV']
    if 'EA_Janak_eV' in props:
        atoms.info['gxtb_EA_Janak_eV'] = props['EA_Janak_eV']
    if 'spin_contamination' in props:
        atoms.info['gxtb_spin_contamination'] = props['spin_contamination']

    # Write output
    if outdir:
        out_path = Path(outdir) / xyz_path.name
    else:
        out_path = xyz_path.parent / f"{xyz_path.stem}_gxtb.xyz"

    write(str(out_path), atoms, format='extxyz')
    print(f"Wrote {out_path}")
    print(f"  Energy: {props.get('total_energy_Eh')} Eh")
    print(f"  Gap:    {props.get('gap_eV')}")
    print(f"  UKS:    {props.get('is_uks')}")
    return atoms


if __name__ == '__main__':
    if len(sys.argv) > 1:
        logs = sys.argv[1:]
    else:
        logs = sorted(glob.glob('*.log'))

    if not logs:
        print("No .log files found.")
        sys.exit(1)

    for log in logs:
        print(f"\n--- Processing {log} ---")
        attach_to_xyz(log)
