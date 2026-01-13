#!/usr/bin/env python3
"""
Split a PDB at the SECOND-LAST TER, then fix ligand bond orders from a SMILES using RDKit.

Behavior:
- protein1.pdb := all lines up to and including the second-last 'TER'
- ligand1.pdb  := ATOM/HETATM lines after that 'TER'
- Then assign bond orders in ligand1.pdb from a SMILES template and write:
    - ligand_with_bo.sdf
    - ligand_with_bo.pdb

Usage:
  python split_at_second_last_ter_and_fix.py --input sample1.pdb --smiles "CC(=O)Oc1ccccc1C(=O)O"
"""

import argparse
import os
import sys

def parse_args():
    p = argparse.ArgumentParser(description="Split PDB at the second-last TER and assign ligand bond orders from SMILES.")
    p.add_argument("--input", "-i", default="sample1.pdb", help="Input PDB file (default: sample1.pdb)")
    p.add_argument("--smiles", "-s", required=True, help="SMILES used as template for bond orders")
    p.add_argument("--protein-out", default="protein1.pdb", help="Output protein PDB (default: protein1.pdb)")
    p.add_argument("--ligand-out", default="ligand1.pdb", help="Output ligand PDB (default: ligand1.pdb)")
    p.add_argument("--sdf-out", default="ligand_with_bo.sdf", help="Output SDF with bond orders (default: ligand_with_bo.sdf)")
    p.add_argument("--pdb-out", default="ligand_with_bo.pdb", help="Output PDB with bond orders (default: ligand_with_bo.pdb)")
    return p.parse_args()

def load_pdb_lines(path):
    if not os.path.isfile(path):
        print(f"ERROR: File not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path, "r") as f:
        return f.readlines()

def find_ter_indices(lines):
    """Return a list of indices where a line starts with 'TER'."""
    return [i for i, ln in enumerate(lines) if ln.startswith("TER")]

def write_lines(path, lines):
    with open(path, "w") as f:
        f.writelines(lines)

def extract_ligand_after_index(lines, idx):
    """Return only ATOM/HETATM lines AFTER the given index."""
    out = []
    for ln in lines[idx+1:]:
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            out.append(ln)
    return out

def main():
    args = parse_args()
    lines = load_pdb_lines(args.input)
    ter_idx = find_ter_indices(lines)

    if len(ter_idx) == 0:
        print("WARNING: No 'TER' records found. Entire file will be treated as protein; ligand will be empty.", file=sys.stderr)
        split_at = -1
    elif len(ter_idx) == 1:
        print("WARNING: Only one 'TER' found. Using that as the split point (cannot use second-last).", file=sys.stderr)
        split_at = ter_idx[0]
    else:
        # Second-last TER
        split_at = ter_idx[-2]

    # Protein: everything up to and including the chosen TER (or entire file if split_at == -1)
    if split_at >= 0:
        protein_lines = lines[: split_at + 1]
        ligand_lines = extract_ligand_after_index(lines, split_at)
    else:
        protein_lines = lines[:]
        ligand_lines = []

    # Write split files
    write_lines(args.protein_out, protein_lines)
    write_lines(args.ligand_out, ligand_lines)

    print(f"Wrote {args.protein_out} ({len(protein_lines)} lines)")
    print(f"Wrote {args.ligand_out} ({len(ligand_lines)} ATOM/HETATM lines)")

    # ---- RDKit step: assign bond orders from SMILES template ----
    if not ligand_lines:
        print("No ligand atoms to process with RDKit. Exiting after split.", file=sys.stderr)
        return

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("ERROR: RDKit is required for bond-order assignment.\n"
              "Try:  conda install -c conda-forge rdkit", file=sys.stderr)
        sys.exit(1)

    template = Chem.MolFromSmiles(args.smiles)
    if template is None:
        print("ERROR: Failed to parse the provided SMILES.", file=sys.stderr)
        sys.exit(1)

    with open(args.ligand_out, "r") as f:
        pdb_block = f.read()

    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
    if mol is None:
        print(f"ERROR: RDKit could not parse {args.ligand_out}. Ensure it contains only the intended ligand atoms.", file=sys.stderr)
        sys.exit(1)

    try:
        mol_bo = AllChem.AssignBondOrdersFromTemplate(template, mol)
    except Exception as e:
        print(f"ERROR during AssignBondOrdersFromTemplate: {e}", file=sys.stderr)
        sys.exit(1)

    # Write outputs with bond orders
    w = Chem.SDWriter(args.sdf_out)
    w.write(mol_bo)
    w.close()
    Chem.MolToPDBFile(mol_bo, args.pdb_out)

    print(f"Wrote {args.sdf_out} and {args.pdb_out} with bond orders from the SMILES template.")

if __name__ == "__main__":
    main()

