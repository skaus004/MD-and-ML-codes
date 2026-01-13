#!/usr/bin/env python3
"""
Generate a single JSON file with your protein template and one ligand SMILES.

Usage examples:
  python make_ligand_json.py "CCO" out.json
  python make_ligand_json.py "CCO" out.json --name protenix_prediction_OCT20_sample
"""

import argparse
import json
from pathlib import Path

PROTEIN_SEQ_A = (
    "AMPLDAGGQNSTQMVLAPGASIFRCRQCGQTISRRDWLLPMGGDHEHVVFNPAGMIFRVWCFSLAQGLRLIGAPSGEFSWFKGYDWTIALCGQCGSHLGWHYEGGSQPQTFFGLIKDRLAEGPAD"
)
PROTEIN_SEQ_B = (
    "QFRHLPMPFHWKQEELKFKTGLRRLQHRVGEIHLLREALQKGAEAGQVSLHSLIETPANGTGPSEALAMLLQETTGELEAAKALVLKRIQIWKRQQQLAGNGAPFEESLAPLQERCESLVDIYSQLQQEVGAAGGELEPKTRASLTGRLDEVLRTLVTSCFLVEKQPPQVLKTQTKFQAGVRFLLGLRFLGAPAKPPLVRADMVTEKQARELSVPQGPGAGAESTGEIINNTVPLENSIPGNCCSALFKNLLLKKIKRCERKGTESVTEEKCAVLFSASFTLGPGKLPIQLQALSLPLVVIVHGNQDNNAKATILWDNAFSEMDRVPFVVAERVPWEKMCETLNLKFMAEVGTNRGLLPEHFLFLAQKIFNDNSLSMEAFQHRSVSWSQFNKEILLGRGFTFWQWFDGVLDLTKRCLRSYWSDRLIIGFISKQYVTSLLLNEPDGTFLLRFSDSEIGGITIAHVIRGQDGSPQIENIQPFSAKDLSIRSLGDRIRDLAQLKNLYPKKPKDEAFRSHYKPEQMGKDGRGYVPATIKMTVERDQPLPT"
)

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a JSON with protein template and a single ligand SMILES."
    )
    p.add_argument(
        "ligand_smiles",
        help="Ligand SMILES string (quote it if it contains special characters)",
    )
    p.add_argument(
        "output_path",
        help="Path to the output .json file",
    )
    p.add_argument(
        "--name",
        default="ligand",
        help="Value for the 'name' field in the JSON (default: %(default)s)",
    )
    return p.parse_args()

def make_entry(ligand_smiles: str, name: str) -> dict:
    return {
        "sequences": [
            {"proteinChain": {"sequence": PROTEIN_SEQ_A, "count": 1}},
            {"proteinChain": {"sequence": PROTEIN_SEQ_B, "count": 1}},
            {"ligand": {"ligand": ligand_smiles, "count": 1}},
        ],
        "name": name,
    }

def main():
    args = parse_args()
    outpath = Path(args.output_path)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    entry = make_entry(args.ligand_smiles, args.name)

    # Keep as [entry] if downstream expects a list; if not, change to just entry
    with open(outpath, "w", encoding="utf-8") as fout:
        json.dump([entry], fout, indent=4, ensure_ascii=False)

    print(f"Wrote {outpath}")

if __name__ == "__main__":
    main()

