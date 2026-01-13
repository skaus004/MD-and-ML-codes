import pandas as pd
import yaml
from pathlib import Path
from collections import OrderedDict

# === Fix PyYAML dumping of OrderedDict ===
def represent_ordereddict(dumper, data):
    return dumper.represent_dict(data.items())

yaml.add_representer(OrderedDict, represent_ordereddict)
yaml.Dumper.ignore_aliases = lambda *args: True

# === Configuration ===
input_csv_file = "WO2025_Kymera_571cmpds.csv"
output_dir = "stat6_20250807"
smiles_column = "Compound Structure"
id_column = "ID"
selected_ids = None  # Optional: restrict to a subset

# === Protein Sequences ===
protein_a_seq = "QFRHLPMPFHWKQEELKFKTGLRRLQHRVGEIHLLREALQKGAEAGQVSLHSLIETPANGTGPSEALAMLLQETTGELEAAKALVLKRIQIWKRQQQLAGNGAPFEESLAPLQERCESLVDIYSQLQQEVGAAGGELEPKTRASLTGRLDEVLRTLVTSCFLVEKQPPQVLKTQTKFQAGVRFLLGLRFLGAPAKPPLVRADMVTEKQARELSVPQGPGAGAESTGEIINNTVPLENSIPGNCCSALFKNLLLKKIKRCERKGTESVTEEKCAVLFSASFTLGPGKLPIQLQALSLPLVVIVHGNQDNNAKATILWDNAFSEMDRVPFVVAERVPWEKMCETLNLKFMAEVGTNRGLLPEHFLFLAQKIFNDNSLSMEAFQHRSVSWSQFNKEILLGRGFTFWQWFDGVLDLTKRCLRSYWSDRLIIGFISKQYVTSLLLNEPDGTFLLRFSDSEIGGITIAHVIRGQDGSPQIENIQPFSAKDLSIRSLGDRIRDLAQLKNLYPKKPKDEAF"
protein_b_seq = "GPSTSLSCKQCQETEITTKNEIFSLSLSGPMAAYVNPHGYVHETLTVYKASNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTI"

# === Output Directory ===
Path(output_dir).mkdir(parents=True, exist_ok=True)

# === Load Input ===
df = pd.read_csv(input_csv_file, index_col=id_column)
if selected_ids:
    df = df[df.index.isin(selected_ids)]

# === Generate YAML Files ===
for compound_id, row in df.iterrows():
    smiles = row[smiles_column]
    if pd.isna(smiles):
        print(f"Skipping {compound_id}: missing SMILES.")
        continue

    # Compose Ordered YAML structure
    data = OrderedDict()
    data['sequences'] = [
        {'protein': {'id': 'A', 'sequence': protein_a_seq}},
        {'protein': {'id': 'B', 'sequence': protein_b_seq}},
        {'ligand': {'id': 'X', 'smiles': str(smiles)}}
    ]
    data['constraints'] = [
        {'pocket': OrderedDict({
            'binder': 'X',
            'contacts': '[[A, 258 ]]',  # Placeholder string
            'max_distance': 5,
            'force': True
        })}
    ]
    data['properties'] = [{'affinity': {'binder': 'X'}}]

    # Dump YAML and fix contact line manually
    yaml_str = yaml.dump(data, default_flow_style=False, sort_keys=False)
    yaml_str = yaml_str.replace("contacts: '[[A, 258 ]]'", "contacts: [[A, 258 ]]")

    # Write output
    with open(Path(output_dir) / f"{compound_id}.yaml", 'w') as f:
        f.write(yaml_str)

    print(f"Generated: {compound_id}.yaml")
