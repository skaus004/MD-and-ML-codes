#usage --> sh run.sh "COO"

# First argument: ligand SMILES (default to "COO" if not provided)
LIGAND_SMILES="${1:-COO}"

# Output JSON file (you can also make this an argument if you want)
OUTPUT_JSON="out.json"

# Name field in the JSON
NAME="protenix_prediction_OCT20_sample"

# Call the Python script
python 01_protenix.py "$LIGAND_SMILES" "$OUTPUT_JSON" --name "$NAME"

protenix predict --input out.json --out_dir ./output --use_msa true

cp output/protenix_prediction_OCT20_sample/seed_101/predictions/protenix_prediction_OCT20_sample_sample_0.cif output/

python cif2pdb.py output/protenix_prediction_OCT20_sample_sample_0.cif output/sample_0.pdb

python split_complex.py -i output/sample_0.pdb -s "$LIGAND_SMILES"

cat protein1.pdb ligand_with_bo.pdb > complex1.pdb
sed -i 's/l01/L01/g' complex1.pdb

gnina -r protein1.pdb -l ligand_with_bo.sdf --minimize -o minimized.sdf.gz
