"""
Author:
    - Danial Gharaie Amirabadi

Usage
-----
python detect_duplicates.py --excel my_molecules.xlsx
"""

import argparse
from collections import defaultdict

import pandas as pd
from rdkit import Chem


def canonicalize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
    except:
        pass
    return None


def main(args):
    # 1) Read Excel ------------------------------------------------------
    df_src = pd.read_excel(args.excel)
    df_src = df_src.rename(columns={c: c.strip() for c in df_src.columns})
    id_col = next(c for c in df_src.columns if c.lower().startswith("molecule"))
    smi_col = next(c for c in df_src.columns if c.lower().startswith("smile"))
    df_src = df_src[[id_col, smi_col]]
    df_src.columns = ["Molecule ChEMBL ID", "Smiles"]

    # 2) Canonicalize SMILES and detect duplicates -----------------------
    canonical_map = {}
    duplicates = defaultdict(list)

    for i, row in df_src.iterrows():
        smi = row["Smiles"]
        mol_id = row["Molecule ChEMBL ID"]
        canon = canonicalize_smiles(smi)

        if canon is None:
            continue

        if canon in canonical_map:
            duplicates[canon].append((mol_id, smi))
        else:
            canonical_map[canon] = (mol_id, smi)

    # 3) Report duplicates ------------------------------------------------
    if duplicates:
        print(f"\n🔍 Found {len(duplicates)} duplicate canonical SMILES entries.")
        duplicate_data = []
        for canon, dups in duplicates.items():
            all_ids = [canonical_map[canon][0]] + [x[0] for x in dups]
            duplicate_data.append({
                "Canonical SMILES": canon,
                "Molecule IDs": ", ".join(all_ids)
            })
        
        df_duplicates = pd.DataFrame(duplicate_data)
        output_csv_path = "duplicates.csv" # Or make this an argument
        df_duplicates.to_csv(output_csv_path, index=False)
        print(f"🧬 Duplicate entries saved to {output_csv_path}")
    else:
        print("\n✅ No duplicate SMILES found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect duplicate canonical SMILES in an Excel file."
    )
    parser.add_argument("--excel", required=True, help="Path to the source XLS/XLSX file.")
    args = parser.parse_args()
    main(args)
