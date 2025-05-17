from neurosnap.conformers import generate

"""
Author:
    - Danial Gharaie Amirabadi
Usage
-----
python conf_gen.py --excel my_molecules.xlsx --num-confs 1500
"""
import os
import argparse
from pathlib import Path
import pandas as pd
from tqdm import tqdm


def _process_row(row, num_confs, min_method, write_multi):
    """
    Worker: run in a subprocess.  Returns (mol_id, success_flag, err_msg)
    """
    mol_id = str(row["Molecule ChEMBL ID"]).strip()
    smiles = str(row["Smiles"]).strip()

    out_dir = Path(mol_id)
    try:
        out_dir.mkdir(exist_ok=True)

        # Tell generate() to drop SDFs directly into the ID folder
        df = generate(
            input_mol=smiles,
            output_name=str(out_dir),  # <- generate writes here
            write_multi=write_multi,
            num_confs=num_confs,
            min_method=min_method,
        )

        # Save conformer stats
        df.to_csv(out_dir / f"{mol_id}_conformers.csv", index=False)
        return mol_id, True, ""
    except Exception as e:
        # Propagate the error message up so the parent can log it
        return mol_id, False, repr(e)


def main(args):
    # 1) Read Excel ------------------------------------------------------
    df_src = pd.read_excel(args.excel)

    # column-name normalization (robust to different capitalisation)
    df_src = df_src.rename(columns={c: c.strip() for c in df_src.columns})
    id_col = next(c for c in df_src.columns if c.lower().startswith("molecule"))
    smi_col = next(c for c in df_src.columns if c.lower().startswith("smile"))
    df_src = df_src[[id_col, smi_col]]
    df_src.columns = ["Molecule ChEMBL ID", "Smiles"]

    # 2) Sequential processing with progress bar -------------------------
    work = [
        (row, args.num_confs, args.min_method, args.write_multi)
        for _, row in df_src.iterrows()
    ]

    successes, failures = [], []
    for args_tuple in tqdm(work, total=len(work), desc="Generating conformers"):
        mol_id, ok, err = _process_row(*args_tuple)
        if ok:
            successes.append(mol_id)
        else:
            failures.append((mol_id, err))

    # 3) Summary ---------------------------------------------------------
    print("\n✅ Finished!")
    print(f"   Successful molecules: {len(successes)} / {len(df_src)}")
    if failures:
        print("   ❌ Failures:")
        for mol_id, err in failures:
            print(f"     - {mol_id}: {err}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate conformers for many molecules in parallel."
    )
    parser.add_argument(
        "--excel", required=True, help="Path to the source XLS/XLSX file."
    )
    parser.add_argument(
        "--num-confs",
        type=int,
        default=250,
        help="How many conformers to generate per molecule.",
    )
    parser.add_argument(
        "--min-method",
        default="auto",
        choices=["auto", "UFF", "MMFF94", "MMFF94s", "none"],
        help='"auto" tries MMFF94→UFF→MMFF94s. Use "none" to skip minimisation.',
    )
    parser.add_argument(
        "--write-multi",
        action="store_true",
        help="Write all unique conformers to a single SDF instead of separate ones.",
    )
    args = parser.parse_args()

    # quick convenience: allow --min-method none
    if args.min_method.lower() == "none":
        args.min_method = None

    main(args)
