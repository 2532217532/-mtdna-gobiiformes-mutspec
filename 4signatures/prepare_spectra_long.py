#!/usr/bin/env python3
"""
Prepare Gobiiformes spectra data in long format for signatures analysis.

Converts wide-format species/gene spectra to long format with family taxonomy,
matching the format expected by the SigProfilerAssignment analysis notebook.

Input:
  - data/processed/species_spectra_192syn.csv  (wide, 295 species x 192 SBS)
  - data/processed/taxonomy.csv                (species -> family mapping)

Output:
  - 4signatures/data/gobiiformes_spectra_long.csv  (long format: family, Species, Mut, MutSpec)
"""

import pandas as pd
from pathlib import Path

PROJECT = Path("/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec")


def main():
    # Load wide-format species spectra
    wide = pd.read_csv(PROJECT / "data" / "processed" / "species_spectra_192syn.csv")
    print(f"Loaded species spectra: {wide.shape} (rows x cols)")

    # Load taxonomy
    tax = pd.read_csv(PROJECT / "data" / "processed" / "taxonomy.csv")
    print(f"Loaded taxonomy: {tax.shape}")

    # Merge spectra with taxonomy
    wide = wide.merge(tax[["species", "family"]], on="species", how="left")
    print(f"Merged with taxonomy: {wide.shape}")

    # Check for missing families
    missing = wide[wide["family"].isna()]["species"].tolist()
    if missing:
        print(f"WARNING: Missing family for {len(missing)} species: {missing}")
        raise ValueError("Missing family assignments")

    # Melt to long format
    sbs_cols = [c for c in wide.columns if c not in ["species", "family"]]
    long_df = wide.melt(
        id_vars=["species", "family"],
        value_vars=sbs_cols,
        var_name="Mut",
        value_name="MutSpec",
    )

    # Rename to match notebook expectations: Species (capitalized)
    long_df = long_df.rename(columns={"species": "Species"})
    # Capitalize family column for consistency with original notebook
    long_df["family"] = long_df["family"].astype(str)
    long_df["Mut"] = long_df["Mut"].astype(str)

    print(f"Long format: {len(long_df)} rows")
    print(f"Species: {long_df['Species'].nunique()}")
    print(f"Families: {sorted(long_df['family'].unique())}")

    # Show family counts
    fam_counts = long_df.groupby("family")["Species"].nunique().sort_values(ascending=False)
    print(f"\nSpecies per family:")
    for fam, cnt in fam_counts.items():
        print(f"  {fam}: {cnt}")

    # Check H-strand notation: first few unique Mut types
    mut_types = long_df["Mut"].unique()
    print(f"\nUnique Mut types: {len(mut_types)}")
    print(f"Sample Mut types: {sorted(mut_types)[:5]} ... {sorted(mut_types)[-5:]}")

    # Confirm they match H-strand notation (A[N>N]N format)
    import re
    h_strand_pattern = re.compile(r"^[ACGT]\[[ACGT]>[ACGT]\][ACGT]$")
    non_h = [m for m in mut_types if not h_strand_pattern.match(m)]
    if non_h:
        print(f"WARNING: {len(non_h)} mut types don't match H-strand notation!")
        print(f"  Examples: {non_h[:10]}")
    else:
        print("All mutation types confirmed in H-strand notation (A[N>N]N)")

    # Save
    outdir = PROJECT / "4signatures" / "data"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "gobiiformes_spectra_long.csv"
    long_df.to_csv(outpath, index=False)
    print(f"\nSaved to: {outpath}")


if __name__ == "__main__":
    main()
