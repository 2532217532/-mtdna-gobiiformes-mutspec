#!/usr/bin/env python3
"""
Extract species-to-family taxonomy mapping from GenBank GBFF files.

Input:  data/raw/gbff/batch_*.gbff (7 batch files)
Output: data/processed/taxonomy.csv (species, accession, family, lineage)

Handles sp. species via genus-level inference from GBFF records.
"""

import re
import pandas as pd
from collections import Counter, defaultdict
from pathlib import Path


def parse_gbff_files(gbff_dir):
    """Parse all GBFF files, return dict: species_name -> {accession, family, lineage}"""
    species_info = {}
    genus_families = defaultdict(list)

    for gbff_file in sorted(gbff_dir.glob("*.gbff")):
        content = gbff_file.read_text()
        for rec in content.split("//\n"):
            if not rec.strip():
                continue

            ver_match = re.search(r"VERSION\s+(\S+)", rec)
            org_match = re.search(r"ORGANISM\s+(.+?)\n", rec)
            lineage_match = re.search(
                r"ORGANISM\s+.+?\n(.+?)(?=\n[A-Z])", rec, re.DOTALL
            )
            if not (ver_match and org_match and lineage_match):
                continue

            acc = ver_match.group(1)
            org_name = org_match.group(1).strip()
            sp = org_name.replace(" ", "_")
            genus = org_name.split()[0]
            lineage = " ".join(lineage_match.group(1).strip().split())

            # Extract family: last rank ending in -idae
            ranks = [r.strip(" ;.") for r in lineage.split(";")]
            family = None
            for r in reversed(ranks):
                if r.lower().endswith("idae"):
                    family = r
                    break
            if family is None:
                fam_match = re.search(r"(\w+idae)", lineage)
                if fam_match:
                    family = fam_match.group(1)

            # Track genus->family mapping for sp. species inference
            if family:
                genus_families[genus].append(family)

            # Keep first record per species
            if sp not in species_info:
                species_info[sp] = {
                    "accession": acc,
                    "family": family,
                    "lineage": lineage,
                }

    return species_info, dict(genus_families)


def resolve_sp_species(species_name, genus_families):
    """Resolve family for sp. (undescribed) species via genus-level majority vote."""
    genus = species_name.split("_")[0]
    if genus == "Gobiidae":
        return "Gobiidae"
    families = genus_families.get(genus, [])
    if not families:
        return None
    return Counter(families).most_common(1)[0][0]


def main():
    project_dir = Path("/home/zengjl/WorkSpace/mtdna-gobiiformes-mutspec")
    gbff_dir = project_dir / "data" / "raw" / "gbff"
    spectra_path = project_dir / "data" / "processed" / "species_spectra_192syn.csv"
    output_path = project_dir / "data" / "processed" / "taxonomy.csv"

    # 1. Parse GBFF
    gbff_species, genus_families = parse_gbff_files(gbff_dir)
    print(f"GBFF unique species: {len(gbff_species)}")

    # 2. Load spectra species list
    spectra_df = pd.read_csv(spectra_path)
    spectra_species = spectra_df["species"].tolist()
    print(f"Spectra species:    {len(spectra_species)}")

    # 3. Match and resolve
    records = []
    direct_match = 0
    inferred = 0
    for sp in spectra_species:
        if sp in gbff_species:
            info = gbff_species[sp]
            direct_match += 1
        else:
            family = resolve_sp_species(sp, genus_families)
            info = {
                "accession": "",
                "family": family,
                "lineage": f"inferred_from_genus:{sp.split('_')[0]}",
            }
            inferred += 1

        records.append(
            {
                "species": sp,
                "accession": info["accession"],
                "family": info["family"],
                "lineage": info["lineage"],
            }
        )

    # 4. Report
    family_counts = Counter(r["family"] for r in records)
    print(f"\nDirect matches:  {direct_match}")
    print(f"Genus-inferred:  {inferred}")
    print(f"Unmatched:       {len(spectra_species) - direct_match - inferred}")
    print(f"\nFamily distribution:")
    for fam, cnt in family_counts.most_common():
        print(f"  {fam}: {cnt}")

    # 5. Save
    taxonomy_df = pd.DataFrame(records)
    taxonomy_df.to_csv(output_path, index=False)
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()
