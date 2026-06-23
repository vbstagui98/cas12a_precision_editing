#!/usr/bin/env python3
"""Export path-independent CSV inputs for the public analysis repository."""

from pathlib import Path
import hashlib
import os

import numpy as np
import pandas as pd


PROJECT = Path(os.environ.get("CAS12A_REPO_ROOT", Path.cwd())).resolve()
ANNOTATIONS = Path(
    os.environ.get(
        "CAS12A_PUBLICATION_ANNOTATIONS",
        PROJECT / "outputs/publication_annotations_20260623",
    )
)
PAPER = os.environ.get("CAS12A_PAPER_NOTEBOOKS_ROOT")
ENAS_PRE_CORRECTION = os.environ.get("CAS12A_ENAS_PRE_CORRECTION_VARIANTS")
OUTPUT = PROJECT / "data/publication_inputs"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def join_unique(values: pd.Series) -> str:
    return "; ".join(sorted(set(values.dropna().astype(str))))


def export_endogenous_inputs() -> list[dict]:
    if PAPER is None:
        raise RuntimeError(
            "Set CAS12A_PAPER_NOTEBOOKS_ROOT to the finalized paper-notebook folder."
        )
    if ENAS_PRE_CORRECTION is None:
        raise RuntimeError(
            "Set CAS12A_ENAS_PRE_CORRECTION_VARIANTS to "
            "editing_window_all_variants_20260122.csv."
        )

    paper_root = Path(PAPER)
    long_source = pd.read_csv(
        paper_root
        / "processed_data/genomic_ampli_enas_fn_harmonized_20260318.csv",
        low_memory=False,
    )
    guide_source = pd.read_csv(
        ANNOTATIONS / "endogenous_panel_guide_annotations.csv",
        low_memory=False,
    )
    enas_pre_correction = pd.read_csv(
        Path(ENAS_PRE_CORRECTION),
        low_memory=False,
    )
    parental_g06 = enas_pre_correction[
        (enas_pre_correction["promoter"] == "RPR1")
        & (enas_pre_correction["Guide"] == "G06")
        & (enas_pre_correction["MUTATION"] == "mnp")
        & (enas_pre_correction["mismatches"] == "TG")
        & (enas_pre_correction["Experiment"] != "E5")
        & ~enas_pre_correction["Replicate"].isin(["C1", "C2"])
    ].copy()

    hdr_rows = long_source[
        (long_source["MUTATION"] == "HDR")
        & long_source["CHROM"].notna()
        & long_source["pos_mismatch"].notna()
        & (long_source["CHROM"].astype(str) != "0")
    ].copy()
    parental_for_assignment = parental_g06.assign(
        cas_variant="enAsCas12a",
        source_batch="editing_window_all_variants_20260122",
        source_file="editing_window_all_variants_20260122.csv",
        source_row=parental_g06["X"],
    )
    hdr_rows = pd.concat(
        [hdr_rows, parental_for_assignment],
        ignore_index=True,
        sort=False,
    )

    intended_loci = (
        hdr_rows.groupby("Guide", as_index=False)
        .agg(
            intended_chromosome=("CHROM", lambda x: x.mode().iloc[0]),
            intended_position=("pos_mismatch", lambda x: int(x.mode().iloc[0])),
        )
        .rename(columns={"Guide": "crRNA_id"})
        .merge(
            guide_source[["crRNA_id", "alternate_allele"]],
            on="crRNA_id",
            how="right",
            validate="one_to_one",
        )
        .sort_values("crRNA_id")
    )
    if intended_loci[
        ["intended_chromosome", "intended_position", "alternate_allele"]
    ].isna().any().any():
        raise ValueError("Incomplete endogenous intended-locus assignments.")

    hdr_rows["candidate_match"] = (
        hdr_rows["promoter"].astype(str)
        + "_"
        + hdr_rows["Guide"].astype(str)
        + "_"
        + hdr_rows["CHROM"].astype(str)
        + "_"
        + hdr_rows["pos_mismatch"].astype(int).astype(str)
        + "_"
        + hdr_rows["mismatches"].astype(str)
    )
    hdr_rows["candidate_variant_id"] = (
        hdr_rows["promoter"].astype(str)
        + "_"
        + hdr_rows["Guide"].astype(str)
        + "_"
        + hdr_rows["CHROM"].astype(str)
        + "_"
        + hdr_rows["POS"].astype(int).astype(str)
        + "_"
        + hdr_rows["REF"].astype(str)
        + "_"
        + hdr_rows["ALT"].astype(str)
        + "_"
        + hdr_rows["TYPE"].astype(str)
    )
    intended_aliases = (
        hdr_rows.groupby(["Guide", "promoter"], as_index=False)
        .agg(
            intended_match=(
                "candidate_match",
                lambda x: "|".join(sorted(set(x.dropna().astype(str)))),
            ),
            intended_variant_id=(
                "candidate_variant_id",
                lambda x: "|".join(sorted(set(x.dropna().astype(str)))),
            ),
        )
        .rename(
            columns={
                "Guide": "crRNA_id",
                "promoter": "crRNA_promoter",
            }
        )
    )
    intended = (
        intended_aliases.merge(
            intended_loci,
            on="crRNA_id",
            how="left",
            validate="many_to_one",
        )
        .sort_values(["crRNA_id", "crRNA_promoter"])
        .reset_index(drop=True)
    )

    guide_columns = [
        "crRNA_id",
        "target",
        "variant_position_annotation",
        "crRNA_sequence",
        "donor_sequence",
        "mutation_type",
        "alternate_allele",
        "pam",
        "pre_pam_nt",
        "distance_from_pam",
        "strand",
        "crRNA_gc_pct",
        "deepcpf1_score",
    ]
    guide_features = guide_source[guide_columns].sort_values("crRNA_id")

    raw_columns = [
        "source_batch",
        "source_file",
        "source_row",
        "cas_variant",
        "DR_cognate",
        "Experiment",
        "promoter",
        "recruit",
        "ys",
        "Guide",
        "Timepoint",
        "Replicate",
        "total_gen_in_liquid",
        "Sample",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "TYPE",
        "AO",
        "RO",
        "DP",
        "total",
        "QUAL",
    ]
    variants = long_source[long_source["MUTATION"] != "REF"][raw_columns].copy()
    sample_metadata_lookup = (
        long_source[
            [
                "Sample",
                "cas_variant",
                "DR_cognate",
                "ys",
                "total_gen_in_liquid",
            ]
        ]
        .drop_duplicates("Sample")
        .set_index("Sample")
    )
    parental_export = parental_g06.copy()
    parental_export["source_batch"] = "editing_window_all_variants_20260122"
    parental_export["source_file"] = "editing_window_all_variants_20260122.csv"
    parental_export["source_row"] = parental_export["X"]
    parental_export["cas_variant"] = "enAsCas12a"
    parental_export["DR_cognate"] = parental_export["Sample"].map(
        sample_metadata_lookup["DR_cognate"]
    )
    parental_export["ys"] = parental_export["Sample"].map(
        sample_metadata_lookup["ys"]
    )
    parental_export["total_gen_in_liquid"] = parental_export["Sample"].map(
        sample_metadata_lookup["total_gen_in_liquid"]
    )
    parental_export = parental_export[raw_columns]
    variants = pd.concat([variants, parental_export], ignore_index=True)
    variants = variants.rename(
        columns={
            "DR_cognate": "direct_repeat",
            "Experiment": "experiment",
            "promoter": "guide_promoter",
            "recruit": "donor_recruitment",
            "ys": "strain",
            "Guide": "crRNA_id",
            "Timepoint": "timepoint",
            "Replicate": "replicate",
            "total_gen_in_liquid": "generations",
            "Sample": "sample",
            "total": "total_count",
        }
    )

    sample_columns = [
        "sample",
        "source_batch",
        "source_file",
        "cas_variant",
        "direct_repeat",
        "experiment",
        "guide_promoter",
        "donor_recruitment",
        "strain",
        "crRNA_id",
        "timepoint",
        "replicate",
        "generations",
    ]
    samples = (
        variants[sample_columns]
        .drop_duplicates()
        .sort_values(["sample", "source_batch"])
    )

    files = {
        "endogenous_intended_loci.csv": intended,
        "endogenous_guide_features.csv": guide_features,
        "endogenous_sample_metadata.csv": samples,
        "endogenous_freebayes_variants_publication.csv": variants,
    }
    manifest = []
    for filename, data in files.items():
        path = OUTPUT / filename
        data.to_csv(path, index=False)
        manifest.append(
            {
                "file": filename,
                "rows": len(data),
                "columns": len(data.columns),
                "workflow": "endogenous amplicon sequencing",
            }
        )
    return manifest


def export_fcs_inputs() -> list[dict]:
    source = pd.read_csv(
        ANNOTATIONS / "fcs_publication_platemap.csv",
        low_memory=False,
    )
    source = source[source["raw_fcs_name"].notna()].copy()

    grouped = (
        source.groupby("raw_fcs_name", as_index=False, dropna=False)
        .agg(
            figure=("figure", join_unique),
            paper_panel=("paper_panel", join_unique),
            experiment_key=("experiment_key", lambda x: x.dropna().iloc[0]),
            plate=("plate", lambda x: x.dropna().iloc[0] if x.notna().any() else None),
            instrument_well=(
                "instrument_well",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            biological_sample=(
                "biological_sample",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            well_position=(
                "well_position",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            assay=("assay", lambda x: x.dropna().iloc[0]),
            cas_variant=("cas_variant", lambda x: x.dropna().iloc[0]),
            cas_promoter=(
                "cas_promoter",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            guide_promoter=(
                "guide_promoter",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            direct_repeat=(
                "direct_repeat",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            donor_recruitment=(
                "donor_recruitment",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            guide_id=(
                "guide_id",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            strain=(
                "strain",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            replicate=("replicate", lambda x: x.dropna().iloc[0]),
            timepoint=("timepoint", lambda x: x.dropna().iloc[0]),
            generations=(
                "generations",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            control=(
                "control",
                lambda x: x.dropna().iloc[0] if x.notna().any() else None,
            ),
            editing_efficiency_pct=(
                "editing_efficiency_pct",
                lambda x: float(pd.to_numeric(x, errors="coerce").dropna().iloc[0]),
            ),
            raw_fcs_available_locally=("raw_fcs_found", "max"),
            source_table=("source_table", join_unique),
        )
        .rename(columns={"raw_fcs_name": "name"})
        .sort_values(["experiment_key", "timepoint", "name"])
    )

    assay = grouped["assay"].astype(str).str.upper()
    efficiency_fraction = grouped["editing_efficiency_pct"] / 100
    grouped["gfp_positive_fraction"] = np.where(
        assay.isin(["GFP_ON", "ON", "RESTORE", "RESTORE_GFP"]),
        efficiency_fraction,
        1 - efficiency_fraction,
    )

    plate_columns = [
        "name",
        "experiment_key",
        "plate",
        "instrument_well",
        "biological_sample",
        "well_position",
        "assay",
        "cas_variant",
        "cas_promoter",
        "guide_promoter",
        "direct_repeat",
        "donor_recruitment",
        "guide_id",
        "strain",
        "replicate",
        "timepoint",
        "generations",
        "control",
        "figure",
        "paper_panel",
        "raw_fcs_available_locally",
        "source_table",
    ]
    platemap = grouped[plate_columns]

    pop_stats = pd.DataFrame(
        {
            "name": grouped["name"],
            "Population": "/GFP",
            "Parent": "root",
            "Frequency": grouped["gfp_positive_fraction"],
            "ParentFrequency": 1.0,
        }
    )

    expected = grouped[
        [
            "name",
            "assay",
            "gfp_positive_fraction",
            "editing_efficiency_pct",
        ]
    ]

    files = {
        "fcs_publication_platemap.csv": platemap,
        "fcs_population_frequencies_publication.csv": pop_stats,
        "fcs_expected_efficiencies_publication.csv": expected,
    }
    manifest = []
    for filename, data in files.items():
        path = OUTPUT / filename
        data.to_csv(path, index=False)
        manifest.append(
            {
                "file": filename,
                "rows": len(data),
                "columns": len(data.columns),
                "workflow": "FCS GFP efficiency",
            }
        )
    return manifest


def write_manifest(rows: list[dict]) -> None:
    descriptions = {
        "endogenous_intended_loci.csv": (
            "Exact intended chromosome, edited nucleotide position, and alternate "
            "allele used to assign HDR."
        ),
        "endogenous_guide_features.csv": (
            "Guide, donor, PAM context, GC content, and DeepCpf1 annotations."
        ),
        "endogenous_sample_metadata.csv": (
            "One row per endogenous amplicon sample included in Figure 4."
        ),
        "endogenous_freebayes_variants_publication.csv": (
            "Count-level concatenated FreeBayes rows. AO, RO, and DP are retained; "
            "allele frequencies are deliberately not supplied."
        ),
        "fcs_publication_platemap.csv": (
            "One row per unique plotted FCS file with path-independent metadata."
        ),
        "fcs_population_frequencies_publication.csv": (
            "GFP-gate frequencies corresponding to manuscript samples. These permit "
            "testing the annotation and efficiency-conversion step without raw FCS files."
        ),
        "fcs_expected_efficiencies_publication.csv": (
            "Expected GFP fraction and editing percentage for regression checking."
        ),
    }
    manifest = pd.DataFrame(rows)
    manifest["description"] = manifest["file"].map(descriptions)
    manifest["sha256"] = manifest["file"].map(lambda name: sha256(OUTPUT / name))
    manifest.to_csv(OUTPUT / "manifest.csv", index=False)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = export_endogenous_inputs() + export_fcs_inputs()
    write_manifest(manifest)
    print(f"Wrote {len(manifest)} reproducibility inputs to {OUTPUT}")


if __name__ == "__main__":
    main()
