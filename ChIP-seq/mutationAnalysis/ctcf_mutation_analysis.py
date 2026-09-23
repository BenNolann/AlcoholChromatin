#!/usr/bin/env python3
import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# CONFIG - edit these to match your layout
# ---------------------------------------------------------------------------

WORKDIR = Path("/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/ENCODE/CTCF/bams/")

# gatk/bcftools/bedtools/bgzip all run inside this singularity container.
SIF_PATH = "/Zulu/bnolan/Projects/Personal/Ethanol/Integrative_analysis/ENCODE/CTCF/bams/gatk_latest.sif"
# Extra singularity bind mounts, if your data isn't visible inside the
# container by default, e.g. "-B /Zulu:/Zulu" - leave empty if /Zulu is
# already visible inside the container without one.
SINGULARITY_BIND = ""

REF = "hg38_noaltchromosomes.fa"
UT_BAM = "UT_merged.rg.bam"
E6_BAM = "E6_merged.rg.bam"
PEAKS = "E6vUT.tsv"
BLACKLIST = "hg38-blacklist.v2.bed"
GNOMAD = "af-only-gnomad.hg38.vcf.gz"
UT_SAMPLE = "UT"
E6_SAMPLE = "E6"

# CTCF gained/lost/unbiased classification, parsed directly from E6vUT.tsv.
# Columns are 1-indexed as you'd view the file: col1-3 = chrom/start/end,
# col5 = log2 fold change (E6 vs UT). Set to 0-indexed pandas positions below.
PEAKS_CHROM_COL = 0
PEAKS_START_COL = 1
PEAKS_END_COL = 2
PEAKS_LOG2FC_COL = 4  # column 5, 1-indexed
PEAKS_HAS_HEADER = False  # set True if E6vUT.tsv has a header row

GAINED_THRESHOLD = 1.0   # log2FC > this  -> gained (more CTCF occupancy in E6)
LOST_THRESHOLD = -1.0    # log2FC < this  -> lost   (less CTCF occupancy in E6)
# anything in between -> unbiased

N_PERMUTATIONS = 10_000
RANDOM_SEED = 0

CATEGORY_COLORS = {
    "gained": "#B2182B",    # dark red  - matches your existing GCA+GGG coloring
    "lost": "#2166AC",      # dark blue - matches your existing CTA+GAT coloring
    "unbiased": "#777777",  # neutral gray
}
CATEGORY_ORDER = ["unbiased", "gained", "lost"]

OUTDIR = WORKDIR / "ctcf_permutation_analysis"

ORIGINAL_CWD = Path.cwd()

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def containerize(cmd):
    """Wrap a shell command to run inside SIF_PATH via singularity exec.

    Everything (gatk, bcftools, bedtools, bgzip, and the plain sort/awk
    steps) is wrapped uniformly rather than picked apart per-tool, so a
    compound command with && / > still runs in one consistent environment.
    `bash -c` + shlex.quote handles cmd strings that already contain their
    own quoting (e.g. the awk step) safely.
    """
    bind = f"{SINGULARITY_BIND} " if SINGULARITY_BIND else ""
    return f"singularity exec {bind}{SIF_PATH} bash -c {shlex.quote(cmd)}"


def run(cmd, **kwargs):
    """Run a shell command inside the singularity container, echoing it,
    and abort loudly on failure."""
    wrapped = containerize(cmd)
    print(f"$ {wrapped}")
    subprocess.run(wrapped, shell=True, check=True, **kwargs)


def exists(path):
    return Path(path).exists() and Path(path).stat().st_size > 0


def skip_or_run(output_path, cmd, force=False):
    """Idempotency helper: skip expensive steps (Mutect2 etc.) if the output
    already exists, unless --force is passed."""
    if exists(output_path) and not force:
        print(f"[skip] {output_path} already exists")
        return
    run(cmd)


# ---------------------------------------------------------------------------
# Stage 1: variant calling (mutect2.sh, reproduced + sample-order fix)
# ---------------------------------------------------------------------------


def stage1_call_variants(force=False):
    os.chdir(ORIGINAL_CWD / WORKDIR)

    # --- 0. Sort peaks, convert to BED ---
    skip_or_run(
        "E6vUT.sorted.bed",
        f"sort -k1,1 -k2,2n {PEAKS} > E6vUT.sorted.tsv && "
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3}}' E6vUT.sorted.tsv "
        f"> E6vUT.sorted.bed",
        force,
    )

    # --- 1. Call somatic-style variants, E6 vs UT ---
    skip_or_run(
        "E6_vs_UT_somatic_unfiltered.vcf.gz",
        f"gatk Mutect2 -R {REF} "
        f"-I {E6_BAM} -tumor {E6_SAMPLE} "
        f"-I {UT_BAM} -normal {UT_SAMPLE} "
        f"-L E6vUT.sorted.bed --germline-resource {GNOMAD} "
        f"-O E6_vs_UT_somatic_unfiltered.vcf.gz --f1r2-tar-gz f1r2.tar.gz",
        force,
    )

    # --- 2. Learn read-orientation model ---
    skip_or_run(
        "read-orientation-model.tar.gz",
        "gatk LearnReadOrientationModel -I f1r2.tar.gz "
        "-O read-orientation-model.tar.gz",
        force,
    )

    # --- 3. Apply statistical filtering ---
    skip_or_run(
        "E6_vs_UT_somatic.filtered.vcf.gz",
        f"gatk FilterMutectCalls -R {REF} -V E6_vs_UT_somatic_unfiltered.vcf.gz "
        f"--ob-priors read-orientation-model.tar.gz "
        f"-O E6_vs_UT_somatic.filtered.vcf.gz",
        force,
    )

    # --- 4. PASS only ---
    skip_or_run(
        "E6_vs_UT_somatic.pass.vcf.gz",
        "bcftools view -f PASS E6_vs_UT_somatic.filtered.vcf.gz -Oz "
        "-o E6_vs_UT_somatic.pass.vcf.gz && "
        "bcftools index -f E6_vs_UT_somatic.pass.vcf.gz",
        force,
    )

    # --- 5. Biallelic SNPs only ---
    skip_or_run(
        "E6_vs_UT_somatic.snps.vcf.gz",
        "bcftools view -v snps -m2 -M2 E6_vs_UT_somatic.pass.vcf.gz -Oz "
        "-o E6_vs_UT_somatic.snps.vcf.gz && "
        "bcftools index -f E6_vs_UT_somatic.snps.vcf.gz",
        force,
    )

    # --- 6. Blacklist subtraction ---
    skip_or_run(
        "E6_vs_UT_somatic.snps.noblacklist.vcf.gz",
        f"bedtools intersect -v -a E6_vs_UT_somatic.snps.vcf.gz -b {BLACKLIST} "
        f"-header > E6_vs_UT_somatic.snps.noblacklist.vcf && "
        f"bgzip -f E6_vs_UT_somatic.snps.noblacklist.vcf && "
        f"bcftools index -f E6_vs_UT_somatic.snps.noblacklist.vcf.gz",
        force,
    )

    # --- Sample-order check (the fix) ---
    vcf = "E6_vs_UT_somatic.snps.noblacklist.vcf.gz"
    sample_list = subprocess.run(
        containerize(f"bcftools query -l {vcf}"), shell=True, check=True,
        capture_output=True, text=True,
    ).stdout.split()
    print(f"VCF sample order: {sample_list}")
    if set(sample_list) != {E6_SAMPLE, UT_SAMPLE}:
        sys.exit(
            f"ERROR: expected samples {{{E6_SAMPLE}, {UT_SAMPLE}}} in {vcf}, "
            f"found {sample_list}. Fix UT_SAMPLE/E6_SAMPLE config."
        )

    # --- 7. Extract full PASS biallelic SNP table, column order pinned ---
    skip_or_run(
        "somatic_candidates_E6vUT.tsv",
        # NOTE: separator goes *before* each per-sample block ([\t%GT...])
        # rather than after (...%AF\t]) - a trailing \t inside the repeated
        # block leaves a stray empty 13th field on every line once you have
        # 2 samples, which silently misaligns fixed-width parsing downstream.
        f"bcftools query -s {E6_SAMPLE},{UT_SAMPLE} "
        f"-f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD\\t%DP\\t%AF]\\n' "
        f"{vcf} > somatic_candidates_E6vUT.tsv",
        force,
    )

    col_names = [
        "CHROM", "POS", "REF", "ALT",
        "E6_GT", "E6_AD", "E6_DP", "E6_AF",
        "UT_GT", "UT_AD", "UT_DP", "UT_AF",
    ]
    df = pd.read_csv("somatic_candidates_E6vUT.tsv", sep="\t", header=None)
    # Robust to files generated by the old (trailing-\t) format string too:
    # those have one extra, empty trailing column.
    if df.shape[1] == len(col_names) + 1:
        print(
            "NOTE: somatic_candidates_E6vUT.tsv has a stray trailing empty "
            "column (old trailing-tab bcftools format) - dropping it."
        )
        df = df.iloc[:, :len(col_names)]
    elif df.shape[1] != len(col_names):
        sys.exit(
            f"ERROR: somatic_candidates_E6vUT.tsv has {df.shape[1]} columns, "
            f"expected {len(col_names)} (or {len(col_names) + 1} with a "
            f"stray trailing column). Delete it and re-run with --force."
        )
    df.columns = col_names
    for c in ("POS", "E6_DP", "E6_AF", "UT_DP", "UT_AF"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    n_bad = df["POS"].isna().sum()
    if n_bad:
        print(f"WARNING: dropping {n_bad} rows with unparseable POS/AF/DP values")
        df = df.dropna(subset=["POS"])
    df["POS"] = df["POS"].astype(int)
    print(f"Full PASS biallelic SNP candidates: {len(df)}")

    # --- AF-diff candidate table (kept for manual review, NOT used in Stage 2) ---
    df["diff"] = df["E6_AF"] - df["UT_AF"]
    af_diff = df[
        (df["diff"].abs() >= 0.10) & (df["E6_DP"] >= 15) & (df["UT_DP"] >= 15)
    ]
    af_diff.to_csv("af_diff_candidates.tsv", sep="\t", index=False)
    print(f"AF-diff candidates (review list only): {len(af_diff)}")

    os.chdir(ORIGINAL_CWD)
    return df


# ---------------------------------------------------------------------------
# Stage 2: classify CTCF peaks as gained / lost / unbiased from E6vUT.tsv
# ---------------------------------------------------------------------------


def load_ctcf_categories():
    peaks_path = WORKDIR / PEAKS
    df = pd.read_csv(
        peaks_path, sep="\t",
        header=0 if PEAKS_HAS_HEADER else None,
    )
    df = df.rename(columns={
        df.columns[PEAKS_CHROM_COL]: "chrom",
        df.columns[PEAKS_START_COL]: "start",
        df.columns[PEAKS_END_COL]: "end",
        df.columns[PEAKS_LOG2FC_COL]: "log2fc",
    })[["chrom", "start", "end", "log2fc"]].copy()

    df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
    n_bad = df["log2fc"].isna().sum()
    if n_bad:
        print(f"WARNING: {n_bad} rows had non-numeric log2FC (col {PEAKS_LOG2FC_COL + 1}) "
              f"and were dropped - check PEAKS_LOG2FC_COL / PEAKS_HAS_HEADER.")
        df = df.dropna(subset=["log2fc"])

    df["category"] = np.select(
        [df["log2fc"] > GAINED_THRESHOLD, df["log2fc"] < LOST_THRESHOLD],
        ["gained", "lost"],
        default="unbiased",
    )

    print("CTCF peak classification from E6vUT.tsv log2FC "
          f"(gained: >{GAINED_THRESHOLD}, lost: <{LOST_THRESHOLD}):")
    print(df["category"].value_counts().reindex(CATEGORY_ORDER))

    return df


# ---------------------------------------------------------------------------
# Stage 3: permutation test
# ---------------------------------------------------------------------------


def run_permutation_test(snv_df, universe, n_perm=N_PERMUTATIONS, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)

    # universe: chrom, start, end, log2fc, category (from load_ctcf_categories)
    universe = universe.copy()
    universe["length"] = universe["end"] - universe["start"]
    universe = universe[universe["length"] > 0].reset_index(drop=True)

    total_bp_by_cat = universe.groupby("category")["length"].sum()
    print("\nCTCF peak universe (bp) by category:")
    print(total_bp_by_cat)

    # Assign each observed SNV to a category via bedtools intersect
    snv_bed_path = OUTDIR / "observed_snvs.bed"
    OUTDIR.mkdir(parents=True, exist_ok=True)
    snv_bed = snv_df[["CHROM", "POS"]].copy()
    snv_bed["start"] = snv_bed["POS"] - 1
    snv_bed["end"] = snv_bed["POS"]
    snv_bed[["CHROM", "start", "end"]].sort_values(["CHROM", "start"]).to_csv(
        snv_bed_path, sep="\t", header=False, index=False
    )

    universe_bed_path = OUTDIR / "ctcf_universe.bed"
    universe.sort_values(["chrom", "start"])[
        ["chrom", "start", "end", "category"]
    ].to_csv(universe_bed_path, sep="\t", header=False, index=False)

    intersect_out = OUTDIR / "snv_x_category.tsv"
    run(
        f"bedtools intersect -wa -wb -a {snv_bed_path} -b {universe_bed_path} "
        f"> {intersect_out}"
    )
    hits = pd.read_csv(
        intersect_out, sep="\t", header=None,
        names=["chrom", "start", "end", "u_chrom", "u_start", "u_end", "category"],
    )
    n_unassigned = len(snv_bed) - hits["start"].nunique()
    if n_unassigned > 0:
        print(
            f"NOTE: {n_unassigned} SNVs fell outside all classified E6vUT.tsv "
            f"peak intervals and are excluded from the permutation test "
            f"(shouldn't happen often since Mutect2 was itself restricted "
            f"to these peaks via -L, but overlapping/duplicate peak rows "
            f"in E6vUT.tsv could cause edge effects - worth a sanity check)."
        )

    observed_counts = hits["category"].value_counts().reindex(CATEGORY_ORDER, fill_value=0)
    observed_density = observed_counts / (total_bp_by_cat.reindex(CATEGORY_ORDER) / 1000.0)
    n_snvs_in_universe = int(observed_counts.sum())
    print(f"\nObserved SNV counts / density (SNVs per kb) [n={n_snvs_in_universe} SNVs in universe]:")
    print(pd.DataFrame({"count": observed_counts, "density_per_kb": observed_density}))

    # Permutation null: place n_snvs_in_universe SNVs uniformly at random
    # across the peak universe, length-weighted by interval size, and
    # recompute per-category density each time.
    weights = universe["length"].to_numpy() / universe["length"].sum()
    cats = universe["category"].to_numpy()
    n_intervals = len(universe)

    null_density = {cat: np.empty(n_perm) for cat in CATEGORY_ORDER}
    for i in range(n_perm):
        idx = rng.choice(n_intervals, size=n_snvs_in_universe, p=weights, replace=True)
        chosen_cats = cats[idx]
        counts = pd.Series(chosen_cats).value_counts().reindex(CATEGORY_ORDER, fill_value=0)
        density = counts / (total_bp_by_cat.reindex(CATEGORY_ORDER) / 1000.0)
        for cat in CATEGORY_ORDER:
            null_density[cat][i] = density[cat]

    results = []
    for cat in CATEGORY_ORDER:
        obs = observed_density[cat]
        null = null_density[cat]
        # one-sided: is observed *more* likely than chance (upper tail)?
        p_upper = (np.sum(null >= obs) + 1) / (n_perm + 1)
        # two-sided: is observed unusual in either direction?
        null_mean = null.mean()
        p_two_sided = (
            np.sum(np.abs(null - null_mean) >= np.abs(obs - null_mean)) + 1
        ) / (n_perm + 1)
        results.append({
            "category": cat,
            "observed_density_per_kb": obs,
            "null_mean_density_per_kb": null_mean,
            "null_sd": null.std(),
            "p_one_sided_enriched": p_upper,
            "p_two_sided": p_two_sided,
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTDIR / "permutation_results.tsv", sep="\t", index=False)
    print("\nPermutation test results:")
    print(results_df.to_string(index=False))

    return observed_density, null_density, results_df


# ---------------------------------------------------------------------------
# Stage 4: plots
# ---------------------------------------------------------------------------


def plot_permutation_boxplot(observed_density, null_density, results_df):
    fig, ax = plt.subplots(figsize=(7, 6))

    box_data = [null_density[cat] for cat in CATEGORY_ORDER]
    positions = np.arange(len(CATEGORY_ORDER))

    bp = ax.boxplot(
        box_data, positions=positions, widths=0.5, showfliers=False,
        patch_artist=True, medianprops={"color": "black"},
        whiskerprops={"color": "black"}, capprops={"color": "black"},
        boxprops={"edgecolor": "black", "linewidth": 0.8},
    )
    for patch, cat in zip(bp["boxes"], CATEGORY_ORDER):
        patch.set_facecolor(CATEGORY_COLORS[cat])
        patch.set_alpha(0.35)

    # Overlay observed value per category
    for i, cat in enumerate(CATEGORY_ORDER):
        ax.scatter(
            [positions[i]], [observed_density[cat]], marker="D", s=90,
            color=CATEGORY_COLORS[cat], edgecolor="black", zorder=5,
            label="Observed" if i == 0 else None,
        )
        p = results_df.loc[results_df["category"] == cat, "p_two_sided"].values[0]
        ax.annotate(
            f"p={p:.3g}", (positions[i], observed_density[cat]),
            textcoords="offset points", xytext=(12, 0), fontsize=10,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [c.capitalize() for c in CATEGORY_ORDER], fontsize=13
    )
    ax.set_ylabel("SNV density (SNVs / kb of CTCF peak)\n[permutation null]", fontsize=12)
    ax.set_title(
        "Observed SNV density vs. permutation null,\nby CTCF occupancy-change category",
        loc="left", fontsize=14, fontweight="500", pad=15,
    )
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_position(("outward", 10))
    ax.spines["bottom"].set_position(("outward", 10))
    # "upper right" visually collides with the Lost box/whisker in this
    # data (its glyph looked like a spurious 4th diamond near y=0.08) -
    # anchor the legend outside the axes instead so it never sits on data.
    ax.legend(
        loc="lower left", bbox_to_anchor=(0, 1.01), frameon=False,
        fontsize=11, ncol=1, borderaxespad=0,
    )

    fig.tight_layout()
    out = OUTDIR / "ctcf_category_snv_permutation_boxplot.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"\nSaved {out}")


def plot_allelic_fraction(snv_df):
    fig, ax = plt.subplots(figsize=(6, 6))

    af = snv_df["E6_AF"].dropna()
    rng = np.random.default_rng(1)
    x_jitter = rng.normal(0, 0.04, size=len(af))

    ax.scatter(
        x_jitter, af, s=10, alpha=0.35, color="#B2182B", edgecolor="none",
    )
    parts = ax.violinplot([af], positions=[0], widths=0.7, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#777777")
        pc.set_alpha(0.25)

    ax.axhline(0.5, color="black", lw=1.0, linestyle="--")
    ax.text(0.42, 0.51, "clonal heterozygous AF (~0.5)", fontsize=10, va="bottom")
    ax.axhspan(0, 0.2, color="#2166AC", alpha=0.08)
    ax.text(0.42, 0.02, "low/subclonal AF (<0.2)", fontsize=10, color="#2166AC")

    ax.set_xlim(-0.6, 0.9)
    ax.set_xticks([])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Allelic fraction in E6 (ethanol)", fontsize=12)
    ax.set_title(
        f"Allelic fraction of candidate SNVs at CTCF sites (n={len(af)})\n"
        "6-day ethanol treatment",
        loc="left", fontsize=13, fontweight="500", pad=15,
    )
    median_af = af.median()
    ax.annotate(
        f"median AF = {median_af:.3f}", xy=(0, median_af),
        xytext=(0.5, 0.9), fontsize=11,
        arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
    )
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    out = OUTDIR / "ctcf_snv_allelic_fraction.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    plt.rcParams['pdf.fonttype'] = 42
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run every stage even if outputs already exist",
    )
    parser.add_argument(
        "--list-categories-only", action="store_true",
        help="Just classify E6vUT.tsv peaks into gained/lost/unbiased and "
             "exit (do this first to sanity-check PEAKS_HAS_HEADER / "
             "PEAKS_LOG2FC_COL / thresholds before running the slow "
             "variant-calling stage)",
    )
    parser.add_argument(
        "--n-perm", type=int, default=N_PERMUTATIONS,
        help="Number of permutations (default: %(default)s)",
    )
    args = parser.parse_args()
    force, list_categories_only, n_perm = args.force, args.list_categories_only, args.n_perm

    universe = load_ctcf_categories()
    if list_categories_only:
        return None

    OUTDIR.mkdir(parents=True, exist_ok=True)

    snv_df = stage1_call_variants(force=force)
    observed_density, null_density, results_df = run_permutation_test(
        snv_df, universe, n_perm=n_perm
    )
    plot_permutation_boxplot(observed_density, null_density, results_df)
    plot_allelic_fraction(snv_df)

    print(f"\nDone. Outputs in {OUTDIR}/")
    return snv_df, observed_density, null_density, results_df


if __name__ == "__main__":
    main()