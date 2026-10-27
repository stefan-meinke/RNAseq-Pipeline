# Pipeline audit (Stage 1–2)

Audit date: 2026-10-27. Read-only: no existing file or result was modified.

## 0. Scope and sources

| Source | Content | Notes |
|---|---|---|
| GitHub `main` (`793e344`) | `pipeline.smk`, `config.yml`, `cluster.yml`, `manifest.scm`, `environment.yml`, `run_pipeline_{guix,conda}.sh`, `colData.xlsx`, `R_scripts/*.R`, `README.md`, `workflow.png` | "Original working implementation" |
| Local `D:\RNAseq_pipeline` | `pipeline.smk`, `config.yml`, `cluster.yml`, `manifest.scm`, `channels.scm`, `run_pipeline_guix_time-machine.sh`, `samples.xlsx`, `meta.xlsx`, `contrasts.xlsx`, `Microsoft Word-Dokument (neu).docx` | Partial modernization attempt |
| Upstream PSI scripts (`github.com/MIAOKUI/PSI`) | `dexseq_count.py`, `exclusion_count.py`, `psi_calculation.py` (Python 2) | Read to understand the algorithm. **Your Python 3 ports (`PSI_scripts/`) are in neither source**, see §6 |
| Guix channels at the commits pinned in `channels.scm` | `guix` `e12455a` (2025-03-21), `guix-bioc` `aa10484` (2025-03-20), `guix-cran` `36cc25c` (2025-03-18) | Fetched and searched to check the manifest against them |

The Word document is an earlier AI review of the local folder. It states that `pipeline.smk` was not available to it, so its rule-level statements are guesses. I reused its valid points and did not rely on the rest.

Legend: **[C]** = confirmed (read in code or verified against the source); **[P]** = potential, needs validation on real data or on the cluster.

## 1. Existing workflow architecture

Both versions are one monolithic Snakefile (`pipeline.smk`, about 500 lines) driven by a flat `config.yml` and an old-style `cluster.yml` (`--cluster-config`).

```
raw FASTQ ──fastp (-D -c)──► trimmed FASTQ ──FastQC (trimmed only)──► MultiQC (FastQC only)
                                   │
                                   ├──► STAR (2-pass, GeneCounts) ──► sorted BAM ──samtools index
                                   │         │      │
                                   │         │      └── SJ.out.tab ──awk──► junction BED ──exclusion_count.py ─┐
                                   │         │                                                                ├─► psi_calculation.py ──► .psi ──► PSIplot.R (optional)
                                   │         └── BAM ──dexseq_count.py (HTSeq, flattened GFF) ──► .inclusion ─┘
                                   │         └── BAM ──featureCounts (per sample) ──► DESeq2 R scripts  [out of scope]
                                   │         └── BAMs ──rMATS-turbo (per contrast) ──► rmats_analysis.R
                                   ├──► Salmon quant (selective alignment, -l A)
                                   └──► SPLASH2 (R1 only)                                              [out of scope]
star_index / salmon_index: built, or "dummy" rules that touch a marker when a prebuilt index is used.
```

PSI depends on: the STAR BAM (inclusion counts), STAR `SJ.out.tab` (exclusion counts), the DEXSeq-flattened GFF (`exonic_part` features), and the read length. It does **not** use featureCounts, Salmon or rMATS. rMATS depends only on the STAR BAMs, the GTF and the group definitions. So the two splicing analyses are already independent, which is correct and will be kept.

### Local version vs GitHub version

The local `pipeline.smk` is the GitHub file plus a half-finished migration pasted into the middle. It **cannot be parsed into a working DAG** as it stands [C]:

- Two `configfile:` directives (`config.yml` and `config/config.yml`) and two different `SAMPLES` definitions.
- `validate(config, "workflow/schemas/config.schema.yaml")`, `samples.schema.yaml` and `workflow/scripts/validate_inputs.py` are referenced but do not exist.
- In the local `config.yml`, `samples:` and `contrasts:` are now *paths* (`samples/samples.tsv`), but the code still calls `config["samples"].keys()` and `config["contrasts"].keys()` → `AttributeError` at parse time.
- Keys used by the rules are missing from the local config: `readLength`, `sjdbOverhang`, `rmats_type`, `libType`, `gff`, `colData`, `psi_gtf`, `target_genes`, `splash`, `organism`.
- `star: sjdb_overhang: 149 (read_length-1)` parses as the *string* `"149 (read_length-1)"`.
- `modules:` (with the typo `trimmung`) and `tools:` are two parallel sets of on/off switches; no rule reads either.
- `get_fastq()` was changed so that read 2 is looked up as `_R3` first, which is specific to one sequencing run.
- The rMATS rule now reads `config["colData"]` as a TSV with columns `Sample`/`group`, but there is no `colData` key and the provided metadata is `.xlsx`.

## 2. File-by-file assessment

Columns: Used? / Correct? / Duplicates? / Hard-coded? / Overhead? / Deps documented? / Tracked by Snakemake? / Verdict.

| File | Purpose | Assessment | Verdict |
|---|---|---|---|
| `pipeline.smk` (GitHub) | All rules | Used. Several scientific/technical bugs (§3). Hard-coded relative output dirs (`rmats/`, `Results/`, `Splash2/`) mixed with absolute config dirs. No `log:` directives (except one), no benchmarks, no threads for fastp/fastqc. Several outputs not declared (fastp reports, STAR `Log.final.out`, Salmon aux dir). | **Replace** with a modular workflow; keep every tool call's validated parameters as defaults (or as a documented "legacy" preset) |
| `pipeline.smk` (local) | Migration attempt | Not runnable (§1). Useful ideas: sample table with `fastq_1/fastq_2`, schema validation, bais as rMATS input, empty-group check, extra Salmon flags | **Replace**; carry the ideas over |
| `config.yml` (GitHub) | Flat config with sample dict | Sample names, groups (implicitly via names) and contrasts encoded in config; three separate layout switches (`rmats_type`, `feature_type`, `psi.paired`) and three strandedness switches (`libType`, `psi.strand_specificity`, Salmon `-l A`), none for featureCounts | **Replace** with one schema-validated `config.yml` |
| `config.yml` (local) | Hybrid of new nested + old flat keys | Absolute `/fast/AG_Gotthardt/...` paths, a `/gnu/store/...rmats.py` path (breaks after `guix gc`), duplicated switches, string-typed numbers | **Replace**; reuse its reference paths in your real config |
| `cluster.yml` (both) | `--cluster-config` resources | `--cluster-config` was deprecated in Snakemake 6 and removed in 8. `__default__` log uses `{wildcards.sample}` and fails for rules without it unless overridden. Memory without units | **Replace** with `resources:` in rules + `profiles/slurm/config.yaml` |
| `manifest.scm` (local) | Guix packages | Pins `r-apeglm@1.34.0`, which **does not exist** in the pinned channels (only 1.28.0) → `guix shell` aborts [C]. rMATS-turbo is not in any pinned channel [C]. Contains many downstream-only R packages, `singularity@2.6.1` (2018) and `python-deeptools`, `wiggletools` that no rule uses | See `environment_strategy.md` |
| `manifest.scm` (GitHub) | Same, older | Duplicate `r-scales`; unpinned `r-apeglm` | Superseded |
| `channels.scm` (local) | Pins Guix channels | Good practice. The `guix` URL is `git.savannah.gnu.org`; Guix moved to Codeberg in 2025 and only guaranteed the Savannah mirror for ≥1 year [P] | See `environment_strategy.md` |
| `environment.yml` (GitHub) | Conda env | Uses `bioconda` and `defaults` (both disallowed on your cluster). Single monolithic env incl. R/DEXSeq | **Remove** (replace by conda-forge driver env) |
| `run_pipeline_guix.sh` (GitHub) | Launcher | Calls `pipeline_guix.smk`, which does not exist [C] | Replace |
| `run_pipeline_conda.sh` (GitHub) | Launcher | Calls `pipeline_deseq_guix_shell.smk`, which does not exist; comment says "Guix" [C] | Replace |
| `run_pipeline_guix_time-machine.sh` (local) | Launcher | **Only `snakemake --version` runs inside Guix.** The actual `snakemake ...` call is a separate command and runs with whatever `snakemake` is on the host PATH, i.e. the cluster's old installation [C]. No `set -euo pipefail` | Replace (one-line launcher or none) |
| `colData.xlsx` (GitHub) | Sample → group | `Sample`, `group`. Excel input needs `readxl`/`openpyxl` | Replace by `meta.tsv` |
| `samples.xlsx` (local) | Draft sample sheet | `sample, fastq_1, fastq_2, group, genotype, treatment, sex, batch` – good shape, but mixes FASTQ paths with biology | Split into `samples.tsv` + `meta.tsv` |
| `meta.xlsx` (local) | Metadata, 63 samples | `group` values are `WTCtrlm` style while `samples.xlsx`/`contrasts.xlsx`/`config.yml` use `WT_Ctrl_m` → no contrast would match [C]. Group `TiNxTiPc_HFpEF_f` has 3 replicates, others 4 (probably real; will be reported, not rejected). Sheet 2 is a file-rename table for another dataset | Convert to `meta.tsv` |
| `contrasts.xlsx` (local) | Contrast table | Second contrast compares `TiN_HFpEF_m` against `TiN_Ctrl_f` (male vs female control) – very likely a typo [C/P] | Move contrasts into `config.yml` |
| `R_scripts/rmats_analysis.R` | Summarise/plot rMATS output | Downstream. Reads `config.yml` from the working directory and the key `rmats_filter`, which the local config renamed → breaks [C] | Keep as optional report step, pass filters as arguments |
| `R_scripts/PSIplot.R` | PSI plots for target genes | Downstream. Reads `colData` as xlsx and `read_yaml("config.yml")` hard-coded | Keep outside the core pipeline (optional) |
| `R_scripts/dge_analysis.R`, `geneExpression_analysis.R`, `enrichment_analysis.R` | DESeq2 / enrichment | Out of scope by your instructions | Leave untouched in a `downstream/` folder |
| SPLASH2 rules | k-mer analysis | Out of scope; external binary at absolute path. Your Word notes also say "remove SPLASH2" | Remove from the core workflow (confirm) |
| `PSI_scripts/*.py` | Exon-level PSI | Referenced, but **not present** in either source | Need your ported copies (§6) |
| Flattened ("reduced") GFF | PSI annotation | Produced manually outside the pipeline with `dexseq_prepare_annotation.py`; options used are undocumented | Add a reproducible rule, or accept a provided file + record checksum |
| `README.md` | Docs | Describes Conda+Guix, scripts that don't exist, manual PSI patching | Rewrite |

## 3. Bugs and fragile assumptions

### 3.1 Confirmed problems that change biological results

| # | Problem | Where | Effect |
|---|---|---|---|
| B1 | **fastp `-D` (deduplication) is enabled.** For RNA-seq without UMIs this removes reads from highly expressed genes and short transcripts and distorts counts, splicing ratios and Salmon estimates. The usual practice is not to deduplicate before alignment (mark duplicates for QC only). | `fastp` | Affects every downstream number. Needs your decision: was it deliberate? |
| B2 | **featureCounts is always run unstranded** (no `-s`), even though the library is configured as reverse-stranded elsewhere. | `featureCounts` | Antisense reads from overlapping genes are counted for the wrong gene |
| B3 | **featureCounts `-p` without `--countReadPairs`.** Since Subread 2.0.2, `-p` only declares paired-end input; fragments are counted only with `--countReadPairs`. With Subread 2.0.3 (your pinned version) each mate is counted separately (~2× counts, mates not independent). | `featureCounts` | Inflated counts; matters for DESeq2 dispersion. Older runs with Subread < 2.0.2 counted fragments, so **results depend on the Subread version** |
| B4 | featureCounts `-O`: a read overlapping two genes is counted for both. Combined with B2 this double-counts reads in antisense-overlapping gene pairs. | `featureCounts` | Probably not intended for gene-level DESeq2 input; needs your decision |
| B5 | **rMATS group assignment by regex on sample names** (GitHub version): `get_group("WT_Ctrl_m_1")` returns `"WT"`. For your HFpEF design every `WT_*` sample lands in group `WT`, so `WT_HFpEF_m` vs `WT_Ctrl_m` gets zero samples (or wrong samples in other naming schemes). | `get_group`, `rmats` | Wrong or empty comparisons |
| B6 | Salmon index built by the pipeline has **no decoys**, while the prebuilt one in the local config is decoy-aware (`..._salmon_index_decoy`). | `salmon_index` | Results differ depending on which path was taken; not recorded |
| B7 | Three inconsistent group naming conventions in the local metadata (`WTCtrlm` vs `WT_Ctrl_m`) and a sex-mismatched contrast. | metadata | Contrasts would fail or compare the wrong groups |

### 3.2 Confirmed technical bugs / silent-failure risks

| # | Problem | Where |
|---|---|---|
| T1 | fastp gets no `-h/-j`, so every job writes `fastp.html`/`fastp.json` into the working directory; concurrent jobs overwrite each other and MultiQC never sees trimming stats. | `fastp` |
| T2 | With `build_star_index: false` the "dummy" rule runs `touch <index>/SA`: if the index is missing or incomplete it **creates an empty `SA` file inside a shared reference directory**, and the error only surfaces later inside STAR. Same pattern writes `salmon_index.built` into the shared Salmon index. | `star_index`, `salmon_index` |
| T3 | No index validation: a non-empty directory is assumed valid; genome/GTF/STAR-version compatibility is never checked. | indices |
| T4 | `psi_calculation` shell starts with `if [ ! -f {output} ]`: works only because Snakemake deletes outputs first; if outputs are protected or the rule is run outside Snakemake it silently skips. | `psi_calculation` |
| T5 | STAR `Log.final.out`, Salmon `aux_info/`, `lib_format_counts.json`, rMATS `summary.txt` are not declared → not tracked, not reported. | several |
| T6 | rMATS depends on *all* BAMs, not just those of the contrast; `--tmp` inside the output directory keeps large temp files forever. | `rmats` |
| T7 | Relative output dirs (`rmats/`, `Results/`, `Splash2/`) depend on the launch directory; config output dirs are absolute. Two runs of different datasets from the same directory overwrite each other's rMATS output. | global |
| T8 | Missing raw-read FastQC; MultiQC only scans the trimmed FastQC dir. | `fastqc`, `multiqc` |
| T9 | Read layout is guessed from file names (`_R1`/`_1`, local: `_R3`); single-end is not supported anywhere (fastp, STAR, Salmon rules hard-code R1+R2). | `get_fastq` |
| T10 | Only `config.yml` read from CWD; `rmats_analysis.R` and `PSIplot.R` also read `config.yml` from CWD, so `--configfile` elsewhere breaks them. | R scripts |
| T11 | Launchers reference non-existent Snakefiles; the time-machine launcher runs Snakemake outside Guix (see §2). | launchers |
| T12 | `manifest.scm` asks for a package version absent from the pinned channels (`r-apeglm@1.34.0`) → environment creation fails. Most likely cause of your Guix problem, together with T11. | `manifest.scm` |

### 3.3 Potential problems (need validation)

| # | Problem | How to check |
|---|---|---|
| P1 | `sj_bed` awk maps STAR strand `0` (undefined, non-canonical motif) to `-`. Exclusion counts for such junctions land on the minus strand. With the annotation filter only annotated boundaries pass, so the effect is probably small. | Count rows with `$4==0` in `SJ.out.tab` that pass the filter |
| P2 | `psi_calculation.py` walks the inclusion and exclusion files in parallel (`zip`) and stops at the first line starting with `_` (HTSeq summary rows). This relies on gene IDs sorting before `_` in ASCII. True for Ensembl `ENS…` IDs; for lowercase IDs it would silently write an empty/truncated PSI file. | Assert line counts in a validation step |
| P3 | `psi_calculation.py -t` (low-coverage threshold) defaults to 0, so the low-inclusion flag is always `1`. `PSIplot.R` reads that column. | Was a threshold intended? (Schafer et al. use the flag to mark low-coverage exons) |
| P4 | PSI normalisation uses one fixed read length (150). After adapter/quality trimming lengths vary; the formula `inc/(L-1+len)`, `exc/(L-1)` is an approximation either way. Keep as is, but make `read_length` explicit and recorded. | — |
| P5 | `dexseq_count.py -r pos` on coordinate-sorted BAMs buffers unmatched mates in memory (64 GB / 12 h in `cluster.yml`). | Benchmarks |
| P6 | STAR without `--limitBAMsortRAM` can fail on deep samples during sorting. | Watch for "not enough memory for BAM sorting" |
| P7 | STAR parameters `--seedPerWindowNmax 15` (default 50), `--alignSJoverhangMin 8`, `--alignSJDBoverhangMin 1` are non-default. Validated in your past runs, so they stay as defaults, but they are documented as deviations. | — |
| P8 | Salmon `-l A` (auto-detect) is fine, but the detected library type is never compared with the configured strandedness. | Parse `lib_format_counts.json` |
| P9 | Which `dexseq_prepare_annotation.py` option (`-r yes/no`, aggregation of overlapping genes) produced the reduced GFF is undocumented; exonic part IDs differ between the two. | Ask / record |
| P10 | Your Python 3 port of `dexseq_count.py` is written against an HTSeq API from ~2012 (`raise SystemError, "..."`, `HTSeq.__version__ < '0.5.4p4'`). HTSeq 2.x still provides the used calls, but behaviour must be regression-tested against old `.inclusion` files. | Regression test |

## 4. Reproducibility problems

1. No working environment definition: `environment.yml` uses forbidden channels; `manifest.scm` fails to resolve; the launcher runs the host Snakemake.
2. rMATS-turbo, the PSI scripts, the flattened GFF and SPLASH2 are outside every environment; rMATS is referenced by a `/gnu/store` path that disappears after garbage collection.
3. No recorded tool versions, reference checksums, config snapshot or pipeline Git commit per run.
4. Parameter changes do not trigger reruns of the right rules (most parameters are interpolated via `config[...]` inside shell strings rather than `params:`).
5. featureCounts output silently depends on the Subread version (B3).
6. Prebuilt indices are not linked to the genome/GTF they came from.

## 5. Computational and I/O bottlenecks

| Area | Observation | Proposal | Changes results? |
|---|---|---|---|
| PSI inclusion counting | `dexseq_count.py` (pure Python, HTSeq, single thread) is by far the slowest step (12 h / 64 GB budget). | Keep as reference implementation. Optional later: an equivalent featureCounts-on-exonic-parts implementation, only after a regression test shows identical counts | Only if the optional replacement is adopted |
| featureCounts | One job per sample re-parses the GTF each time. | One featureCounts call over all BAMs (fast, gives the matrix and one summary directly) | No |
| Trimmed FASTQs | Stored permanently (≈ doubles FASTQ storage). | Configurable `temp()` (default: keep, because Salmon/STAR reruns need them) | No |
| fastp | No `-w`, rule has 1 thread. | Declare threads | No |
| STAR | Genome loaded per job (unavoidable on SLURM); `_STARtmp` written to the shared FS. | `--outTmpDir` on node-local `$TMPDIR`; declare `Log.final.out`; optionally drop `ReadsPerGene` if unused | No |
| Salmon (local) | `--numGibbsSamples 100 --dumpEq --dumpEqWeights` multiply runtime and output size; only needed for fishpond/swish or terminus. | Off by default, configurable | No (quant.sf unchanged) |
| rMATS | Temp files kept inside output dir. | Temp dir on scratch, removed after success | No |
| Shared FS | Many small log files from SLURM + per-rule logs. | Per-rule logs in `logs/`, SLURM logs via executor plugin | No |

## 6. Scientific correctness notes for the splicing modules

- **rMATS** should use the STAR BAMs directly (it does). Correct inputs: `-t paired|single` from the layout, `--libType` derived from the single `strandedness` setting (`reverse` → `fr-firststrand`, `forward` → `fr-secondstrand`, `none` → `fr-unstranded`), `--readLength` + `--variable-read-length` (needed after trimming), groups from `meta.tsv`. Replicates are passed comma-separated via `--b1/--b2` files.
- **Exon-level PSI (Schafer et al. 2015)** is a different quantity from rMATS PSI: per exonic part of the DEXSeq-flattened annotation, inclusion = reads overlapping the part (HTSeq), exclusion = uniquely mapped junction reads (STAR column 7) spanning the part with both ends at annotated exonic-part boundaries, normalised by `L-1+len` and `L-1`. The workflow will keep the two analyses in separate modules and output folders.
- I read the original Python 2 scripts. I need **your Python 3 versions** (the `PSI_scripts/` folder on the cluster) to integrate them unchanged; I will not rewrite the algorithm. They are GPL-3, like this repository, so they can be vendored into `workflow/scripts/psi/` with the original copyright header and a note on your changes.

## 7. Recommended changes, ranked

| Priority | Change |
|---|---|
| Essential | Decide on B1 (fastp `-D`), B2–B4 (featureCounts strand/`-p`/`-O`). Provide a `legacy` parameter preset so old results can be reproduced and compared. |
| Essential | One schema-validated config + `samples.tsv` + `meta.tsv`; groups and comparisons from metadata only (fixes B5, B7, T9). |
| Essential | Working environment (see `environment_strategy.md`), Snakemake run from that environment, not the cluster module (T11, T12). |
| Essential | Safe index handling: validate, reuse, or build into an empty directory; never touch a shared index (T2, T3, B6). |
| Essential | Declare all outputs, per-rule logs, fastp reports (T1, T5). |
| Recommended | Modular rule files, SLURM executor profile, raw + trimmed FastQC, final MultiQC with all modules, provenance folder. |
| Recommended | Single featureCounts call + gene-count matrix; Salmon library-type check; temp dirs on node-local scratch. |
| Optional | Faster PSI inclusion counting (only with regression proof); Salmon Gibbs/EQ output; trimmed FASTQ as `temp()`. |
