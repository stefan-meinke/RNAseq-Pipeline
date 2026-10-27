# Proposed architecture and implementation plan (Stage 3)

Status: **proposal – waiting for approval.** Nothing below has been implemented yet.
Background: [`pipeline_audit.md`](pipeline_audit.md), [`environment_strategy.md`](environment_strategy.md).

Change classes: **E** = essential, **R** = recommended, **O** = optional.

## 1. Environment (E)

As recommended in `environment_strategy.md`: one small driver environment (conda-forge + pinned PyPI: Snakemake 9.x, SLURM executor plugin, MultiQC, HTSeq, pysam, pandas), plus per-rule tools from Apptainer images pinned by digest; fallback to a pinned `tools/` folder or explicitly configured local executables. Snakemake never comes from the cluster module.

```bash
# one-time
conda env create -f envs/driver.yml        # conda-forge only, PyPI packages pinned inside
conda activate rnaseq-driver
# every run
snakemake --configfile config.yml --profile profiles/slurm
```

## 2. Directory structure (E)

Kept deliberately flat; one rule file per module.

```
workflow/
  Snakefile                 # includes rules, defines `all`
  rules/
    common.smk              # load + validate config/samples/meta, helper functions
    references.smk          # STAR/Salmon index check-or-build, gentrome/decoys, flattened GFF
    qc.smk                  # FastQC raw + trimmed, fastp, MultiQC
    align.smk               # STAR, BAM index
    quant.smk               # featureCounts (+ matrix), Salmon
    rmats.smk
    psi.smk
    report.smk              # provenance, final MultiQC, summary report
  scripts/
    validate_inputs.py      # also usable standalone: python workflow/scripts/validate_inputs.py config.yml
    check_index.py          # STAR/Salmon index validation
    merge_counts.py         # featureCounts → gene matrix (TSV, DESeq2-ready)
    summary_report.py       # Markdown/HTML status report
    psi/                    # your Python 3 PSI scripts, vendored unchanged (GPL-3 headers kept)
  schemas/
    config.schema.yaml  samples.schema.yaml  meta.schema.yaml
  envs/driver.yml
config/
  config.example.yml  samples.example.tsv  meta.example.tsv
profiles/
  slurm/config.yaml   local/config.yaml
tests/
  data/ (tiny reference + reads, or a script to create them)   test_validation.py   run_tests.sh
docs/
legacy/                     # old pipeline.smk, cluster.yml, run scripts, manifest (moved, not deleted)
downstream/R_scripts/       # DESeq2/enrichment/PSI-plot R scripts, untouched
```

## 3. The three dataset files (E)

### `samples.tsv` – where the reads are

| column | required | meaning |
|---|---|---|
| `sample` | yes | ID, `[A-Za-z0-9_.-]+`, must exist in `meta.tsv` |
| `unit` | no (default `1`) | lane/run. Several rows per sample → FASTQs are concatenated per sample before fastp (raw FastQC still per unit) |
| `fastq_1` | yes | path (absolute, or relative to `config.yml`) |
| `fastq_2` | no | empty → single-end. Mixed layouts inside one dataset are rejected |

### `meta.tsv` – the biology

`sample`, `group` required; any extra columns allowed (`genotype`, `treatment`, `sex`, `batch`, …). Optional `include` (`true`/`false`) to drop a failed sample without deleting rows. Groups are always read from here, never parsed from names.

### `config.yml` – everything else (sketch)

```yaml
project: HFpEF_TiN_TiPc
samples: samples.tsv
meta: meta.tsv
outdir: /fast/AG_Gotthardt/smeinke/HFpEF/results   # all outputs below this
tmpdir: ${TMPDIR}                                   # node-local scratch for STAR/rMATS temp

library:
  strandedness: reverse      # none | forward | reverse  → mapped to featureCounts -s, rMATS --libType, PSI -s, Salmon check
  read_length: 150           # used by rMATS, PSI, sjdbOverhang (auto = read_length-1)

reference:
  genome_fasta: /fast/.../Mus_musculus.GRCm39.dna.primary_assembly.fa
  gtf: /fast/.../Mus_musculus.GRCm39.115.gtf
  transcriptome_fasta: /fast/.../Mus_musculus.GRCm39.cdna.all.fa   # or "derive" from genome+GTF
  star_index: /fast/AG_Gotthardt/STAR_index/STAR_index_GRCm39.115_sjdb149   # reused if valid, built if absent
  salmon_index: /fast/.../GRCm39.115_salmon_index_decoy
  salmon_decoys: genome      # genome | none
  flattened_gff: auto        # auto = build from GTF with dexseq_prepare_annotation; or a path

modules:                     # true/false
  fastqc: true
  trimming: true
  star: true
  featurecounts: true
  salmon: true
  rmats: true
  psi: true

params:
  fastp: {extra: "-c"}                       # dedup (-D) decision pending, see §9
  star: {twopass: true, extra: "--seedPerWindowNmax 15 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1"}
  featurecounts: {feature_type: exon, attribute: gene_id, count_read_pairs: true, multi_overlap: false, multimapping: false}
  salmon: {libtype: A, extra: "", gibbs_samples: 0}
  rmats: {group_column: group, variable_read_length: true, cutoff: 0.0001}
  psi: {low_coverage_threshold: 0}

comparisons:                 # case vs control, groups from meta.tsv (or: all_pairs: true)
  - {name: WT_HFpEF_m_vs_WT_Ctrl_m, case: WT_HFpEF_m, control: WT_Ctrl_m}

keep_trimmed_fastq: true
tools:                       # optional local executables (recorded as "not controlled")
  rmats: ""                  # e.g. /path/to/rmats.py
  psi_scripts: ""            # default: workflow/scripts/psi
containers: {...}            # pinned image digests, preset in the example config
```

A `legacy: true` preset reproduces the original parameters exactly (fastp `-D`, featureCounts unstranded `-p -O`, no Salmon decoys if that was the case) for regression testing.

### Validation (fails before any job is submitted)

JSON-schema validation of config and both tables, then explicit checks with plain messages: missing/unreadable FASTQs; duplicate `sample`/`unit`; R1 = R2; mixed layouts; samples missing in `meta.tsv` and vice versa; comparison groups absent or with < 2 replicates (warning at 2, error at 0); unbalanced replicate numbers (warning); GTF chromosome names vs genome FASTA; index validity (below); required executables/containers present; strandedness value legal.

## 4. Rule graph (E)

```
[validate] ─► fastqc_raw (per unit)
              merge_units (per sample, temp) ─► fastp ─► fastqc_trimmed
                                                  ├─► star_align ◄── star_index (check-or-build)
                                                  │      └─► bam_index
                                                  │           ├─► featurecounts (all BAMs, one job) ─► gene_counts.tsv
                                                  │           ├─► rmats (per comparison, only its BAMs)
                                                  │           └─► psi_inclusion ─┐
                                                  │   SJ.out.tab ─► junction_bed ─► psi_exclusion ─► psi_calculation
                                                  │   flattened_gff (check-or-build) ──────┘
                                                  └─► salmon_quant ◄── salmon_index (check-or-build)
multiqc_initial  ◄── fastqc_raw, fastp, fastqc_trimmed           (separate sections "Raw"/"Trimmed")
multiqc_final    ◄── + STAR Log.final.out, featureCounts summary, Salmon
summary_report   ◄── multiqc_final, rMATS summaries, PSI status, provenance
```

Outputs (under `outdir`): `qc/{raw,trimmed}/`, `fastp/`, `star/{sample}.bam(.bai)` + logs, `featurecounts/gene_counts.tsv` + `.summary`, `salmon/{sample}/`, `rmats/{comparison}/`, `psi/{sample}.psi`, `multiqc/{initial,final}/`, `report/summary.html`, `provenance/`; per-rule logs in `logs/{rule}/`, benchmarks in `benchmarks/`.

Module switches only change which targets `all` requests; individual modules can be run with `snakemake ... --until featurecounts` or by target file.

## 5. Reference index management (E)

- **STAR:** `check_index.py` checks the configured directory for the full file set (`SA`, `SAindex`, `Genome`, `chrName.txt`, `genomeParameters.txt`, `sjdbList.out.tab`, `exonInfo.tab`, …) and parses `genomeParameters.txt` (STAR genome version, `sjdbOverhang`, `sjdbGTFfile`, `genomeFastaFiles`). It compares these with the config and the chromosome names of the GTF.
  - valid → reused, read-only, never written to;
  - directory absent or empty → built by one `star_index` rule (all alignments depend on it, so there is no race), with a `pipeline_index.json` (tool version, inputs, checksums, parameters) written next to it;
  - present but invalid/incompatible → **stop with a message**; never overwrite.
- **Salmon:** same logic using `versionInfo.json`, `info.json`, `refInfo.json` and the decoy list. Building: decoy-aware "gentrome" (transcripts + genome as decoys), `-k 31`.
- **Flattened GFF for PSI:** built from the GTF by the vendored `dexseq_prepare_annotation.py` with a recorded option, or a provided file whose checksum is recorded.

## 6. Salmon strategy (R)

Keep **mapping-based selective alignment with a decoy-aware index**, on the trimmed reads. That matches your current prebuilt index, avoids the extra transcriptome BAM from STAR (`--quantMode TranscriptomeSAM` would be needed for alignment-based mode, with more storage and no decoy benefit), and is the mode the Salmon authors recommend. `-l A` stays, but the detected library type is compared with `library.strandedness` and a mismatch is flagged in the report. Gibbs samples / EQ classes are optional.

## 7. rMATS and PSI integration (E)

- **rMATS:** one job per comparison, inputs = BAMs of the two groups only, `--b1/--b2` text files written by the rule, `-t` from the layout, `--libType` from strandedness, `--readLength` + `--variable-read-length`, `--tmp` on scratch (deleted on success), `--od` per comparison. Executable: container (default) or `tools.rmats`. Outputs tracked: `summary.txt` + `*.MATS.JC.txt`/`JCEC.txt`. `rmats_analysis.R` stays optional and gets its filters as arguments.
- **PSI:** your Python 3 scripts are called unchanged with the same arguments as before; only file paths, strandedness, layout and read length come from config. The `if [ ! -f ]` guard is removed; a small check confirms that inclusion and exclusion files have the same exonic parts before PSI is computed (guards against P2 in the audit).

## 8. SLURM execution (R)

`profiles/slurm/config.yaml` for the SLURM executor plugin: `executor: slurm`, `jobs: 100`, `latency-wait: 60`, `rerun-incomplete`, `keep-going`, `retries: 1`, `default-resources: slurm_account=…, slurm_partition=…, mem_mb=4000, runtime=60`, and `set-resources` per rule. Rules declare `threads` and `resources: mem_mb, runtime, disk_mb`, so SLURM CPUs always equal the threads the tool uses. Example defaults (from your `cluster.yml`): STAR 16 threads / 64 GB / 6 h, STAR index 8 / 64 GB, featureCounts 16 / 16 GB, Salmon 10 / 24 GB, rMATS 4 / 32 GB, PSI inclusion 1 / 64 GB / 12 h. Memory retries can scale with `attempt`. `profiles/local/` for testing on one machine.

Rerun behaviour: Snakemake 9 reruns on changed code, params, input or software environment by default; parameters will be passed via `params:` so changes trigger exactly the affected rules.

## 9. Testing (E)

1. **Static:** `snakemake --lint`; pytest tests for `validate_inputs.py` (good and deliberately broken configs/tables).
2. **Dry-run matrix:** `snakemake -n` for PE and SE examples, all modules on, rMATS/PSI off, prebuilt vs to-be-built indices.
3. **Integration:** a tiny dataset (one mouse chromosome region + matching GTF slice, ~100k read pairs per sample, 2 groups × 2) built by `tests/make_test_data.sh` from your own FASTQs, including index building. Runs locally in minutes.
4. **Regression:** run the new workflow with `legacy: true` on a subset of an already processed dataset and compare to your existing outputs: STAR `Log.final.out` metrics, featureCounts counts, `quant.sf`, rMATS event tables, `.psi` files. Then switch the corrected parameters on and report the differences separately ("expected from fixes" vs "unexpected").

I can run 1–2 here. 3–4 need tool binaries and your data; where I cannot run something I will document the exact commands instead of claiming it passed.

## 10. Migration (E)

- All new work on a branch; old files moved to `legacy/` (not deleted), downstream R scripts to `downstream/`.
- New runs write to a new `outdir`; existing results on the cluster are never touched.
- A conversion helper turns your `meta.xlsx`/`samples.xlsx` into `meta.tsv`/`samples.tsv` once.

## 11. Implementation order after approval

1. Config/metadata validation + example files → 2. driver env + profiles → 3. FastQC/fastp → 4. MultiQC → 5. STAR index + alignment → 6. featureCounts → 7. Salmon index + quant → 8. rMATS → 9. PSI → 10. final report + provenance → 11. README.
Dry-run + available tests after every module; one PR, commits per module.

## 12. Decisions I need from you

| # | Question | My default if you agree |
|---|---|---|
| D1 | fastp `-D` (deduplication): keep? | Off by default; `legacy` preset keeps it |
| D2 | featureCounts: stranded `-s` from config, `--countReadPairs`, no `-O`? | Yes; `legacy` preset keeps old behaviour |
| D3 | Remove SPLASH2 from the core workflow? | Yes (out of scope) |
| D4 | Your Python 3 PSI scripts: where are they? Which `dexseq_prepare_annotation.py` options made the reduced GFF (`-r yes/no`)? | Need the files |
| D5 | Cluster facts in `environment_strategy.md` §7 (Apptainer, PyPI/registry access, account/partition) | Primary route assumes Apptainer works |
| D6 | Your message ended at "15" – was there a section 15? | — |
