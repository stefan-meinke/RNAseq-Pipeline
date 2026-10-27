# Software environment strategy

Assessment date: 2026-10-27.

## 1. What was actually checked

- **conda-forge and Bioconda**: downloaded the full `linux-64` and `noarch` `repodata.json` of both channels from `conda.anaconda.org` and searched them by package name.
- **PyPI**: queried `pypi.org/pypi/<package>/json`.
- **Guix**: fetched the three channels at the exact commits pinned in your `channels.scm` and searched the package definitions for every entry in `manifest.scm`.
- **Upstream releases**: listed Git tags of STAR, Salmon, fastp and rMATS-turbo.

Not checked (no access from here; must be answered on the cluster, see §7): whether Apptainer/Singularity runs on compute nodes, whether PyPI, quay.io or GitHub are reachable from the login node, and which Snakemake module is installed.

## 2. Package availability (verified)

| Tool | conda-forge | Bioconda | PyPI | Other official source |
|---|---|---|---|---|
| Snakemake | **absent** | 9.27.1 | **9.27.1** (Python ≥ 3.11) | — |
| snakemake-executor-plugin-slurm | **absent** | 2.8.0 | **2.8.0** | — |
| FastQC | **absent** | 0.13.0 | — | Babraham zip (Java; `openjdk` **is** on conda-forge) |
| fastp | **absent** | 1.4.0 | — | Static Linux binary from OpenGene |
| MultiQC | **absent** | 1.35 | **1.35** | — |
| STAR | **absent** | yes | — | Static Linux binaries in the GitHub repo/releases |
| Subread (featureCounts) | **absent** | 2.1.1 | — | Linux binary tarball (SourceForge) |
| Salmon | **absent** | 2.8.0 (1.10.x also available) | — | Linux binary tarball on GitHub releases |
| rMATS-turbo | **absent** | `rmats` 4.4.0 | — | Source tarball on GitHub (v4.3.0, v4.4.0); build deps are all on conda-forge (`gsl`, `openblas`, `cython`, `gfortran`, `numpy`) |
| samtools | **absent** | 1.24 | — | Not needed: `pysam` bundles `samtools index` |
| HTSeq (PSI scripts) | **absent** | 2.1.2 | **2.1.2** (wheels) | — |
| pysam | **absent** | 0.24.1 | **0.24.1** (wheels) | — |
| bedtools | **absent** | 2.31.1 | — | Not needed by the core workflow |
| DEXSeq (`dexseq_prepare_annotation.py`) | **absent** | 1.56.0 | — | The script needs only HTSeq; can be vendored (GPL-3/LGPL) |
| Python, R, pandas, jsonschema, openjdk, pigz, apptainer | **present** | — | — | — |

**Conclusion: a conda-forge-only environment is not possible.** None of the bioinformatics tools, and not even Snakemake itself, are on conda-forge. Any environment spec listing them under `conda-forge` would fail to solve. I therefore will not write such a spec.

Consequences for Snakemake's per-rule `conda:` environments: they would need Bioconda for almost every rule, so they are not an option here. Snakemake can instead use per-rule `container:` images (`--software-deployment-method apptainer`) or plain executables on `PATH`.

Note on Salmon: Bioconda now carries Salmon 2.x. To preserve your validated behaviour, the plan pins **Salmon 1.10.x**; moving to 2.x would be a separate, documented decision.

## 3. Option A – conda-forge + PyPI + official binaries

- Driver environment from conda-forge (Python 3.12, pandas, openjdk, pigz) with Snakemake, the SLURM plugin, MultiQC, HTSeq and pysam installed from PyPI with pinned versions.
- STAR, Salmon, fastp, Subread, FastQC from official upstream binaries, downloaded once by an install script into a project-local `tools/` folder, each pinned by version **and SHA-256 checksum**.
- rMATS-turbo compiled from the v4.3.0 tarball against conda-forge `gsl`/`openblas`/`gfortran`.

| Pro | Con |
|---|---|
| Uses only allowed sources (if PyPI and GitHub are reachable). Snakemake 9 runs independently of the cluster module. | Binaries are pinned by checksum, not built from source; less hermetic than Guix or containers (relies on host glibc for dynamically linked tools). |
| Simple to understand; works without container permissions. | Install script to maintain; rMATS compile can fail on unusual systems. |

## 4. Option B – GNU Guix

**Why it currently fails (confirmed against the pinned channels):**

1. `manifest.scm` asks for `r-apeglm@1.34.0`; the pinned channels only provide 1.28.0 → `guix shell -m manifest.scm` aborts with "package not found".
2. `run_pipeline_guix_time-machine.sh` only runs `snakemake --version` inside Guix. The real `snakemake …` command is a separate line and runs with the host's (old) Snakemake.
3. rMATS-turbo is in none of the pinned channels; `config.yml` points to a `/gnu/store/...rmats-turbo-4.3.0` path from some other profile, which disappears after `guix gc`.
4. Potential: `channels.scm` uses `git.savannah.gnu.org`. Guix moved to `codeberg.org/guix/guix.git` in 2025 and only guaranteed the Savannah mirror for at least one year. If that mirror is gone or stale on your cluster, `guix time-machine` cannot fetch the channel. To check: `git ls-remote https://git.savannah.gnu.org/git/guix.git | head -1`.

All other manifest entries exist at the pinned commits (Snakemake 7.7.0, STAR 2.7.8a, Salmon 1.10.1, Subread 2.0.3, fastp 0.23.2, MultiQC 1.14, HTSeq 2.0.2, pysam 0.20.0, Python 3.10.7).

| Pro | Con |
|---|---|
| Strongest reproducibility (bit-for-bit, source-built). Your past results came from this stack. | Snakemake in the pinned Guix is 7.7.0 (2022): no executor plugins, deprecated `--cluster`. A newer Guix commit would change *all* tool versions at once. |
| | Needs `guix-daemon` and the store on every compute node; `time-machine` is slow on first use; rMATS-turbo needs a custom package definition. |

**Assessment:** Worth keeping as an optional way to provide the tools (one corrected manifest with only the core tools), but not as the primary route; it is the part that broke and it pins Snakemake to an old version.

## 5. Option C – Apptainer/Singularity containers

- One pinned image per tool, referenced **by digest** (`docker://quay.io/biocontainers/star@sha256:…`). BioContainers images are built from Bioconda recipes, but they are pulled from a container registry, not installed from the Bioconda channel. Whether that is allowed depends on *why* Bioconda is blocked on your cluster (network block vs. licence policy) – please check.
- Images are cached once as `.sif` files in a shared directory (`--apptainer-prefix`); Snakemake binds the project, reference and scratch paths (`--apptainer-args "--bind /fast"`).
- Works with SLURM: each job runs `apptainer exec` on the node; no daemon.

| Pro | Con |
|---|---|
| Exact pinned versions for every tool, including rMATS-turbo, with no compiling. Digest pinning is the strongest guarantee short of Guix. | Requires Apptainer on compute nodes and registry access from the login node. |
| Same images work on any cluster/laptop. | ~0.2–1 GB per image; bind-mount configuration for `/fast` and `$TMPDIR`. |

## 6. Option D – Hybrid (recommended)

```
Snakemake 9 + SLURM plugin + MultiQC + HTSeq/pysam   ← conda-forge env + pinned PyPI packages (lockfile)
        │
        ├── per-rule tools ← (primary)  Apptainer images pinned by digest
        │                    (fallback) pinned upstream binaries in tools/ (checksummed), or Guix
        └── rMATS-turbo, PSI scripts ← same as above, or an explicitly configured local executable/script
```

Rules call plain executable names (`STAR`, `salmon`, …). Each rule also declares `container: config["containers"][...]`. The same Snakefile then runs in either mode without edits:

- `--software-deployment-method apptainer` → tools come from the pinned images;
- without it → tools come from `PATH` (the `tools/` folder, Guix, or modules).

Every run writes `results/provenance/`: `tool_versions.tsv` (version string of the executable actually called, path, container digest or "local executable – not controlled"), `config.snapshot.yml`, the three input tables, the pipeline Git commit, `pip freeze`/`conda list --explicit` of the driver env, and checksums of FASTA/GTF/flattened GFF.

### Recommendation

- **Primary:** Option D with Apptainer images pinned by digest for the tools, driver env from conda-forge + PyPI. Most reproducible option that respects the Bioconda restriction, and the simplest to maintain (no compiling, no Guix daemon).
- **Fallback (if containers are not allowed):** the same workflow with tools from a pinned, checksummed `tools/` folder (Option A), plus a configured local rMATS-turbo if compiling fails. The report will then mark these tools as "less controlled".
- **rMATS-turbo:** fully reproducible via a container (pinned digest) in the primary route. In the fallback route it is compiled from the pinned v4.3.0 tarball; the configurable local executable stays as a last resort and is recorded as such.
- **Not recommended:** the conda-forge-only option (impossible), per-rule `conda:` environments (need Bioconda), and Guix as the primary route.

## 7. Open questions for the cluster

1. Does `apptainer --version` (or `singularity --version`) work on a **compute node**, and can the login node pull from `quay.io` / `docker.io`?
2. Can the login node reach `pypi.org` and `github.com`?
3. Is Bioconda blocked by the network or by policy? (Decides whether BioContainers images are acceptable.)
4. Output of `module avail snakemake` / `snakemake --version` on the cluster, plus SLURM account/partition names.
5. Does `git ls-remote https://git.savannah.gnu.org/git/guix.git` still work from the cluster?
