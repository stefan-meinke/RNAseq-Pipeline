import os 
import re

configfile: "config.yml"

def get_fastq(sample, read):
    # Get the base file name from the config mapping.
    base = config["samples"][sample]
    # Build two possible file names using the base.
    candidate1 = f'{config["sample_path"]}/{base}_R{read}.fastq.gz'
    candidate2 = f'{config["sample_path"]}/{base}_{read}.fastq.gz'
    if os.path.exists(candidate1):
        return candidate1
    elif os.path.exists(candidate2):
        return candidate2
    else:
        raise FileNotFoundError(f"Could not find file for sample {sample} (base: {base}) read {read}")


def get_group(sample):
    """
    Extracts a group name from a sample.
    For example, for "WT1" or "WT_1" it returns "WT"; for "KO1" it returns "KO".
    Adjust the regex if your naming convention is different.
    """
    # split by underscore if available:
    if "_" in sample:
        return sample.split("_")[0]
    # Otherwise, assume that the group is all non-digit characters at the start:
    m = re.match(r"(\D+)", sample)
    if m:
        return m.group(1)
    else:
        raise ValueError(f"Could not extract group from sample: {sample}")


SAMPLES = list(config["samples"].keys())

rule all:
    input:
        # fastp outputs:
        expand("{out_dir}/{sample}_R1_trimmed.fastq.gz", out_dir=config["outdir_fastp"], sample=SAMPLES),
        expand("{out_dir}/{sample}_R2_trimmed.fastq.gz", out_dir=config["outdir_fastp"], sample=SAMPLES),
        # fastqc outputs:
        expand("{qc_dir}/{sample}_R{read}_trimmed_fastqc.zip", qc_dir=config["outdir_fastqc"], sample=SAMPLES, read=[1,2]),
	# multiqc output (see next change):
        f'{config["outdir_fastqc"]}/multiqc_report.html',
        # STAR alignment outputs:
        expand("{star_dir}/{sample}_Aligned.sortedByCoord.out.bam", star_dir=config["outdir_star"], sample=SAMPLES),
        # samtools index outputs:
        expand("{star_dir}/{sample}_Aligned.sortedByCoord.out.bam.bai", star_dir=config["outdir_star"], sample=SAMPLES),
	# rmats outputs:
	expand("rmats/{contrast}_b1.txt", contrast=list(config["contrasts"].keys())),
        expand("rmats/{contrast}_b2.txt", contrast=list(config["contrasts"].keys())),
        expand("rmats/{contrast}", contrast=list(config["contrasts"].keys())),
	#"Results/rmats/rmats_results.xlsx",
	# featureCounts output:
	expand("{fc_dir}/{sample}_featureCounts.txt", fc_dir=config["outdir_featureCounts"], sample=SAMPLES),
	# DESeq2 output:
	"Results/DESeq2/DESeq_results.xlsx",
	# geneExpression analysis
	"Results/DESeq2/gene_expression_done.marker",
	# Enrichment analysis:
	"Results/DESeq2/Enrichment/enrichment_done.marker",
	# salmon output:
	expand("{salmon_dir}/{sample}_quant/quant.sf", salmon_dir=config["outdir_salmon"], sample=SAMPLES),
	# PSI output:
	expand("{psi_dir}/{sample}.psi", psi_dir=config["outdir_psi"], sample=SAMPLES),
        expand("{psi_dir}/{sample}.inclusion", psi_dir=config["outdir_psi"], sample=SAMPLES),
        expand("{psi_dir}/{sample}.exclusion", psi_dir=config["outdir_psi"], sample=SAMPLES),
	# Include STAR and Salmon index outputs:
        config["star_index_dir"] + "/SA",
        config["salmon_index_dir"] + "/salmon_index.built",
	# SPLASH2 files:
	"Splash2/input.txt",
        "Splash2/splash_done.marker",
	# PSI plots:
	"Results/PSI_plots/psi_plots_done.marker"



rule fastp:
    input:
        r1 = lambda wc: get_fastq(wc.sample, 1),
        r2 = lambda wc: get_fastq(wc.sample, 2)
    output:
        r1_trim = f'{config["outdir_fastp"]}/{{sample}}_R1_trimmed.fastq.gz',
        r2_trim = f'{config["outdir_fastp"]}/{{sample}}_R2_trimmed.fastq.gz',
        failed  = f'{config["outdir_fastp"]}/{{sample}}_failed_reads.txt'
    shell:
        """
        fastp -i {input.r1} \
              -I {input.r2} \
              -o {output.r1_trim} \
              -O {output.r2_trim} \
              -D -c \
              --failed_out {output.failed}
        """

rule fastqc:
    input:
        trimmed = f'{config["outdir_fastp"]}/{{sample}}_R{{read}}_trimmed.fastq.gz'
    output:
        qc_zip = f'{config["outdir_fastqc"]}/{{sample}}_R{{read}}_trimmed_fastqc.zip'
    shell:
        """
        fastqc --noextract -o {config[outdir_fastqc]} {input.trimmed}
        """

rule multiqc:
    input:
        expand("{qc_dir}/{sample}_R{read}_trimmed_fastqc.zip",
               qc_dir=config["outdir_fastqc"],
               sample=SAMPLES,
               read=[1,2])
    output:
        report = f'{config["outdir_fastqc"]}/multiqc_report.html'
    shell:
        """
        multiqc {config[outdir_fastqc]} -o {config[outdir_fastqc]}
        """


# rule to check for the STAR index output (e.g., STAR creates a file called SA)
if config["build_star_index"]:
    rule star_index:
        output:
            star_sa = config["star_index_dir"] + "/SA"
        params:
            genome_fasta = config["genome_fasta"],
            gtf = config["gtf"],
            sjdbOverhang = config["sjdbOverhang"],
            index_dir = config["star_index_dir"]
        threads: 8
	resources:
            mem_mb=64000
        shell:
            """
            STAR --runThreadN {threads} --runMode genomeGenerate \
                --genomeDir {params.index_dir} \
                --genomeFastaFiles {params.genome_fasta} \
                --sjdbGTFfile {params.gtf} \
                --sjdbOverhang {params.sjdbOverhang}
            """
else:
    # If pre-built index is present, create a dummy rule
    rule star_index:
        output:
            star_sa = config["star_index_dir"] + "/SA"
        shell:
            "echo 'Using pre-built STAR index' && touch {output.star_sa}"



rule star:
    input:
        r1 = f'{config["outdir_fastp"]}/{{sample}}_R1_trimmed.fastq.gz',
        r2 = f'{config["outdir_fastp"]}/{{sample}}_R2_trimmed.fastq.gz',
	star_index = config["star_index_dir"] + "/SA"
    output:
        aligned_bam = f'{config["outdir_star"]}/{{sample}}_Aligned.sortedByCoord.out.bam',
        gene_counts = f'{config["outdir_star"]}/{{sample}}_ReadsPerGene.out.tab',
        sj = f'{config["outdir_star"]}/{{sample}}_SJ.out.tab'
    params:
        star_index_dir = config["star_index_dir"],
        out_prefix = lambda wc: f'{config["outdir_star"]}/{wc.sample}_',
        twopass = "--twopassMode Basic" if config.get("twopass", False) else ""
    threads: 16
    resources:
        mem_mb=64000
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {params.star_index_dir} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --seedPerWindowNmax 15 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             {params.twopass}
        """

rule samtools_index:
    input:
        aligned_bam = f'{config["outdir_star"]}/{{sample}}_Aligned.sortedByCoord.out.bam'
    output:
        bai = f'{config["outdir_star"]}/{{sample}}_Aligned.sortedByCoord.out.bam.bai'
    shell:
        """
        samtools index {input.aligned_bam}
        """

rule rmats:
    input:
        # Assume all STAR BAM files have been created.
        bams = expand("{star_dir}/{sample}_Aligned.sortedByCoord.out.bam", star_dir=config["outdir_star"], sample=SAMPLES)
    output:
        b1 = "rmats/{contrast}_b1.txt",
        b2 = "rmats/{contrast}_b2.txt",
        outdir = directory("rmats/{contrast}")
    params:
        rmats_path = config["rmats_path"],
        gtf = config["gtf"],
        contrast = lambda wc: wc.contrast
    threads: 4
    run:
        # Build a dictionary mapping sample -> group using the helper function.
        groups = {sample: get_group(sample) for sample in SAMPLES}
        
        # For the given contrast (e.g., "KO-vs-WT"), get the two group names.
        # config["contrasts"] should be defined like: { "KO-vs-WT": ["KO", "WT"] }
        group1, group2 = config["contrasts"][wildcards.contrast]
        
        # Filter the BAM files based on group membership.
        b1_bams = [f'{config["outdir_star"]}/{sample}_Aligned.sortedByCoord.out.bam'
                   for sample in SAMPLES if groups[sample] == group1]
        b2_bams = [f'{config["outdir_star"]}/{sample}_Aligned.sortedByCoord.out.bam'
                   for sample in SAMPLES if groups[sample] == group2]
        
        # Write out the lists to the text files.
        with open(output.b1, 'w') as f:
            f.write(",".join(b1_bams))
        with open(output.b2, 'w') as f:
            f.write(",".join(b2_bams))
        
        # run rMATS-turbo.
        shell("""
            python3 {params.rmats_path} \
                --b1 {output.b1} \
                --b2 {output.b2} \
                --gtf {params.gtf} \
                --od {output.outdir} \
                --tmp {output.outdir} \
                -t {config[rmats_type]} \
                --nthread {threads} \
                --readLength {config[readLength]} \
                --variable-read-length \
                --libType {config[libType]}
        """)


rule rmats_analysis:
    input:
        # Wait for all rMATS contrast directories to be ready.
        contrasts = expand("rmats/{contrast}", contrast=list(config["contrasts"].keys()))
    output:
        excel = "Results/rmats/rmats_results.xlsx"
    params:
        rmats_dir = "rmats"
    threads: 1
    shell:
        """
        Rscript R_scripts/rmats_analysis.R \
            --rmats_dir {params.rmats_dir} \
            --output {output.excel}
        """


rule featureCounts:
    input:
        bam = f'{config["outdir_star"]}/{{sample}}_Aligned.sortedByCoord.out.bam'
    output:
        fc = f'{config["outdir_featureCounts"]}/{{sample}}_featureCounts.txt'
    params:
        gtf = config["gtf"],
        paired = "-p" if config.get("feature_type", False) else ""
    threads: 16
    shell:
        """
        mkdir -p {config[outdir_featureCounts]}
        featureCounts {params.paired} -O -t exon -g gene_id -T {threads} -a {params.gtf} -o {output.fc} {input.bam}
        """


# Check for salmon index
if config["build_salmon_index"]:
    rule salmon_index:
        output:
            # Create a marker file to indicate the index was built
            marker = config["salmon_index_dir"] + "/salmon_index.built"
        params:
            transcriptome = config["transcriptome_fasta"],
            index_dir = config["salmon_index_dir"]
        shell:
            """
            salmon index -t {params.transcriptome} -i {params.index_dir} && touch {output.marker}
            """
else:
    rule salmon_index:
        output:
            marker = config["salmon_index_dir"] + "/salmon_index.built"
        shell:
            "echo 'Using pre-built Salmon index' && touch {output.marker}"




rule salmon:
    input:
        r1 = f'{config["outdir_fastp"]}/{{sample}}_R1_trimmed.fastq.gz',
        r2 = f'{config["outdir_fastp"]}/{{sample}}_R2_trimmed.fastq.gz',
	salmon_index = config["salmon_index_dir"] + "/salmon_index.built"
    output:
        quant_sf = f'{config["outdir_salmon"]}/{{sample}}_quant/quant.sf'
    params:
        salmon_index = config["salmon_index_dir"],
        lib = "A"  # Library type; adjust if needed.
    threads: 10
    shell:
        """
        salmon quant -i {params.salmon_index} -l {params.lib} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            --validateMappings \
            -o {config[outdir_salmon]}/{wildcards.sample}_quant
        """


rule sj_bed:
    input:
        sj_tab = f'{config["outdir_star"]}/{{sample}}_SJ.out.tab'
    output:
        bed = f'{config["outdir_psi"]}/{{sample}}_splice_junctions.bed'
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}{{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+" : "-", $2-20-1, $3+20, "255,0,0", 2, "20,20", "0,300"}}' {input.sj_tab} > {output.bed}
        """


rule dexseq_count:
    input:
        bam = f'{config["outdir_star"]}/{{sample}}_Aligned.sortedByCoord.out.bam'
    output:
        inclusion = f'{config["outdir_psi"]}/{{sample}}.inclusion'
    params:
        paired = "yes" if str(config["psi"]["paired"]).lower() in ["yes", "true"] else "no",
        strand = config["psi"]["strand_specificity"]  # e.g. "reverse"
    shell:
        """
        python3 PSI_scripts/dexseq_count.py -p {params.paired} -s {params.strand} -r pos -f bam {config[gff]} {input.bam} {output.inclusion}
        """


rule exclusion_count:
    input:
        bed = rules.sj_bed.output.bed
    output:
        exclusion = f'{config["outdir_psi"]}/{{sample}}.exclusion'
    shell:
        """
        python3 PSI_scripts/exclusion_count.py {config[gff]} {input.bed} {output.exclusion}
        """


rule psi_calculation:
    input:
        inclusion = f'{config["outdir_psi"]}/{{sample}}.inclusion',
        exclusion = f'{config["outdir_psi"]}/{{sample}}.exclusion'
    output:
        psi = f'{config["outdir_psi"]}/{{sample}}.psi'
    params:
        readLength = config["readLength"],
        gff = config["gff"]
    shell:
        """
        if [ ! -f {output.psi} ]; then
            python3 PSI_scripts/psi_calculation.py {params.gff} {input.inclusion} {input.exclusion} -l {params.readLength} {output.psi}
        else
            echo "{output.psi} already exists. Skipping psi_calculation."
        fi
        """

rule psi_plots:
    resources:
        mem_mb=4000
    input:
        psi_files = expand("{psi_dir}/{sample}.psi", psi_dir=config["outdir_psi"], sample=SAMPLES)
    output:
        marker = "Results/PSI_plots/psi_plots_done.marker"
    log:
        "log/PSI_plots/psi_plots.log"
    run:
        if not config.get("generate_psi_plots", False):
            print("Skipping PSI plot generation because generate_psi_plots is set to False")
            shell("mkdir -p Results/PSI_plots && touch {output.marker}")
        else:
            # Run the PSI plot R script and redirect its output to the log file
            shell("""
                Rscript R_scripts/PSIplot.R \
                    --psi_dir {config[outdir_psi]} \
                    --colData {config[colData]} \
                    --gtf {config[psi_gtf]} \
                    --target_genes {config[target_genes]} \
                    --output Results/PSI_plots > {log} 2>&1
            """)
            shell("touch {output.marker}")


rule dge:
    input:
        counts_files = expand("{fc_dir}/{sample}_featureCounts.txt",
                              fc_dir=config["outdir_featureCounts"],
                              sample=SAMPLES),
        colData = config["colData"],
        gtf = config["gtf"]
    output:
        de_results = "Results/DESeq2/DESeq_results.xlsx"
    params:
        script = "R_scripts/dge_analysis.R"
    run:
        shell("mkdir -p Results/DESeq2")
        shell("Rscript {params.script} --counts_dir {config[outdir_featureCounts]} --colData {input.colData} --gtf {input.gtf} --output {output.de_results}")


rule gene_expression_analysis:
    input:
        counts_dir = config["outdir_featureCounts"],
        colData = config["colData"],
        gtf = config["gtf"],
        config_file = "config.yml"
    output:
        marker = "Results/DESeq2/gene_expression_done.marker"
    params:
        script = "R_scripts/geneExpression_analysis.R",
        deseq_dir = "Results/DESeq2"  # future output dir
    shell:
        """
        mkdir -p {params.deseq_dir}
        Rscript {params.script} \
            --deseq_result {params.deseq_dir} \
            --counts_dir {input.counts_dir} \
            --colData {input.colData} \
            --gtf {input.gtf} \
            --config {input.config_file}
        touch {output.marker}
        """


rule enrichment:
    input:
        deseq_result="Results/DESeq2/DESeq_results.xlsx",
        config_file="config.yml"
    output:
        marker="Results/DESeq2/Enrichment/enrichment_done.marker"
    params:
        script="R_scripts/enrichment_analysis.R",
        deseq_dir="Results/DESeq2",
        annotation=config["organism"],
        output_dir="Results/DESeq2/Enrichment"
    shell:
        """
        mkdir -p {params.output_dir}
        Rscript {params.script} \
            --deseq_result {params.deseq_dir} \
            --annotation {params.annotation} \
            --config {input.config_file} \
            --output {params.output_dir}
        touch {output.marker}
        """


# make input.txt file for SPLASH2
rule make_input_txt:
    input:
        fastp_R1 = expand("{out_dir}/{sample}_R1_trimmed.fastq.gz",
                          out_dir=config["outdir_fastp"],
                          sample=SAMPLES)
    output:
        "Splash2/input.txt"
    run:
        import os
        # Create the directory "Splash2" if it doesn't exist
        os.makedirs("Splash2", exist_ok=True)
        outdir = config["outdir_fastp"]
        sample_names = list(config["samples"].keys())
        with open(output[0], "w") as f:
            for sample in sample_names:
                file_path = f"{outdir}/{sample}_R1_trimmed.fastq.gz"
                f.write(f"{sample}\t{file_path}\n")


rule splash2:
    input:
        "Splash2/input.txt"
    output:
        "Splash2/splash_done.marker"
    params:
        splash_executable = config["splash"],
        splash_dir = os.path.abspath("Splash2")  # absolute path to the Splash2 folder
    threads: 4
    shell:
        """
        cd {params.splash_dir} && {params.splash_executable} ../Splash2/input.txt \
            --dump_sample_anchor_target_count_binary \
            --keep_significant_anchors_satc \
            --n_most_freq_targets 4 \
            --bin_path . \
            --n_threads_stage_2 {threads}
        cd .. && touch {output}
        """


