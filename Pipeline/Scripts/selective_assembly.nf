/************************** 
* PARAMETERS
**************************/

// if ( !params.fastq ) {
//     exit 1, "input missing, use [--fastq]"
// }

/************************** 
* CALL
**************************/

// (/scratch/drechselo/nextflow_env) hpc-login02[/scratch/.../selectiveAssembly]
// $ nextflow run ../script_development/selectiveAssembly/selective_assembly.nf \
// --config /scratch/Tausch/drechselo/read_filter_pipeline/samplesheet.csv \
// --reference_gff /scratch/...quastRef/GCF_000006845.1_ASM684v1_genomic.gff \
// --reference_fasta /scratch/.../quastRef/GCF_000006845.1_ASM684v1_genomic.fna \
// --krakendb /scratch/databases/kraken2_20230314/ \
// --taxid 482 \ ### modify here
// --analysesdir /scratch/.../analyses/ \
// -profile rki_slurm,rki_mamba \
// -with-report logs/test.$(date +%T).report.html \
// # use this parameter for an empty test run
// -stub


/************************** 
* DEFINE VARIABLES
**************************/


inputCSV      = params.config
reference_gff   = file(params.reference_gff)
reference_fasta = file(params.reference_fasta)
krakendb      = params.krakendb
taxid         = params.taxid
analysesdir   = params.analysesdir

/************************** 
* DEFINE INCLUDES
**************************/
// necessary, because NextFlow does not allow to use processes twice, if not stored as a module

include { assemble_spades as assemble_spades_filtered } from "./modules/spades.nf"
include { assemble_spades as assemble_spades_unfiltered } from "./modules/spades.nf"

include { rename_spades as rename_spades_filtered } from "./modules/rename_spades.nf"
include { rename_spades as rename_spades_unfiltered } from "./modules/rename_spades.nf"

include { compress_spades as compress_spades_filtered } from "./modules/compress_spades.nf"
include { compress_spades as compress_spades_unfiltered } from "./modules/compress_spades.nf"

include { qc_spades as qc_spades_filtered } from "./modules/qc_spades.nf"
include { qc_spades as qc_spades_unfiltered } from "./modules/qc_spades.nf"

/************************** 
* WORKFLOW
**************************/

// ch_infiles = Channel.fromPath(inputFiles)

// read in csv file (sep=; due to Excel)
// sample;raw_fwd;raw_rev;trimmed_fwd;trimmed_rev;raw_folder;trimmed_folder

ch_infiles = Channel
                .fromPath(inputCSV)
                .splitCsv(header: true, sep: ';')
                .map {row -> tuple(
                                    row.sample, 
                                    file(row.raw_folder + "/" + row.raw_fwd),
                                    file(row.raw_folder + "/" + row.raw_rev),
                                    // file(row.trimmed_folder + "/" + row.trimmed_fwd),
                                    // file(row.trimmed_folder + "/" + row.trimmed_rev),
                                    )}

workflow {

    // ch_infiles.view()

    ch_trimmed = trim(ch_infiles) // run trimming

    // DEBUG
    // ch_trimmed = ch_infiles 
    // END DEBUG

    
    // **********************************
    // Workflow of species filtered reads
    // **********************************
    
    ch_classified = classify(ch_trimmed.fastq, krakendb)

    ch_filtered_reads = filter_reads(ch_classified.fastq, ch_classified.kraken_output, ch_classified.kraken_report, taxid)
    ch_compressed_reads = compress_reads(ch_filtered_reads.fastq)
    ch_assembled_reads = assemble_spades_filtered(ch_compressed_reads.fastq, "filtered")
    ch_renamed_assembly = rename_spades_filtered(ch_assembled_reads, "filtered")
    ch_compressed_assembly = compress_spades_filtered(ch_renamed_assembly, "filtered")
    ch_checked_assembly = qc_spades_filtered(ch_compressed_assembly.contigs, reference_fasta, reference_gff, "filtered")

    // *******************************
    // Workflow without species filter
    // *******************************
    ch_assembled_unfiltered_reads = assemble_spades_unfiltered(ch_trimmed.fastq, "unfiltered")
    ch_renamed_unfiltered_assembly = rename_spades_unfiltered(ch_assembled_unfiltered_reads, "unfiltered")
    ch_compressed_unfiltered_assembly = compress_spades_unfiltered(ch_renamed_unfiltered_assembly, "unfiltered")
    ch_checked_unfiltered_assembly = qc_spades_unfiltered(ch_compressed_unfiltered_assembly.contigs, reference_fasta, reference_gff, "unfiltered")

    // *******************************
    // generate report
    // *******************************
    ch_merged_files = ch_checked_assembly.report_tsv.concat(ch_checked_unfiltered_assembly.report_tsv).collect{ it -> it[1] }
    ch_compile_report = multiqc( ch_merged_files )
}




/************************** 
* PROCESSES
**************************/

process get_references {

    script:
        """
        """
    
    stub:
        """
        """
}

process trim {

    label "fastp"

    conda "${projectDir}/envs/fastp.yaml"

    publishDir "${params.analysesdir}/01_trimming/${sample}", pattern: "*.fastq.gz", mode: "copy"

    input:
        tuple val(sample), path(read1), path(read2)

    output:
        tuple val(sample), path("${sample}.R*.trimmed.fastq.gz"), emit: fastq
        tuple val(sample), path("${sample}.trim_report.json"), emit: json
        tuple val(sample), path("${sample}.trim_report.html"), emit: html

    // SLURM cluster options
    cpus 10
    memory "3 GB"
    time "1h"
    // clusterOptions "--job-name=classify_${sample}"
    tag "trim_${sample}"

    script:
        """
            fastp \
                --in1 ${read1} \
                --in2 ${read2} \
                --out1 ${sample}.R1.trimmed.fastq.gz \
                --out2 ${sample}.R2.trimmed.fastq.gz \
                --trim_poly_g \
                --json ${sample}.trim_report.json \
                --html ${sample}.trim_report.html \
                --thread ${task.cpus}
        """

    stub:
        """
            touch ${sample}.R{1,2}.trimmed.fastq.gz ${sample}.trim_report.json ${sample}.trim_report.html
        """

}

// Kraken2
process classify {

    label "kraken"

    conda "${projectDir}/envs/kraken2.yaml"

    publishDir "${params.analysesdir}/02_classification/${sample}", pattern: "*.txt"

    // SLURM cluster options
    cpus 10
    memory "150 GB"
    time "4h"
    // clusterOptions "--job-name=classify_${sample}"
    tag "classify_${sample}"

    input:
        // tuple val(sample), path(read1), path(read2) // <- if trimmed input from csv is used (trimming run with QCurchin)
        tuple val(sample), path(reads) // <- if fastp trimming process is used
        
        val krakendb

    output:
        tuple val(sample), path("${sample}.classified.R*.fastq"),     emit: fastq
        tuple val(sample), path("${sample}.kraken.out.txt"),          emit: kraken_output
        tuple val(sample), path("${sample}.kraken.report.txt"),       emit: kraken_report

    script:
        """
            kraken2 \
                --threads ${task.cpus} \
                --db ${krakendb} \
                --paired \
                --classified-out ${sample}.classified.R#.fastq \
                --output ${sample}.kraken.out.txt \
                --report ${sample}.kraken.report.txt \
                ${reads[0]} ${reads[1]}
        """

    stub:
        """
            touch ${sample}.classified.R_{1,2}.fastq ${sample}.kraken.out.txt ${sample}.kraken.report.txt
        """
}

// krakenTools
process filter_reads {

    label "krakentools"

    conda "${projectDir}/envs/krakentools.yaml"

    // SLURM cluster options
    cpus 1
    memory "1 GB"
    time "1h"

    tag "filter_${sample}"

    input:
        tuple val(sample), path(reads)
        tuple val(sample), path(kraken_output)
        tuple val(sample), path(kraken_report)
        val(taxid)

    output:
        tuple val(sample), file("${sample}.classified.R*.fastq"),     emit: fastq
    
    script:
        """
            extract_kraken_reads.py \
                -t ${taxid}\
                -k ${kraken_output} \
                -s1 ${reads[0]} \
                -s2 ${reads[1]} \
                -o ${sample}.classified.R1.fastq \
                -o2 ${sample}.classified.R2.fastq \
                --fastq-output \
                --report ${kraken_report} \
                --include-children

        """
    stub:
        """
            touch ${sample}.classified.R1.fastq ${sample}.classified.R2.fastq
        """
}

// compress kraken2 output files, selected files
process compress_reads {

    label "compress"

    publishDir path: "${params.analysesdir}/02_read_filter/${sample}", pattern: "*.fastq.gz", failOnError: true, mode: "copy"

    tag "compress_${sample}"

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}.classified.R*.fastq.gz"), emit: fastq

    script:
        """
            gzip -c \$(realpath ${sample}.classified.R1.fastq) > ${sample}.classified.R1.fastq.gz
            gzip -c \$(realpath ${sample}.classified.R2.fastq) > ${sample}.classified.R2.fastq.gz
        """

    stub:
        """
            touch ${sample}.classified.R1.fastq.gz ${sample}.classified.R2.fastq.gz
        """
}

// Assembly w/ Spades
// included from module


// rename Assembly
// included from module


// compress Assembly
// included from module


// Assembly QC w/ Quast
// included from module

process multiqc {
    label "multiqc"

    conda "${projectDir}/envs/multiqc.yaml"

    publishDir path: "${analysesdir}/06_delivery", pattern: "*html", mode: "copy"

    input:
        path multiqc_files , stageAs: "?/*"
    
    output:
        path "multiqc_report.html", emit: report

    script:
        """
            multiqc .
        """
    stub:
        """
            touch multiqc_report.html
        """
}

