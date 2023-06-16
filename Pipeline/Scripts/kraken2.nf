nextflow.enable.dsl = 2

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


/************************** 
-----------CALL------------
**************************/

// (/scratch/rykalinav/nextflow) hpc-login02[/scratch/.../kraken2]
// $ nextflow run ../scratch/rki_phyloTSI/Pipeline/Scripts/kraken2.nf \
// --krakendb /scratch/databases/kraken2_20230314/ \
// --taxid ### \ ### modify here
// -profile rki_slurm,rki_mamba \
// -with-report logs/test.$(date +%T).report.html \
// # use this parameter for an empty test run
// -stub

/************************** 
DEFINE VARIABLES
**************************/

projectDir = "/home/rykalinav/scratch/rki_phyloTSI/Pipeline"
krakendb      = params.krakendb
//taxid         = params.taxid


/************************** 
---------WORKFLOW----------
**************************/
ch_infiles = channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz")

workflow {

    ch_infiles.view()
    ch_classified = CLASSIFY(ch_infiles, krakendb)

}

/************************** 
PROCESSES
**************************/

// kraken2
process CLASSIFY {

    label "kraken2"
    conda "${projectDir}/Environments/kraken2.yml"
    publishDir "${params.outdir}/01_classification/${sample}", pattern: "*.txt"

    // SLURM cluster options
    cpus 10
    memory "150 GB"
    time "4h"
    // clusterOptions "--job-name=classify_${sample}"
    tag "${sample}_classify"

    input:
        tuple val(sample), path(reads)
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
    conda "${projectDir}/envs/krakentools.yml"

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

