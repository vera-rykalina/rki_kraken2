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
// $ nextflow run ../scratch/rki_kraken2/Pipeline/Scripts/kraken2.nf \
// -profile rki_slurm,rki_mamba \
// -with-report logs/test.$(date +%T).report.html \
// # use this parameter for an empty test run
// -stub

/************************** 
DEFINE VARIABLES
**************************/

projectDir = "/home/rykalinav/scratch/rki_kraken2/Pipeline"


// Parameters for kraken2
params.krakendb = "/scratch/databases/kraken2_20230314/"

/************************** 
---------WORKFLOW----------
**************************/
ch_infiles = channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz")

workflow {

    ch_classified = CLASSIFY(ch_infiles, params.krakendb)
    ch_extracted = EXTRACT(ch_infiles, ch_classified.kraken2_output)
    ch_compressed = COMPRESS(ch_extracted.extracted_fastq)

}

/************************** 
PROCESSES
**************************/

// kraken2
process CLASSIFY {

    label "kraken2"
    conda "${projectDir}/Environments/kraken2.yml"
    publishDir "${params.outdir}/01_classification/${id}", mode: "copy", overwrite: true
    tag "${id}_classify"

    input:
        tuple val(id), path(reads)
        val (db)

    output:
        tuple val(id), path("${id}_classified_R*.fastq"),     emit: classified_fastq
        tuple val(id), path("${id}_unclassified_R*.fastq"),   emit: unclassified_fastq
        tuple val(id), path("${id}_kraken2_out.txt"),         emit: kraken2_output
        tuple val(id), path("${id}_kraken2_report.txt"),      emit: kraken2_report

    script:
        """
            kraken2 \
                --threads 10 \
                --db ${params.krakendb} \
                --paired \
                --classified-out ${id}_classified_R#.fastq \
                --unclassified-out ${id}_unclassified_R#.fastq \
                --output ${id}_kraken2_out.txt \
                --report ${id}_kraken2_report.txt \
                --report-minimizer-data \
                ${reads[0]} ${reads[1]}
        """
     
}

// krakenTools
process EXTRACT {
    label "krakentools"
    conda "${projectDir}/Environments/krakentools.yml"
    tag "${id}_extract"

    input:
        tuple val(id), path(reads)
        tuple val(id), path(kraken2_output)
    
    output:
        tuple val(id), path("${id}_extracted_R*.fastq"), emit: extracted_fastq
    
    script:
        """
            extract_kraken_reads.py \
                -k ${kraken2_output} \
                --taxid 9606 \
                --exclude \
                --fastq-output \
                -s1 ${reads[0]} \
                -s2 ${reads[1]} \
                -o ${id}_extracted_R1.fastq \
                -o2 ${id}_extracted_R2.fastq       
        """
}

// compress extracted output files
process COMPRESS {
    label "compress"
    publishDir "${params.outdir}/02_kraken2_extracted", failOnError: true, mode: "copy", overwrite: true
    tag "${id}_compress"

    input:
        tuple val(id), path(reads)

    output:
        tuple val(id), path("${id}_filtered_R*.fastq.gz"), emit: fastq_gz

    script:
        """
            gzip -c \$(realpath "${id}_extracted_R1.fastq") > ${id}_filtered_R1.fastq.gz
            gzip -c \$(realpath "${id}_extracted_R2.fastq") > ${id}_filtered_R2.fastq.gz
        """
}

