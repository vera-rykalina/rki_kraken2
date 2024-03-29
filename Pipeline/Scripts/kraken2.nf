nextflow.enable.dsl = 2

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


/************************** 
-----------CALL------------
**************************/

// $ nextflow ../scratch/rki_kraken2/Pipeline/Scripts/kraken2.nf \
// -c Scripts/rki_config \
// -profile rki_slurm,rki_mamba \
// -with-report logs/test.$(date +%T).report.html 

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
    ch_extracted = EXTRACT(ch_classified.classified_fastq, ch_classified.kraken2_output)
    ch_merged_ags = ch_classified.unclassified_fastq.combine(ch_extracted, by:0)
    ch_merged_compressed = MERGE(ch_merged_ags)

}

/************************** 
PROCESSES
**************************/

// kraken2
process CLASSIFY {

    label "kraken2"
    conda "/home/rykalinav/.conda/envs/kraken2_v2.1.2"
    //conda "${projectDir}/Environments/kraken2.yml"
    publishDir "${params.outdir}/01_classified_reads/${id}", mode: "copy", overwrite: true

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
                ${reads[0]} ${reads[1]}
        """
     
}

// krakentools
process EXTRACT {
    label "krakentools"
    conda "/home/rykalinav/.conda/envs/krakentools_v1.2"
    //conda "${projectDir}/Environments/krakentools.yml"
    publishDir "${params.outdir}/02_homo_filtered_reads", mode: "copy", overwrite: true

    input:
        tuple val(id), path(reads)
        tuple val(id), path(kraken2_output)
    
    output:
        tuple val(id), path("${id}_filtered_R*.fastq")
    
    script:
        """
            extract_kraken_reads.py \
                -k ${kraken2_output} \
                --taxid 9606 \
                --exclude \
                -s1 ${reads[0]} \
                -s2 ${reads[1]} \
                -o ${id}_filtered_R1.fastq \
                -o2 ${id}_filtered_R2.fastq \
                --fastq-output   
        """
}

// merge unclassified and filtered reads
process MERGE {
    publishDir "${params.outdir}/03_merged_reads", failOnError: true, mode: "copy", overwrite: true

    input:
        tuple val(id), path(unclassified), path(filtered)

    output:
        tuple val("${id}"), path("${id}_R*.fastq.gz")

    script:
        """
        gzip -c ${unclassified[0]} ${filtered[0]} > ${id}_R1.fastq.gz
        gzip -c ${unclassified[1]} ${filtered[1]} > ${id}_R2.fastq.gz
        """
}
