nextflow.enable.dsl = 2

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


/************************** 
-----------CALL------------
**************************/

// $ nextflow ../scratch/rki_kraken2/Pipeline/Scripts/create_hiv_krakendb.nf \
// -c Scripts/rki_config \
// -profile rki_slurm,rki_mamba \
// -with-report logs/test.$(date +%T).report.html 

/************************** 
DEFINE VARIABLES
**************************/

projectDir = "/home/rykalinav/scratch/rki_kraken2/Pipeline"


// Parameters for ncbi
params.taxon = "11676"

/************************** 
---------WORKFLOW----------
**************************/
ch_infiles = channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz")

workflow {

    ch_hiv_genomes = GET_GENOMES(params.taxon)

}

/************************** 
PROCESSES
**************************/

// kraken2
process GET_GENOMES {

    //label "ncbi"
    conda "/home/rykalinav/.conda/envs/ncbi-datasets-cli_v15.25.0"
    //conda "${projectDir}/Environments/ncbi-datasets-cli_v15.25.0.yml"
    publishDir "${params.outdir}/01_complete_hiv_genomes", mode: "copy", overwrite: true

    input:
        val (taxid)

    output:
        path("ncbi_dataset/data/*.fna"),   emit: fna
        path("ncbi_dataset/data/*.md"),   emit: md
        

    script:
        """
        datasets download virus genome taxon ${taxid} \
            --host human   \
            --complete-only \
            --include genome
        unzip ncbi_dataset.zip
        """
     
}
