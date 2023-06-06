nextflow.enable.dsl = 2

projectDir = "/home/beast2/rki_shiver/Pipeline"
params.trimmomatic = "${projectDir}/Scripts/bin/trimmomatic-0.36.jar"
params.gal_primers = "${projectDir}/DataShiverInit/primers_GallEtAl2012.fasta"
params.illumina_adapters = "${projectDir}/DataShiverInit/adapters_Illumina.fasta"
params.alignment = "${projectDir}/DataShiverInit/HIV1_COM_2012_genome_DNA_NoGaplessCols.fasta"
params.config = "${projectDir}/Scripts/bin/config.sh"
params.remove_whitespace = "${projectDir}/Scripts/bin/tools/RemoveTrailingWhitespace.py"
params.maf = "${projectDir}/Scripts/bin/produce_maf.py"

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}

process RENAME_FASTQ {
  //conda "/home/beast2/anaconda3/envs/python3"
  //conda "${projectDir}/Environments/python3.yml"
  publishDir "${params.outdir}/0_renamed_fastq", mode: "copy", overwrite: true

  input:
   path file

  output:
    path "*.fastq.gz"

  script:
   """
   regex-rename '\\d{6}_\\d{2}-\\d{5}_HIV(\\d{2}-\\d{5})_\\w{2,}_?\\d{2,}?_\\w{3}_\\w{4}_(R\\d)_\\d{3}.(\\w{5}.\\w{2})' '\\1_\\2.\\3' --rename
   """
}

process INITIALISATION {
  //conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  conda "/home/beast2/anaconda3/envs/shiver"
  //conda "${projectDir}/Environments/shiver.yml"
  publishDir "${projectDir}/${params.outdir}/1_init_dir", mode: "copy", overwrite: true

  input:
     
  output:
     path "InitDir"
  script:
  
  """
  shiver_init.sh InitDir ${params.config} ${params.alignment} ${params.illumina_adapters} ${params.gal_primers}
  """

}

process IVA_CONTIGS {
  label "iva"
  //errorStrategy 'ignore'
  conda "/home/beast2/anaconda3/envs/iva"
  //conda "/usr/local/Caskroom/miniconda/base/envs/iva"
  //conda "${projectDir}/Environments/iva.yml"
  publishDir "${params.outdir}/2_iva_contigs", mode: "copy", overwrite: true
  
  input:
    tuple val(id), path(reads)

  output:
    tuple val("${id}"), path("${id}/${id}_contigs.fasta")

  script:
    """
    iva -f ${reads[0]} -r ${reads[1]} --pcr_primers ${params.gal_primers} --adapters ${params.illumina_adapters} --trimmomatic ${params.trimmomatic} ${id}
    mv ${id}/contigs.fasta ${id}/${id}_contigs.fasta
    """
}

process ALIGN_CONTIGS {
  //errorStrategy "ignore"
  conda "/home/beast2/anaconda3/envs/shiver"
  //conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  //conda "${projectDir}/Environments/shiver.yml"
  publishDir "${params.outdir}/3_alignments/${id}", mode: "copy", overwrite: true
 
  
  input:
    path initdir
    tuple val(id), path(contigs)

  output:
    tuple val("${id}"), path("${id}*.fasta"), path("${id}.blast")

  script:
    """
    shiver_align_contigs.sh ${initdir} ${params.config} ${contigs} ${id}
    rm temp_*
    rm *_MergedHits.blast*
    """
}


process ID_FASTQ {
  conda "/home/beast2/anaconda3/envs/shiver"
  //conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  //conda "${projectDir}/Environments/shiver.yml"
  publishDir "${params.outdir}/4_id_fastq", mode: "copy", overwrite: true

  input:
    tuple val(id), path(fastq)
  output:
    tuple val("${id}"), path("${fastq[0].getBaseName().split("_R")[0]}_RWS_1.fastq"), path("${fastq[1].getBaseName().split("_R")[0]}_RWS_2.fastq")
  
  script:
   """
   zcat ${fastq[0]} | awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' | sed 's/:N:.*//' > ${fastq[0].getBaseName().split("_R")[0]}_1.fastq
   rm ${fastq[0]}
   python ${params.remove_whitespace} ${fastq[0].getBaseName().split("_R")[0]}_1.fastq > ${fastq[0].getBaseName().split("_R")[0]}_RWS_1.fastq

   zcat ${fastq[1]} | awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' | sed 's/:N:.*//' > ${fastq[1].getBaseName().split("_R")[0]}_2.fastq
   rm ${fastq[1]}
   python ${params.remove_whitespace} ${fastq[0].getBaseName().split("_R")[0]}_2.fastq > ${fastq[0].getBaseName().split("_R")[0]}_RWS_2.fastq
   """
}

process MAP {
  conda "/home/beast2/anaconda3/envs/shiver"
  //conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  //conda "${projectDir}/Environments/shiver.yml"
  publishDir "${params.outdir}/5_mapped/${id}", mode: "copy", overwrite: true
  debug true

  input:
    path initdir
    tuple val(id), path(contigs), path(refs), path(blast), path(read1), path(read2)
    
  output:
    tuple val("${id}"), path("${id}_remap_ref.fasta"), path("${id}*.bam"), path("${id}*.bam.bai"), path("${id}*WithHXB2.csv")
    
  script:
    if (refs instanceof List) {
    """
    shiver_map_reads.sh ${initdir} ${params.config} ${contigs} ${id} ${blast} ${refs[0]} ${read1} ${read2}
    rm temp_* 
    rm *PreDedup.bam
    
    """ 
    } else {
     """
    shiver_map_reads.sh ${initdir} ${params.config} ${contigs} ${id} ${blast} ${refs} ${read1} ${read2}
    rm temp_* 
    rm *PreDedup.bam
    
     """
    }
}
  process MAF {
  conda "/home/beast2/anaconda3/envs/python3"
  //conda "${projectDir}/Environments/python3.yml"
  publishDir "${params.outdir}/6_maf", mode: "copy", overwrite: true
  debug true

  input:
    tuple val(id), path(ref), path(bam), path(bai), path(csv)
    
  output:
    tuple val("${id}"), path("${id}_MAF.csv")
    
  script:
    if (csv instanceof List) {
    """
    produce_maf.py ${csv} ${id}_MAF.csv
    """ 
    } else {
     """
    produce_maf.py ${csv} ${id}_MAF.csv
     """
    }
  }



workflow {
  fastq_pairs = channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz")
  initdir = INITIALISATION()
  iva_contigs = IVA_CONTIGS(fastq_pairs)
  refs = ALIGN_CONTIGS(initdir, iva_contigs)
  id_fastq = ID_FASTQ(fastq_pairs)
  // Combine according to a key that is the first value of every first element, which is a list
  map_args = iva_contigs.combine(refs, by:0).combine(id_fastq, by:0)
  map_out = MAP(initdir, map_args)
  maf_out = MAF(map_out).view()
 
 
}

  // Combine according to a key that is the first value of every first element, which is a list
  //map_args = id_contigs.combine(id_blast, by:0).combine(id_refs, by:0).combine(id_fastq, by:0)
  // Get rid off id
  //no_id_args = map_args.map { id, contigs, blast, ref_cut, read1, read2 -> [contigs, blast, ref_cut, read1, read2]}
  //MAP(initdir,no_id_args)
  


