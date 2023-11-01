# A short desctiption on how to create a custom Kraken2 DB for HIV-1

## Prerequisite
Use NCBI Vrius to download all available complete HIV-1 genomes (NCBI Datasets CLI).

- Install and activate **ncbi-datasets-cli_v15.25.0** environment

```sh
conda create -n ncbi-datasets-cli_v15.25.0 -c conda-forge ncbi-datasets-cli=15.25.0
```

- Prepare multifasta
(NB) The fasta sequence headers must include either NCBI accession numbers or the text kraken:taxid followed by the taxonomy ID for the genome (e.g. >sequence_name|kraken:taxid|11676|).


## Steps
- Install and activate **kraken2_v2.1.2** environment

```sh
conda create -n kraken2_v2.1.2 -c bioconda kraken=2.1.2
```

- Download the NCBI taxonomy

```sh 
kraken2-build --db hiv_krakendb --download-taxonomy
```

- Downloading standard Kraken2 libraries (archeae, bacteria, fungi, human, viral)

```sh 
kraken2-build --db hiv_krakendb --download-library archaea
kraken2-build --db hiv_krakendb --download-library bacteria
kraken2-build --db hiv_krakendb --download-library fungi
kraken2-build --db hiv_krakendb --download-library human
kraken2-build --db hiv_krakendb --download-library viral
```

- Download additional library (all complete HIV-1 genomes)

```sh 
kraken2-build --db hiv_krakendb --add-to-library hiv_complete_genomes.fasta
```

- Build the actual custom database

```sh 
kraken2-build --db hiv_krakendb --build --threads 8
```

- Remove unneeded files from a built database

```sh 
kraken-build --clean --db hiv_krakendb --threads 8
```
