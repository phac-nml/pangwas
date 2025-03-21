# Tutorial 03 - Pyseer Tutorial

This tutorial automates and reproduces the results from the [penicillin resistance GWAS](https://pyseer.readthedocs.io/en/master/tutorial.html) written by the `pyseer` authors.

1. Download the tutorial data.

    ```bash
    git clone --depth 1 https://github.com/phac-nml/pangwas.git
    cd pangwas
    ```

1. Decompress the core genome data.

    ```bash
    gunzip data/tutorial_core/snps.Rtab.gz
    gunzip data/tutorial_core/gff.gz
    gunzip data/tutorial_core/clusters.tsv.gz
    ```

1. Decompress the pan genome data.

    ```bash
    gunzip data/tutorial_pangenome/variants.Rtab.gz
    gunzip data/tutorial_pangenome/clusters.tsv.gz
    ```

1. Run the core genome GWAS.

    ```bash
    nextflow run phac-nml/pangwas -profile tutorial_core --trait penicillin
    ```

1. Run the pangenome GWAS.

    ```bash
    nextflow run phac-nml/pangwas -profile tutorial_pangenome --trait penicillin
    ```

Penicillin resistance is primarily controlled by core genome genes, and we can see that the major genes are identical between a pangenome and core genome GWAS.

![](..//images/core_vs_pangenome.png)
