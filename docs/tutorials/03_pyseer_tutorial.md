# Tutorial 03 - Pyseer Tutorial

This tutorial automates and reproduces the results from the [penicillin resistance GWAS](https://pyseer.readthedocs.io/en/master/tutorial.html) written by the `pyseer` authors.

1. Download the tutorial data.

    ```bash
    git clone --depth 1 https://github.com/phac-nml/pangwas.git
    cd pangwas
    ```

1. Decompress the data.

    ```bash
    gunzip data/tutorial_core/snps.Rtab.gz
    gunzip data/tutorial_pangenome/variants.Rtab.gz
    gunzip data/tutorial_pangenome/clusters.tsv.gz
    ```

1. Run the core genome GWAS on pencillin resistance.

    ```bash
    nextflow run phac-nml/pangwas -profile tutorial_core
    ```

1. Open up the manhattan plot in Edge or Firefox.

    - Path: `results/tutorial_core/manhattan/penicillin/penicillin.plot.svg`
    - Hover your mouse over variants to see contextual information.

    ![](../images/core_manhattan_hovertext.png)

1. Run the pangenome GWAS on penicillin resistance.

    ```bash
    nextflow run phac-nml/pangwas -profile tutorial_pangenome
    ```

Penicillin resistance is primarily controlled by core genome genes, and we can see that the major genes are identical between a pangenome and core genome GWAS.

The following image was made by importing the raw SVG plots into a vector graphics program (ex. [Affinity](https://affinity.serif.com/en-us/designer)) and adding text and box overlays to highlight the significant variants.

![](..//images/core_vs_pangenome.png)
