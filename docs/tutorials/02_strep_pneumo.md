# Tutorial 02 - Streptococcus pneumoniae

This tutorial analyzes 15 _Streptococcus pneumoniae_ genomes from the [penicillin resistance GWAS](https://pyseer.readthedocs.io/en/master/tutorial.html) example.

1. Run the test data through a GWAS on the trait `penicillin`.

    ```bash
    nextflow run phac-nml/pangwas -profile streptococcus_pneumoniae --max_cpus 4
    ```

## Iroki

```bash
csvtk cut -f sample,lineage,penicillin data/streptococcus_pneumoniae/samplesheet.csv \
  | csvtk rename -f sample -n name \
  | csvtk rename -f lineage,penicillin -n bar1_color,bar2_height \
  | csvtk replace -f bar1_color -p '(.*)' -r 'k_$1' \
  | csvtk mutate2 -n bar1_height -e '"1"' \
  | csvtk mutate2 -n bar2_color -e '$bar2_height == 1 ? "k_4" : "k_5"' \
  | csvtk replace -f bar2_height -p 0 -r "-1" \
  | csvtk csv2tab \
  > data/streptococcus_pneumoniae/iroki.tsv
```