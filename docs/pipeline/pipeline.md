# Pipeline

1. [Nextflow](#nextflow)
1. [Command-Line Interface](#command-line-interface)
1. [Python Package](#python-package)

## Nextflow

### Quick Start

An end-to-end pipeline is provided via `nextflow`. The following command runs the `test` profile, which is a small synthetic dataset that runs quickly:

```bash
nextflow run . -resume -profile test --trait resistant
nextflow run . -resume -profile test --trait lineage
```

**panGWAS** provides sensible defaults for each of the tools to get you started. But each step can be extensively customized, such as the following example, which uses custom arguments for clustering (`mmseqs`), alignment (`mafft`), and `iqtree`:

```bash
nextflow run . -resume -profile test --trait lineage \
  --cluster_args "-k 13 --min-seq-id 0.95 -c 0.95" \
  --align_args "--adjustdirection --localpair --maxiterate 1000 --addfragments" \
  --tree_args "-safe -m MFP --ufboot 1000 --alrt 1000 -o sample1"
```

### Samplesheet

Input data is defined in typical nf-core fashion with a `samplesheet.csv`. It must have at minimum the column `sample`. If you wish to annotate samples, paths to genomic fasta assemblies should be specified under `assembly`. If you have your own annotations, provide them under `gff`. Any missing data or paths should be left empty.

<div class="table-responsive">

: {.striped}

|sample |assembly                         |gff                       |lineage|resistant|mic  |
|:------|:--------------------------------|:-------------------------|:------|:--------|-----|
|sample1|data/test/sequences/sample1.fasta|data/test/gff/sample1.gff3|1      |0        | 0.0 |
|sample2|data/test/sequences/sample2.fasta|data/test/gff/sample2.gff3|2      |0        |     |
|sample3|                                 |data/test/gff/sample3.gff3|2      |1        | 2.1 |
|sample4| data/test/sequences/sample4.fasta | | 3      |1        | 2.9 |

</div>

To run the samplesheet with defaults, use:

```bash
nextflow run . -resume --input samplesheet.csv --outdir results
```

If you would like to run the GWAS step, you will also need to have column(s) that represent your trait(s) of interest. In this case, we also have traits like `resistant` and `mic` which can be analyzed with:

```bash
nextflow run . -resume --input samplesheet.csv --outdir results --trait resistant
nextflow run . -resume --input samplesheet.csv --outdir results --trait mic --gwas_args " --continuous"
```

The variable `mic` (minimum inhibitory concentration) is a continuous phenotype, so we specify that with the `--gwas_args` parameter (note the space within the quotations).

### Configuration

To see the full list of configurable parameters, run:

```bash
nextflow run . --help
```

These parameters can be configured either at runtime, or within a `nextflow` config file such as:

```json
params {
  // required paths
  input  = "samplesheet.csv"
  outdir = "results/test"

  // optional: resources
  max_cpus   = 2
  max_memory = '8.GB'

  // optional: step configuration
  extract_args = "--min-len 10"
  snps_args    = "--core 0.75 --indel-window 3 --snp-window 3"
  tree_args    = "-safe -m MFP --ufboot 1000 --alrt 1000 --seed 123456 -o sample1"
  heatmap_args = "--tree-width 100"
}
```

## Command Line Interface

The individual modules of the pipeline can be called as CLI commands. This code block is essentially the entirety of the `nextflow` pipeline.

```bash
# Extract sequence regions, both annotated and unannotated
pangwas extract --gff data/test/gff/sample1.gff3
pangwas extract --gff data/test/gff/sample2.gff3
pangwas extract --gff data/test/gff/sample3.gff3
pangwas extract --gff data/test/gff/sample4.gff3

# Collect sequence regions from all samples
pangwas collect --tsv sample1.tsv sample2.tsv sample3.tsv sample4.tsv --prefix collect

# Cluster sequences
pangwas cluster --fasta collect.sequences.fasta

# Defrag clusters
pangwas defrag --clusters clusters.tsv --representative representative.fasta --prefix defrag

# Summarize clusters
pangwas summarize --clusters clusters.tsv --regions collect.regions.tsv --prefix summarize

# Align clusters
pangwas align --clusters summarize.clusters.tsv --regions collect.regions.tsv

# Call variants
pangwas presence_absence --clusters summarize.clusters.tsv
pangwas structural --clusters summarize.clusters.tsv --alignments alignments
pangwas snps --alignment pangenome.aln --bed pangenome.bed --consensus pangenome.consensus.fasta --structural structural.Rtab
pangwas combine --rtab presence_absence.Rtab structural.Rtab snps.Rtab

# Phylogenetic Tree, using the defaults + some custom iqtree arguments
pangwas tree --alignment snps.core.fasta --constant-sites snps.constant_sites.txt -safe -m MFP --ufboot 1000 --alrt 1000 -o sample1

# GWAS on categorical variable
pangwas gwas --table data/test/samplesheet.csv --column lineage --variants combine.Rtab --tree tree.rooted.nwk --clusters summarize.clusters.tsv
# GWAS on binary variable
pangwas gwas --table data/test/samplesheet.csv --column resistant --variants combine.Rtab --tree tree.rooted.nwk --clusters summarize.clusters.tsv --lineage-column lineage
# GWAS on continuous variable
pangwas gwas --table data/test/samplesheet.csv --column mic --variants combine.Rtab --tree tree.rooted.nwk --clusters summarize.clusters.tsv  --lineage-column lineage --continuous

# Manhattan
pangwas manhattan --gwas resistant.locus_effects.significant.tsv --bed pangenome.bed --prefix resistant
pangwas manhattan --gwas lineage_1.locus_effects.significant.tsv --bed pangenome.bed --prefix lineage
pangwas manhattan --gwas mic.locus_effects.significant.tsv --bed pangenome.bed --prefix mic

# Heatmap
pangwas heatmap --gwas resistant.locus_effects.significant.tsv --tree tree.rooted.nwk --prefix resistant
pangwas heatmap --gwas lineage_1.locus_effects.significant.tsv --tree tree.rooted.nwk --prefix lineage
pangwas heatmap --gwas mic.locus_effects.significant.tsv --tree tree.rooted.nwk --prefix mic
```

## Python Package

The individual modules of the pipeline can be called as `python` functions:

```python
import pangwas

# Extract sequence regions, both annotated and unannotated
samplesheet = "data/test/samplesheet.csv"
sample1 = pangwas.extract(gff="data/test/gff/sample1.gff3")
sample2 = pangwas.extract(gff="data/test/gff/sample2.gff3")
sample3 = pangwas.extract(gff="data/test/gff/sample3.gff3")
sample4 = pangwas.extract(gff="data/test/gff/sample4.gff3")

# Collect sequence regions
(fasta, regions) = pangwas.collect(tsv = [sample1, sample2, sample3, sample4], prefix="collect")

# Cluster sequences
(clusters, representative) = pangwas.cluster(fasta=fasta)

# Defrag clusters
(clusters, representative) = pangwas.defrag(clusters=clusters, representative=representative, prefix="defrag")

# Summarize clusters
summarize = pangwas.summarize(clusters=clusters, regions=regions, prefix="summarize")

# Align clusters
pangwas.align(clusters=summarize, regions=regions)

# Call variants
presence_absence = pangwas.presence_absence(clusters=summarize)
structural       = pangwas.structural(clusters=summarize, alignments="alignments")
snps             = pangwas.snps(alignment="pangenome.aln", bed="pangenome.bed", consensus="pangenome.consensus.fasta", structural=structural)
variants         = pangwas.combine(rtab=[presence_absence, structural, snps])

# Phylogenetic Tree, using the defaults + some custom iqtree arguments
tree  = pangwas.tree(alignment="snps.core.fasta", constant_sites="snps.constant_sites.txt", args=pangwas.TREE_ARGS + ' --ufboot 1000 --alrt 1000 -o sample1')

# GWAS on categorical variable
pangwas.gwas(table=samplesheet, column="lineage",   variants=combine, tree=tree, clusters=summarize)
# GWAS on binary variable
pangwas.gwas(table=samplesheet, column="resistant", variants=combine, tree=tree, clusters=summarize, lineage_column="lineage")
# GWAS on continuous variable
pangwas.gwas(table=samplesheet, column="mic", variants=variants, tree=tree, clusters=summarize, lineage_column="lineage", continuous=True)

# Manhattan
pangwas.manhattan(gwas="resistant.locus_effects.significant.tsv", bed="pangenome.bed", prefix="resistant")
pangwas.manhattan(gwas="lineage_1.locus_effects.significant.tsv", bed="pangenome.bed", prefix="lineage")
pangwas.manhattan(gwas="mic.locus_effects.significant.tsv",       bed="pangenome.bed", prefix="mic")

# Heatmap
pangwas.heatmap(gwas="resistant.locus_effects.significant.tsv", tree=tree, prefix="resistant")
pangwas.heatmap(gwas="lineage_1.locus_effects.significant.tsv", tree=tree, prefix="lineage")
pangwas.heatmap(gwas="mic.locus_effects.significant.tsv",       tree=tree, prefix="mic")
```