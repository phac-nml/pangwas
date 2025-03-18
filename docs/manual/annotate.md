
# annotate

Annotate genomic assemblies with [bakta](https://github.com/oschwengers/bakta).

## Overview

The purpose of the [annotate](../annotate.html) subcommand is to provide a very simple wrapper around the bacterial genome annotator [bakta](https://github.com/oschwengers/bakta). The main added value of using [annotate](../annotate.html) is that it has a post-processing step where if fixes bad quotation marks in the output files (ex. `‘`, `’`), as these particular characters tend to crash downstream applications.

The output GFF file will the be used by the [extract](../extract.html) subcommand to break up the genome into informative regions for clustering.


Genome annotation is the most complicated analysis and therefore most fragile part of the pipeline. If using the `nextflow` runtime, it is recommended to first run the annotation step in isolation with `--only_annotate`, which skips all downstream steps. After the annotation is successful, you can then include your new `gff` file paths in your samplesheet. In the case where `nextflow` fails to resume a previously run properly, this will prevent the pipeline from running annotation all over again.

## Database

If you are not bringing your own annotations, you will need to download a `bakta` annotation database.

**panGWAS** has been verified with `bakta` CLI `v1.9.2` and the `v5.1` database. Specifically, `v5.1` includes the `amrfinderplus-db` release `v2023-11-15.1`, which was the last release before their Database Schema changed. The full datbase can be downloaded with the following command:

>❗ **WARNING**: This download size is +30 GB. If you want to use the light database for testing (1 GB), replace `db` with `db-light`.

```bash
mkdir -p database/bakta/
wget -O database/bakta/db.tar.gz https://zenodo.org/records/10522951/files/db.tar.gz
tar -xvf database/bakta/db.tar.gz -C database/bakta --strip-components 1
```

The database directory can then be tested with the following command:

```bash
pangwas annotate \
  --db database/bakta \
  --fasta data/test/sequences/sample1.fasta \
  --outdir sample1 \
  --threads 2
```

## Usage

The [annotate](../extract.html) subcommand takes as input a FASTA file of genomic sequences and a path to a [bakta](https://github.com/oschwengers/bakta) database. Typically the FASTA file is in the form of contigs from a denovo or reference-guided assembly.

::: {.panel-tabset group="language"}

### Python

```python
pangwas.annotate(fasta="sample1.fasta", db="database/bakta")
```

### CLI

```bash
pangwas annotate --fasta sample1.fasta --db database/bakta
```

:::

While [annotate](../annotate.html) will try to parse the sample identifier from the FASTA file name, it's recommended to manually specify this information:

::: {.panel-tabset group="language"}

### Python

```python
pangwas.annotate(fasta="sample1.fasta", db="database/bakta", sample="sample1", prefix="sample1")
```

### CLI

```bash
pangwas annotate --fasta sample1.fasta --db database/bakta --sample sample1 --prefix sample1
```

:::

If you would like to pass additional arguments to [bakta](https://github.com/oschwengers/bakta), they can be specified like the following:

::: {.panel-tabset group="language"}

### Python

```python
pangwas.annotate(fasta="sample1.fasta", db="database/bakta", args="--genus Streptococcus --species pneumoniae")
```

### CLI

```bash
pangwas annotate --fasta sample1.fasta --db database/bakta --genus Streptococcus --species pneumoniae
```

:::

## Parameters

### Required

- `fasta`: Input FASTA sequences.
- `db`: bakta database directory.

### Optional

Optional output arguments:

- `outdir`: Output directory.
- `prefix`: Output file prefix. If not provided, will be parsed from the fasta file name.
- `tmp`: Temporary directory.

Optional arguments:

- `sample` Sample identifier. If not provided, will be parsed from the fasta file name.
- `threads` CPU threads for bakta.

## Examples

### Minimal

The following is a minimal example of annotating a `Streptococcus pyogenes` contig that contains the genes `smeZ` and `speM` separated by a short intergenic region.

::: {.panel-tabset group="language"}

### Python

```python
pangwas.annotate(
  fasta  = "sample1.fasta",
  db     = "database/bakta",
  locus  = "sample1",
  outdir = "sample1",
  args   = "--genus Streptococcus --species pneumoniae"
)
```

### CLI

```bash
pangwas annotate \
  --fasta sample1.fasta \
  --db database/bakta \
  --locus sample1 \
  --genus Streptococcus \
  --species pneumoniae
```

:::

::: {.panel-tabset}

### Input

```text
>contig
ATGAAAAAAACAAAACTTATTTTTTCTTTTACTTCAATATTCATTGCAATAATTTCTCGTCCTGTGTTTG
GATTAGAAGTAGATAATAATTCCCTTCTAAGGAATATCTATAGTACGATTGTGTATGAATATTCAGATAC
AGTAATTGAGTTTAAAACCAGTCATAACTTAGTGACTAAGAAACTTGATGTTAGAGATGCTAGAGATTTT
TTTATTAACTCCGAAATGGACGAATATGCAGCCAATGATTTTAAAGATGGAGATAAAATAGCTATGTTCT
CCGTCCCATTTGATTGGAACTACTTGTCAGAAGGAAAAGTCATAGCATATACGTATGGTGGAATGACGCC
TTATCAAGAAGAACCAATATCTAAAAATATCCCTGTTAATTTATGGATTAATGGAAAGCAGATCTCTGTT
CCTTACAACGAAATATCAACTAACAAAACAACAGTTACAGCTCAAGAAATTGATCTAAAGGTTAGAAAAT
TTTTAATATCACAACATCAATTATATTCTTCTGGTTCTAGCTACAAAAGTGGTAAATTAGTTTTTCATAC
AAATGATAATTCAGATAAATATTCTCTCGATCTTTTCTATGTAGGATATAGAGATAAAGAAAGTATTTTT
AAAGTATACAAAGACAATAAATCTTTCAATATAGATAAAATTGGGCATTTAGATATAGAAATTGACTCCT
AA
AGTGACGACATGCTGACTGTGCTCTGACGTTGACTGACTGA
CTAATTTTTAGAAAAATCTTCGTTTAAGTAAATATCAAAGTGACTTACTTTACTCATATCAATCGTTTCATTATCTGT
ATAGTTAGGATGAGTGAATAAATCGGTAAACTTTGTTGTATTATCTTTATAATGAATTCCCCAATACCCTTTTTTACA
AATTGAGTTATGTTCATATAACTTTATTCTATTATCGCTCATCAAACTTTTCCTAAGTCTAACATCAATTTCTTGAAA
AGTTACAAACTTATTCTGAAATGTGATTTTATATTTTGATTGCTCTTTTAGCGGTATCTGTTCCCCAAAAATATTCAT
ATATATTGTTGAATCAAGTTTCTCTCTGTCACTTGTTCTTATCAAACCTCCATCAACATAATTATATTGTTCCTTACA
TATCACACTGTAGGATTTTATTAGAGCATAAATATCAACTTCTTCTTCCTTAAAGCGTCTTTCTTGCGCTGGAGAAAC
ATTGCTAGATATAACTTTATTATAATTATCATCATCCCAGACTCTAGTTTTTTCATTCGTGTTAAATATTAACTGGGT
GCCAATTTTCTTTGTTATCTTCATATTGGTTCTATTAATTACATCTTTCGTATAGATATTTTTTAATTCGCTATTAAC
CAACACAGCATCTGAAAAGACACTCTCAGTAGTGTATAGAGCAAGCGATGCACACACAAGGAATAACAAAGTCAAGGT
ATTTTTTTTCAT
```

> Line wrapping is optional, and is used for illustrative purposes in this example to show the two separate genes present and the intergenic sequence between them.


### Output

```text
##gff-version 3
##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo
# organism Streptococcus pyogenes
# Annotated with Bakta
# Software: v1.9.2
# Database: v5.1, light
# DOI: 10.1099/mgen.0.000685
# URL: github.com/oschwengers/bakta
##sequence-region sample1_1 1 1457
sample1_1	Bakta	region	1	1457	.	+	.	ID=sample1_1;Name=sample1_1
sample1_1	Prodigal	CDS	1	702	.	+	0	ID=AFLIOK_00005;Name=streptococcal mitogenic exotoxin SmeZ;locus_tag=AFLIOK_00005;product=streptococcal mitogenic exotoxin SmeZ;Dbxref=BlastRules:WP_010922705,SO:0001217,UniRef:UniRef50_Q9RQQ5;gene=smeZ
sample1_1	Prodigal	CDS	744	1457	.	-	0	ID=AFLIOK_00010;Name=streptococcal pyrogenic exotoxin SpeM;locus_tag=AFLIOK_00010;product=streptococcal pyrogenic exotoxin SpeM;Dbxref=BlastRules:WP_011017838,SO:0001217,UniRef:UniRef50_Q8L3E1;gene=speM
##FASTA
>sample1_1
ATGAAAAAAACAAAACTTATTTTTTCTTTTACTTCAATATTCATTGCAATAATTTCTCGT
CCTGTGTTTGGATTAGAAGTAGATAATAATTCCCTTCTAAGGAATATCTATAGTACGATT
GTGTATGAATATTCAGATACAGTAATTGAGTTTAAAACCAGTCATAACTTAGTGACTAAG
AAACTTGATGTTAGAGATGCTAGAGATTTTTTTATTAACTCCGAAATGGACGAATATGCA
GCCAATGATTTTAAAGATGGAGATAAAATAGCTATGTTCTCCGTCCCATTTGATTGGAAC
TACTTGTCAGAAGGAAAAGTCATAGCATATACGTATGGTGGAATGACGCCTTATCAAGAA
GAACCAATATCTAAAAATATCCCTGTTAATTTATGGATTAATGGAAAGCAGATCTCTGTT
CCTTACAACGAAATATCAACTAACAAAACAACAGTTACAGCTCAAGAAATTGATCTAAAG
GTTAGAAAATTTTTAATATCACAACATCAATTATATTCTTCTGGTTCTAGCTACAAAAGT
GGTAAATTAGTTTTTCATACAAATGATAATTCAGATAAATATTCTCTCGATCTTTTCTAT
GTAGGATATAGAGATAAAGAAAGTATTTTTAAAGTATACAAAGACAATAAATCTTTCAAT
ATAGATAAAATTGGGCATTTAGATATAGAAATTGACTCCTAAAGTGACGACATGCTGACT
GTGCTCTGACGTTGACTGACTGACTAATTTTTAGAAAAATCTTCGTTTAAGTAAATATCA
AAGTGACTTACTTTACTCATATCAATCGTTTCATTATCTGTATAGTTAGGATGAGTGAAT
AAATCGGTAAACTTTGTTGTATTATCTTTATAATGAATTCCCCAATACCCTTTTTTACAA
ATTGAGTTATGTTCATATAACTTTATTCTATTATCGCTCATCAAACTTTTCCTAAGTCTA
ACATCAATTTCTTGAAAAGTTACAAACTTATTCTGAAATGTGATTTTATATTTTGATTGC
TCTTTTAGCGGTATCTGTTCCCCAAAAATATTCATATATATTGTTGAATCAAGTTTCTCT
CTGTCACTTGTTCTTATCAAACCTCCATCAACATAATTATATTGTTCCTTACATATCACA
CTGTAGGATTTTATTAGAGCATAAATATCAACTTCTTCTTCCTTAAAGCGTCTTTCTTGC
GCTGGAGAAACATTGCTAGATATAACTTTATTATAATTATCATCATCCCAGACTCTAGTT
TTTTCATTCGTGTTAAATATTAACTGGGTGCCAATTTTCTTTGTTATCTTCATATTGGTT
CTATTAATTACATCTTTCGTATAGATATTTTTTAATTCGCTATTAACCAACACAGCATCT
GAAAAGACACTCTCAGTAGTGTATAGAGCAAGCGATGCACACACAAGGAATAACAAAGTC
AAGGTATTTTTTTTCAT
```

::: 