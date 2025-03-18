# extract

Extract sequences and annotations from GFF files.

## Overview

The purpose of the [extract](../extract.html) subcommand is to break up the genome into regions according to the annotations. While many tools offer this functionality, [extract](../extract.html) is unique in that it will output both the annotated _and_ the unannotated regions.

These regions will later be used in the [cluster](../cluster.html) command to identify which parts of the genome should be aligned together for variant calling. If a region is annotated, that contextual information will also be used in the [summarize](../summarize.html) step, to give the cluster informative identifiers and gene/product names.

## Usage

The [extract](../extract.html) subcommand takes as input a single GFF file of sequence annotations. If the GFF file does not contain the sequences, a separate FASTA file should be additionally be provided.

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample1.gff3")
pangwas.extract(gff="sample1.gff3", fasta="sample1.fasta")
```

### CLI

```bash
pangwas extract --gff sample1.gff3
pangwas extract --gff sample1.gff3 --fasta sample1.fasta
```

:::

While [extract](../extract.html) will try to parse the sample identifier from the GFF file name, it's recommended to manually specify this information:

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample1.gff3", sample="sample1", prefix="sample1")
```

### CLI

```bash
pangwas extract --gff sample1.gff3 --sample sample1 --prefix sample1
```

:::

## Parameters

### Required

- `gff`: Input GFF annotations.

### Optional

Optional output parameters:

- `outdir`: Output directory.
- `prefix`: Output file prefix. If not provided, will be parsed from the gff file name.

Optional parameters:

- `fasta`: Input FASTA sequences, if not provided at the end of the GFF.
- `max-len`: Maximum length of sequences to extract (default: 100000).
- `min-len`: Minimum length of sequences to extract (default: 20).
- `sample`: Sample identifier to use. If not provided, is parsed from the gff file name.
- `regex`: Only extract gff lines that match this regular expression.

## Examples

### Minimal

The following is a minimal example of an input GFF file. The only attributes that are required is `ID`, all others are optional (ex. `gene`, `product`, etc.)

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample2.gff3")
```

### CLI

```bash
pangwas extract --gff sample2.gff3
```

:::

::: {.panel-tabset}

### Input

```text
sample2_multi	Bakta	region	1	130	.	+	.	ID=sample2_multi
sample2_multi	Prodigal	CDS	1	36	.	+	0	ID=multi_CDS2
sample2_multi	Prodigal	CDS	66	101	.	-	0	ID=multi_CDS1;gene=CDS1

##FASTA
>sample2_multi
ATGCGACGTAGCATGCAGCGCAGCTGAGCATCATAA
GGATCGATGCATCGGCGATTCACTGCATC
TTATGATGCTCAGCTGCGCTGCATGCTACGTCGCAT
GATGCAGTGAATCGCCGATGCATCGATCC
```

### Output

<div class="table-responsive">

: {.striped}

|sample |contig       |locus                                     |feature    |start|end|length|strand|upstream                              |downstream                                |attributes                                                                               |sequence_id                               |sequence                            |
|:------|:------------|:-----------------------------------------|:----------|:----|:--|:-----|:-----|:-------------------------------------|:-----------------------------------------|:----------------------------------------------------------------------------------------|:-----------------------------------------|:-----------------------------------|
|sample2|sample2_multi|sample2_multi_CDS2                        |CDS        |1    |36 |36    |+     |sample2_multi_TERMINAL                |sample2_multi_CDS2__sample2_multi_CDS1    |ID=multi_CDS2;Name=Fake protein;locus_tag=multi_CDS2;product=Fake protein;gene=multi_CDS2|sample2_multi_CDS2                        |ATGCGACGTAGCATGCAGCGCAGCTGAGCATCATAA|
|sample2|sample2_multi|sample2_multi_CDS2__sample2_multi_CDS1    |unannotated|37   |65 |29    |+     |sample2_multi_CDS2                    |sample2_multi_CDS1                        |                                                                                         |sample2_multi_CDS2__sample2_multi_CDS1    |GGATCGATGCATCGGCGATTCACTGCATC       |
|sample2|sample2_multi|sample2_multi_CDS1                        |CDS        |66   |101|36    |-     |sample2_multi_CDS2__sample2_multi_CDS1|sample2_multi_CDS1__sample2_multi_TERMINAL|ID=multi_CDS1;Name=Fake protein;locus_tag=multi_CDS1;product=Fake protein;gene=multi_CDS1|sample2_multi_CDS1                        |ATGCGACGTAGCATGCAGCGCAGCTGAGCATCATAA|
|sample2|sample2_multi|sample2_multi_CDS1__sample2_multi_TERMINAL|unannotated|102  |130|29    |-     |sample2_multi_CDS1                    |sample2_multi_TERMINAL                    |                                                                                         |sample2_multi_CDS1__sample2_multi_TERMINAL|GGATCGATGCATCGGCGATTCACTGCATC       |

</div>
::: 

### Real

These are real biological sequences that were annotated in the [annotate example](../annotate.html#minimal).

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample1.gff3")
```

### CLI

```bash
pangwas extract --gff sample1.gff3
```

:::

::: {.panel-tabset}

### Input

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

### Output

<div class="table-responsive">

: {.striped}

|sample |contig   |locus   |feature    |start|end |length|strand|upstream  |downstream  |attributes  |sequence_id  |sequence   |
|:------|:--------|:-------|:----------|:----|:---|:-----|:-----|:-----------------------------------------|:-----------------------------------------|:------|:-----------------------------------------|:-----------|
|sample1|sample1_1|sample1_AFLIOK_00005                      |CDS        |1    |702 |702   |+     |sample1_1_TERMINAL                        |sample1_AFLIOK_00005__sample1_AFLIOK_00010|ID=AFLIOK_00005;Name=streptococcal mitogenic exotoxin SmeZ;locus_tag=AFLIOK_00005;product=streptococcal mitogenic exotoxin SmeZ;Dbxref=BlastRules:WP_010922705,SO:0001217,UniRef:UniRef50_Q9RQQ5;gene=smeZ|sample1_AFLIOK_00005                      |ATGAAAAAAACAAAACTTATT... |
|sample1|sample1_1|sample1_AFLIOK_00005__sample1_AFLIOK_00010|unannotated|703  |743 |41    |+     |sample1_AFLIOK_00005  |sample1_AFLIOK_00010                      |  |sample1_AFLIOK_00005__sample1_AFLIOK_00010|AGTGACGACATGCTGACTGTGCTCTGACGTTGACTGACTGA |
|sample1|sample1_1|sample1_AFLIOK_00010                      |CDS        |744  |1457|714   |-     |sample1_AFLIOK_00005__sample1_AFLIOK_00010|sample1_1_TERMINAL                        |ID=AFLIOK_00010;Name=streptococcal pyrogenic exotoxin SpeM;locus_tag=AFLIOK_00010;product=streptococcal pyrogenic exotoxin SpeM;Dbxref=BlastRules:WP_011017838,SO:0001217,UniRef:UniRef50_Q8L3E1;gene=speM|sample1_AFLIOK_00010                      |ATGAAAAAAAATACCTTGACT...|

</div>

:::

Some important things to note:

1. **The locus identifiers have been given a prefix based on the sample identifier.**
    - `AFLIOK_00005` -> `sample1_AFLIOK_00005`
    - This will be a fail-safe against duplicate identifiers from other samples in later steps.
2. **Intergenic sequences are identified based on their upstream/downstream loci.**
    - The unannotated region between loci `AFLIOK_00005` and  `AFLIOK_00010` is given the identifier `sample1_AFLIOK_00005__sample1_AFLIOK_00010`, where the `__` character is the delimiter.

### Contigs

Contig coordinates will be identified by either the `sequence-region` comment line or a feature with the `region` type:

::: {.panel-tabset}

### Comment Line

```text
##sequence-region sample2_multi 1 130
```

### Feature

```text
sample2_multi   Bakta   region  1       130     .       +       .       ID=sample2_multi;sample2_multi
```

:::

### Duplicates

If duplicate identifiers are found, the sequence regions will be given a numeric suffix (ex. `loci.1`, `loci.2`).

::: {.panel-tabset}

### Input

```text
sample2_dup_gene        Bakta   region  1       120     .       +       .       ID=sample2_dup_gene;sample2_dup_gene
sample2_dup_gene        Prodigal        CDS     1       30      .       +       0       ID=dup_gene;Name=Duplicate gene name 1;locus_tag=dup_gene;Duplicate gene name 1;gene=dup_gene
sample2_dup_gene        Prodigal        CDS     61      90      .       +       0       ID=dup_gene;Name=Duplicate gene name 2.1;locus_tag=dup_gene;Duplicate gene name 2.1;gene=dup_gene
sample2_dup_gene        Prodigal        CDS     91      120     .       +       0       ID=dup_gene;Name=Duplicate gene name 2.2;locus_tag=dup_gene;Duplicate gene name 2.2;gene=dup_gene
```

### Output

<div class="table-responsive">

: {.striped}

|sample |contig          |locus                                 |feature    |start|end|length|strand|upstream                              |downstream                            |attributes               |sequence_id                           |sequence                      |
|:------|:---------------|:-------------------------------------|:----------|:----|:--|:-----|:-----|:-------------------------------------|:-------------------------------------|:------------------------|:-------------------------------------|:-----------------------------|
|sample2|sample2_dup_gene|sample2_dup_gene.1                    |CDS        |1    |30 |30    |+     |sample2_dup_gene_TERMINAL             |sample2_dup_gene.1__sample2_dup_gene.2|ID=dup_gene;gene=dup_gene|sample2_dup_gene.1                    |GACTATGCACTGCTCGCGCAGATCGTAGCG|
|sample2|sample2_dup_gene|sample2_dup_gene.1__sample2_dup_gene.2|unannotated|31   |60 |30    |+     |sample2_dup_gene.1                    |sample2_dup_gene.2                    |                         |sample2_dup_gene.1__sample2_dup_gene.2|ATATTATTCGGCATCTGATGTCTGCATGTG|
|sample2|sample2_dup_gene|sample2_dup_gene.2                    |CDS        |61   |90 |30    |+     |sample2_dup_gene.1__sample2_dup_gene.2|sample2_dup_gene.3                    |ID=dup_gene;gene=dup_gene|sample2_dup_gene.2                    |TTTAGAGTGCATGTGCACGGATGCATGCAC|
|sample2|sample2_dup_gene|sample2_dup_gene.3                    |CDS        |91   |120|30    |+     |sample2_dup_gene.2                    |sample2_dup_gene_TERMINAL             |ID=dup_gene;gene=dup_gene|sample2_dup_gene.3                    |TTTAGAGTGCATGTGCACGGATGCATGCAC|

</div>

:::

### Regex

> ❗ Using the `regex` parameter disables the extraction of unnannotated regions.

You can use the `regex` parameter to control which sequence regions are extracted. You might want to only extract annotated regions (which contain an `ID` attribute):

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample1.gff3", regex=".*ID=.*")
```

### CLI

```bash
pangwas extract --gff sample1.gff3 --regex ".*ID=.*"
```

:::

Or perhaps just a particular contig of interest:

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample1.gff3", regex="^sample1_pim")
```

### CLI

```bash
pangwas extract --gff sample1.gff3 --regex "^sample1_pim"
```

:::

Or a selection of genes:

::: {.panel-tabset group="language"}

### Python

```python
pangwas.extract(gff="sample1.gff3", regex=".*gene=(speG|speB).*")
```

### CLI

```bash
pangwas extract --gff sample1.gff3 --regex ".*gene=(speG|speB).*"
```

:::