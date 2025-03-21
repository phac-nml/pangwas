#!/usr/bin/env python3

# Base python packages
from collections import OrderedDict, Counter
import copy
import itertools
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap

from tqdm import tqdm

NUCLEOTIDES         = ["A", "C", "G", "T"]

# Default arguments for programs
CLUSTER_ARGS = "-k 13 --min-seq-id 0.90 -c 0.90 --cluster-mode 2 --max-seqs 300"
DEFRAG_ARGS  = "-k 13 --min-seq-id 0.90 -c 0.90 --cov-mode 1"
ALIGN_ARGS   = "--adjustdirection --localpair --maxiterate 1000 --addfragments"
TREE_ARGS    = "-safe -m MFP"
GWAS_ARGS    = "--lmm"

LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
logging.basicConfig(level=LOGLEVEL, stream=sys.stdout, format='%(asctime)s %(funcName)20s %(levelname)8s: %(message)s')

# Notify users when they use an algorithm from ppanggolin
PPANGGOLIN_NOTICE =  """
    The defrag algorithm is adapted from PPanGGOLiN (v2.2.0) which is  
    distributed under the CeCILL FREE SOFTWARE LICENSE AGREEMENT LABGeM.
    - Please cite: Gautreau G et al. (2020) PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph.
                    PLOS Computational Biology 16(3): e1007732. https://doi.org/10.1371/journal.pcbi.1007732
    - PPanGGOLiN license: https://github.com/labgem/PPanGGOLiN/blob/2.2.0/LICENSE.txt
    - Defrag algorithm source: https://github.com/labgem/PPanGGOLiN/blob/2.2.0/ppanggolin/cluster/cluster.py#L317
"""

pangwas_description          = "pangenome wide association study (panGWAS) pipeline."
annotate_description         = "Annotate genomic assemblies with bakta."
extract_description          = "Extract sequences and annotations from GFF files."
collect_description          = "Collect extracted sequences from multiple samples into one file."
cluster_description          = "Cluster nucleotide sequences with mmseqs."
defrag_description           = "Defrag clusters by associating fragments with their parent cluster."
summarize_description        = "Summarize clusters according to their annotations."
align_description            = "Align clusters using mafft and create a pangenome alignment."
snps_description             = "Extract SNPs from a pangenome alignment."
structural_description       = "Extract structural variants from cluster alignments."
presence_absence_description = "Extract presence absence of clusters."
combine_description          = "Combine variants from multiple Rtab files."
table_to_rtab_description    = "Convert a TSV/CSV table to an Rtab file based on regex filters."
tree_description             = "Estimate a maximum-likelihood tree with IQ-TREE."
root_tree_description        = "Root tree on outgroup taxa."
binarize_description         = "Convert a categorical column to multiple binary (0/1) columns."
vcf_to_rtab_description      = "Convert a VCF file to an Rtab file."
gwas_description             = "Run genome-wide association study (GWAS) tests with pyseer."
heatmap_description          = "Plot a heatmap of variants alongside a tree."
manhattan_description        = "Plot the distribution of variant p-values across the genome."

def get_options(args:str=None):
    import argparse

    sys_argv_original = sys.argv

    if args != None:
        sys.argv = args.split(" ")

    description = textwrap.dedent(
    f"""\
    {pangwas_description}

    ANNOTATE
      annotate:         {annotate_description}

    CLUSTER
      extract:          {extract_description}
      collect:          {collect_description}
      cluster:          {cluster_description}
      defrag:           {defrag_description}
      summarize:        {summarize_description}

    ALIGN
      align:            {align_description}

    VARIANTS
      structural:       {structural_description}    
      snps:             {snps_description}
      presence_absence: {presence_absence_description}

    TREE
      tree:             {tree_description}

    GWAS
      gwas:             {gwas_description}

    PLOT
      manhattan:        {manhattan_description}
      heatmap:          {heatmap_description}

    UTILITY
      root_tree:        {root_tree_description}
      binarize:         {binarize_description}
      table_to_rtab:    {table_to_rtab_description}
      vcf_to_rtab:      {vcf_to_rtab_description}
    """)
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', help='Display program version.', action="store_true")
    subbcommands = parser.add_subparsers(dest="subcommand")

    # -------------------------------------------------------------------------
    # Annotate

    description = textwrap.dedent(
    f"""\
    {annotate_description}

    Takes as input a FASTA file of genomic assemblies. Outputs a GFF file
    of annotations, among many other formats from bakta.

    All additional arguments with be passed to the `bakta` CLI.

    Examples:
    > pangwas annotate --fasta sample1.fasta --db database/bakta
    > pangwas annotate --fasta sample1.fasta --db database/bakta --sample sample1 --threads 2 --genus Streptococcus
    """)

    annotate_parser = subbcommands.add_parser('annotate', description = description, formatter_class=argparse.RawTextHelpFormatter)

    annotate_req_parser = annotate_parser.add_argument_group("required arguments")  
    annotate_req_parser.add_argument("--fasta",  required=True, help='Input FASTA sequences.')
    annotate_req_parser.add_argument("--db",     required=True, help='bakta database directory.')

    annotate_out_opt_parser = annotate_parser.add_argument_group("optional output arguments")
    annotate_out_opt_parser.add_argument("--outdir",  help='Output directory. (default: .)', default=".")
    annotate_out_opt_parser.add_argument('--prefix',  help='Output file prefix. If not provided, will be parsed from the fasta file name.')
    annotate_out_opt_parser.add_argument("--tmp",     help='Temporary directory.')

    annotate_opt_parser = annotate_parser.add_argument_group("optional arguments")
    annotate_opt_parser.add_argument("--sample",  help='Sample identifier. If not provided, will be parsed from the fasta file name.')
    annotate_opt_parser.add_argument("--threads", help='CPU threads for bakta. (default: 1)', default=1)

    # -------------------------------------------------------------------------
    # Extract Sequences

    description = textwrap.dedent(
    f"""\
    {extract_description}

    Takes as input a GFF annotations file. If sequences are not included, a FASTA
    of genomic contigs must also be provided. Both annotated and unannotated regions 
    will be extracted. Outputs a TSV table of extracted sequence regions.

    Examples:
    > pangwas extract --gff sample1.gff3
    > pangwas extract --gff sample1.gff3 --fasta sample1.fasta --min-len 10
    """)

    extract_parser = subbcommands.add_parser('extract', description = description, formatter_class=argparse.RawTextHelpFormatter)

    extract_req_parser = extract_parser.add_argument_group("required arguments")  
    extract_req_parser.add_argument('--gff', required=True, help='Input GFF annotations.')

    extract_out_opt_parser = extract_parser.add_argument_group("optional output arguments")
    extract_out_opt_parser.add_argument("--outdir",  help='Output directory. (default: .)', default=".")
    extract_out_opt_parser.add_argument('--prefix',  help='Output file prefix. If not provided, will be parsed from the gff file name.')

    extract_opt_parser = extract_parser.add_argument_group("optional arguments")
    extract_opt_parser.add_argument('--fasta',   help='Input FASTA sequences, if not provided at the end of the GFF.')
    extract_opt_parser.add_argument('--max-len', help='Maximum length of sequences to extract (default: 100000).', type=int, default=100000)
    extract_opt_parser.add_argument('--min-len', help='Minimum length of sequences to extract (default: 20).', type=int, default=20)
    extract_opt_parser.add_argument('--sample',  help='Sample identifier to use. If not provided, is parsed from the gff file name.')
    extract_opt_parser.add_argument('--regex',  help='Only extract gff lines that match this regular expression.')

    # -------------------------------------------------------------------------
    # Collect Sequences

    description = textwrap.dedent(
    f"""\
    {collect_description}

    Takes as input multiple TSV files from extract, which can be supplied 
    as either space separate paths, or a text file containing paths. 
    Duplicate sequence IDs will be identified and given the suffix '.#'.
    Outputs concatenated FASTA and TSV files.

    Examples:
    > pangwas collect --tsv sample1.tsv sample2.tsv sample3.tsv sample4.tsv
    > pangwas collect --tsv-paths tsv_paths.txt
    """)
    
    collect_parser = subbcommands.add_parser('collect', description = description, formatter_class=argparse.RawTextHelpFormatter)

    collect_req_parser = collect_parser.add_argument_group("required arguments (mutually-exclusive)")
    collect_input = collect_req_parser.add_mutually_exclusive_group(required=True)
    collect_input.add_argument('--tsv',       help='TSV files from the extract subcommand.', nargs='+')
    collect_input.add_argument('--tsv-paths', help='TXT file containing paths to TSV files.')

    collect_out_opt_parser = collect_parser.add_argument_group("optional output arguments")
    collect_out_opt_parser.add_argument("--outdir",  help='Output directory. (default: .)', default=".")
    collect_out_opt_parser.add_argument('--prefix', help='Prefix for output files.')

    # -------------------------------------------------------------------------
    # Cluster Sequences (mmseqs)

    description = textwrap.dedent(
    f"""\
    {cluster_description}

    Takes as input a FASTA file of sequences for clustering. Calls MMSeqs2 
    to cluster sequences and identify a representative sequence. Outputs a
    TSV table of sequence clusters and a FASTA of representative sequences.

    Any additional arguments will be passed to `mmseqs cluster`. If no additional
    arguments are used, the following default args will apply:
      {CLUSTER_ARGS}

    Examples:
    > pangwas cluster --fasta collect.fasta
    > pangwas cluster --fasta collect.fasta --threads 2 -k 13 --min-seq-id 0.90 -c 0.90
    """)

    cluster_parser = subbcommands.add_parser('cluster', description = description, formatter_class=argparse.RawTextHelpFormatter)

    cluster_req_parser = cluster_parser.add_argument_group("required arguments")  
    cluster_req_parser.add_argument('-f', '--fasta', required=True, help='FASTA file of input sequences to cluster.')

    cluster_out_opt_parser = cluster_parser.add_argument_group("optional output arguments")
    cluster_out_opt_parser.add_argument("--outdir",  help='Output directory. (default: .)', default=".")
    cluster_out_opt_parser.add_argument('--prefix',  help='Prefix for output files.')
    cluster_out_opt_parser.add_argument('--tmp',     help='Tmp directory (default: tmp).', default="tmp")

    cluster_opt_parser = cluster_parser.add_argument_group("optional arguments")
    cluster_opt_parser.add_argument('--memory',   help='Memory limit for mmseqs split. (default: 1G)', default="1G")
    cluster_opt_parser.add_argument('--no-clean', help="Don't clean up intermediate files.", action="store_false", dest="clean")
    cluster_opt_parser.add_argument('--threads',  help='CPU threads for mmseqs. (default: 1)', default=1)

    # -------------------------------------------------------------------------
    # Defrag clusters

    description = textwrap.dedent(
    f"""\
    {defrag_description}

    Takes as input the TSV clusters and FASTA representatives from cluster.
    Outputs a new cluster table and representative sequences fasta.

    {PPANGGOLIN_NOTICE}

    Any additional arguments will be passed to `mmseqs search`. If no additional
    arguments are used, the following default args will apply:
      {DEFRAG_ARGS}

    Examples:
    > pangwas defrag --clusters clusters.tsv --representative representative.fasta --prefix defrag
    > pangwas defrag --clusters clusters.tsv --representative representative.fasta --prefix defrag --threads 2 -k 13 --min-seq-id 0.90 -c 0.90
    """)

    defrag_parser = subbcommands.add_parser('defrag', description=description, formatter_class=argparse.RawTextHelpFormatter)
    defrag_req_parser = defrag_parser.add_argument_group("required arguments")  
    defrag_req_parser.add_argument('--clusters',             required=True, help='TSV file of clusters from mmseqs.')
    defrag_req_parser.add_argument('--representative',       required=True, help='FASTA file of cluster representative sequences.')

    defrag_out_opt_parser = defrag_parser.add_argument_group("optional output arguments")
    defrag_out_opt_parser.add_argument("--outdir",  help='Output directory. (default: .)', default=".")
    defrag_out_opt_parser.add_argument('--prefix',  help='Prefix for output files.')
    defrag_out_opt_parser.add_argument('--tmp',     help='Tmp directory (default: tmp).', default="tmp")

    defrag_opt_parser = defrag_parser.add_argument_group("optional arguments")
    defrag_opt_parser.add_argument('--memory',   help='Memory limit for mmseqs split. (default: 2G)', default="2G")
    defrag_opt_parser.add_argument('--no-clean', help="Don't clean up intermediate files.", action="store_false", dest="clean")
    defrag_opt_parser.add_argument('--threads',  help='CPU threads for mmseqs. (default: 1)', default=1)

    # -------------------------------------------------------------------------
    # Summarize clusters

    description = textwrap.dedent(
    f"""\
    {summarize_description}

    Takes as input the TSV table from collect, and the clusters table from 
    either cluster or defrag. Outputs a TSV table of summarized clusters 
    with their annotations.

    Examples:
    > pangwas summarize --clusters defrag.clusters.tsv --regions regions.tsv
    """)

    summarize_parser = subbcommands.add_parser('summarize', description=description, formatter_class=argparse.RawTextHelpFormatter) 
    summarize_req_parser = summarize_parser.add_argument_group("required arguments")    
    summarize_req_parser.add_argument('--clusters',  required=True, help='TSV file of clusters from cluster or defrag.')
    summarize_req_parser.add_argument('--regions',   required=True, help='TSV file of sequence regions from collect.')

    summarize_opt_parser = summarize_parser.add_argument_group("optional arguments")
    summarize_opt_parser.add_argument("--max-product-len", help='Truncate the product description to this length if used as an identifie. (default: 50)', type=int, default=50)
    summarize_opt_parser.add_argument("--min-samples",     help='Cluster must be observed in at least this many samples to be summarized.', type=int, default=1)
    summarize_opt_parser.add_argument("--outdir",          help='Output directory. (default: . )', default=".")
    summarize_opt_parser.add_argument('--prefix',          help='Prefix for output files.')
    summarize_opt_parser.add_argument('--threshold',       help='Required this proportion of samples to have annotations in agreement. (default: 0.5)', type=float, default=0.5)

    # -------------------------------------------------------------------------
    # Align clusters

    description = textwrap.dedent(
    f"""\
    {align_description}

    Takes as input the clusters from summarize and the sequence regions
    from collect. Outputs multiple sequence alignments per cluster 
    as well as a pangenome alignment of concatenated clusters.\n

    Any additional arguments will be passed to `mafft`. If no additional
    arguments are used, the following default args will apply:
      {ALIGN_ARGS}

    Examples:
    > pangwas align --clusters clusters.tsv --regions regions.tsv
    > pangwas align --clusters clusters.tsv --regions regions.tsv --threads 2 --exclude-singletons --localpair --maxiterate 100
    """)

    align_parser = subbcommands.add_parser('align', description=description, formatter_class=argparse.RawTextHelpFormatter) 

    align_req_parser = align_parser.add_argument_group("required arguments")
    align_req_parser.add_argument('--clusters',  help='TSV of clusters from summarize.', required=True)
    align_req_parser.add_argument('--regions', help='TSV of sequence regions from collect.', required=True)

    align_out_opt_parser = align_parser.add_argument_group("optional output arguments")
    align_out_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    align_out_opt_parser.add_argument('--prefix',   help='Prefix for output files.')

    align_opt_parser = align_parser.add_argument_group("optional arguments")
    align_opt_parser.add_argument('--threads',            help='CPU threads for running mafft in parallel. (default: 1)', type=int, default=1)
    align_opt_parser.add_argument('--exclude-singletons', help='Exclude clusters found in only 1 sample.', action="store_true")

    # -------------------------------------------------------------------------
    # Variants: Structural

    description = textwrap.dedent(
    f"""\
    {structural_description}

    Takes as input the summarized clusters and their individual alignments.
    Outputs an Rtab file of structural variants.

    Examples:
    > pangwas structural --clusters clusters.tsv --alignments alignments
    > pangwas structural --clusters clusters.tsv --alignments alignments --min-len 100 --min-indel-len 10
    """)

    structural_parser = subbcommands.add_parser('structural', description=description, formatter_class=argparse.RawTextHelpFormatter)

    structural_req_parser = structural_parser.add_argument_group("required arguments")
    structural_req_parser.add_argument('--clusters',  required=True, help='TSV of clusters from summarize.')
    structural_req_parser.add_argument('--alignments', required=True, help='Directory of cluster alignments (not consensus alignments!).') 
 
    structural_out_opt_parser = structural_parser.add_argument_group("optional output arguments")
    structural_out_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    structural_out_opt_parser.add_argument('--prefix',   help='Prefix for output files.')

    structural_opt_parser = structural_parser.add_argument_group("optional arguments")  
    structural_opt_parser.add_argument('--min-len', help='Minimum variant length. (default: 10)', type=int,   default=10)
    structural_opt_parser.add_argument('--min-indel-len', help='Minimum indel length. (default: 3)', type=int,   default=3)


    # -------------------------------------------------------------------------
    # Variants: SNPs

    description = textwrap.dedent(
    f"""\
    {snps_description}

    Takes as input the pangenome alignment fasta, bed, and consensus file from align.
    Outputs an Rtab file of SNPs.

    Examples:
    > pangwas snps --alignment pangenome.aln --bed pangenome.bed --consensus pangenome.consensus.fasta
    > pangwas snps --alignment pangenome.aln --bed pangenome.bed --consensus pangenome.consensus.fasta --structural structural.Rtab --core 0.90 --indel-window 3 --snp-window 10
    """)

    snps_parser = subbcommands.add_parser('snps', description=description, formatter_class=argparse.RawTextHelpFormatter) 

    snps_req_parser = snps_parser.add_argument_group("required arguments")
    snps_req_parser.add_argument('--alignment', required=True, help='Fasta sequences alignment.')
    snps_req_parser.add_argument('--bed',       required=True, help='Bed file of coordinates.')
    snps_req_parser.add_argument('--consensus', required=True, help='Fasta consensus/representative sequence.')

    snps_rec_parser = snps_parser.add_argument_group("optional but recommended arguments")
    snps_rec_parser.add_argument('--structural',   help='Rtab from the structural subcommand, used to avoid treating terminal ends as indels.')     

    snps_out_opt_parser = snps_parser.add_argument_group("optional output arguments")
    snps_out_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    snps_out_opt_parser.add_argument('--prefix',   help='Prefix for output files.')

    snps_opt_parser = snps_parser.add_argument_group("optional arguments")
    snps_opt_parser.add_argument('--core',         help='Core threshold for calling core SNPs. (default: 0.95)', type=float, default=0.95)
    snps_opt_parser.add_argument('--indel-window', help='Exclude SNPs that are within this proximity to indels. (default: 0)', type=int, default=0)
    snps_opt_parser.add_argument('--snp-window',   help='Exclude SNPs that are within this proximity to another SNP. (default: 0)', type=int, default=0)

    # -------------------------------------------------------------------------
    # Variants: Presence Absence

    description = textwrap.dedent(
    f"""\
    {presence_absence_description}

    Takes as input a TSV of summarized clusters from summarize.
    Outputs an Rtab file of cluster presence/absence.

    Examples:
    > pangwas presence_absence --clusters clusters.tsv
    """)

    pres_parser = subbcommands.add_parser('presence_absence', description=description, formatter_class=argparse.RawTextHelpFormatter) 
    pres_req_parser = pres_parser.add_argument_group("required arguments")
    pres_req_parser.add_argument('--clusters', help='TSV of clusters from summarize.', required=True)

    pres_out_opt_parser = pres_parser.add_argument_group("optional output arguments")
    pres_out_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    pres_out_opt_parser.add_argument('--prefix',   help='Prefix for output files.')

    # -------------------------------------------------------------------------
    # Variants: Combine

    description = textwrap.dedent(
    f"""\
    {combine_description}

    Takes as input a list of file paths to Rtab files. Outputs an Rtab file with
    the variants concatenated, ensuring consistent ordering of the sample columns.

    Examples:
    > pangwas combine --rtab snps.Rtab presence_absence.Rtab
    > pangwas combine --rtab snps.Rtab structural.Rtab presence_absence.Rtab
    """)
    combine_parser = subbcommands.add_parser('combine', description=description, formatter_class=argparse.RawTextHelpFormatter)

    combine_req_parser = combine_parser.add_argument_group("required arguments")    
    combine_req_parser.add_argument('--rtab', required=True, help="Rtab variants files.", nargs='+')

    combine_opt_parser = combine_parser.add_argument_group("optional arguments")
    combine_opt_parser.add_argument("--outdir", help='Output directory. (default: . )', default=".")
    combine_opt_parser.add_argument('--prefix', help='Prefix for output files.')
    
    # -------------------------------------------------------------------------
    # Variants: Table to Rtab

    description = textwrap.dedent(
        f"""\
        {table_to_rtab_description}

        Takes as input a TSV/CSV table to convert, and a TSV/CSV of regex filters.
        The filter table should have the header: column, regex, name. Where column
        is the 'column' to search, 'regex' is the regular expression pattern, and
        'name' is how the output variant should be named in the Rtab.

        An example `filter.tsv` might look like this:

        column    regex         name
        assembly  .*sample2.*   sample2
        lineage   .*2.*         lineage_2

        Where the goal is to filter the assembly and lineage columns for particular values.

        Examples:
        > pangwas table_to_rtab --table samplesheet.csv --filter filter.tsv
        """
    )

    table_to_rtab_parser = subbcommands.add_parser('table_to_rtab', description=description, formatter_class=argparse.RawTextHelpFormatter)
    table_to_rtab_req_parser = table_to_rtab_parser.add_argument_group("required arguments")
    table_to_rtab_req_parser.add_argument('--table',  required=True, help='TSV or CSV table.')
    table_to_rtab_req_parser.add_argument('--filter', required=True, help='TSV or CSV filter table.')

    table_to_rtab_opt_parser = table_to_rtab_parser.add_argument_group("optional arguments")
    table_to_rtab_opt_parser.add_argument("--outdir", help='Output directory. (default: . )', default=".")
    table_to_rtab_opt_parser.add_argument('--prefix', help='Prefix for output files.')

    # -------------------------------------------------------------------------
    # Variants: VCF to Rtab

    description = textwrap.dedent(
        f"""\
        {vcf_to_rtab_description}

        Takes as input a VCF file to convert to a SNPs Rtab file.

        Examples:
        > pangwas vcf_to_rtab --vcf snps.vcf
        """
    )

    vcf_to_rtab_parser = subbcommands.add_parser('vcf_to_rtab', description=description, formatter_class=argparse.RawTextHelpFormatter)
    vcf_to_rtab_req_parser = vcf_to_rtab_parser.add_argument_group("required arguments")
    vcf_to_rtab_req_parser.add_argument('--vcf',  required=True, help='VCF file.')

    vcf_to_rtab_opt_parser = vcf_to_rtab_parser.add_argument_group("optional arguments")
    vcf_to_rtab_opt_parser.add_argument("--bed", help='BED file with names by coordinates.')
    vcf_to_rtab_opt_parser.add_argument("--outdir", help='Output directory. (default: . )', default=".")
    vcf_to_rtab_opt_parser.add_argument('--prefix', help='Prefix for output files.')

    # -------------------------------------------------------------------------
    # Tree

    description = textwrap.dedent(
        f"""\
        {tree_description}

        Takes as input a multiple sequence alignment in FASTA format. If a SNP
        alignment is provided, an optional text file of constant sites can be 
        included for correction. Outputs a maximum-likelihood tree, as well as 
        additional rooted trees if an outgroup is specified in the iqtree args.

        Any additional arguments will be passed to `iqtree`. If no additional
        arguments are used, the following default args will apply:
        {TREE_ARGS}

        Examples:
        > pangwas tree --alignment snps.core.fasta --constant-sites snps.constant_sites.txt
        > pangwas tree --alignment pangenome.aln --threads 4 -o sample1 --ufboot 1000
        """
    )
    tree_parser = subbcommands.add_parser('tree', description=description, formatter_class=argparse.RawTextHelpFormatter)

    tree_req_parser = tree_parser.add_argument_group("required arguments")
    tree_req_parser.add_argument('--alignment', required=True, help='Multiple sequence alignment.')

    tree_opt_parser = tree_parser.add_argument_group("optional arguments")
    tree_opt_parser.add_argument('--constant-sites', help='Text file containing constant sites correction for SNP alignment.')
    tree_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    tree_opt_parser.add_argument('--prefix',   help='Prefix for output files.')     
    tree_opt_parser.add_argument('--threads',  help='CPU threads for IQTREE. (default: 1)', default=1)

   # -------------------------------------------------------------------------
   # Root tree

    description = textwrap.dedent(
        f"""\
        {root_tree_description}

        Takes as input a path to a phylogenetic tree and outgroup taxa. This is
        a utility script that is meant to fix IQ-TREE's creation of a multifurcated
        root node. It will position the root node along the midpoint of the branch
        between the outgroup taxa and all other samples. If no outgroup is selected,
        the tree will be rooted using the first taxa. Outputs a new tree in the specified
        tree format.

        Note: This functionality is already included in the tree subcommand.

        Examples:
        > pangwas root_tree --tree tree.treefile
        > pangwas root_tree --tree tree.treefile --outgroup sample1
        > pangwas root_tree --tree tree.treefile --outgroup sample1,sample4
        > pangwas root_tree --tree tree.nex --outgroup sample1 --tree-format nexus
        """
    )
    root_tree_parser = subbcommands.add_parser('root_tree', description=description, formatter_class=argparse.RawTextHelpFormatter)

    root_tree_req_parser = root_tree_parser.add_argument_group("required arguments")
    root_tree_req_parser.add_argument('--tree', help='Path to phylogenetic tree.')

    root_tree_opt_parser = root_tree_parser.add_argument_group("optional arguments")
    root_tree_opt_parser.add_argument("--outdir", help='Output directory. (default: . )', default=".")
    root_tree_opt_parser.add_argument("--outgroup", help='Outgroup taxa as CSV string. If not specified, roots on first taxon.')
    root_tree_opt_parser.add_argument('--prefix', help='Prefix for output files.')
    root_tree_opt_parser.add_argument('--tree-format',  help='Tree format. (default: newick)', type=str, default="newick")

    # -------------------------------------------------------------------------
    # GWAS: Binarize

    description = textwrap.dedent(
        f"""\
        {binarize_description}

        Takes as input a table (TSV or CSV) and the name of a column to binarize.
        Outputs a new table with separate columns for each categorical value.

        Any additional arguments will be passed to `pyseer`.

        Examples:
          pangwas binarize --table samplesheet.csv --column lineage --output lineage.binarize.csv
          pangwas binarize --table samplesheet.tsv --column outcome --output outcome.binarize.tsv
        """)

    binarize_parser = subbcommands.add_parser('binarize', description=description, formatter_class=argparse.RawTextHelpFormatter)
    binarize_req_parser = binarize_parser.add_argument_group("required arguments")
    binarize_req_parser.add_argument('--table',  required=True, help='TSV or CSV table.')
    binarize_req_parser.add_argument('--column', required=True, help='Column name to binarize.')

    binarize_opt_parser = binarize_parser.add_argument_group("optional arguments")
    binarize_opt_parser.add_argument('--column-prefix', help='Prefix to add to column names.')
    binarize_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    binarize_opt_parser.add_argument("--output-delim",   help='Output delimiter. (default: \t )', default="\t")
    binarize_opt_parser.add_argument('--prefix',   help='Prefix for output files')
    binarize_opt_parser.add_argument('--transpose', help='Tranpose output table.', action='store_true')

    # -------------------------------------------------------------------------
    # GWAS: pyseer

    description = textwrap.dedent(
        f"""\
        {gwas_description}

        Takes as input the TSV file of summarized clusters from summarize, an Rtab file of variants,
        a TSV/CSV table of phenotypes, and a column name representing the trait of interest.
        Outputs tables of locus effects, and optionally lineage effects (bugwas) if specified.

        Any additional arguments will be passed to `pyseer`. If no additional
        arguments are used, the following default args will apply:
        {GWAS_ARGS}

        Examples:
        > pangwas gwas --variants combine.Rtab --table samplesheet.csv --column lineage --no-distances
        > pangwas gwas --variants combine.Rtab --table samplesheet.csv --column resistant --lmm --tree tree.rooted.nwk --clusters clusters.tsv --lineage-column lineage '
        """)
    gwas_parser = subbcommands.add_parser('gwas', description=description, formatter_class=argparse.RawTextHelpFormatter)
    # Required options
    gwas_req_parser = gwas_parser.add_argument_group("required arguments")
    gwas_req_parser.add_argument('--variants', required=True,  help='Rtab file of variants.')
    gwas_req_parser.add_argument('--table',    required=True,  help='TSV or CSV table of phenotypes (variables).')
    gwas_req_parser.add_argument('--column',   required=True,  help='Column name of trait (variable) in table.')
    # Lineage effects options
    gwas_lin_parser = gwas_parser.add_argument_group("optional arguments (bugwas)")
    gwas_lin_parser.add_argument('--lineage-column', help='Column name of lineages in table. Enables bugwas lineage effects.')
    # LMM options
    gwas_lmm_parser = gwas_parser.add_argument_group("optional arguments (lmm)")
    gwas_lmm_parser.add_argument('--no-midpoint', help='Disable midpoint rooting of the tree for the kinship matrix.', dest="midpoint", action="store_false")
    gwas_lmm_parser.add_argument('--tree', help='Newick phylogenetic tree. Not required if pyseer args includes --no-distances.')
    # Data Options
    gwas_data_parser = gwas_parser.add_argument_group("optional arguments (data)")
    gwas_data_parser.add_argument('--continuous', help='Treat the column trait as a continuous variable.', action="store_true")
    gwas_data_parser.add_argument('--exclude-missing', help='Exclude samples missing phenotype data.', action="store_true")

    gwas_opt_parser = gwas_parser.add_argument_group("optional arguments")
    gwas_opt_parser.add_argument('--clusters', help='TSV of clusters from summarize.')
    gwas_opt_parser.add_argument("--outdir",   help='Output directory. (default: . )', default=".")
    gwas_opt_parser.add_argument('--prefix',   help='Prefix for output files.')
    gwas_opt_parser.add_argument('--threads',  help='CPU threads for pyseer. (default: 1)', default=1)

    # -------------------------------------------------------------------------
    # Plot

    description = textwrap.dedent(
        f"""\
        {heatmap_description}

        Takes as input a table of variants and/or a newick tree. The table can be either
        an Rtab file, or the locus effects TSV output from the gwas subcommand.
        If both a tree and a table are provided, the tree will determine the sample order 
        and arrangement. If just a table is provided, sample order will follow the 
        order of the sample columns. A TXT of focal sample IDs can also be supplied
        with one sample ID per line. Outputs a plot in SVG and PNG format.

        Examples:
        > pangwas heatmap --tree tree.rooted.nwk
        > pangwas heatmap --rtab combine.Rtab
        > pangwas heatmap --gwas resistant.locus_effects.significant.tsv
        > pangwas heatmap --tree tree.rooted.nwk --rtab combine.Rtab --focal focal.txt
        > pangwas heatmap --tree tree.rooted.nwk --gwas resistant.locus_effects.significant.tsv
        > pangwas heatmap --tree tree.rooted.nwk --tree-width 500 --png-scale 2.0
        """
    )

    heatmap_parser = subbcommands.add_parser('heatmap', description=description, formatter_class=argparse.RawTextHelpFormatter) 

    heatmap_req_parser = heatmap_parser.add_argument_group("optional variant arguments (mutually-exclusive)")
    heatmap_variants_input = heatmap_req_parser.add_mutually_exclusive_group(required=False)
    heatmap_variants_input.add_argument('--gwas', help='TSV table of variants from gwas subcommand.')
    heatmap_variants_input.add_argument('--rtab', help='Rtab table of variants.')

    heatmap_tree_parser = heatmap_parser.add_argument_group("optional tree arguments")
    heatmap_tree_parser.add_argument('--tree',  help='Tree file.')
    heatmap_tree_parser.add_argument('--tree-format', help='Tree format. (default: newick)', type=str, default="newick")
    heatmap_tree_parser.add_argument('--tree-width',  help='Width of the tree in pixels. (default: 200)', type=int, default=200)
    heatmap_tree_parser.add_argument('--root-branch', help='Root branch length in pixels. (default: 10)', type=int, default=10)

    heatmap_opt_parser = heatmap_parser.add_argument_group("optional arguments")
    heatmap_opt_parser.add_argument('--focal',         help='TXT file of focal samples.')
    heatmap_opt_parser.add_argument('--font-size',     help='Font size of the tree labels. (default: 16)', type=int, default=16)
    heatmap_opt_parser.add_argument('--font-family',   help='Font family of the tree labels. (default: Roboto)', type=str, default="Roboto")    
    heatmap_opt_parser.add_argument('--heatmap-scale', help='Scaling factor of heatmap boxes relative to the text. (default: 1.5)', type=float, default=1.5)    
    heatmap_opt_parser.add_argument('--margin',        help='Margin sizes in pixels. (default: 20)', type=int, default=20)
    heatmap_opt_parser.add_argument('--min-score',     help='Filter GWAS variants by a minimum score (range: -1.0 to 1.0).', type=float)
    heatmap_opt_parser.add_argument("--outdir",        help='Output directory. (default: . )', default=".")
    heatmap_opt_parser.add_argument('--png-scale',     help='Scaling factor of the PNG file. (default: 2.0)', type=float, default=2.0)
    heatmap_opt_parser.add_argument('--prefix',        help='Prefix for output files.')
    heatmap_opt_parser.add_argument('--tip-pad',       help='Tip padding. (default: 10)', type=int, default=10)

    # -------------------------------------------------------------------------
    # Manhattan

    description = textwrap.dedent(
        f"""\
        {manhattan_description}

        Takes as input a table of locus effects from the subcommand and a bed
        file such as the one producted by the align subcommand. Outputs a 
        manhattan plot in SVG and PNG format.

        Examples:
        > pangwas manhattan --gwas locus_effects.tsv --bed pangenome.bed
        > pangwas manhattan --gwas locus_effects.tsv --bed pangenome.bed --syntenies chromosome --clusters pbpX --variant-types snp presence_absence
        """
    )

    manhattan_parser = subbcommands.add_parser('manhattan', description=description, formatter_class=argparse.RawTextHelpFormatter) 

    manhattan_req_parser = manhattan_parser.add_argument_group("required arguments")
    manhattan_req_parser.add_argument('--bed',  required=True, help='BED file with region coordinates and names.')    
    manhattan_req_parser.add_argument('--gwas', required=True, help='TSV table of variants from gwas subcommand.')

    manhattan_opt_parser = manhattan_parser.add_argument_group("optional filter arguments")
    manhattan_opt_parser.add_argument('--clusters',      help='Only plot these clusters. (default: all)', nargs='+', default=["all"])
    manhattan_opt_parser.add_argument('--syntenies',     help='Only plot these synteny blocks. (default: all)', nargs='+', default=["all"])
    manhattan_opt_parser.add_argument('--variant-types', help='Only plot these variant types. (default: all)', nargs='+', default=["all"])

    manhattan_opt_parser = manhattan_parser.add_argument_group("optional arguments")
    manhattan_opt_parser.add_argument('--font-size',   help='Font size of the tree labels. (default: 16)', type=int, default=16)
    manhattan_opt_parser.add_argument('--font-family', help='Font family of the tree labels. (default: Roboto)', type=str, default="Roboto")        
    manhattan_opt_parser.add_argument('--height',      help='Plot height in pixels.', type=int, default=500)    
    manhattan_opt_parser.add_argument('--margin',      help='Margin sizes in pixels. (default: 20)', type=int, default=20)
    manhattan_opt_parser.add_argument('--max-blocks',  help='Maximum number of blocks to draw before switching to pangenome coordinates. (default: 20)', type=int, default=20)
    manhattan_opt_parser.add_argument("--outdir",      help='Output directory. (default: . )', default=".")
    manhattan_opt_parser.add_argument('--png-scale',   help='Scaling factor of the PNG file. (default: 2.0)', type=float, default=2.0)    
    manhattan_opt_parser.add_argument('--prefix',      help='Prefix for output files.')
    manhattan_opt_parser.add_argument('--prop-x-axis', help='Make x-axis proporional to genomic length.', action="store_true")
    manhattan_opt_parser.add_argument('--width',       help='Plot width in pixels.', type=int, default=1000)
    manhattan_opt_parser.add_argument('--ymax',        help='A -log10 value to use as the y axis max, to synchronize multiple plots.', type=float)

    # -------------------------------------------------------------------------
    # Finalize

    defined, undefined = parser.parse_known_args()
    if len(undefined) > 0:
        defined.args = " ".join(undefined)
        logging.warning(f"Using undefined args: {defined.args}")

    sys.argv = sys_argv_original
    return defined


def check_output_dir(output_dir):
    if output_dir != None and output_dir != "." and output_dir != "" and not os.path.exists(output_dir):
        logging.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)


def reverse_complement(sequence):
    lookup = {"A": "T", "C": "G", "G": "C", "T": "A"}
    split = []
    for nuc in sequence[::-1]:
        nuc = nuc.upper()
        compl = lookup[nuc] if nuc in lookup else nuc
        split.append(compl)
    return "".join(split)


def extract_cli_param(args: str, target: str):
    """
    Extract a value flag from a string of CLI arguments.

    :param args: String of CLI arguments (ex. 'iqtree -T 3 --seed 123').
    :param args: String of target param (ex. '-T', '--seed')
    """
    if not target.endswith(" "):
        target += " "
    val = None
    if target in args:
        start = args.index(target) + len(target)
        try:
            end = args[start:].index(' ')
            val = args[start:][:end]
        except ValueError:
            val = args[start:]
    return val


def run_cmd(cmd: str, output=None, err=None, quiet=False, display_cmd=True):
    """
    Run a subprocess command.
    """
    if display_cmd:
        logging.info(f"{cmd}")
    cmd_run = [str(c) for c in cmd.replace("\n", " ").split(" ") if c != ""]
    try:
        result = subprocess.run(cmd_run, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as error:
        logging.error(f"stdout:\n{error.stdout}")
        logging.error(f"stderr:\n{error.stderr}")
        raise Exception(f"Error in command: {cmd}\n" + error.stderr)
    else:
        if not quiet:
            logging.info(result.stdout)
        if output:
            with open(output, 'w') as outfile:
                outfile.write(result.stdout)
        if err:
            with open(err, 'w') as errfile:
                errfile.write(result.stderr)
        return (result.stdout, result.stderr)

def annotate(
    fasta:   str,
    db:      str,
    outdir:  str = ".",
    prefix:  str = None,
    sample:  str = None,
    threads: int = 1,
    tmp:     str = None,
    args:    str = None
    ):
    """
    Annotate genomic assemblies with bakta.

    Takes as input a FASTA file of genomic assemblies. Outputs a GFF file
    of annotations, among many other formats from bakta.

    Any additional arguments in `args` will be passed to `bakta`.

    >>> annotate(fasta="sample1.fasta", db="database/bakta")
    >>> annotate(fasta="sample2.fasta", db="database/bakta", threads=2, args="--genus Streptococcus")

    :param fasta:   File path to the input FASTA file.
    :param db:      Directory path of the bakta database.
    :param outdir:  Output directory.
    :param prefix:  Prefix for output files.   
    :param sample:  Sample identifier for output files.
    :param threads: CPU threads for bakta.
    :param tmp:     Temporary directory.
    :param args:    Str of additional arguments to pass to bakta

    :return:        Path to the output GFF annotations.
    """
    from bakta.main import main as bakta

    output_dir, tmp_dir = outdir, tmp

    # Argument Checking
    args = args if args != None else ""
    sample = sample if sample else os.path.basename(fasta).split(".")[0]
    logging.info(f"Using sample identifier: {sample}")
    prefix = prefix if prefix != None else os.path.basename(fasta).split(".")[0]
    logging.info(f"Using prefix: {prefix}")
    if tmp != None:
        args += f" --tmp-dir {tmp}"

    # bakta needs the temporary directory to already exist
    # bakta explicity wants the output directory to NOT exist
    # so that is why we don't create it
    check_output_dir(tmp_dir)

    # Set up cmd to execute, we're not going to use the run_cmd, because
    # bakta is a python package, and we call it directly! This also has the 
    # benefit of better/immediate logging.
    cmd = f"bakta --force --prefix {prefix} --db {db} --threads {threads} --output {output_dir} {args} {fasta}"
    logging.info(f"{cmd}")

    # Split cmd string and handle problematic whitespace
    cmd_run = [str(c) for c in cmd.replace("\n", " ").split(" ") if c != ""]
    # Manually set sys.argv for bakta
    sys_argv_original = sys.argv
    sys.argv = cmd_run
    bakta()

    # Fix bad quotations
    for ext in ["embl", "faa", "ffn", "gbff", "gff3", "tsv"]:
        file_path = os.path.join(output_dir, f"{prefix}.{ext}")
        logging.info(f"Fixing bad quotations in file: {file_path}")
        with open(file_path) as infile:
            content = infile.read()
        with open(file_path, 'w') as outfile:
            outfile.write(content.replace('‘', '').replace('’', ''))
    
    # Restore original sys.argv
    sys.argv = sys_argv_original

    gff_path = os.path.join(output_dir, f"{sample}.gff3")
    return gff_path


def extract(
        gff:     str,
        sample:  str = None,
        fasta:   str = None,
        outdir:  str = ".",
        prefix:  str = None,
        min_len: int = 20,
        max_len: int = 100000,
        regex: str = None,
        args:    str = None
    ) -> OrderedDict:
    """
    Extract sequence records from GFF annotations.

    Takes as input a GFF annotations file. If sequences are not included, a FASTA
    of genomic contigs must also be provided. Both annotated and unannotated regions 
    will be extracted. Outputs a TSV table of extracted sequence regions.

    >>> extract(gff='sample1.gff3')
    >>> extract(gff='sample2.gff3', fasta="sample2.fasta", min_len=10)

    :param gff:      File path to the input GFF file.
    :param sample:   A sample identifier that will be used in output files.
    :param fasta:    File path to the input FASTA file (if GFF doesn't have sequences).
    :param outdir:   Output directory.
    :param prefix:   Output file path prefix.
    :param min_len:  The minimum length of sequences/annotations to extract.
    :param regex:    Only extract lines that match this regex.
    :param args:     Str of additional arguments [not implemented]

    :return: Sequence records and annotations as an OrderedDict.
    """

    from Bio import SeqIO
    from io import StringIO
    import re

    gff_path, fasta_path, output_dir = gff, fasta, outdir

    # Argument checking
    args = args if args != None else ""    
    sample = sample if sample != None else os.path.basename(gff).split(".")[0]
    logging.info(f"Using sample identifier: {sample}")
    prefix = f"{prefix}" if prefix != None else os.path.basename(gff).split(".")[0]
    logging.info(f"Using prefix: {prefix}")
    if "/" in prefix:
        msg = "Prefix cannot contain slashes (/)"
        logging.error(msg)
        raise Exception(msg)
    min_len = int(min_len)
    max_len = int(max_len) if max_len != None else None

    # Check output directory
    check_output_dir(output_dir)

    # A visual separator for creating locus names for unannotated regions
    # ex. 'SAMPL01_00035___SAMPL01_00040' indicating the region
    # falls between two annotated regions: SAMPL01_00035 and SAMPL01_00040
    delim = "__"

    # -------------------------------------------------------------------------
    # Extract contig sequences

    logging.info(f"Reading GFF: {gff_path}")

    contigs = OrderedDict()
    with open(gff) as infile:
        gff = infile.read()
    if "##FASTA" in gff:
        logging.info(f"Extracting sequences from GFF.")
        fasta_i = gff.index("##FASTA")
        fasta_string = gff[fasta_i + len("##FASTA"):].strip()
        fasta_io = StringIO(fasta_string) 
        records = SeqIO.parse(fasta_io, "fasta")
    elif fasta_path == None:
        msg = f"A FASTA file must be provided if no sequences are in the GFF: {gff_path}"
        logging.error(msg)
        raise Exception(msg)

    if fasta_path != None:
        logging.info(f"Extracting sequences from FASTA: {fasta_path}")
        if "##FASTA" in gff:
            logging.warning(f"Sequences found in GFF, fasta will be ignored: {fasta}")
        else: 
            records = SeqIO.parse(fasta_path, "fasta")
            fasta_i = None

    sequences_seen = set()
    for record in records:
        if record.id in sequences_seen:
            msg = f"Duplicate sequence ID found: {record.id}"
            logging.error(msg)
            raise Exception(msg)
        sequences_seen.add(record.id)
        contigs[record.id] = record.seq

    # -------------------------------------------------------------------------
    # Extract annotations
    
    logging.info(f"Extracting annotations from gff: {gff_path}")

    # Parsing GFF is not yet part of biopython, so this is done manually:
    # https://biopython.org/wiki/GFF_Parsing

    annotations = OrderedDict()
    contig = None
    sequence = ""

    # Keep track of loci, in case we need to flag duplicates
    locus_counts = {}
    gene_counts = {}
    locus_to_contig = {}

    comment_contigs = set()

    for i,line in enumerate(gff.split("\n")):
        # If we hit the fasta, we've seen all annotations
        if line.startswith("##FASTA"): break
        if line.startswith("##sequence-region"):
            _comment, contig, start, end = line.split(" ")
            start, end = int(start), int(end)
            if contig not in annotations:
                annotations[contig] = {"start": start, "end": end, "length": 1 + end - start, "loci": OrderedDict()}
                comment_contigs.add(contig)
            continue
        # Skip over all other comments
        if line.startswith("#"): continue
        # Skip over empty lines
        if line == "": continue

        # Parse standardized annotation fields
        line = line.replace("\"", "")
        line_split = line.split("\t")
        if len(line_split) < 9:
            msg = f"GFF record does not contain 9 fields: {line}"
            logging.error(msg)
            raise Exception(msg)
        contig, _source, feature, start, end, _score, strand, _frame = [line_split[i] for i in range(0,8)]
        start, end = int(start), int(end)
        if feature == "region" or feature == "databank_entry":
            if contig in comment_contigs:
                continue
            elif contig in annotations:
                msg = f"Duplicate contig ID found: {contig}"
                logging.error(msg)
                raise Exception(msg)
            annotations[contig] = {"start": start, "end": end, "length": 1 + end - start, "loci": OrderedDict()}
            continue
        elif regex != None and not re.match(regex, line, re.IGNORECASE):
            continue

        attributes = line_split[8]
        attributes_dict = {a.split("=")[0].replace(" ", ""):a.split("=")[1] for a in attributes.split(";") if "=" in a}
        locus = f"{attributes_dict['ID']}"
        if "gene" in attributes_dict:
            gene = attributes_dict["gene"]
            if gene not in gene_counts:
                gene_counts[gene] = {"count": 0, "current": 1}
            gene_counts[gene]["count"] += 1

        # Add sample prefix to help duplicate IDs later
        if not locus.startswith(sample):
            locus = f"{sample}_{locus}"

        # Check for duplicate loci, how does this happen in NCBI annotations?
        if locus not in locus_counts:
            locus_counts[locus] = 1
        else:
            locus_counts[locus] += 1
            dup_i = locus_counts[locus]
            locus = f"{locus}.{dup_i}"
            logging.debug(f"Duplicate locus ID found, flagging as: {locus}")                        

        if contig not in contigs:
            msg = f"Contig {contig} in {sample} annotations is not present in the sequences."
            logging.error(msg)
            raise Exception(msg)
        sequence = contigs[contig][start - 1:end]
        if strand == "-":
            sequence = reverse_complement(sequence)

        data = OrderedDict({
            "sample"      : sample,
            "contig"      : contig,                        
            "locus"       : locus,
            "feature"     : feature,
            "start"       : start, 
            "end"         : end,
            "length"      : 1 + end - start,
            "strand"      : strand,
            "upstream"    : "",
            "downstream"  : "",
            "attributes"  : attributes,
            "sequence_id" : f"{locus}",
            "sequence"    : sequence,
        })
        logging.debug(f"\tcontig={contig}, locus={locus}")
        annotations[contig]["loci"][locus] = data
        locus_to_contig[locus] = contig

    # Check for duplicates, rename the first duplicate loci
    logging.info(f"Checking for duplicate locus IDs.")
    for locus,count in locus_counts.items():
        if count <= 1: continue
        contig = locus_to_contig[locus]
        data = annotations[contig]["loci"][locus]
        new_locus = f"{locus}.1"
        logging.debug(f"Duplicate locus ID found, flagging as: {new_locus}")
        data["locus"], data["sequence_id"] = new_locus, new_locus
        annotations[contig]["loci"] = OrderedDict({
            new_locus if k == locus else k:v
            for k,v in annotations[contig]["loci"].items()
        })

    # -------------------------------------------------------------------------
    # Extract unannotated regions

    if regex != None:
        logging.info("Skipping extraction of unannotate regions due to regex.")
    else:
        logging.info(f"Extracting unannotated regions.")
        for contig,contig_data in annotations.items():
            if contig not in contigs:
                msg = f"Contig {contig} in {sample} annotations is not present in the sequences."
                logging.error(msg)
                raise Exception(msg)
            contig_sequence = contigs[contig]
            num_annotations = len([k for k,v in contig_data["loci"].items() if v["feature"] != "region"])
            contig_start    = contig_data["start"]
            contig_end      = contig_data["end"]

            # Contig minimum length check
            contig_len = (1 + contig_end - contig_start)
            if contig_len < min_len: continue

            logging.debug(f"\tcontig={contig}")

            # If there were no annotations on the contig, extract entire contig as sequence
            if num_annotations == 0:
                logging.debug(f"\t\tlocus={contig}")
                locus = f"{contig}"
                if not locus.startswith(sample):
                    locus = f"{sample}_{locus}"
                l_data = {
                    "sample"     : sample,
                    "contig"     : contig,
                    "locus"      : locus,
                    "feature"    : "unannotated",
                    "start"      : contig_start, 
                    "end"        : contig_end,
                    "length"     : 1 + contig_end - contig_start,
                    "strand"     : "+",
                    "upstream"   : f"{contig}_TERMINAL",
                    "downstream" : f"{contig}_TERMINAL",
                    "attributes" : "",
                    "sequence_id" : locus,
                    "sequence"   : contig_sequence,
                }
                annotations[contig]["loci"][contig] = l_data
                logging.debug(f"\t\t\tfull contig unannotated, locus={locus}")
                continue

            contig_annotations = list(contig_data["loci"].keys())

            for i,locus in enumerate(contig_annotations):
                logging.debug(f"\t\tlocus={locus}")
                l_data = contig_data["loci"][locus]
                l_start, l_end = l_data["start"], l_data["end"]

                # Make a template for an unannotated region
                l_data = {
                        "sample"      : sample,
                        "contig"      : contig,
                        "locus"       : None,
                        "feature"     : "unannotated",
                        "start"       : None, 
                        "end"         : None,
                        "length"      : 0,
                        "strand"      : "+",
                        "upstream"    : "",
                        "downstream"  : "",
                        "attributes"  : "",
                        "sequence_id" : None,
                        "sequence"    : None,
                }
                # Find the inter-genic regions and upstream/downstream loci
                start, end, sequence, upstream, downstream = None, None, None, None, None

                # Case 1. Unannotated at the start of the contig
                if i == 0 and l_start != contig_start:
                    start = 1
                    end = l_start - 1
                    length = 1 + end - start
                    if length >= min_len and (max_len == None or length <= max_len):
                        upstream = f"{contig}_TERMINAL"
                        downstream = locus
                        # Base strand on the downstream loci
                        strand = annotations[contig]["loci"][downstream]["strand"]      
                        sequence_id = f"{upstream}{delim}{downstream}"
                        sequence = contig_sequence[start - 1:end]
                        # If the downstream is reverse, treat unannotated as reversed
                        if strand == "-":
                            sequence = reverse_complement(sequence)
                        l_data_start = copy.deepcopy(l_data)
                        for k,v in zip(["start", "end", "length", "sequence_id", "sequence", "locus", "strand"], [start, end, length, sequence_id, sequence, sequence_id, strand]):
                            l_data_start[k] = v
                        l_data_start["upstream"], l_data_start["downstream"] = upstream, downstream
                        logging.debug(f"\t\t\tunannotated at start, contig={contig}, sequence_id={sequence_id}, start: {start}, end: {end}")
                        annotations[contig]["loci"][sequence_id] = l_data_start

                        # Update upstream for annotation
                        annotations[contig]["loci"][locus]["upstream"] = sequence_id

                # Case 2. Unannotated at the end of the contig
                if i == (num_annotations - 1) and l_end != contig_end:
                    start = l_end + 1
                    end = contig_end
                    length = 1 + end - start
                    if length >= min_len and (max_len == None or length <= max_len):
                        upstream = locus
                        # base strand on the upstream loci
                        strand = annotations[contig]["loci"][upstream]["strand"]
                        downstream = f"{contig}_TERMINAL"
                        sequence_id = f"{upstream}{delim}{downstream}"
                        sequence = contig_sequence[start - 1:end]
                        # If the upstream is reversed, treat unannotated as reversed
                        if strand == "-":
                            sequence = reverse_complement(sequence)
                        l_data_end = copy.deepcopy(l_data) 
                        for k,v in zip(["start", "end", "length", "sequence_id", "sequence", "locus", "strand"], [start, end, length, sequence_id, sequence, sequence_id, strand]):
                            l_data_end[k] = v
                        l_data_end["upstream"], l_data_end["downstream"] = upstream, downstream
                        logging.debug(f"\t\t\tunannotated at end, contig={contig}, sequence_id={sequence_id}, start: {start}, end: {end}")          
                        annotations[contig]["loci"][sequence_id] = l_data_end

                        # Update downstream for annotation
                        annotations[contig]["loci"][locus]["downstream"] = sequence_id                        

                # Case 3. Unannotated in between annotations
                if num_annotations > 1 and i != (num_annotations - 1):

                    upstream = locus
                    downstream = contig_annotations[i+1]                 
                    start = l_end + 1
                    end = contig_data["loci"][downstream]["start"] - 1
                    length = 1 + end - start
                    # Set upstream downstream based on order in GFF
                    upstream_strand = annotations[contig]["loci"][upstream]["strand"]
                    downstream_strand = annotations[contig]["loci"][downstream]["strand"]      

                    # Check that the region is long enough
                    if length >= min_len and (max_len == None or length <= max_len):
                        sequence_id = f"{upstream}{delim}{downstream}"
                        sequence = contig_sequence[start - 1:end]
                        # Identify strand and orientation                      
                        if upstream_strand == "-" and downstream_strand == "-":
                            strand = "-"
                            sequence = reverse_complement(sequence)
                        else:
                            strand = "+"
                        l_data_middle = copy.deepcopy(l_data)
                        for k,v in zip(["start", "end", "length", "sequence_id", "sequence", "locus", "strand"], [start, end, length, sequence_id, sequence, sequence_id, strand]):
                            l_data_middle[k] = v
                        l_data_middle["upstream"], l_data_middle["downstream"] = upstream, downstream
                        logging.debug(f"\t\t\tunannotated in middle, contig={contig}, sequence_id={sequence_id}, start: {start}, end: {end}")              
                        annotations[contig]["loci"][sequence_id] = l_data_middle

    # -------------------------------------------------------------------------
    # Order and Filter

    logging.info(f"Ordering records by contig and coordinate.")
    contig_order = list(annotations.keys())
    for contig in contig_order:
        loci = annotations[contig]["loci"]
        # sort by start and end position
        loci_sorted = OrderedDict(sorted(loci.items(), key=lambda item: [item[1]["start"], item[1]["end"]]) )
        annotations[contig]["loci"] = OrderedDict()
        for locus,locus_data in loci_sorted.items():          
            # final checks, minimum length and not entirely ambiguous characters
            non_ambig = list(set([n for n in locus_data["sequence"] if n in NUCLEOTIDES]))
            length = locus_data["length"]
            if (length >= min_len) and (len(non_ambig) > 0):
                if max_len == None or length <= max_len:
                    annotations[contig]["loci"][locus] = locus_data


    # -------------------------------------------------------------------------
    # Upstream/downtream loci

    for contig in annotations:
        contig_data = annotations[contig]
        c_start, c_end = contig_data["start"], contig_data["end"]
        loci = list(contig_data["loci"].keys())
        for i,locus in enumerate(contig_data["loci"]):
            upstream, downstream = "", ""
            locus_data = contig_data["loci"][locus]
            l_start, l_end = locus_data["start"], locus_data["end"]
            # Annotation at very start
            if i > 0:
                upstream = loci[i-1]
            elif l_start == c_start:
                upstream = f"{contig}_TERMINAL"
            # Annotation at very end
            if i < (len(loci) - 1):
                downstream = loci[i+1]
            elif l_end == c_end:
                downstream = f"{contig}_TERMINAL"
            annotations[contig]["loci"][locus]["upstream"] = upstream
            annotations[contig]["loci"][locus]["downstream"] = downstream

    # -------------------------------------------------------------------------
    # Write Output

    tsv_path = os.path.join(output_dir, prefix + ".tsv")
    logging.info(f"Writing output tsv: {tsv_path}")
    with open(tsv_path, 'w') as outfile:
        header = None
        for contig,contig_data in annotations.items():
            for locus,locus_data in contig_data["loci"].items():
                if locus_data["feature"] == "region": continue                   
                if header == None:
                    header = list(locus_data.keys())
                    outfile.write("\t".join(header) + "\n")
                row = [str(v) for v in locus_data.values()]
                # Restore commas from there %2C encoding
                line = "\t".join(row).replace("%2C", ",")
                outfile.write(line + "\n")   

    return tsv_path


def collect(
        tsv:         list = None,
        tsv_paths:   str  = None,
        outdir :     str  = ".",
        prefix :     str  = None,
        args:        str  = None,
    ) -> str:
    """
    Collect sequences from multiple samples into one file.

    Takes as input multiple TSV files from the extract subcommand, which can 
    be supplied as either space separate paths, or a text file containing paths. 
    Duplicate sequence IDs will be identified and given the suffix '.#'.
    Outputs concatenated FASTA and TSV files.

    >>> collect(tsv=["sample1.tsv", "sample2.tsv"])
    >>> collect(tsv_paths='paths.txt')

    :param tsv:       List of TSV file paths from the output of extract.
    :param tsv_paths: TXT file with each line containing a TSV file path.
    :param output:    Path to the output directory.
    :param prefix:    Output file prefix.
    :param args:      Additional arguments [not implemented]

    :return: Tuple of output FASTA and TSV paths.
    """

    from Bio import SeqIO

    tsv_paths, tsv_txt_path, output_dir = tsv, tsv_paths, outdir
    prefix = f"{prefix}." if prefix != None else ""

    # Check output directory
    check_output_dir(output_dir)

    # TSV file paths
    tsv_file_paths = []
    if tsv_txt_path != None:
        logging.info(f"Reading tsv paths from file: {tsv_txt_path}")
        with open(tsv_txt_path) as infile:
            for line in infile:
                file_path = [l.strip() for l in line.split("\t")][0]
                tsv_file_paths.append(file_path.strip())
    if tsv_paths != None:
        for file_path in tsv_paths:
            tsv_file_paths.append(file_path)

    logging.info(f"Checking for duplicate samples and sequence IDs.")
    all_samples = []
    sequence_id_counts = {}
    for file_path in tqdm(tsv_file_paths):
        with open(file_path) as infile:
            header = [l.strip() for l in infile.readline().split("\t")]
            for i,line in enumerate(infile):
                row = [l.strip() for l in line.split("\t")]
                data = {k:v for k,v in zip(header, row)}
                
                # Check for duplicate samples
                if i == 0:
                    sample = data["sample"]
                    if sample in all_samples:
                        msg = f"Duplicate sample ID found: {sample}"
                        logging.error(msg)
                        raise Exception(msg)
                    all_samples.append(sample)
    
                # Check for duplicate IDs
                sequence_id = data["sequence_id"]
                if sequence_id not in sequence_id_counts:
                    sequence_id_counts[sequence_id] = {"count": 0, "current": 1}
                sequence_id_counts[sequence_id]["count"] += 1

    output_fasta = os.path.join(outdir, f"{prefix}sequences.fasta")
    logging.info(f"Writing output fasta file to: {output_fasta}")
    output_tsv = os.path.join(outdir, f"{prefix}regions.tsv")
    logging.info(f"Writing output tsv file to: {output_tsv}")

    fasta_outfile = open(output_fasta, 'w')
    tsv_outfile = open(output_tsv, "w")
    header = None

    for file_path in tqdm(tsv_file_paths):
        with open(file_path) as infile:
            header_line = infile.readline()
            if header == None:
                header = [l.strip() for l in header_line.split("\t")]
                tsv_outfile.write(header_line)
            for line in infile:
                row = [l.strip() for l in line.split("\t")]
                data = {k:v for k,v in zip(header, row)}

                # Handle duplicate sequence IDs with a .# suffix        
                sequence_id =           data["sequence_id"]
                count = sequence_id_counts[sequence_id]["count"]
                if count > 1:
                    dup_i = sequence_id_counts[sequence_id]["current"]
                    sequence_id_counts[sequence_id]["current"] += 1
                    data["sequence_id"] = sequence_id = f"{sequence_id}.{dup_i}"
                    logging.debug(f"Duplicate sequence ID found in TSV, flagging as: {sequence_id}")

                line = "\t".join([str(data[col]) for col in header])
                tsv_outfile.write(line + "\n")

                sequence = data["sequence"]
                line = f">{sequence_id}\n{sequence}"
                fasta_outfile.write(line + "\n")

    fasta_outfile.close()
    tsv_outfile.close()

    return (output_fasta, output_tsv)


def cluster(
        fasta: str,
        outdir:  str  = ".",
        prefix:  str  = None,
        threads: int  = 1,
        memory:  str  = "1G",
        tmp:     str  = "tmp",
        clean:   bool = True,
        args:    str  = CLUSTER_ARGS,
    ):
    """
    Cluster nucleotide sequences with mmseqs.

    Takes as input a FASTA file of sequences for clustering from collect.
    Calls MMSeqs2 to cluster sequences and identify a representative sequence. 
    Outputs a TSV table of sequence clusters and a FASTA of representative sequences.

    Note: The default kmer size (15) requires at least 9G of memory. To use less
    memory, please set the kmer size to '-k 13'.

    >>> cluster(fasta='sequences.fasta')
    >>> cluster(fasta='sequences.fasta', threads=2, memory='2G', args='-k 13 --min-seq-id 0.90 -c 0.90')
    >>> cluster(fasta='sequences.fasta', threads=4, args='-k 13 --min-seq-id 0.90 -c 0.90')

    :param fasta:   Path fo FASTA sequences.
    :param outdir:  Output directory.
    :param prefix:  Prefix for output files.
    :param threads: CPU threads for MMSeqs2.
    :param memory:  Memory for MMSeqs2.
    :param tmp:     Path to a temporary directory.
    :param clean:   True if intermediate files should be cleaned up.
    :param args:    Additional parameters for MMSeqs2 cluster command. 
    """

    from Bio import SeqIO

    fasta_path, tmp_dir, output_dir = fasta, tmp, outdir
    prefix = f"{prefix}." if prefix != None else ""

    # Wrangle the output directory
    check_output_dir(output_dir)
    check_output_dir(tmp_dir)

    args = args if args != None else ""

    # fix memory formatting (ex. '6 GB' -> '6G')
    memory = memory.replace(" ", "").replace("B", "")

    # 1. Cluster Sequences
    seq_db = os.path.join(output_dir, f"{prefix}seqDB")
    clust_db = os.path.join(output_dir, f"{prefix}clustDB")
    tsv_path = os.path.join(output_dir, f"{prefix}clusters.tsv")
    run_cmd(f"mmseqs createdb {fasta_path} {seq_db}")
    try:
        run_cmd(f"mmseqs cluster {seq_db} {clust_db} {tmp_dir} --threads {threads} --split-memory-limit {memory} {args}")
    except Exception as e:
        if "Segmentation fault" in f"{e}":
            logging.error(f"Segmentation fault. Try changing your threads and/or memory. Memory consumption is {threads} x {memory}.")
        raise Exception(e)
    run_cmd(f"mmseqs createtsv {seq_db} {seq_db} {clust_db} {tsv_path} --threads {threads} --full-header")

    # Sort and remove quotations in clusters
    with open(tsv_path) as infile:
        lines = sorted([l.strip().replace("\"", "") for l in infile.readlines()])
    with open(tsv_path, 'w') as outfile:
        outfile.write("\n".join(lines) + "\n")

    # 2. Identify representative sequences
    rep_db = os.path.join(output_dir, f"{prefix}repDB")
    rep_fasta_path = os.path.join(output_dir, f"{prefix}representative.fasta")
    run_cmd(f"mmseqs result2repseq {seq_db} {clust_db} {rep_db} --threads {threads}")
    run_cmd(f"mmseqs result2flat {seq_db} {seq_db} {rep_db} {rep_fasta_path} --use-fasta-header")

    # Sort sequences, and remove trailing whitespace
    sequences = {}
    for record in SeqIO.parse(rep_fasta_path, "fasta"):
        locus = record.id.strip()
        seq = record.seq.strip()
        sequences[locus] = seq
    with open(rep_fasta_path, 'w') as outfile:
        for locus in sorted(list(sequences.keys())):
            seq = sequences[locus]
            outfile.write(f">{locus}\n{seq}\n")

    # Sort clusters
    with open(tsv_path) as infile:
        lines = sorted([line.strip() for line in infile])
    with open(tsv_path, 'w') as outfile:
        outfile.write("\n".join(lines))

    # Cleanup
    if clean == True:    
        for file_name in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file_name)
            db_prefix = os.path.join(output_dir, prefix)
            if (
                file_path.startswith(f"{db_prefix}seqDB") or 
                file_path.startswith(f"{db_prefix}clustDB") or
                file_path.startswith(f"{db_prefix}repDB")
            ):
                logging.info(f"Cleaning up file: {file_path}")
                os.remove(file_path)

    return(tsv_path, rep_fasta_path)


def defrag(
        clusters:       str,
        representative: str,
        outdir:         str = ".",
        prefix:         str = None,
        tmp:            str = "tmp",
        threads:        str = 1,
        memory:         str = "2G",
        clean:          str = True,
        args:           str = DEFRAG_ARGS,
    ) -> OrderedDict:
    """
    Defrag clusters by associating fragments with their parent cluster.

    Takes as input the TSV clusters and FASTA representatives from the cluster subcommand.
    Outputs a new cluster table and representative sequences fasta.

    This is a modification of ppanggolin's refine_clustering function:
    https://github.com/labgem/PPanGGOLiN/blob/2.2.0/ppanggolin/cluster/cluster.py#L317

    >>> defrag(clusters='clusters.tsv', representative='representative.fasta', prefix="defrag")
    >>> defrag(clusters='clusters.tsv', representative='representative.fasta', threads=2, memory='2G', args="-k 13 --min-seq-id 0.90 -c 0.90 --cov-mode 1")

    :param clusters:       TSV file of clusters.
    :param representative: FASTA file of representative sequences from cluster (mmseqs)
    :param outdir:         Path to the output directory.
    :param prefix:         Prefix for output files.
    :param tmp:            Path to a temporary directory.
    :param threads:        CPU threads for mmseqs.
    :param memory:         Memory for mmseqs.
    :param clean:          True if intermediate DB files should be cleaned up.
    :param args:           Additional arguments for `mmseqs search`.


    :return: Ordered Dictionary of new defragmented clusters
    """

    from Bio import SeqIO

    clusters_path, representative_path, output_dir, tmp_dir = clusters, representative, outdir, tmp
    prefix = f"{prefix}." if prefix != None else ""

    # Check output directory
    check_output_dir(output_dir)
    check_output_dir(tmp_dir)

    args = args if args != None else ""

    # fix memory formatting (ex. '6 GB' -> '6G')
    memory = memory.replace(" ", "").replace("B", "")

    # -------------------------------------------------------------------------
    # Align representative sequences against each other

    logging.info("Aligning representative sequences.")

    seq_db = os.path.join(output_dir, f"{prefix}seqDB")
    aln_db = os.path.join(output_dir, f"{prefix}alnDB")
    tsv_path = os.path.join(output_dir, f"{prefix}align.tsv")
    run_cmd(f"mmseqs createdb {representative_path} {seq_db}")
    run_cmd(f"mmseqs search {seq_db} {seq_db} {aln_db} {tmp_dir} --threads {threads} --split-memory-limit {memory} --search-type 3 {args}")
    columns="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits"
    run_cmd(f"mmseqs convertalis {seq_db} {seq_db} {aln_db} {tsv_path} --search-type 3 --format-output {columns}")

    # Sort align rep stats for reproducibility
    with open(tsv_path, 'r') as infile:
        lines = sorted([l.strip() for l in infile.readlines()])
    with open(tsv_path, 'w') as outfile:
        header = "\t".join(columns.split(","))
        outfile.write(header + "\n")
        outfile.write("\n".join(lines) + "\n")


    # -------------------------------------------------------------------------
    logging.info(f"Reading clusters: {clusters_path}")
    loci = OrderedDict()
    clusters = OrderedDict()
    with open(clusters_path) as infile:
        lines = infile.readlines()
        for line in tqdm(lines):
            cluster, locus = [l.strip() for l in line.split("\t")]
            loci[locus] = cluster
            if cluster not in clusters:
                clusters[cluster] = []
            if locus not in clusters[cluster]:
                clusters[cluster].append(locus)

    # -------------------------------------------------------------------------
    logging.info(f"Reading representative: {representative_path}")
    representative = OrderedDict()
    for record in SeqIO.parse(representative_path, "fasta"):
        representative[record.id] = record.seq

    # -------------------------------------------------------------------------
    # Similarity Graph

    # Create a graph of cluster relationships, based on the alignment of their
    # representative sequences. The edges between clusters will represent the
    # pairwise alignment score (bits).

    # This is a modification of ppanggolin's refine_clustering function:
    # https://github.com/labgem/PPanGGOLiN/blob/2.2.0/ppanggolin/cluster/cluster.py#L317
    # The major difference is that the original function imposes a constraint
    # that the fragmented cluster cannot contain more loci than the new parent. 
    # This function does not use the number of loci, which allows a cluster
    # with many small fragments to be reassigned to a longer parent that might
    # be represented by only one intact sequence. This function also prefers
    # a OrderedDict over a formal graph object, so we can avoid the networkx
    # dependency.

    logging.info(f"Creating similarity graph from alignment: {tsv_path}")

    graph = OrderedDict()

    with open(tsv_path) as infile:
        header = infile.readline().split()
        cluster_i, locus_i, qlen_i, tlen_i, bits_i = [header.index(c) for c in ["query", "target", "qlen", "tlen", "bits"]]
        lines = infile.readlines()
        for line in tqdm(lines):
            row =  [r.strip() for r in line.replace('"', '').split()]
            if row == []: continue
            query, target, qlen, tlen, bits = [row[i] for i in [cluster_i, locus_i, qlen_i, tlen_i, bits_i]]
            if query != target:
                if query not in graph:
                    graph[query] = OrderedDict()
                if target not in graph:
                    graph[target] = OrderedDict()
                graph[query][target] = {"score": float(bits), "length": int(qlen)}
                graph[target][query] = {"score": float(bits), "length": int(tlen)}

    # -------------------------------------------------------------------------
    logging.info(f"Identifying fragmented loci.")

    # Reassign fragmented loci to their new 'parent' cluster which must:
    # 1. Be longer (length)
    # 2. Have a higher score.

    defrag_clusters = OrderedDict({c:{
        "loci": clusters[c],
        "sequence": representative[c],
        "fragments": []}
        for c in clusters
    })

    reassigned = {}

    nodes = list(graph.keys())

    for node in tqdm(nodes):
        # Get the current parent node
        candidate_node = None
        candidate_score = 0
        # Iterate through the candidate targets (clusters it could be aligned against)
        logging.debug(f"query: {node}")
        for target in graph[node]:
            ndata = graph[node][target]
            tdata = graph[target][node]
            nlen, tlen, tscore = ndata["length"], tdata["length"], tdata["score"]
            # Compare lengths and scores
            logging.debug(f"\ttarget: {target}, tlen: {tlen}, nlen: {nlen}, tscore: {tscore}, cscore: {candidate_score}")
            if tlen > nlen and candidate_score < tscore:
                candidate_node = target
                candidate_score = tscore

        # Check candidate
        if candidate_node is not None:
            # candidate node might have also been a fragment that got reassigned
            while candidate_node not in defrag_clusters and candidate_node in reassigned:
                new_candidate_node = reassigned[candidate_node]
                logging.debug(f"Following fragments: {node}-->{candidate_node}-->{new_candidate_node}")
                candidate_node = new_candidate_node
            for locus in clusters[node]:
                defrag_clusters[candidate_node]["loci"].append(locus)
                defrag_clusters[candidate_node]["fragments"].append(locus)
            del defrag_clusters[node]
            reassigned[node] = candidate_node

    # Sort order for reproducibility
    defrag_cluster_order = sorted(list(defrag_clusters.keys()))

    defrag_clusters_path = os.path.join(output_dir, f"{prefix}clusters.tsv")
    defrag_rep_path  = os.path.join(output_dir, f"{prefix}representative.fasta")

    with open(defrag_clusters_path, 'w') as clust_file:
        with open(defrag_rep_path, 'w') as seq_file:
            for cluster in defrag_cluster_order:
                info = defrag_clusters[cluster]
                # Write sequence
                seq_file.write(f">{cluster}\n{info['sequence']}" + "\n")
                # Write cluster loci, we sort for reproducibility
                for locus in sorted(info["loci"]):
                    fragment = "F" if locus in info["fragments"] else ""
                    line = f"{cluster}\t{locus}\t{fragment}"
                    clust_file.write(line + '\n')

    # Cleanup
    if clean == True:
        for file_name in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file_name)
            db_prefix = os.path.join(output_dir, prefix)
            if (
                file_path.startswith(f"{db_prefix}seqDB") or 
                file_path.startswith(f"{db_prefix}alnDB")
            ):
                logging.info(f"Cleaning up file: {file_path}")
                os.remove(file_path)

    logging.info(f"IMPORTANT!\n{PPANGGOLIN_NOTICE}")

    return (defrag_clusters_path, defrag_rep_path)


def summarize(
        clusters:  str,
        regions: str,
        outdir:    str = ".",
        product_clean:  dict = {
            " ": "_",
            "putative" : "",
            "hypothetical": "",
            "(": "",
            ")": "",
            ",": "",
        },
        max_product_len: int = 50,
        min_samples:     int = 1,
        threshold: float = 0.5,
        prefix:    str = None,
        args:      str = None,
    ):
    """
    Summarize clusters according to their annotations.

    Takes as input the clusters TSV from either the cluster or defrag subcommand,
    and the TSV table of annotations from the collect subcommand.
    Outputs a TSV table of clusters and their summarized annotations.

    >>> summarize(clusters="clusters.tsv", sequences="sequences.tsv", prefix="summarize")

    :param clusters:  TSV file of clusters from cluster or defrag (mmseqs).
    :param sequences: TSV file of sequences from collect.
    :param outdir:    Output directory.
    :param product:   Remove these words from the product description when it is the identifier.
    :param max_product_len: Truncate product to this length if it's being used as an identifier.
    :param prefix:    Prefix for output files.
    :param args:      Str of additional arguments for mmseqs search.  

    :return: Ordered Dictionary of summarized clusters
    """

    import networkx
    from networkx.exception import NetworkXNoCycle

    sequences_path, clusters_path, output_dir = regions, clusters, outdir
    prefix = f"{prefix}." if prefix != None else ""
    check_output_dir(output_dir)

    args = args if args != None else ""

    # Type conversions as fallback
    threshold = float(threshold)
    max_product_len, min_samples = int(max_product_len), int(min_samples)

    # -------------------------------------------------------------------------
    # Read Sequence records

    logging.info(f"Reading sequence regions: {sequences_path}")
    sequences = OrderedDict()
    all_samples = []
    with open(sequences_path) as infile:
        header = [line.strip() for line in infile.readline().split("\t")]
        lines = infile.readlines()
        for line in tqdm(lines):
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            sample, sequence_id = data["sample"], data["sequence_id"]
            if sample not in all_samples:
                all_samples.append(sample)
            if sequence_id not in sequences:
                sequences[sequence_id] = data
            # At this point, we should not have duplicate sequence IDs
            # That would have been handled in the extract/collect command
            else:
                msg = f"Duplicate sequence ID found: {sequence_id}"
                logging.error(msg)
                raise Exception(msg)

    # -------------------------------------------------------------------------
    # Read Clusters

    logging.info(f"Reading clusters: {clusters_path}")
    seen = set()
    clusters = OrderedDict()
    representative_to_cluster = {}
    i = 0

    with open(clusters_path) as infile:
        lines = infile.readlines()
        for line in tqdm(lines):

            # Extract the 3 columns: cluster representative, sequence Id, and (optional) fragment
            row = [l.strip() for l in line.split("\t")]
            representative, sequence_id = row[0], row[1]
            fragment = True if len(row) > 2 and row[2] == "F" else False
            if representative not in seen:
               i += 1
            cluster = f"Cluster_{i}"
            seen.add(representative)
            representative_to_cluster[representative] = cluster

            if sequence_id not in sequences:
                msg = f"Sequence is present in clusters but not in regions: {sequence_id}"
                logging.error(msg)
                raise Exception(msg)
            sequences[sequence_id]["cluster"] = cluster

            # Start the cluster data
            if cluster not in clusters:
                clusters[cluster] = OrderedDict({
                    "cluster": cluster,
                    "representative": representative, 
                    "sequences": OrderedDict()
                })
            if sequence_id not in clusters[cluster]["sequences"]:
                clusters[cluster]["sequences"][sequence_id] = sequences[sequence_id]
                clusters[cluster]["sequences"][sequence_id]["fragment"] = fragment
            else:
                msg = f"Duplicate cluster sequence ID found: {sequence_id}"
                logging.error(msg)
                raise Exception(msg)
            
    logging.info(f"Found {len(clusters)} clusters.")

    # -------------------------------------------------------------------------
    # Summarize
    # -------------------------------------------------------------------------    

    logging.info(f"Summarizing clusters.")
    summarized = OrderedDict()
    gene_counts = {}
    product_counts = {}

    for cluster,cluster_data in tqdm(clusters.items()):

        # Identify number of samples
        cluster_samples = set()
        for s in cluster_data["sequences"].values():
            cluster_samples.add(s["sample"])
        num_samples = len(list(cluster_samples))
        # If user requested a minimum number of samples
        if num_samples < min_samples:
            continue        
        samples_non_fragmented = set()
        sequences_per_sample = {}
        summarized[cluster] = OrderedDict({
            "cluster": cluster,
            "cluster_id": cluster,
            "synteny": "",
            "synteny_pos": "",
            "num_samples": num_samples,
            "num_samples_non_fragmented": 0,
            "num_sequences": 0,
            "mean_sequences_per_sample": 0,
            "representative": cluster_data["representative"],
            "upstream": "",
            "downstream": "",
            "upstream_alt": "",
            "downstream_alt": ""
            }
        )
        cluster_sequences = cluster_data["sequences"]
        features = Counter()
        genes = Counter()
        products = Counter()
        names = Counter()
        dbxrefs = Counter()
        strands = Counter()
        contexts = Counter()
        contexts_uniq = Counter()

        # -------------------------------------------------------------------------
        # Collect sequence attributes

        for sequence_id,seq_data in cluster_sequences.items():

            summarized[cluster]["num_sequences"] += 1

            sample = seq_data["sample"]
            if sample not in sequences_per_sample:
                sequences_per_sample[sample] = 0
            sequences_per_sample[sample] += 1

            # Ignore annotations from fragments starting here
            if seq_data["fragment"] == True: continue

            samples_non_fragmented.add(sample)
            features[seq_data["feature"]] += 1
            strands[seq_data["strand"]] += 1            

            # collect the annotation attributes
            attributes = {a.split("=")[0]:a.split("=")[1] for a in seq_data["attributes"].split(";") if "=" in a}
            if "gene" in attributes:
                genes[attributes["gene"]] += 1
            if "product" in attributes:
                products[attributes["product"]] += 1
            if "Name" in attributes:
                names[attributes["Name"]] += 1
            if "Dbxref" in attributes:
                dbxrefs[attributes["Dbxref"]] += 1

            # Get the upstream/downstream locus IDs
            upstream, downstream = seq_data["upstream"], seq_data["downstream"]

            # Simply upstream/downstream loci if they represent the start/end of a contig
            upstream   = "TERMINAL" if upstream.endswith("_TERMINAL") and "__" not in upstream else upstream
            downstream = "TERMINAL" if downstream.endswith("_TERMINAL") and "__" not in downstream else downstream

            # Get the upstream/downstream cluster IDs
            # There is a possibility that the upstream/downstream sequences
            # didn't actually get classified into a cluster
            if upstream in sequences and "cluster" in sequences[upstream]:
                upstream = sequences[upstream]["cluster"]
            if downstream in sequences and "cluster" in sequences[downstream]:
                downstream = sequences[downstream]["cluster"]

            # In what cases would it be itself?
            if upstream != cluster and downstream != cluster:
                context = [upstream, downstream]
                contexts["__".join(context)] += 1
                contexts_uniq["__".join(sorted(context))] += 1


        num_samples_non_fragmented = len(samples_non_fragmented)
        summarized[cluster]["num_samples_non_fragmented"] = num_samples_non_fragmented

        mean_sequences_per_sample = sum(sequences_per_sample.values()) / len(sequences_per_sample)
        summarized[cluster]["mean_sequences_per_sample"] = round(mean_sequences_per_sample, 1)

        # -------------------------------------------------------------------------
        # Summarize upstream/downstream

        # Part 1. Are there neighboring loci (regardless of directionality)?
        neighbors = None
        if len(contexts_uniq) > 1:
            most_common, count = contexts_uniq.most_common(1)[0]
            # Are these loci observed frequently enough?
            prop = count / num_samples_non_fragmented
            if prop >= threshold:
                # Handle ties, by simply picking the first option alphabetically
                candidates = sorted([c for c,v in contexts_uniq.items() if v == count])
                neighbors = candidates[0]
                if len(candidates) > 1:
                    logging.debug(f"{cluster} tie broken alphabetically: neighbors={neighbors}, {contexts_uniq}")

        # Part 2. Filter the neighbors to our top matches, now allowing either direction
        # We will summarize in the next step
        if neighbors != None:
            c1, c2 = most_common.split("__")
            forward, reverse = f"{c1}__{c2}", f"{c2}__{c1}"
            contexts = Counter({k:v for k,v in contexts.items() if k == forward or k == reverse})
            
         
        # -------------------------------------------------------------------------
        # Summarize sequence attributes

        for key,values in zip(["feature", "strand", "gene", "product", "name", "dbxref", "contexts"], [features, strands, genes, products, names, dbxrefs, contexts]):
            value = ""
            value_alt = ""

            if len(values) > 0:
                # Check if the most common value passes the threshold
                most_common_count = values.most_common(1)[0][1]
                most_common_prop = most_common_count / num_samples_non_fragmented
                # We don't do a second threshold filter for contexts
                if key == "contexts" or most_common_prop >= threshold :
                    # Handle ties, by simply picking the first option alphabetically
                    candidates = sorted([c for c,v in values.items() if v == most_common_count])
                    value = candidates[0]
                    if len(candidates) > 1:
                        logging.debug(f"{cluster} tie broken alphabetically: {key}={value}, {values}")
                    value_alt = ";".join([v for v in values.keys() if v != value])
                else:
                    value_alt = ";".join([v for v in values.keys()])

            if key == "contexts":
                if value != "":
                    upstream, downstream = value.split("__")
                    summarized[cluster]["upstream"] = upstream
                    summarized[cluster]["downstream"] = downstream
                continue

            summarized[cluster][key] = value
            summarized[cluster][f"{key}_alt"] = value_alt

            # gene/product identifiers need to be checked for case!
            
            if key == "gene" and value != "":
                value_lower = value.lower()
                if value_lower not in gene_counts:
                    gene_counts[value_lower] = {"count": 0, "current": 1}
                gene_counts[value_lower]["count"] += 1

            if key == "product":
                summarized[cluster]["product_clean"] = ""
                if value != "":
                    # Clean up product as a potential identifier
                    # These values are definitely not allowed
                    clean_value = value.replace(" ", "_").replace("/", "_").replace(",", "_")
                    for k,v in product_clean.items():
                        clean_value = clean_value.replace(k, v)
                    # Restrict the length
                    if max_product_len != None:
                        clean_value = clean_value[:max_product_len]
                    while "__" in clean_value:
                        clean_value = clean_value.replace("__", "_")
                    while clean_value.startswith("_") or clean_value.endswith("_"):
                        clean_value = clean_value.lstrip("_").rstrip("_")
                    if clean_value != "":
                        clean_value_lower = clean_value.lower()
                        if clean_value_lower not in product_counts:
                            product_counts[clean_value_lower] = {"count": 0, "current": 1}
                        product_counts[clean_value_lower]["count"] += 1
                        summarized[cluster]["product_clean"] = clean_value

            # Handle if feature was 'unannotated'
            # Can happen in tie breaking situations
            if (
                (key == "gene" or key == "product") and 
                value != "" and 
                summarized[cluster]["feature"] == "unannotated"
                ):
                candidates = [f for f in summarized[cluster]["feature_alt"].split(";") if f != "unannotated"]
                feature = candidates[0]
                summarized[cluster]["feature"] = feature
                features_alt = ["unannotated"] + [c for c in candidates if c != feature]
                summarized[cluster]["feature_alt"] = ";".join(features_alt)
                logging.debug(f"Updating {cluster} from an unannotated feature to {feature} based on {key}={value}.")

        summarized[cluster]["sequences"] = cluster_sequences

    # -------------------------------------------------------------------------
    # Identifiers: Part 1
    # -------------------------------------------------------------------------    

    # Give identifiers to clusters based on their gene or product.
    # We do this now, so that the synteny graph has some nice 
    # helpful names for clusters. We will give names to the 
    # unnanotated clusters based on their upstream/downstream
    # loci after the synteny graph is made.

    logging.info(f"Assigning identifiers to annotated clusters.")

    identifiers = OrderedDict()

    # Pass 1. Identifiers based on gene/product
    for cluster,cluster_data in tqdm(summarized.items()):
        identifier = None       
        # Option 1. Try to use gene as cluster identifier
        if cluster_data["gene"] != "":
            gene = cluster_data["gene"]
            # we use the lowercase to figure out the duplicate
            # number, but use the original name in the table
            gene_lower = gene.lower()
            if gene_counts[gene_lower]["count"] == 1:
                identifier = gene
            else:
                dup_i = gene_counts[gene_lower]["current"]
                gene_counts[gene_lower]["current"] += 1
                new_gene = f"{gene}.{dup_i}"
                logging.debug(f"Duplicate gene identifer found for {cluster}: {new_gene}")
                identifier = new_gene

        # Option 2. Try to use product as cluster identifier
        #           Useful for non-gene annotations (ex. tRNA)
        elif cluster_data["product_clean"] != "":
            product = cluster_data["product_clean"]
            product_lower = product.lower()
            if product_counts[product_lower]["count"] == 1:
                identifier = product
            else:
                dup_i = product_counts[product_lower]["current"]
                product_counts[product_lower]["current"] += 1
                new_product = f"{product}.{dup_i}"
                logging.debug(f"Duplicate product identifer found for {cluster}: {new_product}")
                identifier = new_product

        if identifier != None:
            identifiers[cluster] = identifier
            summarized[cluster]["cluster"] = identifier

    # Pass 2: Update the upstream/downstream identifiers
    for cluster,cluster_data in tqdm(summarized.items()):
        upstream, downstream = cluster_data["upstream"], cluster_data["downstream"]
        if upstream in identifiers:
            summarized[cluster]["upstream"] = identifiers[upstream]
        if downstream in identifiers:
            summarized[cluster]["downstream"] = identifiers[downstream]

    # Update the cluster keys
    summarized = OrderedDict({v["cluster"]:v for k,v in summarized.items()})

    # -------------------------------------------------------------------------
    # Synteny
    # -------------------------------------------------------------------------

    logging.info(f"Computing initial synteny graph.")

    # ------------------------------------------------------------------------- 
    # Create initial graph 

    # Create simple graph based on what we know about the
    # upstream/downstream loci so far. This will be our "full"
    # graph which retains all the cycles and multifurcations
    synteny_full = networkx.DiGraph()
    seen_nodes = set()

    for cluster, cluster_data in summarized.items():
        if cluster not in seen_nodes:
            synteny_full.add_node(cluster)
        upstream, downstream = cluster_data["upstream"], cluster_data["downstream"]
        if upstream != "" and upstream != "TERMINAL":
            synteny_full.add_edge(upstream, cluster)
        if downstream != "" and downstream != "TERMINAL":
            synteny_full.add_edge(cluster, downstream)

    # ------------------------------------------------------------------------- 
    # Filter nodes

    # Remove clusters from the graph. This can happen if the 
    # user requested min_samples>0, so we need to remove the 
    # low prevalence clusters. It can also happen if a 
    # cluster had an upstream/downstream loci that didn't
    # make it to the final clustering list
    logging.info(f"Filtering graph for missing clusters.")
    for node in list(networkx.dfs_tree(synteny_full)):
        if node in summarized: continue
        in_nodes = [e[0] for e in synteny_full.in_edges(node)]
        out_nodes = [e[1] for e in synteny_full.out_edges(node)]
        # Remove this node from the graph
        logging.debug(f"Removing {node} from the graph.")
        synteny_full.remove_node(node)
        # Connect new edges between in -> out
        for n1 in [n for n in in_nodes if n not in out_nodes]:
            for n2 in [n for n in out_nodes if n not in in_nodes]:
                # Not sure if we need this check?
                if not synteny_full.has_edge(n1, n2):
                    synteny_full.add_edge(n1, n2)

    # ------------------------------------------------------------------------- 
    # Break up multifurcations

    synteny_linear = copy.deepcopy(synteny_full)

    logging.info(f"Breaking up multifurcations.")
    for node in tqdm(list(networkx.dfs_tree(synteny_linear))):
        in_nodes = [e[0] for e in synteny_linear.in_edges(node)]
        out_nodes = [e[1] for e in synteny_linear.out_edges(node)]
        
        if len(in_nodes) > 1:
            for n in in_nodes:
                logging.debug(f"Removing multifurcation in_edge: {n} -> {node}")
                synteny_linear.remove_edge(n, node)
        if len(out_nodes) > 1:
            for n in out_nodes:
                logging.debug(f"Removing multifurcation out_edge: {node} -> {n}")
                synteny_linear.remove_edge(node, n)

    # ------------------------------------------------------------------------- 
    # Isolate and linear cycles

    # I'm not sure in what cases this is still needed
    # after switching the graph to directed
    logging.info(f"Isolating cycles and linearizing.")
    cycles = True
    while cycles == True:
        try:
            cycle_raw = [x for xs in networkx.find_cycle(synteny_linear) for x in xs]
            seen = set()
            cycle = []
            for n in cycle_raw:
                if n not in seen:
                    cycle.append(n)
                    seen.add(n)
            logging.debug(f"Cycle found: {cycle}")
            cycle_set = set(cycle)
            # Separate this cycle from the rest of the graph
            for n1 in synteny_linear.nodes():
                if n1 in cycle: continue
                neighbors = set(synteny_linear.neighbors(n1))
                for n2 in cycle_set.intersection(neighbors):
                    logging.debug(f"Isolating cycle by removing edge: {n1} <-> {n2}")
                    synteny_linear.remove_edge(n1, n2)

            # Break the final cycle between the 'first' and 'last' nodes
            first, last = cycle[0], cycle[-1]            
            logging.debug(f"Breaking cycle by removing edge: {last} -> {first}")
            synteny_linear.remove_edge(last, first)

        except NetworkXNoCycle:
            cycles = False

    # ------------------------------------------------------------------------- 
    # Identify synteny blocks

    # We need to use the function connected_components, but that only 
    # works on undirected graphs
    synteny_linear_undirected = synteny_linear.to_undirected()

    # Get synteny blocks, sorted largest to smallest
    logging.info(f"Identifying synteny blocks.")
    synteny_blocks = sorted([
        synteny_linear_undirected.subgraph(c).copy()
        for c in networkx.connected_components(synteny_linear_undirected)
    ], key=len, reverse=True)

    # Now we have enough information to finalize the order
    summarized_order = OrderedDict()

    for i_b, block in enumerate(tqdm(synteny_blocks)):
        i_b += 1
        clusters = list(block)

        if len(clusters) > 1:
            terminals = [n for n in block.nodes if len(list(block.neighbors(n))) == 1]
            # Legacy error from troubleshooting, unclear if still relevant                
            if len(terminals) != 2:
                # Check if it's a cycle we somehow didn't handle
                try:
                    cycle = networkx.find_cycle(block)
                    cycle = True
                except NetworkXNoCycle:
                    cycle = False
                msg = f"Synteny block {i_b} has an unhandled error, cycle={cycle}, terminals={len(terminals)}: {terminals}"
                logging.error(msg)
                print(networkx.write_network_text(block))
                raise Exception(msg)

            # Figure out which terminal is the 5' end
            first_upstream = summarized[terminals[0]]["upstream"]
            last_upstream = summarized[terminals[-1]]["upstream"]

            if first_upstream == "TERMINAL":
                first, last = terminals[0], terminals[1]
            elif last_upstream == "TERMINAL":
                first, last = terminals[1], terminals[0]
            else:
                # If it's fully ambiguous, we'll sort the terminals for reproducibility
                terminals = sorted(terminals)
                first, last = terminals[0], terminals[1]

            # Manually walk through the graph, this is not ideal
            # but networkx's bfs/dfs was giving me odd values in testing
            clusters, neighbors = [first], []
            curr_node, next_node = first, None
            while next_node != last:
                neighbors = [n for n in block.neighbors(curr_node) if n not in clusters]
                # Legacy error from troubleshooting, unclear if still relevant
                if len(neighbors) != 1:
                    msg = f"Synteny error, unhandled multifurcation in {curr_node}: {neighbors}"
                    logging.error(msg)
                    raise Exception(msg)
                next_node = neighbors[0]
                clusters.append(next_node)
                curr_node = next_node

        for i_c, cluster in enumerate(clusters):
            summarized_order[cluster] = summarized[cluster]
            summarized_order[cluster]["synteny"] = str(i_b)
            summarized_order[cluster]["synteny_pos"] = str(i_c + 1)

            upstream = summarized_order[cluster]["upstream"]
            downstream = summarized_order[cluster]["downstream"]

            # Use the synteny block to finalize upstream/downstream
            upstream = summarized_order[cluster]["upstream"]
            downstream = summarized_order[cluster]["downstream"]

            # We'll save a copy of how it was before as the alt
            upstream_orig, downstream_orig = copy.deepcopy(upstream), copy.deepcopy(downstream)

            # 'TERMINAL' will now refer to the ends of the synteny block
            if i_c > 0:
                upstream = clusters[i_c - 1]
            else:
                upstream = "TERMINAL"

            if i_c < (len(clusters) - 1):
                downstream = clusters[i_c + 1]
            else:
                downstream = "TERMINAL"

            if upstream != upstream_orig:
                summarized_order[cluster]["upstream_alt"] = upstream_orig
            if downstream != downstream_orig:
                summarized_order[cluster]["downstream_alt"] = downstream_orig              

            summarized_order[cluster]["upstream"] = upstream
            summarized_order[cluster]["downstream"] = downstream

    summarized = summarized_order

    # -------------------------------------------------------------------------
    # Create directed

    logging.info(f"Converting synteny blocks to directed graph.")

    # Now we need to go back to a directed graph
    synteny_linear_directed = networkx.DiGraph()
    synteny_seen = set()

    for cluster,cluster_data in summarized.items():
        upstream, downstream = cluster_data["upstream"], cluster_data["downstream"]
        if cluster not in synteny_seen:
            synteny_seen.add(cluster)
            synteny_linear_directed.add_node(cluster)

        if upstream != "TERMINAL":
            if upstream not in synteny_seen:
                synteny_seen.add(upstream)
                synteny_linear_directed.add_edge(upstream, cluster)
                                    
        if downstream != "TERMINAL":
            if downstream not in synteny_seen:
                synteny_seen.add(downstream)
                synteny_linear_directed.add_edge(cluster, downstream)

    # -------------------------------------------------------------------------
    # Identifiers: Part 2
    # -------------------------------------------------------------------------

    # Give unannotated clusters identifiers based on their upstream/downstream
    # loci in the synteny graph.

    # Note: the synteny reconstruction ensures that we're not going to have
    #       any duplicate cluster IDs for the unannotated regions
    #       because any loops/multifurcations have been removed.

    logging.info(f"Assigning identifiers to unannotated clusters.")

    # -------------------------------------------------------------------------
    # Pass #1: Give identifiers to unannotated clusters based on upstream/downstream

    for cluster,cluster_data in tqdm(summarized.items()):
        # Skip this cluster if it already has a new identifier
        # Example, annotated clusters based on gene/product
        if cluster in identifiers and identifiers[cluster] != cluster:
            continue
        # Skip this cluster if it has gene/product info, just a safety fallback
        if cluster_data["gene"] != "" or cluster_data["product"] != "":
            continue
        upstream, downstream = cluster_data["upstream"], cluster_data["downstream"]
        # Option 1. No known neighbors
        if upstream == "TERMINAL" and downstream == "TERMINAL":
            identifier = cluster
        # Option 2. Upstream/downstream already has the notation
        elif "__" in upstream or "__" in downstream:
            identifier = cluster
        # Option 3. At least one side is known
        else:
            identifier = f"{upstream}__{downstream}"

        identifiers[cluster] = identifier
        summarized[cluster]["cluster"] = identifier

    # -------------------------------------------------------------------------
    # Pass #2: Finalize upstream/downstream

    # Now that we know the identifiers, we need to update the 
    # following fields: cluster, upstream, and downstream

    logging.info(f"Finalizing upstream/downstream identifiers.")

    for cluster, cluster_data in tqdm(summarized.items()):
        upstream, downstream = cluster_data["upstream"], cluster_data["downstream"]
        if upstream in identifiers:
            summarized[cluster]["upstream"] = upstream
        if downstream in identifiers:
            summarized[cluster]["downstream"] = downstream

    # Update the keys in the graph
    summarized = OrderedDict({v["cluster"]:v for k,v in summarized.items()})

    # -------------------------------------------------------------------------
    # Update synteny graph with new identifiers
    # -------------------------------------------------------------------------

    logging.info(f"Updating cluster identifiers in the synteny graphs.")

    # We will both update the original "full" graph, as well as our new
    # "linear" graph

    # -------------------------------------------------------------------------
    # Full Graph

    networkx.relabel_nodes(synteny_full, mapping=identifiers)
    edges = list(synteny_full.out_edges())
    for c1, c2 in edges:
        synteny_full.remove_edge(c1, c2)
        c1 = identifiers[c1] if c1 in identifiers else c1
        c2 = identifiers[c2] if c2 in identifiers else c2
        synteny_full.add_edge(c1, c2)

    synteny_full_path = os.path.join(output_dir, f"{prefix}synteny.full.graphml")
    logging.info(f"Writing full synteny GraphML: {synteny_full_path}")
    networkx.write_graphml(synteny_full, synteny_full_path)

    gfa_path = os.path.join(output_dir, f"{prefix}synteny.full.gfa")
    logging.info(f"Writing full synteny GFA: {gfa_path}")

    with open(gfa_path, 'w') as outfile:
        outfile.write("H\tVN:Z:1.0\n")
        for cluster in summarized:
            outfile.write(f"S\t{cluster}\t*\tLN:i:1\n")
        for c1, c2 in synteny_full.out_edges():
            c1_strand, c2_strand = summarized[c1]["strand"], summarized[c2]["strand"]
            outfile.write(f"L\t{c1}\t{c1_strand}\t{c2}\t{c2_strand}\t0M\n")

    # -------------------------------------------------------------------------
    # Linearized Graph

    networkx.relabel_nodes(synteny_linear_directed, mapping=identifiers)
    edges = list(synteny_linear_directed.out_edges())

    for c1, c2 in edges:
        synteny_linear_directed.remove_edge(c1, c2)
        c1 = identifiers[c1] if c1 in identifiers else c1
        c2 = identifiers[c2] if c2 in identifiers else c2
        synteny_linear_directed.add_edge(c1, c2)

    synteny_linear_path = os.path.join(output_dir, f"{prefix}synteny.linear.graphml")
    logging.info(f"Writing linear synteny GraphML: {synteny_linear_path}")
    networkx.write_graphml(synteny_linear_directed, synteny_linear_path)

    gfa_path = os.path.join(output_dir, f"{prefix}synteny.linear.gfa")
    logging.info(f"Writing linear synteny GFA: {gfa_path}")

    with open(gfa_path, 'w') as outfile:
        outfile.write("H\tVN:Z:1.0\n")
        for cluster in summarized:
            outfile.write(f"S\t{cluster}\t*\tLN:i:1\n")
        for c1, c2 in synteny_linear_directed.out_edges():
            c1_strand, c2_strand = summarized[c1]["strand"], summarized[c2]["strand"]
            outfile.write(f"L\t{c1}\t{c1_strand}\t{c2}\t{c2_strand}\t0M\n")

    # -------------------------------------------------------------------------
    # Write Output tsv

    tsv_path = os.path.join(output_dir, f"{prefix}clusters.tsv")
    logging.info(f"Writing summarized clusters tsv: {tsv_path}")

    with open(tsv_path, 'w') as outfile:
        header = None
        for i, cluster in enumerate(tqdm(summarized)):
            cluster_data = summarized[cluster]
            if not header:
                header = [k for k in cluster_data if k != "sequences" and k != "product_clean"] + all_samples
                outfile.write("\t".join(header) + "\n")
            row = [str(v) for k,v in cluster_data.items() if k != "sequences" and k != "product_clean"]
            # Add info about which sequences map to each sample in the cluster
            sample_to_seq_id = OrderedDict({sample:[] for sample in all_samples})
            for seq_id,seq_data in cluster_data["sequences"].items():
                sample = seq_data["sample"]
                sample_to_seq_id[sample].append(seq_id)
            for sample in all_samples:
                row += [",".join(sample_to_seq_id[sample])]

            outfile.write("\t".join(row) + "\n")

    # -------------------------------------------------------------------------
    # Write table (for phandango)

    phandango_path = os.path.join(output_dir, f"{prefix}phandango.csv")
    logging.info(f"Writing table for phandango: {phandango_path}")

    logging.info(f"Sorting synteny blocks according to number of samples: {phandango_path}")
    syntenies = OrderedDict()
    for cluster,c_data in summarized.items():
        synteny = c_data["synteny"]
        if synteny not in syntenies:
            syntenies[synteny] = {"max_samples": 0, "clusters": []}
        num_samples = c_data["num_samples"]
        syntenies[synteny]["max_samples"] = max(num_samples, syntenies[synteny]["max_samples"])
        syntenies[synteny]["clusters"].append(cluster)
    
    syntenies = OrderedDict(
        sorted(
            syntenies.items(), 
            key=lambda item: item[1]["max_samples"], reverse=True
        )
    )

    # This it the roary gene_presence_absence.csv format
    with open(phandango_path, 'w') as outfile:
        header = [
            "Gene","Non-unique Gene name","Annotation","No. isolates",
            "No. sequences","Avg sequences per isolate","Genome Fragment",
            "Order within Fragment","Accessory Fragment",
            "Accessory Order with Fragment","QC"
        ] + all_samples
        outfile.write(",".join(header) + "\n")
        for synteny,s_data in tqdm(syntenies.items()):
            for cluster in s_data["clusters"]:
                c_data = summarized[cluster]
                data = OrderedDict({k:"" for k in header})
                # This is deliberately reversed
                data["Gene"] = c_data["cluster"]
                data["Non-unique Gene name"] = c_data["cluster_id"]
                # We're going to use the cleaned product, because
                # we need at least some bad characters to removed 
                # (like commas)
                data["Annotation"] = c_data["product_clean"]
                data["No. isolates"] = c_data["num_samples"]
                data["No. sequences"] = cluster_data["num_sequences"]
                data["Avg sequences per isolate"] = c_data["mean_sequences_per_sample"]
                data["Genome Fragment"] = c_data["synteny"]
                data["Order within Fragment"] = c_data["synteny_pos"]
                # If a sample has a sequence then it will be given a "1"
                # If not, is recorded as the empty string ""
                # This is based on the file minimizing back from phandango:
                #   https://github.com/jameshadfield/phandango/blob/master/scripts/minimiseROARY.py            
                for s_data in c_data["sequences"].values():
                    sample = s_data["sample"]
                    data[sample] = "1"
                line = ",".join([str(v) for v in data.values()])
                outfile.write(line + "\n")

    return tsv_path


# Parallel
def run_mafft(kwargs: dict):
    """A wrapper function to pass multi-threaded pool args to mafft."""
    cmd, output, quiet, display_cmd = [kwargs[k] for k in ["cmd", "output", "quiet", "display_cmd"]]
    run_cmd(cmd=cmd, output=output, quiet=quiet, display_cmd=display_cmd)

def align(
        clusters:           str,
        regions:            str,
        outdir:             str  = ".",
        prefix:             str  = None,
        exclude_singletons: bool = False,
        threads:            int  = 1,
        args:               str  = ALIGN_ARGS,
    ):
    """
    Align clusters using mafft and create a pangenome alignment.

    Takes as input the summarized clusters from summarize and sequence regions from collect.
    Outputs multiple sequence alignments per cluster as well as a pangenome alignment of 
    concatenated clusters.

    >>> align(clusters="summarize.clusters.tsv", sequences="sequences.tsv")
    >>> align(clusters="summarize.clusters.tsv", sequences="sequences.tsv", exclude_singletons=True, args="--localpair")

    :param clusters:           TSV file of clusters from summarize.
    :param sequences:          TSV file of sequence regions from collect.
    :param outdir:             Output directory.
    :param prefix:             Prefix for output files.
    :param exclude_singletons: True is clusters found in only one sample should be excluded.
    :param threads:            Number of cpu threads to parallelize mafft across.
    :param args:               Additional arguments for MAFFT.
    """

    from multiprocessing import get_context

    clusters_path, sequences_path, output_dir = clusters, regions, outdir

    # Check output directory
    check_output_dir(output_dir)

    args = args if args != None else ""
    prefix = f"{prefix}." if prefix != None else ""
    threads = int(threads)

    # -------------------------------------------------------------------------
    # Read Sequence Regions

    all_samples = []
    sequences = {}

    logging.info(f"Reading sequence regions: {sequences_path}")
    with open(sequences_path) as infile:
        header = [line.strip() for line in infile.readline().split("\t")]
        for line in infile:
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            sample, sequence_id = data["sample"], data["sequence_id"]
            if sample not in all_samples:
                all_samples.append(sample)
            sequences[sequence_id] = data

    # -------------------------------------------------------------------------
    # Read Summarized Clusters

    logging.info(f"Reading summarized clusters: {clusters_path}")
    clusters = OrderedDict()
    all_samples = []
    with open(clusters_path) as infile:
        header = [line.strip() for line in infile.readline().split("\t")]
        # sample columns begin after dbxref_alt
        all_samples = header[header.index("dbxref_alt")+1:]
        lines = infile.readlines()
        for line in tqdm(lines):
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            cluster = data["cluster"]
            clusters[cluster] = data
            clusters[cluster]["sequences"] = OrderedDict()
            for sample in all_samples:
                sequence_ids = data[sample].split(",") if data[sample] != "" else []
                # Add the sequence IDs, we will associate it with the actual
                # sequence from the collect TSV
                for sequence_id in sequence_ids:
                    clusters[cluster]["sequences"][sequence_id] = {
                        "sequence": sequences[sequence_id]["sequence"],
                        "sample": sample
                    }

    # -------------------------------------------------------------------------
    # Write Sequences

    representative_output = os.path.join(output_dir, prefix + "representative")
    sequences_output  = os.path.join(output_dir, prefix + "sequences")
    alignments_output = os.path.join(output_dir, prefix + "alignments")
    consensus_output  = os.path.join(output_dir, prefix + "consensus")

    check_output_dir(representative_output)
    check_output_dir(sequences_output)
    check_output_dir(alignments_output)
    check_output_dir(consensus_output)

    # -------------------------------------------------------------------------
    # Write Representative sequences for each cluster to file

    skip_align = set()
    clusters_exclude_singletons = OrderedDict()

    logging.info(f"Writing representative sequences: {representative_output}")
    for cluster,cluster_data in tqdm(clusters.items()):
        rep_seq_id = cluster_data["representative"]
        rep_seq = cluster_data["sequences"][rep_seq_id]["sequence"]
        rep_path = os.path.join(representative_output, f"{cluster}.fasta")
        samples = list(set([s["sample"] for s in cluster_data["sequences"].values()]))

        # If only found in 1 sample, and the user requested to exclude these, skip
        if len(samples) == 1 and exclude_singletons:
            logging.debug(f"Skipping singleton cluster: {cluster}")
            continue

        # If this is cluster only has one sequence, write to final output
        if len(cluster_data["sequences"]) == 1:
            skip_align.add(cluster)
            file_path = os.path.join(alignments_output, f"{cluster}.aln")
            with open(file_path, 'w') as outfile:
                outfile.write(f">{rep_seq_id}\n{rep_seq}\n")
        # Otherwise, save rep seq to separate folder
        # We might need it if user has requeseted the --add* mafft args
        else:
            with open(rep_path, 'w') as outfile:
                outfile.write(f">{rep_seq_id}\n{rep_seq}\n")

        clusters_exclude_singletons[cluster] = cluster_data
    
    clusters = clusters_exclude_singletons

    # -------------------------------------------------------------------------
    # Write DNA sequences for each cluster to file

    # A queue of commands to submit to mafft in parallel 
    mafft_queue = []

    logging.info(f"Writing cluster sequences: {sequences_output}")
    for i,cluster in enumerate(tqdm(clusters)):
        # Skip singleton clusters, we already handled them in previous block
        if cluster in skip_align: continue
        cluster_data = clusters[cluster]
        rep_seq_id = cluster_data["representative"]
        representative_path = os.path.join(representative_output,  f"{cluster}.fasta")
        sequences_path = os.path.join(sequences_output,  f"{cluster}.fasta")
        with open(sequences_path, "w") as outfile:
            for sequence_id,seq_data in cluster_data["sequences"].items():
                # Skip representative if we're using addfragments
                if "--add" in args and sequence_id == rep_seq_id:
                    continue
                sequence = seq_data["sequence"]
                line = f">{sequence_id}\n{sequence}\n"
                outfile.write(line)
        alignment_path = os.path.join(alignments_output,  f"{cluster}.aln")
        cmd = f"mafft --thread 1 {args} {sequences_path}"
        # If we're using addfragments, we're aligning against representative/reference
        if "--add" in args:
            cmd += f" {representative_path}"
        mafft_queue.append({"cmd": cmd, "output": alignment_path, "quiet": True, "display_cmd": False})

    # Display first command
    if len(mafft_queue) > 0:
        logging.info(f"Command to run in parallel: {mafft_queue[0]['cmd']}")

    # -------------------------------------------------------------------------
    # Align DNA sequences with MAFFT

    # Parallel: MAFFT is very CPU-bound, so parallel is appropriate
    logging.info(f"Aligning cluster sequences in parallel with {threads} threads: {alignments_output}")
    with get_context('fork').Pool(threads) as p:
        with tqdm(total=len(mafft_queue), unit="cluster") as bar:
            for _ in p.imap_unordered(run_mafft, mafft_queue):
               bar.update()

    # -------------------------------------------------------------------------
    # Unwrap and dedup cluster alignments

    # Duplicates occur due to multi-copy/fragments. We'll try to dedup the 
    # sequence by reconstructing consensus bases where possible.

    logging.info(f"Unwrapping and dedupping alignments: {consensus_output}")
    consensus_alignments = OrderedDict()

    for cluster,cluster_data in tqdm(clusters.items()):
        original_path = os.path.join(alignments_output, cluster + ".aln")
        tmp_path = original_path + ".tmp"
        consensus_path = os.path.join(consensus_output, cluster + ".aln")

        # Unwrap sequences and convert to uppercase
        alignment = {}
        with open(original_path, 'r') as infile:
            with open(tmp_path, 'w') as outfile:
                records = infile.read().split(">")[1:]
                for record in records:
                    record_split = record.split("\n")
                    sequence_id = record_split[0]
                    # Check for indicator that it was reverse complemented
                    if sequence_id not in cluster_data["sequences"] and sequence_id.startswith("_R_"):
                        sequence_id = sequence_id[3:]

                    sample = cluster_data["sequences"][sequence_id]["sample"]
                    sequence = "".join(record_split[1:]).replace("\n", "").upper()
                    if sample not in alignment:
                        alignment[sample] = []
                    alignment[sample].append(sequence)

                    # Write the sequence to the output file, using the original header
                    outfile.write(f">{sequence_id}\n{sequence}\n")

        # Replace the original file, now that unwrapped and uppercased
        # We can use this in a different script for structural variant detection
        shutil.move(tmp_path, original_path)

        if len(alignment) == 0:
            logging.info(f"WARNING: No sequences written for cluster: {cluster}")
            continue

        # Create consensus sequence from fragments and multi-copies
        duplicate_samples = [sample for sample,seqs in alignment.items() if len(seqs) > 1]
        # Create consensus alignment, first with the non-duplicated samples
        alignment_consensus = {sample:seqs[0] for sample,seqs in alignment.items() if len(seqs) == 1}

        for sample in duplicate_samples:
            seqs = alignment[sample]
            length = len(seqs[0])
            consensus = []

            # Reconstruct Consensus sequence
            for i in range(0, length):
                nuc_raw = set([s[i] for s in seqs])
                nuc = list(set([n for n in nuc_raw if n in NUCLEOTIDES] ))
                if len(nuc) == 1:      nuc_consensus = nuc[0]
                elif nuc_raw == set("-"): nuc_consensus = "-"
                else:                  nuc_consensus = "N"
                consensus.append(nuc_consensus)

            consensus_sequence = "".join(consensus)
            alignment_consensus[sample] = consensus_sequence

        consensus_alignments[cluster] = alignment_consensus

        # Write unwrapped, dedupped, defragged, consensus alignment
        with open(consensus_path, 'w') as outfile:
            for (sample, sequence) in alignment_consensus.items():
                outfile.write(f">{sample}\n{sequence}\n")

    # -------------------------------------------------------------------------
    # Concatenate into pangenome alignment

    logging.info(f"Creating pangenome alignment.")

    pangenome = {
        "bed" : {},
        "alignment" : {s: [] for s in all_samples},
    }

    curr_pos = 0

    for cluster,alignment in tqdm(consensus_alignments.items()):
        # identify samples missing from cluster
        observed_samples = list(alignment.keys())
        missing_samples  = [s for s in all_samples if s not in alignment]

        # write gaps for missing samples
        seq_len = len(alignment[observed_samples[0]])
        for sample in missing_samples:
            alignment[sample] = "-" * seq_len

        # concatenate cluster sequence to phylo alignment
        for sample, seq in alignment.items():
            pangenome["alignment"][sample].append(seq)

        # update bed coordinates
        prev_pos = curr_pos
        curr_pos = curr_pos + seq_len

        pangenome["bed"][prev_pos] = {
            "start" : prev_pos,
            "end" : curr_pos,
            "cluster" : cluster,
            "synteny": clusters[cluster]["synteny"]
        }

    # Final concatenation of alignment
    logging.info("Performing final concatenation.")
    for sample in tqdm(all_samples):
        pangenome["alignment"][sample] = "".join( pangenome["alignment"][sample])

    # -------------------------------------------------------------------------
    # Write pangenome Consensus Sequence

    consensus_file_path = os.path.join(output_dir, prefix + "pangenome.consensus.fasta")
    logging.info(f"Writing pangenome consensus: {consensus_file_path}")
    with open(consensus_file_path, 'w') as out_file:
        out_file.write(">consensus\n")
        # get total length from first sample in alignment
        first_sample = list(pangenome["alignment"].keys())[0]
        consensus_len = len(pangenome["alignment"][first_sample])
        consensus = []
        logging.info(f"Consensus Length: {consensus_len}")
        for cluster_start,bed_info in tqdm(pangenome["bed"].items()):
            cluster_end, cluster = bed_info["end"], bed_info["cluster"]
            for i in range(cluster_start, cluster_end):
                nuc = [pangenome["alignment"][s][i] for s in all_samples]
                nuc_non_ambig = [n for n in nuc if n in NUCLEOTIDES]
                # All samples excluded due to ambiguity
                if len(set(nuc_non_ambig)) == 0:
                    consensus_nuc = "N"
                # Invariant position, all the same
                elif len(set(nuc_non_ambig)) == 1:
                    consensus_nuc = nuc[0]
                # Variant position, get consensus of nonambiguous nucleotides
                else:
                    # choose which ever nucleotide has highest count as the consensus
                    counts = {n:nuc.count(n) for n in NUCLEOTIDES}
                    consensus_nuc = max(counts, key=counts.get)
                consensus.append(consensus_nuc)
        # write final end of line
        out_file.write("".join(consensus) + "\n")

    # -------------------------------------------------------------------------
    # Write Bed File

    bed_file_path = os.path.join(output_dir, prefix + "pangenome.bed")
    logging.info(f"Writing bed file: {bed_file_path}")
    with open(bed_file_path, "w") as bed_file:
        for start,info in tqdm(pangenome["bed"].items()):
            cluster, synteny = info["cluster"], info["synteny"]
            name = f"cluster={cluster};synteny={synteny}"
            row = ["pangenome", start, info["end"], name]
            line = "\t".join([str(val) for val in row])
            bed_file.write(line + "\n")

    # -------------------------------------------------------------------------
    # Write Alignment

    aln_file_path = os.path.join(output_dir, prefix + "pangenome.aln")
    logging.info(f"Writing alignment file: {aln_file_path}")
    with open(aln_file_path, "w") as aln_file:
        for sample in tqdm(all_samples):
            seq = pangenome["alignment"][sample].upper()
            aln_file.write(">" + sample + "\n" + seq + "\n")  

    # -------------------------------------------------------------------------
    # Write GFF

    gff_file_path = os.path.join(output_dir, prefix + "pangenome.gff3")
    logging.info(f"Writing gff file: {gff_file_path}")
    with open(gff_file_path, "w") as outfile:
        outfile.write("##gff-version 3\n")
        outfile.write(f"##sequence-region pangenome 1 {consensus_len}\n")
        row = ["pangenome", ".", "region", 1, consensus_len, ".", "+", ".", "ID=pangenome;Name=pangenome"]
        line = "\t".join([str(v) for v in row])
        outfile.write(line + "\n")

        for start in pangenome["bed"]:
            start, end, cluster, synteny = pangenome["bed"][start].values()
            cluster_data = clusters[cluster]
            representative, feature  = cluster_data["representative"], cluster_data["feature"]
            strand = sequences[representative]["strand"]
            attributes = OrderedDict({"ID": cluster, "locus_tag": cluster})
            for col,k in zip(["gene", "name", "product", "dbxref"], ["gene", "Name", "product", "Dbxref"]):
                if cluster_data[col] != "":
                    attributes[k] = cluster_data[col]
            attributes[synteny] = synteny
            attributes = ";".join([f"{k}={v}" for k,v in attributes.items()])
            row = ["pangenome", ".", feature, start, end, ".", strand, ".", attributes ]
            line = "\t".join([str(v) for v in row])
            outfile.write(line + "\n")


def get_ranges(numbers):
    """
    Author: user97370, bossylobster
    Source: https://stackoverflow.com/a/4629241
    """
    for _a, b in itertools.groupby(enumerate(numbers), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def collapse_ranges(ranges, min_len: int, min_gap: int):
    """Collapse a list of ranges into overlaps."""
    collapse = []
    curr_range = None
    # Collapse based on min gap
    for i, coord in enumerate(ranges):
        start, stop = coord
        if curr_range:
            gap_size = start - curr_range[1] - 1
            # Collapse tiny gap into previous fragment
            if gap_size < min_gap:
                curr_range = (curr_range[0], stop)
            # Otherwise, start new range
            else:
                collapse.append(curr_range)
                curr_range = (start, stop)
        else:
            curr_range = (start, stop)
        # Last one
        if curr_range not in collapse and i == len(ranges) - 1:
            collapse.append(curr_range)
    # Collapse based on min length size
    collapse = [(start, stop) for start,stop in collapse if (stop - start + 1) >= min_len]
    return collapse


def structural(
        clusters:     str,
        alignments:    str,
        outdir:        str = ".",
        prefix:        str = None,
        min_len:       int = 10,
        min_indel_len: int = 3,
        args:          str = None,
    ):
    """
    Extract structural variants from cluster alignments.

    Takes as input the summarized clusters TSV and their individual alignments.
    Outputs an Rtab file of structural variants.

    >>> structural(clusters="summarize.clusters.tsv", alignments="alignments")
    >>> structural(clusters="summarize.tsv", alignments="alignments", min_len=100, min_indel_len=10)

    :param clusters:      TSV file of clusters from summarize.
    :param alignments:    Directory of alignments from align.
    :param outdir:        Output directory.
    :param prefix:        Prefix for output files.
    :param min_len:       Minimum length of structural variants to extract.
    :param min_indel_len: Minimum length of gaps that should separate variants.
    :param args:          Str of additional arguments [not implemented]

    :return: Ordered Dictionary of structural variants.
    """

    from Bio import SeqIO

    clusters_path, alignments_dir, output_dir = clusters, alignments, outdir
    args = args if args != None else ""
    prefix = f"{prefix}." if prefix != None else ""

    # Check output directory
    check_output_dir(output_dir)

    variants = OrderedDict()

    # -------------------------------------------------------------------------
    # Read Clusters

    logging.info(f"Reading summarized clusters: {clusters_path}")
    all_samples = []

    with open(clusters_path) as infile:
        header = [line.strip() for line in infile.readline().split("\t")]
        # sample columns begin after dbxref_alt
        all_samples = header[header.index("dbxref_alt")+1:]
        lines = infile.readlines()
        for line in tqdm(lines):
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            cluster = data["cluster"]

            # Get mapping of sequence IDs to samples
            seq_to_sample = {}
            for sample in all_samples:
                sequence_ids = data[sample].split(",") if data[sample] != "" else []
                for sequence_id in sequence_ids:
                    seq_to_sample[sequence_id] = sample

            # Read cluster alignment from the file
            alignment_path = os.path.join(alignments_dir, f"{cluster}.aln")
            if not os.path.exists(alignment_path):
                # Alignment might not exist if user excluded singletons in align
                num_samples = len(list(set(seq_to_sample.values())))
                if num_samples == 1:
                    logging.debug(f"Singleton {cluster} alignment was not found: {alignment_path}")
                    continue
                # otherwise, this alignment should exist, stop here
                else:
                    msg = f"{cluster} alignment was not found: {alignment_path}"
                    logging.error(msg)
                    raise Exception(msg)

            alignment = OrderedDict()

            for record in SeqIO.parse(alignment_path, "fasta"):
                sample = seq_to_sample[record.id]
                if sample not in alignment:
                    alignment[sample] = OrderedDict()
                alignment[sample][record.id] = record.seq

            # Parse out the fragment pieces based on the pattern of '-'
            cluster_variants = OrderedDict()

            present_samples = set()
            for sample in alignment:
                sample_frag_ranges = []
                for seq in alignment[sample].values():
                    non_ambig_nuc = set([nuc for nuc in seq if nuc in NUCLEOTIDES])
                    if len(non_ambig_nuc) > 0:
                        present_samples.add(sample)
                    frag_ranges = [(start+1, stop+1) for start,stop in list(get_ranges([i for i,nuc in enumerate(seq) if nuc != "-"]))]
                    frag_ranges_collapse = collapse_ranges(ranges=frag_ranges, min_len=min_len, min_gap=min_indel_len)
                    sample_frag_ranges.append(frag_ranges_collapse)

                # Sort them for output consistency
                sample_frag_ranges = sorted(sample_frag_ranges)

                for frag in sample_frag_ranges:
                    # Simple presence                    
                    frag_text = f"{cluster}|structural:" + "_".join([f"{start}-{stop}" for start,stop in frag])
                    if frag_text not in cluster_variants: 
                        cluster_variants[frag_text] = []
                    if sample not in cluster_variants[frag_text]:
                        cluster_variants[frag_text].append(sample)
                    # Copy number (ex. 2X)
                    count = sample_frag_ranges.count(frag)
                    if count > 1:
                        copy_number_text = f"{frag_text}|{count}X"
                        if copy_number_text not in cluster_variants: 
                            cluster_variants[copy_number_text] = []
                        if sample not in cluster_variants[copy_number_text]:
                            cluster_variants[copy_number_text].append(sample)

            # If there is only one structure, it's not variant
            if len(cluster_variants) <= 1: continue
            # Samples that are truly missing will be given a "."
            missing_samples = [s for s in all_samples if s not in present_samples]
            for variant,samples in cluster_variants.items():
                if len(samples) == len(all_samples): continue
                variants[variant] = {"present": samples, "missing": missing_samples}

    # -------------------------------------------------------------------------
    # Write Structural Variants Rtab

    rtab_path = os.path.join(output_dir, f"{prefix}structural.Rtab")
    logging.info(f"Writing variants: {rtab_path}")

    with open(rtab_path, 'w') as outfile:
        header = ["Variant"] + all_samples
        outfile.write("\t".join(header) + "\n")     
        for variant,data in tqdm(variants.items()):
            observations = ["1" if s in data["present"] else "." if s in data["missing"] else "0" for s in all_samples]
            line = "\t".join([variant] + observations)
            outfile.write(line + "\n")

    return rtab_path


def snps(
        alignment:    str,
        bed:          str,
        consensus:    str,
        outdir:       str   = ".",
        prefix:       str   = None,
        structural:   str   = None,
        core:         float = 0.95,
        indel_window: int   = 0,
        snp_window:   int   = 0,
        args:         str   = None,
    ):
    """
    Extract SNPs from a pangenome alignment.

    Takes as input the pangenome alignment fasta, bed, and consensus file.
    Outputs an Rtab file of SNPs.

    >>> snps(alignment="pangenome.aln", bed="pangenome.bed", consensus="pangenome.consensus.fasta")
    >>> snps(alignment="pangenome.aln", bed="pangenome.bed", consensus="pangenome.consensus.fasta", structural="structural.Rtab", indel_window=10, snp_window=3)

    :param alignment:    FASTA file of the pangenome alignment from align.
    :param bed:          BED file of the pangenome coordinates from align.
    :param consensus:    FASTA file of the pangenome consensus from align.
    :param outdir:       Output directory.
    :param prefix:       Prefix for output files.
    :param structural:   Rtab file of structural variants from structural.
    :param core:         Core genome threshold for calling core SNPs.
    :param indel_window: Exclude SNPs that are within this proximity to indels.
    :param snp_window:   Exclude SNPs that are within this proximity to another SNP.
    :param args:         Str of additional arguments [not implemented]
    """

    alignment_path, bed_path, consensus_path, structural_path, output_dir = alignment, bed, consensus, structural, outdir
    args = args if args != None else ""
    prefix = f"{prefix}." if prefix != None else ""

    # Check output directory
    check_output_dir(output_dir)

    # -------------------------------------------------------------------------
    # Read Pangenome Bed

    clusters = {}
    logging.info(f"Reading bed: {bed_path}")
    with open(bed_path) as infile:
        lines = infile.readlines()
        for line in tqdm(lines):
            line = [l.strip() for l in line.split("\t")]
            start, end, name = (int(line[1]), int(line[2]), line[3])
            info = {n.split("=")[0]:n.split("=")[1] for n in name.split(";")}
            cluster = info["cluster"]
            synteny = info["synteny"]
            clusters[cluster] = {"start": start, "end": end, "synteny": synteny }

    # -------------------------------------------------------------------------
    # Read Pangenome Alignment

    alignment = {}
    all_samples = []
    logging.info(f"Reading alignment: {alignment_path}")
    with open(alignment_path) as infile: 
        records = infile.read().split(">")[1:]
        for record in tqdm(records):
            record_split = record.split("\n")
            sample = record_split[0]
            if sample not in all_samples:
                all_samples.append(sample)
            sequence = "".join(record_split[1:]).replace("\n", "").upper()
            alignment[sample] = sequence

    # -------------------------------------------------------------------------
    # Read Consensus Sequence

    representative = ""
    logging.info(f"Reading consensus sequence: {consensus_path}")
    with open(consensus_path) as infile:
        record = infile.read().split(">")[1:][0]
        record_split = record.split("\n")
        representative = "".join(record_split[1:]).replace("\n", "").upper()       

    alignment_len = len(alignment[all_samples[0]])

    # -------------------------------------------------------------------------
    # Read Optional Structural Rtab

    # Use the structural variants to locate the terminal ends of sequences

    structural = OrderedDict()
    if structural_path:
        logging.info(f"Reading structural: {structural_path}")        
        with open(structural_path) as infile:
            header = infile.readline().strip().split("\t")
            samples = header[1:]
            lines = infile.readlines()
            for line in tqdm(lines):
                row = [v.strip() for v in line.split("\t")]
                variant = row[0].split("|")
                cluster = variant[0]
                if cluster not in structural:
                    structural[cluster] = OrderedDict()
                # Cluster_1|structural|1-195_202-240 --> Terminal=(1,240)
                coords = variant[1].split(":")[1].split("_")
                start = int(coords[0].split("-")[0])
                end   = int(coords[len(coords) - 1].split("-")[1])
                for sample,observation in zip(samples, row[1:]):
                    if observation == "1":
                        if sample not in structural[cluster]:
                            structural[cluster][sample] = []
                        structural[cluster][sample].append((start, end))

    # -------------------------------------------------------------------------
    # Extract SNPs from Alignment
    # -------------------------------------------------------------------------

    logging.info("Extracting SNPs.")

    constant_sites = {n:0 for n in NUCLEOTIDES}

    snps_data = OrderedDict()

    # Iterate over alignment according to cluster
    for cluster,cluster_data in tqdm(clusters.items()):
        synteny = cluster_data["synteny"]
        cluster_start, cluster_stop = cluster_data["start"], cluster_data["end"]
        for i in range(cluster_start - 1, cluster_stop):
            # Extract nucleotides for all samples
            nuc = {s:alignment[s][i] for s in all_samples}
            nuc_non_ambig = [n for n in nuc.values() if n in NUCLEOTIDES]
            nuc_non_ambig_set = set(nuc_non_ambig)
            prop_non_ambig = len(nuc_non_ambig) / len(nuc)

            # Option #1. If all missing/ambiguous, skip over            
            if len(nuc_non_ambig_set) == 0:
                continue
            # Option #2. All the same/invariant
            elif len(nuc_non_ambig_set) == 1:
                # record for constant sites if it's a core site
                if prop_non_ambig >= core:
                    n = list(nuc_non_ambig)[0]
                    constant_sites[n] += 1
                continue
            # Option #3. Variant, process it further
            # Make snp positions 1-based (like VCF)
            pangenome_pos = i + 1
            cluster_pos = pangenome_pos - cluster_start
            # Treat the representative nucleotide as the reference
            ref = representative[i]
            alt = []
            genotypes = []
            sample_genotypes = {}

            # Filter on multi-allelic and indel proximity
            for s,n in nuc.items():

                # -------------------------------------------------------------
                # Indel proximity checking

                if indel_window > 0:
                    # Check window around SNP for indels
                    upstream_i = i - indel_window
                    if upstream_i < cluster_start:
                        upstream_i = cluster_start
                    downstream_i = i + indel_window
                    if downstream_i > (cluster_stop - 1):
                        downstream_i = (cluster_stop - 1)

                    # Check if we should truncate the downstream/upstream i by terminals
                    if cluster in structural and s in structural[cluster]:
                        for start,stop in structural[cluster][s]:
                            # Adjust start,stop coordinates from cluster to whole genome
                            start_i, stop_i = (cluster_start + start) - 1, (cluster_start + stop) - 1
                            # Original: 226-232, Terminal: 0,230 --> 226-230
                            if stop_i > upstream_i and stop_i < downstream_i:
                                downstream_i = stop_i          
                            # Original: 302-308, Terminal: 303,389 --> 303-308
                            if start_i > upstream_i and start_i < downstream_i:
                                upstream_i = start_i

                    context = alignment[s][upstream_i:downstream_i + 1]

                    # If indel is found nearby, mark this as missing/ambiguous
                    if "-" in context:
                        logging.debug(f"{cluster} {ref}{cluster_pos}{n} ({ref}{pangenome_pos}{n}) in {s} was filtered out due to indel proximity: {context}")                
                        genotypes.append(".")
                        nuc[s] = "."
                        sample_genotypes[s] = "./."
                        continue

                # -------------------------------------------------------------                    
                # sample nuc is different than ref
                if n != ref:
                    # handle ambiguous/missing
                    if n not in NUCLEOTIDES:
                        genotypes.append(".")
                        nuc[s] = "."
                        sample_genotypes[s] = "./."
                        continue

                    # add this as a new alt genotype
                    if n not in alt:
                        alt.append(n)
                    sample_genotypes[s] = "1/1"
                else:
                    sample_genotypes[s] = "0/0"

                genotypes.append(([ref] + alt).index(n))

            # Update our non-ambiguous nucleotides
            nuc_non_ambig = [n for n in nuc.values() if n in NUCLEOTIDES]
            nuc_non_ambig_set = set(nuc_non_ambig)

            genotypes_non_ambig = [g for g in genotypes if g != "."]

            # Check if it's still a variant position after indel filtering
            if len(alt) == 0 or len(genotypes_non_ambig) == 1:
                logging.debug(f"{cluster} {ref}{cluster_pos}{n} ({ref}{pangenome_pos}{n}) was filtered out as mono-allelic: {ref}")
                constant_sites[ref] += 1
                continue
            # If more than 1 ALT alleles, this is non-biallelic (multi-allelic) so we'll filter it out
            elif len(alt) > 1:
                logging.debug(f"{cluster} {ref}{cluster_pos}{n} ({ref}{pangenome_pos}{n}) was filtered out as multi-allelic: {','.join([ref] + alt)}")
                continue

            alt = list(alt)[0]
            # Use "." to indicate missing values (ex. Rtab format)
            observations = ["1" if n == alt else "." if n != ref else "0" for n in nuc.values()]
            # One more check if it's no longer variant...
            observations_non_ambig = set([o for o in observations if o == "1" or o == "0"])
            if len(observations_non_ambig) <= 1:
                observations_nuc = alt[0] if "1" in observations_non_ambig else ref if "0" in observations_non_ambig else "N"
                logging.debug(f"{cluster} {ref}{cluster_pos}{n} ({ref}{pangenome_pos}{n}) was filtered out as mono-allelic: {observations_nuc}")
                continue

            # This is a valid SNP! Calculate extra stats
            allele_frequency = len([n for n in nuc_non_ambig if n != ref]) / len(nuc)

            # Update data
            pangenome_snp = f"{ref}{pangenome_pos}{alt}"
            cluster_snp = f"{cluster}|snp:{ref}{cluster_pos}{alt}"
            snps_data[cluster_snp] = OrderedDict({
                "pangenome_snp": pangenome_snp,
                "prop_non_ambig"  : prop_non_ambig,
                "allele_frequency": allele_frequency,
                "cluster" : cluster,
                "synteny":  synteny,
                "observations": observations,
                "pangenome_pos": pangenome_pos,
                "cluster_pos": cluster_pos,
                "nuc": nuc,
                "ref": ref,
                "alt": alt,
                "sample_genotypes": sample_genotypes,
            })

    if snp_window > 0:
        snps_exclude = set()
        logging.info(f"Filtering out SNPs within {snp_window} bp of each other.")
        snps_order = list(snps_data.keys())
        for i,snp in enumerate(snps_order):
            cluster, coord = snps_data[snp]["cluster"], snps_data[snp]["cluster_pos"]
            if i > 0:
                prev_snp = snps_order[i-1]
                prev_cluster, prev_coord = snps_data[prev_snp]["cluster"], snps_data[prev_snp]["cluster_pos"],
                if (cluster == prev_cluster) and (coord - prev_coord <= snp_window):
                    snps_exclude.add(snp)
                    snps_exclude.add(prev_snp)
                    logging.debug(f"SNP {snp} and {prev_snp} are filtered out due to proximity <= {snp_window}.")
        snps_data = OrderedDict({snp:data for snp,data in snps_data.items() if snp not in snps_exclude})

    # -------------------------------------------------------------------------
    # Prepare Outputs

    snp_all_alignment = { s:[] for s in all_samples}
    snp_core_alignment = { s:[] for s in all_samples}
    snp_all_table   = open(os.path.join(output_dir, f"{prefix}snps.all.tsv"), 'w')
    snp_core_table  = open(os.path.join(output_dir, f"{prefix}snps.core.tsv"), 'w')
    snp_all_vcf   = open(os.path.join(output_dir, f"{prefix}snps.all.vcf"), 'w')
    snp_core_vcf  = open(os.path.join(output_dir, f"{prefix}snps.core.vcf"), 'w')
    snp_rtab_path   = os.path.join(output_dir, f"{prefix}snps.Rtab")
    snp_rtab        = open(snp_rtab_path, 'w')

    # TSV Table Header
    header = ["snp", "pangenome_snp", "prop_non_ambig", "allele_frequency", "cluster"] + all_samples
    snp_all_table.write("\t".join(header) + "\n")
    snp_core_table.write("\t".join(header) + "\n")
    
    # Rtab Header
    header = ["Variant"] + all_samples
    snp_rtab.write("\t".join(header) + "\n")

    # VCF Header
    all_samples_header = "\t".join(all_samples)
    header = textwrap.dedent(
        f"""\
        ##fileformat=VCFv4.2
        ##contig=<ID=pangenome,length={alignment_len}>
        ##INFO=<ID=CR,Number=0,Type=Flag,Description="Consensus reference allele, not based on real reference genome">
        ##INFO=<ID=C,Number=1,Type=String,Description="Cluster">
        ##INFO=<ID=CP,Number=1,Type=Int,Description="Cluster Position">
        ##INFO=<ID=S,Number=1,Type=String,Description="Synteny">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FILTER=<ID=IW,Description="Indel window proximity filter {indel_window}">
        ##FILTER=<ID=SW,Description="SNP window proximity filter {snp_window}">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{all_samples_header}
        """
    )
    snp_all_vcf.write(header)
    snp_core_vcf.write(header)

    logging.info(f"Finalizing all and core SNPs ({core}).")
    for snp,d in snps_data.items():
        # Table
        table_row = [snp] + [d["pangenome_snp"], d["prop_non_ambig"], d["allele_frequency"], d["cluster"]] + d["observations"]
        snp_all_table.write("\t".join(str(v) for v in table_row) + "\n") 
        # Alignment
        for s in all_samples:
            snp_all_alignment[s].append(d["nuc"][s])
        # Rtab
        rtab_row = [snp] + d["observations"]
        snp_rtab.write("\t".join(rtab_row) + "\n")
        # VCF
        pangenome_pos, cluster_pos = d["pangenome_pos"], d["cluster_pos"]
        cluster, synteny = d["cluster"], d["synteny"]
        ref, alt = d["ref"], d["alt"]
        info = f"CR;C={cluster};CP={cluster_pos},S={synteny}"
        genotypes = []
        for sample in all_samples:
            genotypes.append(d["sample_genotypes"][sample])
        vcf_row = ["pangenome", pangenome_pos, snp, ref, alt, ".", "PASS", info, "GT"] + genotypes
        vcf_line = "\t".join([str(v) for v in vcf_row])
        snp_all_vcf.write(vcf_line + "\n")

        # core SNPs
        if d["prop_non_ambig"] >= core:
            # Table
            snp_core_table.write("\t".join(str(v) for v in table_row) + "\n")
            for s in all_samples:
                snp_core_alignment[s].append(d["nuc"][s])
            # VCF
            snp_core_vcf.write(vcf_line + "\n")

    # -------------------------------------------------------------------------
    # Write SNP fasta alignments

    snp_all_alignment_path = os.path.join(output_dir, f"{prefix}snps.all.fasta")
    logging.info(f"Writing all SNP alignment: {snp_all_alignment_path}")
    with open(snp_all_alignment_path, 'w') as outfile:
        for sample in tqdm(all_samples):
            sequence = "".join(snp_all_alignment[sample])
            outfile.write(f">{sample}\n{sequence}\n")

    snp_core_alignment_path = os.path.join(output_dir, f"{prefix}snps.core.fasta")
    logging.info(f"Writing {core} core SNP alignment: {snp_core_alignment_path}")
    with open(snp_core_alignment_path, 'w') as outfile:
        for sample in tqdm(all_samples):
            sequence = "".join(snp_core_alignment[sample])
            outfile.write(f">{sample}\n{sequence}\n")            

    # -------------------------------------------------------------------------
    # Write constant sites

    constant_sites_path = os.path.join(output_dir, f"{prefix}snps.constant_sites.txt")
    logging.info(f"Writing constant sites: {constant_sites_path}")
    with open(constant_sites_path, 'w') as outfile:
        line = ",".join([str(v) for v in list(constant_sites.values())])
        outfile.write(line + "\n")

    # -------------------------------------------------------------------------
    # Cleanup

    snp_all_table.close()
    snp_core_table.close()
    snp_rtab.close()
    snp_all_vcf.close()
    snp_core_vcf.close()

    return snp_rtab_path


def presence_absence(
        clusters: str,
        outdir:   str = ".",
        prefix:   str = None,
        args:     str = None
    ):
    """
    Extract presence absence of summarized clusters.

    Takes as input the TSV file of summarized clusters from summarize.
    Outputs an Rtab file of cluster presence/absence.

    Examples:
    >>> presence_absence(clusters="summarize.clusters.tsv")

    :param clusters: Path to TSV of summarized clusters from summarize.
    :param outdir:   Output directory.
    :param prefix:   Prefix for output files.
    :param args:     Str of addtional arguments [not implemented].
    """

    clusters_path, output_dir = clusters, outdir
    args = args if args != None else ""
    prefix = f"{prefix}." if prefix != None else ""

    # Check output directory
    check_output_dir(output_dir)

    all_samples = []

    # -------------------------------------------------------------------------
    # Read Clusters

    logging.info(f"Reading summarized clusters: {clusters_path}")
    variants = OrderedDict()

    with open(clusters_path) as infile:
        header = [line.strip() for line in infile.readline().split("\t")]
        # sample columns begin after dbxref_alt
        all_samples = header[header.index("dbxref_alt")+1:]
        lines = infile.readlines()
        for line in tqdm(lines):
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            cluster = data["cluster"]
            variants[cluster] = []

            # Get samples with sequences for this cluster
            for sample in all_samples:
                sequence_ids = data[sample].split(",") if data[sample] != "" else []
                if len(sequence_ids) > 0:
                    variants[cluster].append(sample)
            # If the cluster is all samples, it's not "variable"
            if len(variants[cluster]) == len(all_samples):
                del variants[cluster]

    # -------------------------------------------------------------------------
    # Write Presence Absence Rtab

    rtab_path = os.path.join(output_dir, f"{prefix}presence_absence.Rtab")
    logging.info(f"Writing variants: {rtab_path}")

    with open(rtab_path, 'w') as outfile:
        header = ["Variant"] + all_samples
        outfile.write("\t".join(header) + "\n")
        for variant,samples in tqdm(variants.items()):
            variant = f"{variant}|presence_absence"
            observations = ["1" if s in samples else "0" for s in all_samples]
            line = "\t".join([variant] + observations)
            outfile.write(line + "\n")

    return rtab_path


def combine(
        rtab:   list,
        outdir: str = ".",
        prefix: str = None,
        args:   str = None,
    ):
    """
    Combine variants from multiple Rtab files.

    Takes as input a list of file paths to Rtab files. Outputs an Rtab file with
    the variants concatenated, ensuring consistent ordering of the sample columns.

    >>> combine(rtab=["snps.Rtab", "structural.Rtab", "presence_absence.Rtab"])

    :param rtab:   List of paths to input Rtab files to combine.
    :param outdir: Output directory.
    :param prefix: Prefix for output files.
    :param args:   Str of additional arguments [not implemented].
    """

    output_dir = outdir
    prefix = f"{prefix}." if prefix != None else ""
    # Check output directory
    check_output_dir(output_dir)
    
    all_samples = []
    rtab_path = os.path.join(output_dir, f"{prefix}combine.Rtab")
    logging.info(f"Writing combined Rtab: {rtab_path}")
    with open(rtab_path, 'w') as outfile:
        for file_path in rtab:
            if not file_path: continue
            logging.info(f"Reading variants: {file_path}")
            with open(file_path) as infile:
                header = infile.readline().strip().split("\t")
                samples = header[1:]
                if len(all_samples) == 0:
                    all_samples = samples
                    outfile.write("\t".join(["Variant"] + all_samples) + "\n")
                # Organize the sample order consistently
                lines = infile.readlines()
                for line in tqdm(lines):
                    row = [v.strip() for v in line.split("\t")]
                    variant = row[0]
                    observations = {k:v for k,v in zip(samples, row[1:])}
                    out_row = [variant] + [observations[s] for s in all_samples]
                    outfile.write("\t".join(out_row) + "\n")

    return rtab_path


def root_tree(
        tree:     str,
        outgroup: list = [],
        tree_format: str  = "newick",
        outdir:   str  = ".",
        prefix:   str  = None,
        args:     str  = None,
    ):
    """
    Root a tree on outgroup taxa.
    """

    tree_path, output_dir = tree, outdir
    prefix = f"{prefix}." if prefix != None else ""
    # Check output directory
    check_output_dir(output_dir)

    if type(outgroup) == str:
        outgroup = outgroup.split(",")
    elif type(outgroup) != list:
        msg = f"The outgroup must be either a list or CSV string: {outgroup}"
        logging.error(msg)
        raise Exception(msg)

    # Prioritize tree path
    logging.info(f"Reading tree: {tree_path}")
    tree = read_tree(tree_path, tree_format=tree_format)
    tree.is_rooted = True
    # Cleanup node/taxon labels
    for node in list(tree.preorder_node_iter()):
        if node.taxon:
            node.label = str(node.taxon)
            node.label = node.label.replace("\"", "").replace("'", "") if node.label != None else None

    logging.info(f"Rooting tree with outgroup: {outgroup}")

    tree_tips = [n.label for n in tree.leaf_node_iter()]
    if sorted(tree_tips) == sorted(outgroup):
        msg = "Your outgroup contains all taxa in the tree."
        logging.error(msg)
        raise Exception(msg)

    # sort the outgroup labels for easier checks later
    outgroup = sorted(outgroup)

    # Case 1: No outgroup requested, use first tip
    if len(outgroup) == 0:
        outgroup_node = [n for n in tree.leaf_node_iter()][0]
        logging.info(f"No outgroup requested, rooting on first tip: {outgroup_node.label}")
        edge_length = outgroup_node.edge_length / 2
    # Case 2: 1 or more outgroups requested
    else:
        # Case 2a: find a clade that is only the outgroup
        logging.info(f"Searching for outgroup clade: {outgroup}")
        outgroup_node = None
        for node in tree.nodes():
            node_tips =  sorted([n.label for n in node.leaf_iter()])
            if node_tips == outgroup:
                logging.info(f"Outgroup clade found.")
                outgroup_node = node
                break
        # Case 2b: find a clade that is only the ingroup
        if outgroup_node == None:
            logging.info(f"Searching for ingroup clade that does not include: {outgroup}")
            for node in tree.nodes():
                node_tips =  sorted([n.label for n in node.leaf_iter() if n.label in outgroup])
                if len(node_tips) == 0:
                    # If this is a taxon, make sure the outgroup is complete
                    if node.taxon != None:
                        remaining_taxa = sorted([t for t in tree_tips if t != node.label])
                        if remaining_taxa != outgroup:
                            continue
                    logging.info(f"Ingroup clade found.")
                    outgroup_node = node
                    break
        # If we couldn't find an outgroup, it's probably not monophyletic
        # If we created the tree with the same outgroup in IQTREE, this 
        # shouldn't be possible, but we'll handle just in case
        if outgroup_node == None:
            msg = "Failed to find the outgroup clade. Are you sure your outgroup is monophyletic?"
            logging.error(msg)
            raise Exception(msg)
        edge_length = outgroup_node.edge_length / 2

    tree.reroot_at_edge(outgroup_node.edge, update_bipartitions=True, length1=edge_length, length2=edge_length )    

    rooted_path = os.path.join(output_dir, f"{prefix}rooted.{tree_format}")
    logging.info(f"Writing rooted tree: {rooted_path}")
    # Suppress the [&R] prefix, which causes problems in downstream programs 
    tree.write(path=rooted_path, schema=tree_format, suppress_rooting=True)

    return rooted_path

def root_tree_midpoint(tree):
        """
        Reroots the tree at the midpoint of the longest distance between
        two taxa in a tree.

        This is a modification of dendropy's reroot_at_midpoint function:
        https://github.com/jeetsukumaran/DendroPy/blob/49cf2cc2/src/dendropy/datamodel/treemodel/_tree.py#L2616

        This is a utility function that is used by gwas for the kinship
        similarity matrix.

        :param tree: dendropy.datamodel.treemodel._tree.Tree
        :return:     dendropy.datamodel.treemodel._tree.Tree
        """
        from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
        pdm = PhylogeneticDistanceMatrix.from_tree(tree)

        tax1, tax2 = pdm.max_pairwise_distance_taxa()
        plen = float(pdm.patristic_distance(tax1, tax2)) / 2

        n1 = tree.find_node_with_taxon_label(tax1.label)
        n2 = tree.find_node_with_taxon_label(tax2.label)

        # Find taxa furthest from the root
        mrca = pdm.mrca(tax1, tax2)
        n1_dist, n2_dist = n1.distance_from_root(), n2.distance_from_root()

        # Start at furthest taxa, walk up to root, look for possible edge 
        # representing midpoint
        node = n1 if n1_dist >= n2_dist else n2
        remaining_len = plen

        while node != mrca:
            edge_length = node.edge.length
            remaining_len -= edge_length
            # Found midpoint node
            if remaining_len == 0:
                node = node._parent_node
                logging.info(f"Rerooting on midpoint node: {node}")
                tree.reroot_at_node(node, update_bipartitions=True)
                break
            # Found midpoint edge
            elif remaining_len < 0:
                l2 = edge_length + remaining_len
                l1 = edge_length - l2
                logging.info(f"Rerooting on midpoint edge, l1={l1}, l2={l2}")
                tree.reroot_at_edge(node.edge, update_bipartitions=True, length1=l1, length2=l2 ) 
                break            
            node = node._parent_node

        return tree


def read_tree(tree, tree_format:str="newick"):
    """Read a tree from file into an object, and clean up labels."""
    from dendropy import Tree

    tree = Tree.get(path=tree, schema=tree_format, tree_offset=0, preserve_underscores=True)
    tree.is_rooted = True

    # Give tip nodes the name of their taxon, cleanup tip labels
    for node in list(tree.preorder_node_iter()):
        if node.taxon:
           node.label = str(node.taxon)
        node.label = node.label.replace("\"", "").replace("'", "") if node.label != None else None

    return tree


def tree(
        alignment:      str,
        outdir:         str = ".",
        prefix:         str = None,
        threads:        int = 1,
        constant_sites: str = None,
        args:           str = TREE_ARGS,
    ):
    """
    Estimate a maximum-likelihood tree with IQ-TREE.

    Takes as input a multiple sequence alignment in FASTA format. If a SNP
    alignment is provided, an optional text file of constant sites can be 
    included for correction. Outputs a maximum-likelihood tree, as well as 
    additional rooted trees if an outgroup is specified in the iqtree args.

    tree(alignment="snps.core.fasta", constant_sites="snps.constant_sites.txt")
    tree(alignment="pangenome.aln", threads=4, args='--ufboot 1000 -o sample1')

    :param alignment:      Path to multiple sequence alignment for IQ-TREE.
    :param outdir:         Output directory.
    :param prefix:         Prefix for output files.
    :param threads:        CPU threads for IQ-TREE.
    :param constant_sites: Path to text file of constant site corrections.
    :param args:           Str of additional arguments for IQ-TREE.
    """

    alignment_path, constant_sites_path, output_dir = alignment, constant_sites, outdir
    args = args if args != None else ""
    prefix = f"{prefix}." if prefix != None else ""

    # Check output directory
    check_output_dir(output_dir)

    # Check for seed and threads conflict
    if int(threads) > 1 and "--seed" in args:
        logging.warning(
        f"""
        Using multiple threads is not guaranteed to give reproducible
        results, even when using the --seed argument:
        See this issue for more information: https://github.com/iqtree/iqtree2/discussions/233
        """)
    args += f" -T {threads}"

    # Check for outgroup rooting
    if "-o " in args:
        outgroup = extract_cli_param(args, "-o")
        logging.info(f"Tree will be rooted on outgroup: {outgroup}")
        outgroup_labels = outgroup.split(",")
    else:
        outgroup = None
        outgroup_labels = []

    # Check for constant sites input for SNP
    if constant_sites_path:
        with open(constant_sites_path) as infile:
            constant_sites = infile.read().strip()
        constant_sites_arg = f"-fconst {constant_sites}"
    else:
        constant_sites_arg = ""

    # -------------------------------------------------------------------------
    # Estimate maximum likelihood tree
    iqtree_prefix = os.path.join(output_dir, f"{prefix}tree")
    run_cmd(f"iqtree -s {alignment_path} --prefix {iqtree_prefix} {constant_sites_arg} {args}")
    tree_path = f"{iqtree_prefix}.treefile"

    # Check for outgroup warnings in log
    log_path = f"{iqtree_prefix}.log"
    with open(log_path) as infile:
        for line in infile:
            outgroup_msg = "Branch separating outgroup is not found"
            if outgroup_msg in line:
                msg = f"Branch separating outgroup is not found. Please check if your ({outgroup}) is monophyletic in: {tree_path}"
                logging.error(msg)
                raise Exception(msg)

    # -------------------------------------------------------------------------
    # Fix root position 

    tree_path = root_tree(tree=tree_path, outgroup=outgroup_labels, outdir=outdir, prefix=f"{prefix}tree")
    tree = read_tree(tree_path, tree_format="newick")
    tips = sorted([n.label for n in tree.leaf_node_iter()])
    # Sort child nodes, this extra sort should help us 
    # organize tip nodes by their branch length?
    for node in tree.preorder_node_iter():
        node._child_nodes.sort(
            key=lambda node: getattr(getattr(node, "taxon", None), "label", ""),
            reverse=True
        )

    tree_path = f"{iqtree_prefix}.rooted.nwk"
    logging.info(f"Writing rooted tree: {tree_path}")
    # Suppress the [&R] prefix, which causes problems in downstream programs 
    tree.write(path=tree_path, schema="newick", suppress_rooting=True)

    # -------------------------------------------------------------------------
    # Tidy branch supports

    node_i = 0

    branch_support_path = tree_path.replace(".treefile", ".branch_support.tsv").replace(".nwk", ".branch_support.tsv")
    logging.info(f"Extracting branch supports to table: {branch_support_path}")
    with open(branch_support_path, 'w') as outfile:
        header = ["node", "original", "ufboot", "alrt", "significant"]
        outfile.write('\t'.join(header) + '\n')
        for n in tree.preorder_node_iter():
            # Sanity check quotations
            if n.label:
                n.label = n.label.replace("'", "")
            if not n.is_leaf():
                if n.label:
                    original = n.label
                else:
                    original = "NA"
                support = n.label
                significant = ""
                ufboot = "NA"
                alrt = "NA"
                if support:
                    if "/" in support:
                        ufboot = float(support.split("/")[0])
                        alrt = float(support.split("/")[1])
                        if ufboot > 95 and alrt > 80:
                            significant = "*"
                    else:
                        ufboot = float(support.split("/")[0])
                        alrt = "NA"
                        if ufboot > 95:
                            significant = "*"
                n.label = f"NODE_{node_i}"
                row = [n.label, str(original), str(ufboot), str(alrt), significant]
                outfile.write('\t'.join(row) + '\n')
                node_i += 1

    # Tree V1: Branch support is replaced by node labels
    node_labels_path = tree_path.replace(".treefile", ".labelled_nodes.nwk").replace(".nwk", ".labelled_nodes.nwk")
    logging.info(f"Writing tree with node labels: {node_labels_path}")
    tree.write(path=node_labels_path, schema="newick", suppress_rooting=True) 
    
    # Tree V2: No branch support or node labels (plain)
    plain_path = node_labels_path.replace(".labelled_nodes.nwk", ".plain.nwk")
    logging.info(f"Writing plain tree: {plain_path}")
    for n in tree.preorder_node_iter():
        n.label = None
    tree.write(path=plain_path, schema="newick", suppress_rooting=True) 
    return tree_path


def get_delim(table:str):
    if table.endswith(".tsv") or table.endswith(".Rtab"):
        return "\t"
    elif table.endswith(".csv"):
        return ","
    else:
        msg = f"Unknown file extension of table: {table}"
        logging.error(msg)
        raise Exception(msg)

def table_to_rtab(
        table:  str,
        filter: str,
        outdir: str = ".",
        prefix: str = None,
        args:   str = None,
    ):
    """
    Convert a TSV/CSV table to an Rtab file based on regex filters.

    Takes as input a TSV/CSV table to convert, and a TSV/CSV of regex filters.
    The filter table should have a header with the fields: column, regex, name. Where column
    is the 'column' to search, 'regex' is the regular expression pattern, and
    'name' is how the output variant should be named in the Rtab.

    An example `filter.tsv` might look like this:

    column    regex         name        extra
    assembly  .*sample2.*   sample2     data that will be ignored
    lineage   .*2.*         lineage_2   more data

    Where the goal is to filter the assembly and lineage columns for particular values.

    >>> table_to_rtab(table="samplesheet.csv", filter="filter.tsv")
    """

    table_path, filter_path, output_dir = table, filter, outdir
    prefix = f"{prefix}." if prefix != None else ""

    check_output_dir(output_dir)

    logging.info(f"Checking delimiters.")
    table_delim  = get_delim(table_path)
    filter_delim = get_delim(filter_path)

    logging.info(f"Reading filters: {filter_path}")
    filters = OrderedDict()
    with open(filter_path) as infile:
        header = [h.strip() for h in infile.readline().split(filter_delim)]
        col_i, regex_i, name_i = [header.index(c) for c in ["column", "regex", "name"]]
        for line in infile:
            row = [l.strip() for l in line.split(filter_delim)]
            column,regex,name = row[col_i], row[regex_i], row[name_i]

            logging.debug(f"name: {name}, column: {column}, regex: {regex}")
            if name in filters:
                msg = f"Filter name {name} is not unique."
                logging.error(msg)
                raise Exception(msg)
            filters[name] = {k:v for k,v in zip(header, row)}

    logging.info(f"Searching input table: {table_path}")
    rtab_data = OrderedDict({name:{} for name in filters})
    all_samples = []
    with open(table_path) as infile:
        header = [h.strip() for h in infile.readline().split(table_delim)]      
        # Locate filter columns in the file
        filter_names = list(filters.keys())
        for name in filter_names:
            column = filters[name]["column"]
            if column in header:
                filters[name]["i"] = header.index(column)
            else:
                logging.warning(f"Filter column {column} is not present in the table.")
                del filters[name]
        if len(filters) == 0:
            msg = "No filters were found matching columns in the input table."
            logging.error(msg)
            raise Exception(msg)
        # Parse the table lines
        for line in infile:
            row = [l.strip() for l in line.split(table_delim)]
            sample = row[0]
            if sample not in all_samples:
                all_samples.append(sample)
            for name,data in filters.items():
                regex = data["regex"]
                text = row[data["i"]]
                val = 1 if re.match(regex, text) else 0
                rtab_data[name][sample] = val

    rtab_path = os.path.join(output_dir, f"{prefix}output.Rtab")
    extended_path = os.path.join(output_dir, f"{prefix}output.tsv")
    output_delim = "\t"
    logging.info(f"Writing Rtab output: {rtab_path}")
    logging.info(f"Writing extended output: {extended_path}")
    with open(rtab_path, 'w') as rtab_outfile:
        header = ["Variant"] + [str(s) for s in all_samples]
        rtab_outfile.write(output_delim.join(header) + "\n")
        with open(extended_path, 'w') as extended_outfile:
            first_filter = list(filters.keys())[0]
            header = ["sample"] + [col for col in filters[first_filter].keys() if col != "i"]
            extended_outfile.write(output_delim.join(header) + "\n")
    
            for name,data in rtab_data.items():
                row = [name] + [str(data[s]) for s in all_samples]
                line = output_delim.join([str(val) for val in row])
                rtab_outfile.write(line + "\n")

                positive_samples = [s for s,o in data.items() if str(o) == "1"]
                for sample in positive_samples:
                    row = [sample] + [v for k,v in filters[name].items() if k != "i"]
                    line = output_delim.join([str(val) for val in row])
                    extended_outfile.write(line + "\n")

def vcf_to_rtab(
        vcf:    str,
        bed:    str = None,
        outdir: str = ".",
        prefix: str = None,
        args:   str = None,
    ):
    """
    Convert a VCF file to an Rtab file.

    >>> vcf_to_rtab(vcf="snps.csv")
    """

    vcf_path, bed_path, output_dir = vcf, bed, outdir
    prefix = f"{prefix}." if prefix != None else ""
    check_output_dir(output_dir)

    bed = OrderedDict()
    if bed_path != None:
        logging.info(f"Reading BED: {bed_path}")
        with open(bed_path) as infile:
            lines = infile.readlines()
            for line in tqdm(lines):
                row = [l.strip() for l in line.split("\t")]
                _chrom, start, end, name = row[0], row[1], row[2], row[3]
                bed[name] = {"start": int(start) + 1, "end": int(end) + 1}

    logging.info(f"Reading VCF: {vcf_path}")

    all_samples = []

    rtab_path = os.path.join(output_dir, f"{prefix}output.Rtab")

    with open(vcf_path) as infile:
        with open(rtab_path, 'w') as outfile:
            logging.info(f"Counting lines in VCF.")
            lines = infile.readlines()

            logging.info(f"Parsing variants.")
            for line in tqdm(lines):
                if line.startswith("##"): continue
                # Write the Rtab header
                if line.startswith("#"):
                    header = [l.strip().replace("#", "") for l in line.split("\t")]
                    # Samples start after 'FORMAT'
                    samples_i = header.index("FORMAT") + 1
                    all_samples = header[samples_i:]
                    line = "\t".join(["Variant"] + [str(s) for s in all_samples])
                    outfile.write(line + "\n")
                    continue

                row = [l.strip() for l in line.split("\t")]
                data = {k:v for k,v in zip(header, row)}
                chrom, pos = f"{data['CHROM']}", int(data['POS'])
                ref = data['REF']
                alt = [nuc for nuc in data['ALT'].split(",") if nuc != "."]

                 # Skip multiallelic/missing
                if len(alt) != 1: continue

                alt = alt[0]
                observations = []
                for sample in all_samples:
                    genotype = data[sample]
                    if genotype == "0/0":
                        value = "0"
                    elif genotype == "1/1":
                        value = "1"
                    else:
                        value = "."

                    observations.append(value)

                for name,region in bed.items():
                    start, end = region["start"], region["end"]
                    if pos >= start and pos < end:
                        chrom = name
                        pos = (pos - start) + 1
                        break

                variant = f"{chrom}|snp:{ref}{pos}{alt}"

                line = "\t".join([variant] + observations)
                outfile.write(line + "\n")

def binarize(
        table:         str,
        column:        str,
        outdir:        str  = ".",
        prefix:        str  = None,
        output_delim:  str  = "\t",
        column_prefix: str  = None,
        transpose:     bool = False,
    ):
    """
    Convert a categorical column to multiple binary (0/1) columns.

    Takes as input a table (tsv or csv) as well as a column to binarize.

    >>> binarize(table="samplesheet.csv", column="lineage", output_delim=",")
    >>> binarize(table="samplesheet.csv", column="resistant", output_delim="\t", transpose=True)

    :param table:  Path to table (TSV or CSV)
    :param column: Name of column in table.
    :param outdir: Path to the output directory.
    """

    table_path, output_dir = table, outdir
    column_prefix = column_prefix if column_prefix else ""
    prefix = f"{prefix}." if prefix != None else ""
    output_ext = "tsv" if output_delim == "\t" else "csv" if output_delim == "," else "txt"
    
    # Check output directory
    check_output_dir(output_dir)
    input_delim = get_delim(table_path)

    logging.info(f"Reading table: {table_path}")

    all_samples = []
    values = {}

    with open(table_path) as infile:
        header = infile.readline().strip().split(input_delim)
        sample_i = 0
        column_i = header.index(column)
        for line in infile:
            line = [v.strip() for v in line.split(input_delim)]
            sample, val = line[sample_i], line[column_i]
            if sample not in all_samples:
                all_samples.append(sample)
            if val == "": continue
            if val not in values:
                values[val] = []
            values[val].append(sample)

    # Sort values
    values = OrderedDict({v:values[v] for v in sorted(values.keys())})

    binarize_path = os.path.join(output_dir, f"{prefix}binarize.{output_ext}")
    logging.info(f"Writing binarized output: {binarize_path}")
    with open(binarize_path, 'w') as outfile:
        first_col = column if transpose else header[sample_i]
        if transpose:
            out_header = [first_col] + all_samples
            outfile.write(output_delim.join(out_header) + "\n")
            for val,samples in values.items():
                observations = ["1" if sample in samples else "0" for sample in all_samples]
                row = [f"{column_prefix}{val}"] + observations
                outfile.write(output_delim.join(row) + "\n")
        else:
            out_header = [first_col] + [f"{column_prefix}{k}" for k in values.keys()]
            outfile.write(output_delim.join(out_header) + "\n")            
            for sample in all_samples:
                row = [sample]
                for val,val_samples in values.items():
                    row.append("1") if sample in val_samples else row.append("0")
                outfile.write(output_delim.join(row) + "\n")


def qq_plot(
        locus_effects:str,
        outdir: str = ".",
        prefix: str = None,
    ):
    """
    Modified version of pyseer's qq_plot function.

    Source: https://github.com/mgalardini/pyseer/blob/master/scripts/qq_plot.py
    """

    import matplotlib.pyplot as plt
    import numpy as np
    import statsmodels
    import statsmodels.api as sm
    import pandas as pd

    locus_effects_path, output_dir = locus_effects, outdir
    prefix = f"{prefix}." if prefix != None else ""

    logging.info(f"Reading locus effects: {locus_effects_path}")
    lrt_pvalues = []
    with open(locus_effects_path) as infile:
        header = [l.strip() for l in infile.readline().split("\t")]
        for line in infile:
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            lrt_pvalue = data["lrt-pvalue"]
            if lrt_pvalue != "":
                lrt_pvalues.append(float(lrt_pvalue))

    # Exclude p-value=0, which has no log
    y = -np.log10([val for val in lrt_pvalues if val != 0])
    x = -np.log10(np.random.uniform(0, 1, len(y)))

    # check for statsmodels version (issue #212)
    old_stats = False
    try:
        vmajor, vminor = statsmodels.__version__.split('.')[:2]
        if int(vmajor) == 0 and int(vminor) < 13:
            old_stats = True
        else:
            old_stats = False
    except:
        msg = "Failed to identify the statsmodel version, QQ plot will not be created."
        logging.warning(msg)
        return None

    if old_stats:
        xx = y
        yy = x
    else:
        xx = x
        yy = y

    # Plot
    logging.info(f"Creating QQ plot.")
    plt.figure(figsize=(4, 3.75))
    ax = plt.subplot(111)
 
    fig = sm.qqplot_2samples(xx,
                             yy,
                             xlabel='Expected $-log_{10}(pvalue)$',
                             ylabel='Observed $-log_{10}(pvalue)$',
                             line='45',
                             ax=ax)

    ax = fig.axes[0]
    ax.lines[0].set_color('k')
    ax.lines[0].set_alpha(0.3)

    # Handle inf values if a p-value was 0
    x_max = x.max() if not np.isinf(x.max()) else max([val for val in x if not np.isinf(val)])
    y_max = y.max() if not np.isinf(y.max()) else max([val for val in y if not np.isinf(val)])
    ax.set_xlim(-0.5, x_max + 0.5)
    ax.set_ylim(-0.5, y_max + 0.5)

    plt.tight_layout()

    plot_path = os.path.join(output_dir, f"{prefix}qq_plot.png")
    logging.info(f"Writing plot: {plot_path}")
    plt.savefig(plot_path, dpi=150)
    plt.close('all')

    return plot_path


def gwas(
        variants:       str,
        table:          str,
        column:         str,
        clusters:       str = None,        
        continuous:     str = False,
        tree:           str = None,
        midpoint:       bool= True,
        outdir:         str = ".",
        prefix:         str = None,
        lineage_column: str = None,
        exclude_missing: bool = False,
        threads:        int = 1,
        args:           str = GWAS_ARGS,
    ):
    """
    Run genome-wide association study (GWAS) tests with pyseer.

    Takes as input the TSV file of summarized clusters, an Rtab file of variants,
    a TSV/CSV table of phenotypes, and a column name representing the trait of interest.
    Outputs tables of locus effects, and optionally lineage effects (bugwas) if specified.

    >>> gwas(clusters='summarize.clusters.tsv', variants='combine.Rtab', table='samplesheet.csv', column='resistant')
    >>> gwas(clusters='summarize.clusters.tsv', variants='combine.Rtab', table='samplesheet.csv', column='lineage', args='--no-distances')
    >>> gwas(clusters='summarize.clusters.tsv', variants='combine.Rtab', table='samplesheet.csv', column='resistant', args='--lmm --min-af 0.05')

    :param variants:       Path to Rtab file of variants.
    :param table:          Path to TSV/CSV table of traits.
    :param column:         Column name of variable in table.
    :param clusters:       Path to TSV file of summarized clusters.    
    :param continuous:     Treat column as a continuous variable.
    :param tree:           Path to newick tree.
    :param midpoint:       True if the tree should be rerooted at the midpoint.
    :param outdir:         Output directory.
    :param prefix:         Prefix for output files.
    :param lineage_column: Name of lineage column in table (enables bugwas lineage effects).
    :param exclude_missing: Exclude samples missing phenotype data.
    :param threads:        CPU threads for pyseer.
    :param args:           Str of additional arguments for pyseer.
    """

    clusters_path, variants_path, table_path, tree_path, output_dir = clusters, variants, table, tree, outdir

    # Check output directory
    check_output_dir(output_dir)

    # Read in the table that contains the phenotypes (optional lineage column)
    logging.info(f"Reading table: {table_path}")
    with open(table_path) as infile:
        table = infile.read().replace(",", "\t")
    table_split = table.split("\n")
    table_header = table_split[0].split("\t")
    table_rows = [[l.strip() for l in line.split("\t") ] for line in table_split[1:] if line.strip() != ""]

    # Is the user didn't provide any custom args initialize as empty
    args = args if args != None else ""
    args += f" --cpu {threads}"
    prefix = f"{prefix}." if prefix != None else ""

    # Check for conflicts between input arguments
    if "--no-distances" not in args:
        if tree_path == None and "--lmm" in args:
            msg = "You must supply a phylogeny if you have not requested --no-distances."
            logging.error(msg)
            raise Exception(msg)
    else:
        if lineage_column:
            msg = "The param --no-distances cannot be used if a lineage column was specified."
            logging.error(msg)
            raise Exception(msg)

    if "--wg enet" in args:
        msg = "The whole genome elastic net model is not implemented yet."
        logging.error(msg)
        raise Exception(msg)

    # If a lineage column was specified, extract it from the table into 
    # it's own file for pyseer.
    lineage_path = None
    sample_i = 0
    
    if lineage_column:
        lineage_path = os.path.join(output_dir, f"{prefix}{column}.lineages.tsv")
        logging.info(f"Extracting lineage column {lineage_column}: {lineage_path}")
        with open(lineage_path, 'w') as lineage_file:
            lineage_i = table_header.index(lineage_column)
            for row in table_rows:
                sample,lineage = row[sample_i], row[lineage_i]
                lineage_file.write(f"{sample}\t{lineage}\n")

    # -------------------------------------------------------------------------
    # Table / Phenotype

    if column not in table_header:
        msg = f"Column {column} is not present in the table: {table_path}"
        logging.error(msg)
        raise Exception(msg)

    column_i = table_header.index(column)

    # Exclude samples that are missing phenotype if requested
    exclude_samples = set()
    if exclude_missing:
        logging.info("Excluding samples missing phenotype data.")
        table_filter = []
        for row in table_rows:
            sample = row[sample_i]
            value = row[column_i]
            if value == "NA" or value == "":
                exclude_samples.add(sample)
                continue
            table_filter.append(row)
        table_rows = table_filter
        logging.info(f"Number of samples excluded: {len(exclude_samples)}")

    # Convert the phenotypes table from csv to tsv for pyseer
    if table_path.endswith(".csv"):
        logging.info(f"Converting csv to tsv for pyseer: {table_path}")
        file_name = os.path.basename(table_path)
        table_tsv_path = os.path.join(output_dir, f"{prefix}{column}." + file_name.replace(".csv", ".tsv"))
        with open(table_tsv_path, 'w') as outfile:
            outfile.write("\t".join(table_header) + "\n")
            for row in table_rows:
                outfile.write("\t".join(row) + "\n")

    # Check for categorical, binary, continuous trait
    logging.info(f"Checking type of column: {column}")

    observations = OrderedDict()
    for row in table_rows:
        sample, val = row[sample_i], row[column_i]
        observations[sample] = val

    all_samples = list(observations.keys())
    unique_observations = sorted(list(set(observations.values())))

    if len(unique_observations) <=1 :
        msg = f"Column {column} must have at least two different values."
        logging.error(msg)
        raise Exception(msg)
    
    # Now remove missing data
    unique_observations = [o for o in unique_observations if o != ""]
    if unique_observations == ["0", "1"] or unique_observations == ["0"] or unique_observations == ["1"]:
        logging.info(f"Column {column} is a binary variable with values 0/1.")
        unique_observations = [column]
    elif continuous == True:
        logging.info(f"Treating column {column} as continuous.")
        try:
            _check = [float(o) for o in unique_observations]
        except:
            msg = f"Failed to convert all values of {column} as numeric."
            logging.error(msg)
            raise Exception(msg)
        unique_observations = [column]
        args += " --continuous"
    else:
        logging.info(f"Binarizing categorical column {column} to multiple binary columns.")
        binarize_prefix = f"{prefix}{column}"
        binarize(table = table_tsv_path, column = column, column_prefix=f"{column}_", outdir=output_dir, prefix=binarize_prefix)
        table_tsv_path = os.path.join(output_dir, f"{binarize_prefix}.binarize.tsv")
        unique_observations = [f"{column}_{o}" for o in unique_observations]

    # -------------------------------------------------------------------------
    # (optional) Cluster annotations

    clusters = OrderedDict()
    clusters_header = []

    if clusters_path != None:
        logging.info(f"Reading summarized clusters: {clusters_path}")
        with open(clusters_path) as infile:
            clusters_header = [line.strip() for line in infile.readline().split("\t")]

            lines = infile.readlines()
            for line in tqdm(lines):
                row = [l.strip() for l in line.split("\t")]
                # Store the cluster data, exclude all sample columns, they will
                # be provided by the variant input
                data = {
                    k:v for k,v in zip(clusters_header, row) 
                    if k not in exclude_samples and k not in all_samples
                }
                cluster = data["cluster"]
                clusters[cluster] = data

            # Exclude samples from clusters_header
            clusters_header = [
                col for col in clusters_header
                if col not in exclude_samples and col not in all_samples
            ]

    # -------------------------------------------------------------------------
    # Filter variants

    # Exclude missing samples and invariants
    # TBD: Might want to convert 'missing' chars all to '.'

    logging.info(f"Reading variants: {variants_path}")
    variants_filter_path = os.path.join(output_dir, f"{prefix}{column}.filter.Rtab")
    logging.info(f"Writing filtered variants: {variants_filter_path}")

    variants = {}
    variants_header = []
    exclude_variants = set()

    with open(variants_filter_path, "w") as outfile:
        with open(variants_path) as infile:
            header = infile.readline()
            # Save the header before filtering, to match up values to samples
            original_header = header.strip().split("\t")[1:]
            # Filter samples in the variant files
            variants_header = [s for s in original_header if s not in exclude_samples]
            outfile.write("\t".join(["Variant"] + variants_header) + "\n")

            # Filter invariants
            lines = infile.readlines()
            for line in tqdm(lines):
                row = [r.strip() for r in line.split("\t")]
                variant = row[0]
                observations = [
                    o 
                    for s,o in zip(original_header, row[1:]) 
                    if s not in exclude_samples
                ]
                unique_obs = set([o for o in observations])
                # Exclude simple invariants (less two possibilities) and invariants
                # with missing, ex. [".", "0"], [".", "1"]
                if (len(unique_obs) < 2) or (len(unique_obs) == 2 and ("." in unique_obs)):
                    exclude_variants.add(variant)
                else:
                    variants[variant] = observations
                    filter_row = [variant] + observations
                    outfile.write("\t".join(filter_row) + "\n")

    logging.info(f"Number of invariants excluded: {len(exclude_variants)}")
    logging.info(f"Number of variants included: {len(variants)}")
    variants_path = variants_filter_path

    # -------------------------------------------------------------------------
    # Distance Matrices

    patristic_path = os.path.join(output_dir, f"{prefix}{column}.patristic.tsv")
    kinship_path   = os.path.join(output_dir, f"{prefix}{column}.kinship.tsv")

    if "--lmm" in args:
        logging.info("Running pyseer with similarity matrix due to: --lmm")
        args += f" --similarity {kinship_path}"
    if lineage_path:
        logging.info("Running pyseer with patristic distances due to: --lineage-column")
        args += f" --distances {patristic_path}"

    # The user might want to run with --no-distances if population structure (lineage)
    # is the variable of interest.
    if ("--distances" not in args and 
        "--similarity" not in args):
        logging.info("Running pyseer with --no-distances")
        #args += " --no-distances"
    # Otherwise, we setup the distance matrices as needed
    else:
        logging.info(f"Creating distance matrices: {tree_path}")
        # A modified version of pyseer's phylogeny_distance function
        # Source: https://github.com/mgalardini/pyseer/blob/master/scripts/phylogeny_distance.py
        tree = read_tree(tree_path, tree_format="newick")

        if exclude_missing:
            logging.info(f"Excluding missing samples from the tree: {tree_path}")
            tree.prune_taxa_with_labels(labels=exclude_samples)
            tree.prune_leaves_without_taxa()

        tree_path = os.path.join(output_dir, f"{prefix}{column}.filter.nwk")
        tree.write(path=tree_path, schema="newick", suppress_rooting=True)

        if midpoint == True:
            tree_path = os.path.join(output_dir, f"{prefix}{column}.midpoint.nwk")
            logging.info(f"Rerooting tree at midpoint: {tree_path}")
            tree = root_tree_midpoint(tree)
            tree.write(path=tree_path, schema="newick", suppress_rooting=True)
        
        # Reread in tree after all that adjustments
        tree = read_tree(tree_path, tree_format="newick")

        patristic = OrderedDict()
        kinship   = OrderedDict()
        distance_matrix = tree.phylogenetic_distance_matrix()

        for n1 in tree.taxon_namespace:
            if n1.label not in all_samples:
                logging.warning(f"Sample {n1.label} is present in the tree but not in the phenotypes, excluding from distance matrices.")
                continue
            kinship[n1.label]  = kinship.get(n1.label, OrderedDict())
            patristic[n1.label] = patristic.get(n1.label, OrderedDict())
            for n2 in tree.taxon_namespace:
                if n2.label not in all_samples:
                    continue
                # Measure kinship/similarity as distance to root of MRCA
                # TBD: Investigating impact of precision/significant digits
                # https://github.com/mgalardini/pyseer/issues/286
                if n2.label not in kinship[n1.label].keys():
                    mrca = distance_matrix.mrca(n1, n2)
                    distance = mrca.distance_from_root()
                    kinship[n1.label][n2.label] = round(distance, 8)
                # Measure patristic distance between the nodes
                if n2.label not in  patristic[n1.label].keys():
                    distance = distance_matrix.patristic_distance(n1, n2)
                    patristic[n1.label][n2.label] = round(distance, 8)

        if len(patristic) == 0 or len(kinship) == 0:
            msg = "No samples are present in the distance matrices! Please check that your table sample labels match the tree."
            logging.error(msg)
            raise Exception(msg)

        logging.info(f"Saving patristic distances to: {patristic_path}")
        logging.info(f"Saving similarity kinship to: {kinship_path}")

        with open(patristic_path, 'w') as patristic_file:
            with open(kinship_path, 'w') as kinship_file:
                distance_samples = [s for s in all_samples if s in kinship]
                distance_header = "\t" + "\t".join(distance_samples)
                patristic_file.write(distance_header + "\n")
                kinship_file.write(distance_header + "\n")
                # Row order based on all_samples for consistency
                for s1 in distance_samples:
                    # Kinship/similarity
                    row = [s1] + [str(kinship[s1][s2]) for s2 in distance_samples]
                    kinship_file.write("\t".join(row) + "\n")
                    # Patristic
                    row = [s1] + [str(patristic[s1][s2]) for s2 in distance_samples]
                    patristic_file.write("\t".join(row) + "\n")

    # -------------------------------------------------------------------------
    #  Pyseer GWAS

    for new_column in unique_observations:
        logging.info(f"Running GWAS on column: {new_column}")

        log_path                 = os.path.join(output_dir, f"{prefix}{new_column}.pyseer.log")
        focal_path               = os.path.join(output_dir, f"{prefix}{new_column}.focal.txt")
        patterns_path            = os.path.join(output_dir, f"{prefix}{new_column}.locus_effects.patterns.txt")
        locus_effects_path       = os.path.join(output_dir, f"{prefix}{new_column}.locus_effects.tsv")
        locus_significant_path   = os.path.join(output_dir, f"{prefix}{new_column}.locus_effects.significant.tsv")
        locus_filter_path        = os.path.join(output_dir, f"{prefix}{new_column}.locus_effects.significant.filter.tsv")
        lineage_effects_path     = os.path.join(output_dir, f"{prefix}{new_column}.lineage_effects.tsv")
        lineage_significant_path = os.path.join(output_dir, f"{prefix}{new_column}.lineage_effects.significant.tsv")
        bonferroni_path          = os.path.join(output_dir, f"{prefix}{new_column}.locus_effects.bonferroni.txt")
        cmd = textwrap.dedent(
        f"""\
        pyseer {args}
          --pres {variants_path}
          --phenotypes {table_tsv_path}
          --phenotype-column {new_column}
          --output-patterns {patterns_path}
        """)
        if lineage_path:
            cmd += f"  --lineage --lineage-clusters {lineage_path} --lineage-file {lineage_effects_path}"
        cmd = cmd.replace("\n", " ")
        cmd_pretty = cmd + f" 1> {locus_effects_path} 2> {log_path}"
        logging.info(f"pyseer command: {cmd_pretty}")
        run_cmd(cmd, output=locus_effects_path, err=log_path, quiet=True, display_cmd=False)

        # Record which samples are 'positive' for this trait
        positive_samples = []
        negative_samples = []
        missing_samples = []
        with open(table_tsv_path) as infile:
            header = [l.strip() for l in infile.readline().split("\t")]
            for line in infile:
                row = [l.strip() for l in line.split("\t")]
                data = {k:v for k,v in zip(header, row)}
                sample, observation = data["sample"], data[new_column]
                if observation != ".":
                    if float(observation) != 0:
                        positive_samples.append(sample)
                    else:
                        negative_samples.append(sample)
                else:
                    missing_samples.append(sample)

        logging.info(f"{len(positive_samples)}/{len(all_samples)} samples are positive (>0) for {new_column}.")
        logging.info(f"{len(negative_samples)}/{len(all_samples)} samples are negative (=0) for {new_column}.")
        logging.info(f"{len(missing_samples)}/{len(all_samples)} samples are missing (.) for {new_column}.")

        with open(focal_path, 'w') as outfile:
            outfile.write("\n".join(positive_samples))
        # -------------------------------------------------------------------------
        # Significant threshold

        logging.info(f"Determining significance threshold.")

        with open(patterns_path) as infile:
            patterns = list(set([l.strip() for l in infile.readlines()]))
            num_patterns = len(patterns)
            logging.info(f"Number of patterns: {num_patterns}")
            bonferroni = 0.05 / float(num_patterns) if num_patterns > 0 else 0
            logging.info(f"Bonferroni threshold (0.05 / num_patterns): {bonferroni}")
            with open(bonferroni_path, 'w') as outfile:
                outfile.write(str(bonferroni) + "\n")

        # -------------------------------------------------------------------------
        # Extract variants into dict

        logging.info(f"Extracting variants.")

        gwas_variants = OrderedDict()
        # Keep track of the smallest non-zero pvalue for the log10 transformation
        # of 0 values
        min_pvalue = 1
        all_pvalues = set()

        with open(locus_effects_path) as infile:
            locus_effects_header = infile.readline().strip().strip().split("\t")
            for line in infile:
                row = [l.strip() for l in line.split("\t")]
                data = {k:v for k,v in zip(locus_effects_header, row)}
                data["notes"] = ",".join(sorted(data["notes"].split(",")))
                # Convert pvalue
                pvalue = data["lrt-pvalue"]
                if pvalue != "":
                    pvalue = float(pvalue)
                    # Keep track of the smallest non-zero pvalue
                    if pvalue != 0:
                        min_pvalue = min(pvalue, min_pvalue)
                    all_pvalues.add(pvalue)
                data["-log10(p)"] = ""
                data["bonferroni"] = bonferroni
                variant = data["variant"]
                gwas_variants[variant] = data
            locus_effects_header += ["-log10(p)", "bonferroni"]

        logging.info(f"Minimum pvalue observed: {min(all_pvalues)}")
        logging.info(f"Minimum pvalue observed (non-zero): {min_pvalue}")

        # -------------------------------------------------------------------------
        # -log10 transformation
        logging.info("Applying -log10 transformation.")

        # For p-values that are 0, we will use this tiny value instead
        # Remember, the min_pvalue is the minimum pvalue that is NOT 0
        small_val = float("1.00E-100")
        if small_val < min_pvalue:
            if 0 in all_pvalues:
                logging.warning(f"Using small float for pvalue=0 transformations: {small_val}")
        else:
            logging.warning(f"Using pvalue min for pvalue=0 transformations: {min_pvalue}")
            small_val = min_pvalue

        for variant, v_data in gwas_variants.items():
            pvalue = v_data["lrt-pvalue"]
            if pvalue == "": continue
            pvalue = float(pvalue)
            if pvalue == 0:
                logging.warning(f"{variant} has a pvalue of 0. log10 transformation will use: {small_val}")
                pvalue = small_val
            gwas_variants[variant]["-log10(p)"] = -math.log10(pvalue)

        # -------------------------------------------------------------------------
        # Annotate and sort the output files

        # Extract the cluster information for each variant
        logging.info(f"Extracting cluster identifiers.")
        locus_effects = OrderedDict()

        for variant,v_data in gwas_variants.items():
            variant_split = variant.split("|")
            cluster = variant_split[0]
            row = list(v_data.values())
            if cluster not in locus_effects:
                locus_effects[cluster] = []
            locus_effects[cluster].append(row)

        logging.info(f"Adding cluster annotations and variant observations.")
        locus_effects_rows = []
        match_header = ["match", "mismatch", "score"]
        for cluster in tqdm(locus_effects):
            cluster_row = [clusters[cluster][col] for col in clusters_header]
            for effect_row in locus_effects[cluster]:
                effect_data = {k:v for k,v in zip(locus_effects_header, effect_row)}
                try:
                    beta = float(effect_data["beta"])
                except ValueError:
                    beta = ""
                variant = effect_row[0]
                variant_row = variants[variant]
                variant_data = {k:v for k,v in zip(variants_header, variant_row)}
                # A simple summary statistic of how well the variant observations match the phenotype
                denom = len(all_samples)                
                # Presence of variant should 'increase' or 'decrease' phenotype
                if beta == "" or beta > 0:
                    v1, v2 = ["1", "0"]
                else:
                    v1, v2 = ["0", "1"]
                match_num = len([
                    k for k,v in variant_data.items() 
                    if (v == v1 and k in positive_samples) or (v == v2 and k in negative_samples)
                ])
                mismatch_num = len([
                    k for k,v in variant_data.items()
                    if (v == v1 and k not in positive_samples) or (v == v2 and k not in negative_samples)
                ])

                match_value = f"{match_num}/{denom}|{round(match_num / denom, 2)}"
                mismatch_value = f"{mismatch_num}/{denom}|{round(mismatch_num / denom, 2)}"
                match_score = (match_num/denom) - (mismatch_num/denom)
                # TBD What if we are on operating on lrt-filtering-fail p values?
                # We have no way to know the beta direction, and which way the score should go...
                if beta == "":
                    match_score = abs(match_score)
                match_row = [match_value, mismatch_value, round(match_score, 2)]
                locus_effects_rows += [effect_row + match_row + cluster_row + variant_row]

        logging.info(f"Sorting locus effects by lrt-pvalue: {locus_effects_path}")
        locus_pvalues = OrderedDict()
        locus_null = OrderedDict()
        for i,row in enumerate(locus_effects_rows):
            data = {k:v for k,v in zip(locus_effects_header, row)}
            variant, lrt_pvalue = data["variant"], data["lrt-pvalue"]
            if lrt_pvalue == "":
                if lrt_pvalue not in locus_null:
                    locus_null[lrt_pvalue] = []
                locus_null[lrt_pvalue].append((variant, i))
            else:
                lrt_pvalue = float(lrt_pvalue)
                if lrt_pvalue not in locus_pvalues:
                    locus_pvalues[lrt_pvalue] = []
                locus_pvalues[lrt_pvalue].append((variant, i))
        locus_pvalues = sorted(locus_pvalues.items(), key=lambda kv: (kv[0], kv[1]))
        locus_null = sorted(locus_null.items(), key=lambda kv: (kv[0], kv[1]))
        locus_order = []
        for data in locus_pvalues + locus_null:
            for variant, i in data[1]:
                locus_order.append(i)
        locus_effects_rows = [locus_effects_rows[i] for i in locus_order]

        # Save the results
        with open(locus_effects_path, 'w') as outfile:
            header = "\t".join(locus_effects_header + match_header + clusters_header + variants_header)
            outfile.write(header + "\n")
            for row in locus_effects_rows:
                line = "\t".join([str(r) for r in row])
                outfile.write(line + "\n")

        logging.info(f"Sorting patterns: {patterns_path}")
        with open(patterns_path) as infile:
            lines = sorted([l for l in infile.readlines()])
        with open(patterns_path, 'w') as outfile:
            outfile.write("".join(lines))

        if lineage_path:
            logging.info(f"Sorting lineage effects: {lineage_effects_path}")
            with open(lineage_effects_path) as infile:
                header = infile.readline()
                lines = sorted([l for l in infile.readlines()])
            with open(lineage_effects_path, 'w') as outfile:
                outfile.write(header + "".join(lines))

        # -------------------------------------------------------------------------
        # Locus Effects

        logging.info("Identifying significant locus effects.")
        with open(locus_effects_path) as infile:
            header_line = infile.readline()
            header = header_line.strip().split("\t")     
            with open(locus_significant_path, 'w') as significant_outfile:
                significant_outfile.write(header_line)    
                with open(locus_filter_path, 'w') as filter_outfile:
                    filter_outfile.write(header_line)
                    for line in infile:
                        row = [l.strip() for l in line.split("\t")]
                        data = {k:v for k,v in zip(header, row)}
                        lrt_pvalue, notes = data["lrt-pvalue"], data["notes"]
                        if lrt_pvalue != "" and float(lrt_pvalue) < bonferroni:
                            significant_outfile.write(line)
                            if "bad-chisq" not in notes and "high-bse" not in notes:
                                filter_outfile.write(line)

        # -------------------------------------------------------------------------
        # Lineage effects

        if lineage_path:
            logging.info("Identifying significant lineage effects.")
            with open(lineage_effects_path) as infile:
                header_line = infile.readline()
                header = header_line.strip().split("\t")
                pvalue_i = header.index("p-value")         
                with open(lineage_significant_path, 'w') as significant_outfile:
                    significant_outfile.write(header_line)
                    for line in infile:
                        row = [l.strip() for l in line.split("\t")]
                        pvalue = float(row[pvalue_i])
                        if pvalue < bonferroni:
                            significant_outfile.write(line)

        # -------------------------------------------------------------------------
        # QQ Plot

        if len(locus_effects_rows) == 0:
            logging.info("Skipping QQ Plot as no variants were observed.")
        else:
            logging.info("Creating QQ Plot.")
            _plot_path = qq_plot(
                locus_effects = locus_effects_path,
                outdir        = output_dir,
                prefix        = f"{prefix}{new_column}"
            )


def text_to_path(text, tmp_path="tmp.svg", family="Roboto", size=16, clean=True):
    """
    Convert a string of text to SVG paths.

    :param text:     String of text.
    :param tmp_path: File path to temporary file to render svg.
    :param family:   Font family.
    :param size:     Font size.
    :param clean:    True if temporary file should be deleted upon completion.
    """
    
    import cairo
    from xml.dom import minidom
    from svgpathtools import parse_path


    # Approximate and appropriate height and width for the tmp canvas
    # This just needs to be at least big enough to hold it
    w = size * len(text)
    h = size * 2
    # Render text as individual glyphs
    with cairo.SVGSurface(tmp_path, w, h) as surface:
        context = cairo.Context(surface)
        context.move_to(0, h/2)
        context.set_font_size(size)
        context.select_font_face(family)
        context.show_text(text)

    # Parse text data and positions from svg DOM
    doc = minidom.parse(tmp_path)

    # Keep track of overall text bbox
    t_xmin, t_xmax, t_ymin, t_ymax = None, None, None, None
    glyphs = []

    # Absolute positions of each glyph are under <use>
    # ex. <use xlink:href="#glyph-0-1" x="11.96875" y="16"/>
    for u in doc.getElementsByTagName('use'):
        glyph_id = u.getAttribute("xlink:href").replace("#", "")
        # Get absolute position of the glyph
        x = float(u.getAttribute("x"))
        y = float(u.getAttribute("y"))
        # Get relative path data of the glyph
        for g in doc.getElementsByTagName('g'):
            g_id = g.getAttribute("id")
            if g_id != glyph_id: continue
            for path in g.getElementsByTagName('path'):
                d = path.getAttribute("d")
                p = parse_path(d)
                # Get glyph relative bbox
                xmin, xmax, ymin, ymax = [c for c in p.bbox()]
                # Convert bbox to absolute coordinates
                xmin, xmax, ymin, ymax = xmin + x, xmax + x, ymin + y, ymax + y
                # Update the text bbox
                if t_xmin == None:
                    t_xmin, t_xmax, t_ymin, t_ymax = xmin, xmax, ymin, ymax
                else:
                    t_xmin = min(t_xmin, xmin)
                    t_xmax = max(t_xmax, xmax)
                    t_ymin = min(t_ymin, ymin)
                    t_ymax = max(t_ymax, ymax)
                bbox = xmin, xmax, ymin, ymax

                w, h = xmax - xmin, ymax - ymin
                glyph = {"id": glyph_id, "x": x, "y": y, "w": w, "h": h, "d": d, "bbox": bbox}
                glyphs.append(glyph)

    if not t_xmax:
        logging.error(f"Failed to render text: {text}")
        raise Exception(f"Failed to render text: {text}")
    # Calculate the final dimensions of the entire text
    w, h = t_xmax - t_xmin, t_ymax - t_ymin
    bbox = t_xmin, t_xmax, t_ymin, t_ymax

    result = {"text": text, "glyphs": glyphs, "w": w, "h": h, "bbox": bbox }

    if clean:
        os.remove(tmp_path)
    
    return result


def linear_scale(value, value_min, value_max, target_min=0.0, target_max=1.0):
    if value_min == value_max:
        return value
    else:
        return (value - value_min) * (target_max - target_min) / (value_max - value_min) + target_min

def manhattan(
        gwas:         str,
        bed:          str,
        outdir:       str   = ".",
        prefix:       str   = None,
        font_size:    int   = 16,
        font_family:  str   ="Roboto",
        margin:       int   = 20,
        width:        int   = 1000,
        height:       int   = 500,
        png_scale:    float = 2.0,
        prop_x_axis:  bool = False,
        ymax:        float = None,
        max_blocks:   int = 20,
        syntenies:    bool = ["all"],
        clusters:     bool = ["all"],
        variant_types: list = ["all"],
        args:         str = None,
    ):
    """
    Plot the distribution of variant p-values across the genome.

    Takes as input a table of locus effects from the subcommand and a bed
    file such as the one producted by the align subcommand. Outputs a 
    manhattan plot in SVG and PNG format.

    >>> manhattan(gwas="locus_effects.tsv", bed="pangenome.bed")
    >>> manhattan(gwas="locus_effects.tsv", bed="pangenome.bed", syntenies=["chromosome"], clusters=["pbpX"], variant_types=["snp", "presence_absence"])
            
    :param gwas:          Path to tsv table of locus effects from gwas subcommand.
    :param bed:           Path to BED file of coordinates.
    :param outdir:        Output directory.
    :param prefix:        Output prefix.
    :param width:         Width in pixels of plot.
    :param height:        Height in pixels of plot.
    :param margin:        Size in pixels of plot margins.
    :param font_size:     Font size.
    :param font_family:   Font family.
    :param png_scale:     Float that adjusts png scale relative to svg.
    :param prop_x_axis:   Scale x-axis based on the length of the synteny block.
    :param ymax:          Maximum value for the y-axis -log10(p).
    :param sytenies:      Names of synteny blocks to plot.
    :param clusters:      Names of clusters to plot.
    :params variant_types: Names of variant types to plot.
    """

    import numpy as np

    logging.info(f"Importing cairo.")
    import cairosvg

    gwas_path, bed_path, output_dir = gwas, bed, outdir
    prefix = f"{prefix}." if prefix != None else ""
    check_output_dir(outdir)

    # Type conversions fallback
    width, height, margin = int(width), int(height), int(margin)
    png_scale = float(png_scale)
    if ymax != None: ymax = float(ymax)

    # handle space values if given
    syntenies = syntenies if type(syntenies) != str else syntenies.split(" ")
    syntenies = [str(s) for s in syntenies]

    clusters = clusters if type(clusters) != str else clusters.split(" ")
    clusters = [str(c) for c in clusters]
    
    variant_types = variant_types if type(variant_types) != str else variant_types.split(" ")
    variant_types = [str(v) for v in variant_types]

    plot_data = OrderedDict()

    # -------------------------------------------------------------------------
    # Parse BED (Synteny Blocks)

    pangenome_length = 0

    logging.info(f"Reading BED coordinates: {bed_path}")
    with open(bed_path) as infile:
        for line in infile:
            row = [l.strip() for l in line.split("\t")]
            # Adjust 0-based bed coordinates back to 1 base
            c_start, c_end, name = int(row[1]) + 1, int(row[2]), row[3]
            pangenome_length = c_end
            c_length = (c_end - c_start) + 1
            # ex. cluster=geneA;synteny=443
            info = {n.split("=")[0]:n.split("=")[1] for n in name.split(";")}
            cluster, synteny = info["cluster"], info["synteny"]

            if synteny not in plot_data:
                plot_data[synteny] = {
                    "pangenome_pos": [c_start, c_end],
                    "length": c_length,
                    "variants": OrderedDict(),
                    "clusters": OrderedDict(),
                    "ticks": OrderedDict(),
                    "prop": 0.0,
                }

            # update the end position and length of the current synteny block
            plot_data[synteny]["pangenome_pos"][1] = c_end
            plot_data[synteny]["length"] = (
                c_end - plot_data[synteny]["pangenome_pos"][0]
            ) + 1

            s_start = plot_data[synteny]["pangenome_pos"][0]
            c_data = {
                "pangenome_pos": [c_start, c_end],
                "synteny_pos": [1 + c_start - s_start, 1+ c_end - s_start],
                "cluster_pos": [1, c_length],
                "length": c_length,
            }
            plot_data[synteny]["clusters"][cluster] = c_data

    # -------------------------------------------------------------------------
    # Parse GWAS (Variants)

    logging.info(f"Reading GWAS table: {gwas_path}")

    alpha = None
    log10_pvalues = set()

    with open(gwas_path) as infile:

        header = [l.strip() 
        for l in infile.readline().split("\t")]
        if "cluster" not in header:
            msg = "GWAS table must contain cluster annotations for manhattan plot."
            logging.error(msg)
            raise Exception(msg)

        lines = infile.readlines()

        if len(lines) == 0:
            msg = "GWAS table contains no variants."
            logging.warning(msg)
            return 0

        for line in tqdm(lines):
            row = [l.strip() for l in line.split("\t")]
            data = {k:v for k,v in zip(header, row)}
            variant, cluster, synteny = data["variant"], data["cluster"], data["synteny"]

            alpha = float(data["bonferroni"]) if alpha == None else alpha

            if data["lrt-pvalue"] == "": continue
            pvalue = float(data["lrt-pvalue"])
            log10_pvalue = float(data["-log10(p)"])
            log10_pvalues.add(log10_pvalue)

            # try to find matching synteny block by name
            c_data = None
            if synteny in plot_data:
                c_data = plot_data[synteny]["clusters"][cluster]
            # Otherwise, find by cluster
            else:
                for s,s_data in plot_data.items():
                    if cluster in s_data["clusters"]:
                        synteny = s
                        c_data = s_data["clusters"][cluster]
                        
            if c_data == None:
                msg = f"Synteny block {synteny} for {cluster} is not present in the BED file."
                logging.error(msg)
                raise Exception(msg)
    
            s_data = plot_data[synteny]
            s_start, s_end = s_data["pangenome_pos"]
            c_start, c_end = c_data["synteny_pos"]

            # 3 different coordinates systems
            cluster_coords = []
            synteny_coords = []
            pangenome_coords = []

            if "|snp" in variant:
                variant_type = "snp"
                snp = variant.split("|")[1].split(":")[1]
                pos = int("".join([c for c in snp if c.isnumeric()]))
                cluster_pos = [pos, pos]
                cluster_coords.append(cluster_pos)
            elif "|presence_absence" in variant:
                variant_type = "presence_absence"
                cluster_pos = c_data["cluster_pos"]
                cluster_coords.append(cluster_pos)
            elif "|structural" in variant:
                variant_type = "structural"
                variant_split = variant.split("|")[1].split(":")[1]
                for interval in variant_split.split("_"):
                    cluster_pos = [int(v) for v in interval.split("-")]
                    cluster_coords.append(cluster_pos)
            else:
                logging.warning(f"Skipping unknown variant type: {variant}")
                continue

            # Apply filter
            if "all" not in variant_types and variant_type not in variant_types: continue
            if "all" not in syntenies and synteny not in syntenies: continue
            if "all" not in clusters and cluster not in clusters: continue

            # Convert cluster coords to synteny and pangenome coords
            for pos in cluster_coords:
                synteny_pos = [c_start + pos[0], c_start + pos[1]]
                synteny_coords.append(synteny_pos)
                pangenome_pos = [s_start + synteny_pos[0], s_start + synteny_pos[1]]
                pangenome_coords.append(pangenome_pos)

            plot_data[synteny]["variants"][variant] = {
                "variant": variant,
                "synteny": synteny,
                "synteny_pos": int(data["synteny_pos"]),
                "cluster": cluster,
                "pvalue": pvalue,
                "-log10(p)": log10_pvalue,
                "variant_h2": float(data["variant_h2"]),
                "type": variant_type,                
                "cluster_coords": cluster_coords,
                "synteny_coords": synteny_coords,
                "pangenome_coords": pangenome_coords,
            }

    alpha_log10 = -math.log10(alpha)
    if ymax == None:
        max_log10 = math.ceil(max(log10_pvalues) if max(log10_pvalues) > alpha_log10 else alpha_log10)
    else:
        max_log10 = ymax

    min_log10 = 0
    log10_vals= np.arange(min_log10, max_log10 + 1, max_log10 / 4)

    # -------------------------------------------------------------------------
    # Exclude synteny blocks with no variants
    plot_data = OrderedDict({k:v for k,v in plot_data.items() if len(v["variants"]) > 0})

    # Adjust the total length based on just the synteny blocks we've observed
    total_length = sum(
        [sum([c["length"] for c in s["clusters"].values()])
        for s in plot_data.values()]
    )

    observed_clusters = set()
    for s in plot_data.values():
        for v in s["variants"].values():
            observed_clusters.add(v["cluster"])
    for s,s_data in plot_data.items():
        plot_data[s]["clusters"] = OrderedDict({
            k:v for k,v in plot_data[s]["clusters"].items()
            if k in observed_clusters
        })

    # If there's just a single cluster, adjust the total length
    if len(observed_clusters) == 1:
        only_synteny = list(plot_data.keys())[0]
        only_cluster = list(observed_clusters)[0]
        total_length = plot_data[only_synteny]["clusters"][only_cluster]["length"]
    elif "all" not in clusters:
        total_length = sum(
            [sum([c["length"] for c in s["clusters"].values()])
            for s in plot_data.values()]
        )

    if len(plot_data) == 0:
        msg = "No variants remain after filtering."
        logging.warning(msg)
        return None

    for s,s_data in plot_data.items():
        s_length = sum([c["length"] for c in s_data["clusters"].values()])
        s_prop = s_length / total_length
        plot_data[s]["prop"] = s_prop

    # -------------------------------------------------------------------------
    # Phandango data

    phandango_path = os.path.join(output_dir, f"{prefix}phandango.plot")
    logging.info(f"Creating phandango input: {phandango_path}")

    with open(phandango_path, 'w') as outfile:
        header = ["#CHR", "SNP", "BP", "minLOG10(P)", "log10(p)", "r^2"]
        outfile.write("\t".join(header) + "\n")

        for synteny,s_data in plot_data.items():
            pangenome_pos = s_data["pangenome_pos"]            
            for variant, v_data in s_data["variants"].items():
                log10_pvalue = v_data["-log10(p)"]
                variant_h2 = v_data["variant_h2"]
                cluster = v_data["cluster"]
                c_data = s_data["clusters"][cluster]

                # If a snp, we take the actual position
                # for presence/absence or structural, we take the start
                if "|snp" in variant:
                    pos = v_data["pangenome_coords"][0][0]
                else:
                    pos = c_data["pangenome_pos"][0]
                row = ["pangenome", variant, pos, log10_pvalue, log10_pvalue, variant_h2]
                line = "\t".join([str(v) for v in row])
                outfile.write(line + "\n")

    # -------------------------------------------------------------------------
    # Y-Axis Label Dimensions

    logging.info(f"Calculating y-axis label dimensions.")

    # Main y-axis label
    y_axis_text = "-log10(p)"
    y_axis_label = text_to_path(text=y_axis_text, size=font_size, family=font_family)
    y_axis_label_h = y_axis_label["h"]
    y_axis_label_x = margin + y_axis_label_h
    y_axis_label["rx"] = y_axis_label_x
    # To set the y position, we'll need to wait until after x axis labels are done

    y_tick_hmax = 0
    y_tick_wmax = 0
    y_tick_labels = OrderedDict()

    for val in log10_vals:
        label = text_to_path(text=str(val), size=font_size * 0.75, family=font_family)
        w, h = label["w"], label["h"]
        y_tick_hmax, y_tick_wmax = max(y_tick_hmax, h), max(y_tick_wmax, w)        
        y_tick_labels[val] = label

    # We now know some initial dimensions
    tick_len = y_tick_hmax / 2
    tick_pad = tick_len
    y_axis_x = y_axis_label_x + y_axis_label_h + y_tick_wmax + tick_pad + tick_len
    y_axis_y1 = margin
    # To set the y2 position, we'll need to wait until after x axis labels are done

    # -------------------------------------------------------------------------
    # X-Axis Label Dimensions

    logging.info(f"Calculating x-axis label dimensions.")

    # X-axis main label
    if len(plot_data) > max_blocks:
        x_axis_text = "Pangenome"
    elif len(plot_data) == 1:
        x_axis_text = str(list(plot_data.keys())[0])
        if len(observed_clusters) == 1:
            x_axis_text = str(list(observed_clusters)[0])
    else:
        x_axis_text = "Synteny Block"

    x_axis_label = text_to_path(text=x_axis_text, size=font_size, family=font_family)
    x_axis_label_w, x_axis_label_h = x_axis_label["w"], x_axis_label["h"]
    x_axis_label_y = height - margin

    x_axis_x1 = y_axis_x
    x_axis_x2 = width - margin
    x_axis_w  = x_axis_x2 - x_axis_x1
    x_axis_step = x_axis_w / len(plot_data)

    x_axis_label_x = x_axis_x1 + (x_axis_w / 2) - (x_axis_label_w / 2)

    x_tick_hmax = 0
    x_tick_wmax = 0
    
    start_coord, end_coord = 0, total_length


    if len(plot_data) > max_blocks or len(plot_data) == 1:

        if len(plot_data) > max_blocks:
            logging.info(f"Plot data exceeds max_blocks of {max_blocks}, enabling pangenome coordinates.")
            start_coord, end_coord = 0, pangenome_length
        elif len(plot_data) == 1:
            only_synteny = list(plot_data.keys())[0]
            logging.info(f"Single synteny blocked detected, enabling synteny coordinates.")
            # For one cluster, we might have variable start positions
            if len(observed_clusters) == 1:
                only_cluster = list(observed_clusters)[0]
                c_data = plot_data[only_synteny]["clusters"][only_cluster]
                start_coord, end_coord = c_data["cluster_pos"]
            # Round down/up to nearest 10?
            start_coord = start_coord - (start_coord % 10)
            end_coord = end_coord + (end_coord % 10)

        label = text_to_path(text="0123456789", size=font_size * 0.75, family=font_family)
        h = label["h"]
        # Put 1/2 h spacing in between
        max_labels = x_axis_w / (2 * h)
        # Prefer either 10 or 4 ticks
        if max_labels >= 10:
            max_labels = 10
        elif max_labels >= 4:
            max_labels = 4

        step_size = 1
        num_labels = max_labels + 1
        # 10, 50, 100, 500, 1000...
        while num_labels > max_labels:
            if str(step_size).startswith("5"):
                step_size = step_size * 2
            else:
                step_size = step_size * 5
            num_labels = total_length / step_size

        x_tick_vals = list(range(start_coord, end_coord + 1, step_size))

        x_axis_step = x_axis_w / num_labels
        for val in x_tick_vals:
            text = str(val)
            label = text_to_path(text=text, size=font_size * 0.75, family=font_family)
            w, h = label["w"], label["h"]
            x_tick_hmax, x_tick_wmax = max(x_tick_hmax, h), max(x_tick_wmax, w)
            plot_data[synteny]["ticks"][text] = {"label": label}
    else:
        for synteny in plot_data:
            label = text_to_path(text=str(synteny), size=font_size * 0.75, family=font_family)
            w, h = label["w"], label["h"]
            x_tick_hmax, x_tick_wmax = max(x_tick_hmax, h), max(x_tick_wmax, w)
            plot_data[synteny]["ticks"][synteny] = {"label": label}


    # We now have enough information to finalize the axes coordinates
    x_axis_y = x_axis_label_y - (x_axis_label_h * 2) - x_tick_wmax - tick_pad - tick_len    
    y_axis_y2 = x_axis_y
    y_axis_h = y_axis_y2 - y_axis_y1
    y_axis_label_y = y_axis_y1 + (y_axis_h / 2)
    y_axis_label["ry"] = y_axis_label_y

    # -------------------------------------------------------------------------
    logging.info(f"Positioning labels.")

    for i_g,g in enumerate(x_axis_label["glyphs"]):
        x_axis_label["glyphs"][i_g]["x"] = x_axis_label_x + g["x"]
        x_axis_label["glyphs"][i_g]["y"] = x_axis_label_y

    for i_g,g in enumerate(y_axis_label["glyphs"]):
        y_axis_label["glyphs"][i_g]["x"] = y_axis_label_x + g["x"]
        y_axis_label["glyphs"][i_g]["y"] = y_axis_label_y

    # -------------------------------------------------------------------------
    # X-Axis

    logging.info("Creating x-axis.")

    prev_x = x_axis_x1
    for synteny in plot_data:
        s_data = plot_data[synteny]

        if len(plot_data) == 1:
            plot_data[synteny]["x1"] = x_axis_x1
            plot_data[synteny]["x2"] = x_axis_x2
    
        for tick,t_data in plot_data[synteny]["ticks"].items():

            if prop_x_axis == True:
                xw = x_axis_w * s_data["prop"]
            else:
                xw = x_axis_step

            if len(plot_data) > 1 and len(plot_data) <= max_blocks:
                plot_data[synteny]["x1"] = prev_x
                plot_data[synteny]["x2"] = prev_x + xw

            # Pixel coordinates of this synteny block
            if len(plot_data) == 1 or len(plot_data) > max_blocks:
                x = linear_scale(int(tick), start_coord, end_coord, x_axis_x1, x_axis_x2)
            else:
                x = prev_x + (xw / 2)

            # Tick Label
            label = t_data["label"]
            y = x_axis_y + tick_len + label["w"] + (tick_pad)
            t_data["label"]["rx"] = x
            t_data["label"]["ry"] = y

            # Tick Line
            t_data["line"] = {"x": x, "y1": x_axis_y, "y2": x_axis_y + tick_len}
            # Center label on the very last glyph
            center_glyph = label["glyphs"][-1]
            label_y = y + (center_glyph["h"] / 2)
            for i_g,g in enumerate(label["glyphs"]):
                t_data["label"]["glyphs"][i_g]["x"] = x + g["x"]
                t_data["label"]["glyphs"][i_g]["y"] = label_y
            
            plot_data[synteny]["ticks"][tick] = t_data

            prev_x += xw
    
    # -------------------------------------------------------------------------
    # Y-Axis

    logging.info("Creating y-axis.")

    y_axis = OrderedDict({k:{"label": v, "line": None} for k,v in y_tick_labels.items()})

    for val,label in y_tick_labels.items():
        y = linear_scale(val, max_log10, min_log10, y_axis_y1, y_axis_y2)
        x1, x2 = y_axis_x - tick_len, y_axis_x

        # Tick Line
        y_axis[val]["line"] = { "x1": x1, "x2": x2, "y": y }

        # Tick Label
        y_axis[val]["label"] = label   
        center_glyph = label["glyphs"][-1]
        label_y = y + (center_glyph["h"] / 2)
        label_x = x1 - tick_pad - label["w"]
    
        for i_g,g in enumerate(label["glyphs"]):
            label["glyphs"][i_g]["x"] = label_x + g["x"]
            label["glyphs"][i_g]["y"] = label_y

        y_axis[val]["label"] = label

    # -------------------------------------------------------------------------
    # Render

    radius = 2

    svg_path = os.path.join(output_dir, f"{prefix}plot.svg")
    logging.info(f"Rendering output svg ({width}x{height}): {svg_path}")
    with open(svg_path, 'w') as outfile:
        header = textwrap.dedent(
        f"""\
        <svg 
            version="1.1"
            xmlns="http://www.w3.org/2000/svg"
            xmlns:xlink="http://www.w3.org/1999/xlink"
            preserveAspectRatio="xMidYMid meet"
            width="{width}"
            height="{height}"
            viewbox="0 0 {width} {height}">
        """)
        outfile.write(header)

        indent = "    "

        # White canvas background
        outfile.write(f"{indent * 1}<g id='Background'>\n")
        background = f"{indent * 2}<rect width='{width}' height='{height}' x='0' y='0' style='fill:white;stroke-width:1;stroke:white'/>"
        outfile.write(background + "\n")
        outfile.write(f"{indent * 1}</g>\n")

        # Axes
        outfile.write(f"{indent * 1}<g id='axis'>\n")

        # ---------------------------------------------------------------------
        # X Axis

        outfile.write(f"{indent * 2}<g id='x'>\n")

        # Background panels
        if len(plot_data) <= max_blocks:
            outfile.write(f"{indent * 3}<g id='Background Panels'>\n")
            for i_s,synteny in enumerate(plot_data):
                if i_s % 2 == 0: continue
                x1, x2 = plot_data[synteny]["x1"], plot_data[synteny]["x2"]
                y = margin
                w = x2 - x1
                h = y_axis_h
                rect = f"{indent * 2}<rect width='{w}' height='{h}' x='{x1}' y='{y}' style='fill:grey;fill-opacity:0.10'/>"
                outfile.write(rect + "\n")
            outfile.write(f"{indent * 3}</g>\n")

        # X axis line
        outfile.write(f"{indent * 3}<g id='line'>\n")
        line = f"{indent * 4}<path d='M {x_axis_x1} {x_axis_y} H {x_axis_x2}' style='stroke:black;stroke-width:2;fill:none'/>"
        outfile.write(line + "\n")
        outfile.write(f"{indent * 3}</g>\n")

        # X axis label
        outfile.write(f"{indent * 3}<g id='label'>\n")
        for g in x_axis_label["glyphs"]:
            line = f"{indent * 4}<path transform='translate({g['x']},{g['y']})' d='{g['d']}'/>"
            outfile.write(line + "\n")
        outfile.write(f"{indent * 3}</g>\n")

        # X axis ticks
        outfile.write(f"{indent * 3}<g id='ticks'>\n")
        for synteny in plot_data:
            for tick,t_data in plot_data[synteny]["ticks"].items():
                # Tick line
                t = t_data["line"]
                line = f"{indent * 4}<path d='M {t['x']} {t['y1']} V {t['y2']}' style='stroke:black;stroke-width:1;fill:none'/>"
                outfile.write(line + "\n")
                # Tick label
                label = t_data["label"]
                rx, ry = label["rx"], label["ry"]
                for g in label["glyphs"]:
                    line = f"{indent * 4}<path transform='rotate(-90, {rx}, {ry}) translate({g['x']},{g['y']})' d='{g['d']}'/>"
                    outfile.write(line + "\n")

        outfile.write(f"{indent * 3}</g>\n")

        # Close x-axis
        outfile.write(f"{indent * 2}</g>\n")

        # ---------------------------------------------------------------------
        # Y Axis

        outfile.write(f"{indent * 2}<g id='y-axis'>\n")

        # Y-axis line
        outfile.write(f"{indent * 3}<g id='line'>\n")
        line = f"{indent * 4}<path d='M {y_axis_x} {y_axis_y1} V {y_axis_y2}' style='stroke:black;stroke-width:2;fill:none'/>"
        outfile.write(line + "\n")
        outfile.write(f"{indent * 3}</g>\n")

        # Y axis label
        outfile.write(f"{indent * 3}<g id='label'>\n")
        for g in y_axis_label["glyphs"]:
            rx, ry = y_axis_label["rx"], y_axis_label["ry"]
            line = f"{indent * 4}<path transform='rotate(-90, {rx}, {ry}) translate({g['x']},{g['y']})' d='{g['d']}'/>"
            outfile.write(line + "\n")
        outfile.write(f"{indent * 3}</g>\n")        

        outfile.write(f"{indent * 3}<g id='ticks'>\n")
        # Ticks
        for t,t_data in y_axis.items():
            # Tick line
            t = t_data["line"]
            line = f"{indent * 3}<path d='M {t['x1']} {t['y']} H {t['x2']}' style='stroke:black;stroke-width:1;fill:none'/>"
            outfile.write(line + "\n")
            # Tick Label
            label = t_data["label"]
            for g in label["glyphs"]:
                line = f"{indent * 4}<path transform='translate({g['x']},{g['y']})' d='{g['d']}'/>"
                outfile.write(line + "\n")
        outfile.write(f"{indent * 3}</g>\n")

        # Close y-axis
        outfile.write(f"{indent * 2}</g>\n")

        # Close all Axes
        outfile.write(f"{indent * 1}</g>\n")


        # ------------------------------------------------------------- 
        # Data
        outfile.write(f"{indent * 1}<g id='variants'>\n")

        for i_s, synteny in enumerate(plot_data):
            s_data = plot_data[synteny]
            if len(plot_data) <= max_blocks:
                s_x1, s_x2 = s_data["x1"], s_data["x2"]
            else:
                s_x1, s_x2 = x_axis_x1, x_axis_x2
            s_len = s_data["length"]

            variants = plot_data[synteny]["variants"].values()
            variants_order = [v for v in variants if v["type"] == "presence_absence"]
            variants_order += [v for v in variants if v["type"] == "structural"]
            variants_order += [v for v in variants if v["type"] == "snp"]
            opacity = 0.50

            for v_data in variants_order:
                # If we're plotting according to pangenome coordinates
                # we use the panGWAS green color
                if len(plot_data) >= max_blocks:
                    f = "#356920"
                # Otherwise we alternate between the classic blue and orange
                else:
                    f = "#242bbd" if i_s % 2 == 0 else "#e37610"
                vy = linear_scale(v_data["-log10(p)"], min_log10, max_log10, y_axis_y2, y_axis_y1)
                variant = v_data["variant"].replace("'", "")

                # In many browsers, the title will become hover text!
                hover_text = []
                for k,v in v_data.items():
                    if "coords" not in k:
                        text = f"{k}: {v}"
                    else:
                        coords = ", ".join(["-".join([str(xs) for xs in x]) for x in v])
                        text = f"{k}: {coords}"
                    hover_text.append(text)
                hover_text = "\n".join(hover_text)
                title = f"{indent * 4}<title>{hover_text}</title>"

                outfile.write(f"{indent * 2}<g id='{variant}'>\n")

                if len(observed_clusters) == 1:
                    s_start = start_coord
                    s_end = end_coord
                else:
                    s_start = 1
                    s_end = s_len

                # Decide on the coordinate system
                if len(plot_data) > max_blocks:
                    coordinate_system = "pangenome"
                    s_start, s_end = 0, pangenome_length
                else:
                    coordinate_system = "synteny"

                pos = v_data[f"{coordinate_system}_coords"][0][0]
                vx = linear_scale(pos, s_start, s_end, s_x1, s_x2)

                # Circle
                if v_data["type"] == "snp":
                    circle = f"{indent * 3}<circle cx='{vx}' cy='{vy}' r='{radius}' style='fill:{f};fill-opacity:{opacity}'>"
                    outfile.write(circle + "\n")
                    outfile.write(title + "\n")
                    outfile.write(f"{indent * 3}</circle>" + "\n")

                # Square
                elif v_data["type"] == "presence_absence":
                    pos = v_data[f"{coordinate_system}_coords"][0][0]
                    vx = linear_scale(pos, s_start, s_end, s_x1, s_x2)
                    rect = f"{indent * 3}<rect width='{radius*2}' height='{radius*2}' x='{vx}' y='{vy - radius}' style='fill:{f};fill-opacity:{opacity}'>"
                    outfile.write(rect + "\n")
                    outfile.write(title + "\n")
                    outfile.write(f"{indent * 3}</rect>" + "\n")

                # Diamond
                elif v_data["type"] == "structural":
                    rx, ry = vx + radius, vy + radius
                    rect = f"{indent * 3}<rect transform='rotate(-45, {rx}, {ry})' width='{radius*2}' height='{radius*2}' x='{vx}' y='{vy - radius}' style='fill:{f};fill-opacity:{opacity}'>"
                    outfile.write(rect + "\n")
                    outfile.write(title + "\n")
                    outfile.write(f"{indent * 3}</rect>" + "\n")

                outfile.write(f"{indent * 2}</g>\n")
        outfile.write(f"{indent * 1}</g>\n")

        # p-value threshold
        outfile.write(f"{indent * 1}<g id='alpha'>\n")
        x1, x2 = x_axis_x1, x_axis_x2
        y = linear_scale(alpha_log10, min_log10, max_log10, y_axis_y2, y_axis_y1)
        line = f"{indent * 2}<path d='M {x1} {y} H {x2}' style='stroke:grey;stroke-width:1;fill:none' stroke-dasharray='4 4'/>"
        outfile.write(line + "\n")
        outfile.write(f"{indent * 1}</g>\n")

        outfile.write("</svg>" + "\n")

    png_path = os.path.join(output_dir, f"{prefix}plot.png")
    logging.info(f"Rendering output png ({width}x{height}): {png_path}")
    cairosvg.svg2png(url=svg_path, write_to=png_path, output_width=width, output_height=height, scale=png_scale)

    return svg_path


def heatmap(
        tree: str=None,
        tree_format: str="newick",
        gwas: str=None,
        rtab: str=None,
        outdir: str = ".",
        prefix: str=None,
        focal: str=None,
        min_score: float = None,
        tree_width=100,
        margin=20,
        root_branch=10,
        tip_pad=10,
        font_size=16,
        font_family="Roboto",
        png_scale=2.0,
        heatmap_scale=1.5,
        palette={"presence_absence": "#140d91", "structural": "#7a0505", "snp": "#0e6b07", "unknown": "#636362"},
        args: str = None,
    ):
    """
    Plot a tree and/or a heatmap of variants.

    Takes as input a newick tree and/or a table of variants. The table can be either
    an Rtab file, or the locus effects TSV output from the gwas subcommand.
    If both a tree and a table are provided, the tree will determine the sample order 
    and arrangement. If just a table is provided, sample order will follow the 
    order of the sample columns. A TXT of focal sample IDs can also be supplied
    with one sample ID per line. Outputs a plot in SVG and PNG format.

    >>> plot(tree="tree.rooted.nwk")
    >>> plot(rtab="combine.Rtab")
    >>> plot(gwas="resistant.locus_effects.significant.tsv")
    >>> plot(tree="tree.rooted.nwk", rtab="combine.Rtab", focal="focal.txt")
    >>> plot(tree="tree.rooted.nwk", gwas="resistant.locus_effects.significant.tsv")
    >>> plot(tree="tree.rooted.nwk, tree_width=500)

    :param tree:          Path to newick tree.
    :param gwas:          Path to tsv table of locus effects from gwas subcommand.
    :param rtab:          Path to Rtab file of variants.
    :param prefix:        Output prefix.
    :param focal:         Path to text file of focal sample IDs to highlight.
    :param min_score:     Filter GWAS variants for a minimum score.
    :param tree_width:    Width in pixels of tree.
    :param margin:        Size in pixels of plot margins.
    :param root_branch:   Width in pixels of root branch.
    :param tip_pad:       Pad in pixels between tip labels and branches.
    :param font_size:     Font size.
    :param font_family:   Font family.
    :param png_scale:     Float that adjusts png scale relative to svg.
    :param heatmap_scale: Float that adjusts heatmap box scale relative to text.
    :param palette:       Dict of variant types to colors.
    """

    tree_path, gwas_path, focal_path, rtab_path, output_dir = tree, gwas, focal, rtab, outdir
    prefix = f"{prefix}." if prefix != None else ""

    if not tree_path and not gwas_path and not rtab_path:
        msg = "Either a tree (--tree), a GWAS table (--gwas) or an Rtab --rtab) must be supplied."
        logging.error(msg)
        raise Exception(msg)

    elif gwas_path and rtab_path:
        msg = "A GWAS table (--gwas) is mutually exclusive with an Rtab (--rtab) file."
        logging.error(msg)
        raise Exception(msg)

    # Check output directory
    check_output_dir(output_dir)

    logging.info(f"Importing cairo.")
    import cairosvg

    # Todo: If we want to color by beta/p-value
    # logging.info(f"Importing matplotlib color palettes.")
    # import matplotlib as mpl
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import rgb2hex

    # -------------------------------------------------------------------------
    # Parse Focal Samples

    focal = []
    if focal_path:
        logging.info(f"Parsing focal samples: {focal_path}")
        with open(focal_path) as infile:
            for line in infile:
                sample = line.strip()
                if sample != "" and sample not in focal:
                    focal.append(sample)

    # -------------------------------------------------------------------------
    # Parse Input Tree

    if tree_path:
        logging.info(f"Parsing {tree_format} tree: {tree_path}")
        tree = read_tree(tree=tree_path, tree_format=tree_format)
        tree.is_rooted = True
        tree.ladderize()
    else:
        root_branch = 0
        tree_width = 0

    # -------------------------------------------------------------------------
    # Node labels

    node_labels = {}
    node_labels_wmax = 0
    node_labels_wmax_text = None
    node_labels_hmax = 0
    node_labels_hmax_text = None
    all_samples = []
    tree_tips = []
    tree_nodes = []

    # Option 1: Tips from the tree
    if tree_path:
        # Cleanup up labels, dendropy sometimes inserts quotations
        logging.info(f"Parsing tree labels.")
        node_i = 0
        for node in list(tree.preorder_node_iter()):
            # Tip labels
            if node.taxon:
                node.label = str(node.taxon)
                logging.debug(f"Tree Tip: {node.label}")
            # Internal node labels
            else:
                if not node.label or "NODE_" not in node.label:
                    if node_i == 0:
                        logging.warning(f"Using 'NODE_<i> nomenclature for internal node labels.")
                    node.label = f"NODE_{node_i}"
                    node_i += 1
            node.label = node.label.replace("\"", "").replace("'", "")

        # Keep a list of tip labels
        tree_tips = [n.label for n in tree.leaf_node_iter()]
        tree_nodes = [n.label for n in tree.preorder_node_iter()]
        all_samples = copy.deepcopy(tree_tips)

    # Option 2: Tips from GWAS columns
    if gwas_path:
        logging.info(f"Parsing gwas labels: {gwas_path}")
        with open(gwas_path) as infile:
            header = infile.readline().strip().split("\t")
            # Sample IDs come after the dbxref_alt column (if we used summarized clusters)
            # otherwise they come after the "score" column
            if "dbxref_alt" in header:
                samples = header[header.index("dbxref_alt")+1:]
            else:
                samples = header[header.index("score")+1:]
            for sample in samples:
                if sample not in all_samples:
                    all_samples.append(sample)

    # Option 3: Tips from Rtab columns
    if rtab_path:
        logging.info(f"Parsing rtab labels: {rtab_path}")
        with open(rtab_path) as infile:
            header = infile.readline().strip().split("\t")
            # Sample IDs are all columns after the first
            samples = header[1:]
            for sample in samples:
                if sample not in all_samples:
                    all_samples.append(sample)

    # Figure out the largest sample label, we'll need to know
    # this to position elements around it.
    logging.info(f"Calculating sample label dimensions.")

    for node in tqdm(all_samples):
        label = text_to_path(text=node, size=font_size, family=font_family)
        node_labels[node] = label
        w, h = label["w"], label["h"]
        if w > node_labels_wmax:
            node_labels_wmax = w
            node_labels_wmax_text = node
        if h > node_labels_hmax:
            node_labels_hmax = h
            node_labels_hmax_text = node

    node_labels_wmax = math.ceil(node_labels_wmax)
    node_labels_hmax = math.ceil(node_labels_hmax)

    logging.info(f"The widest sample label is {math.ceil(node_labels_wmax)}px: {node_labels_wmax_text}")
    logging.info(f"The tallest sample label is {math.ceil(node_labels_hmax)}px: {node_labels_hmax_text}")

    # -------------------------------------------------------------------------
    # Variant Labels

    # Figure out the largest variant label for the heatmap

    variant_labels_wmax = 0
    variant_labels_wmax_text = ""
    variants = OrderedDict()

    # Option 1. Variants from GWAS rows

    # The score can range from -1 to +1
    # We will make a gradient palette based on the score for the fill
    score_min = None
    score_max = None
    if gwas_path:
        logging.info(f"Parsing GWAS variants: {gwas_path}")
        with open(gwas_path) as infile:
            header = [line.strip() for line in infile.readline().split("\t")]
            for line in infile:
                row = [l.strip() for l in line.split("\t")]
                data = {k:v for k,v in zip(header, row)}
                variant = f"{data['variant']}"
                if "score" in data:
                    score = float(data["score"])
                    if min_score != None and score < min_score:
                        logging.info(f"Excluding {variant} due to score {score} < {min_score}")
                        continue
                    if score < 0:
                        logging.warning(f"Converting {variant} score of {score} to 0.00")
                        score = 0.0
                    data["score"] = score
                    score_min = score if score_min == None else min(score_min, score)
                    score_max = score if score_max == None else max(score_max, score)
                variants[variant] = {"data": data }

    # Option 2. Variants from Rtab rows
    elif rtab_path:
        logging.info(f"Parsing Rtab variants: {rtab_path}")
        with open(rtab_path) as infile:
            header = infile.readline().strip().split("\t")
            lines = infile.readlines()
            for line in tqdm(lines):
                row = [l.strip() for l in line.split("\t")]
                variant = f"{row[0]}"
                data = {col:val for col,val in zip(header,row)}
                variants[variant] = {"data": data }

    # Render variant label text to svg paths
    variant_types = set()
    if len(variants) > 0:
        logging.info(f"Calculating variant label dimensions.")
        for variant in tqdm(variants):
            try:
                variant_type = variant.split("|")[1].split(":")[0]
            except IndexError:
                variant_type = "unknown"
            # Space out content around the pipe delim
            text = variant.replace("|", " | ").replace("presence_absence", "present")
            label = text_to_path(text, size=font_size, family=font_family)
            w = label["w"]
            if w > variant_labels_wmax:
                variant_labels_wmax = w
                variant_labels_wmax_text = text
            variants[variant]["text"] = text
            variants[variant]["type"] = variant_type
            variants[variant]["label"] = label
            variant_types.add(variant_type)
            variants[variant]["fill"] = palette[variant_type]

        logging.info(f"Grouping variants by type.")
        variants_order = OrderedDict()
        for variant_type in palette:
            for variant,variant_data in variants.items():
                if variant_type != variant_data["type"]: continue
                variants_order[variant] = variant_data

        variants = variants_order

        if variant_labels_wmax_text:
            logging.info(f"The widest variant label is {math.ceil(variant_labels_wmax)}px: {variant_labels_wmax_text}")
    # Order variant types by palette
    variant_types = [v for v in palette if v in variant_types]

    # -------------------------------------------------------------------------
    # Initialize plot data

    plot_data = OrderedDict()

    # If tree provided, orient plot around it
    if tree_path:
        for node in tree.preorder_node_iter():
            plot_data[node.label] = {
                "x": None,
                "y": None,
                "parent": node.parent_node.label if node.parent_node else None,
                "children": [c.label for c in node.child_nodes()],
                "label": node_labels[node.label] if node.label in node_labels else OrderedDict(),
                "heatmap": OrderedDict()
            }

    # Check if any variant samples are missing from the tree tips
    for node in all_samples:
        if node not in plot_data:
            plot_data[node] = {
                "x": None,
                "y": None,
                "parent": None,
                "children": [],
                "label": node_labels[node] if node in node_labels else OrderedDict(),
                "heatmap": OrderedDict()
            }

    # -------------------------------------------------------------------------
    # Node coordinates

    logging.info(f"Calculating node coordinates.")

    # Identify the branch that sticks out the most, this will
    # determine tip label placement
    node_xmax = 0
    node_i = 0
    node_xmax_text = node_labels_wmax_text

    # X Coordinates: Distance to the root
    if tree_path:
        for node,dist in zip(
            tree.preorder_node_iter(),
            tree.calc_node_root_distances(return_leaf_distances_only=False)
        ):
            plot_data[node.label]["x"] = dist
            if  dist > node_xmax:
                node_xmax = dist
                node_xmax_text = node.label

        logging.info(f"The most distant tree node is: {node_xmax_text}")

        # Y Coordinates: Start with tips
        for node in list(tree.leaf_nodes()):
            plot_data[node.label]["y"] = node_i
            node_i += 1

        # Y Coordinates: Place internal nodes at the midpoint of (immediate?) children
        for node in tree.postorder_internal_node_iter():
            children = [c.label for c in node.child_nodes()]
            children_y = [plot_data[c]["y"] for c in children]
            y = sum(children_y) / len(children)
            plot_data[node.label]["y"] = y
    
    # Set coords for non-tree samples
    for node in plot_data:
        if plot_data[node]["x"] == None:
            plot_data[node]["x"] = 0
            plot_data[node]["y"] = node_i
            node_i += 1

    # -------------------------------------------------------------------------
    # Coordinate scaling

    # Rescale coordinates to pixels
    # X Scaling:  Based on the user's requested tree width
    # Y Scaling:  Based on the number of nodes in the tree

    logging.info(f"Rescaling branch lengths to pixels.")
    node_xmax = max([c["x"] for c in plot_data.values()])
    node_ymax = max([c["y"] for c in plot_data.values()])

    # Heatmap dimensions for scaling
    heatmap_s   = (node_labels_hmax * heatmap_scale)
    heatmap_pad = heatmap_s / 4
    heatmap_r   = 2
    heatmap_h   = (heatmap_s * len(all_samples)) + (heatmap_pad * (len(all_samples) - 1))
    logging.info(f"Heatmap boxes will be {heatmap_s}px wide with {heatmap_pad}px of padding.")

    # Rescale the x coordinates based on pre-defined maximum width
    if node_xmax > 0:
        x_scale = math.floor(tree_width / node_xmax)
    else:
        x_scale = 1
    # Rescale y coordinates based on the heatmap and text label dimensions
    if len(all_samples) > 1:
        tree_h = math.ceil((heatmap_s * len(all_samples)) + (heatmap_pad * (len(all_samples) - 1)) - heatmap_s)
        y_scale = math.floor(tree_h / node_ymax)
    else:
        tree_h = math.ceil((heatmap_s * len(all_samples)))
        y_scale = math.floor(tree_h)

    # Rescale the node coordinates, and position them
    logging.info(f"Positioning nodes.")
    for node, data in plot_data.items():
        # Locate the node's x coordinate, based on the distance to root
        plot_data[node]["x"] = math.floor(margin + root_branch + (data["x"] * x_scale))
        # Locate the node's y coordinate, if heatmap is available
        if variant_labels_wmax > 0:
            y = math.floor(margin + variant_labels_wmax + heatmap_s + (data["y"] * y_scale))
        else:
            y = math.floor(margin + data["y"] * y_scale)
        plot_data[node]["y"] = y

    # Get the new X coord where tips start
    node_xmax = plot_data[node_xmax_text]["x"]
    node_labels_x = node_xmax + tip_pad

    # Rescale the node label coordinates
    logging.info(f"Positioning tip labels.")
    for node in all_samples:
        data = plot_data[node]
        # Adjust the label y position to center on the first glyph
        first_glyph = data["label"]["glyphs"][0]
        first_glyph_h = first_glyph["h"]
        label_y = data["y"] + (first_glyph_h / 2)
        plot_data[node]["label"]["y"] = label_y
        label_w = plot_data[node]["label"]["w"]
        for i,glyph in enumerate(data["label"]["glyphs"]):
            # Left alignment
            #plot_data[node]["label"]["glyphs"][i]["x"] = node_labels_x + glyph["x"]
            # Right alignment
            plot_data[node]["label"]["glyphs"][i]["x"] = node_labels_x + (node_labels_wmax - label_w) + glyph["x"]
            plot_data[node]["label"]["glyphs"][i]["y"] = label_y

    # Heatmap labels
    logging.info(f"Positioning variant labels.")
    heatmap_x = node_labels_x + node_labels_wmax + tip_pad
    heatmap_w = 0

    # Add spaces between the different variant types
    prev_var_type = None
    prev_x = None

    for variant in variants:
        # Adjust the position to center on the first glyph
        data               = variants[variant]["data"]
        variant_type       = variants[variant]["type"]
        glyphs             = variants[variant]["label"]["glyphs"]
        first_glyph        = glyphs[0]
        first_glyph_h      = first_glyph["h"]
        first_glyph_offset = first_glyph_h + ((heatmap_s - first_glyph_h)/ 2)

        # Set the absolute position of the label
        if prev_x == None:
            x = heatmap_x + first_glyph_offset
        elif prev_var_type != None and variant_type != prev_var_type:
            x += heatmap_s + (heatmap_pad * 3)
        else:
            x = prev_x + (heatmap_s + heatmap_pad)

        # Set the absolute position of the box
        variants[variant]["box"] = {"x": x - first_glyph_offset}
        prev_x = x
        prev_var_type = variant_type
        # Add an extra 2 pixels to make it reach just beyond the final box
        heatmap_w = (2 + heatmap_s + x - first_glyph_offset) -  heatmap_x

        y = margin + variant_labels_wmax
        variants[variant]["label"]["x"] = x
        variants[variant]["label"]["y"] = y
        for i_g,g in enumerate(glyphs):
            variants[variant]["label"]["glyphs"][i_g]["x"] = x + g["x"]
            variants[variant]["label"]["glyphs"][i_g]["y"] = y

    # -------------------------------------------------------------------------
    # Heatmap Table Data

    score_min_alpha, score_max_alpha = 0.10, 1.0

    logging.info(f"Positioning heatmap boxes.")
    for i,variant in enumerate(variants):
        variant_type = variants[variant]["type"]
        data = variants[variant]["data"]
        label = variants[variant]["label"]
        x = variants[variant]["box"]["x"]
        for node in all_samples:
            # Ex. a sample in tree but not heatmap
            if node not in data: continue
            v = data[node]
            v = int(v) if v.isdigit() else str(v)
            fill_opacity = 1.0
            # Light grey for missing values (".")
            score = data["score"] if "score" in data else None
            if v == 1:
                box_fill = variants[variant]["fill"]
                # Use score range if possible
                if score_min != None and "score" in data:
                    score = float(data["score"])
                    # Option 1: Fixed range from 0 to 1
                    fill_opacity = linear_scale(score, 0, 1.0, score_min_alpha, score_max_alpha)
                    # Option 2. Relative range from score_min to score_max
                    # fill_opacity = linear_scale(score, score_min, score_max, score_min_alpha, score_max_alpha)
                box_stroke = "black"
            elif v == 0:
                box_fill = "white"
                box_stroke = "black"
            else:
                box_fill = "none"
                box_stroke = "none"
            y = plot_data[node]["y"] - (heatmap_s / 2)
            d = {
                "v": v, "x": x, "y": y, "w": heatmap_s, "h": heatmap_s,
                "r": heatmap_r, "fill": box_fill, "fill_opacity": fill_opacity, "stroke": box_stroke,
                "hovertext": data,
            }
            plot_data[node]["heatmap"][variant] = d

    # -----------------------------------------------------------------------------
    # Legend

    legend_x = (heatmap_x + heatmap_w + heatmap_s)
    legend_y = margin + variant_labels_wmax + tip_pad
    legend_w = 0
    legend_h = 0
    legend = OrderedDict()
    legend_labels_hmax = 0
    legend_labels_wmax = 0
    legend_title_wmax = 0
    legend_title_hmax = 0 

    if score_min != None and score_max != None:    
        logging.info(f"Drawing legend at: {legend_x}, {legend_y}")
        legend_h = (heatmap_s * 3) + (heatmap_pad * 2)
        # TBD: Think about negative score handling?
        # Option 1. Fixed Values
        legend_values = ["0.00", "0.25", "0.50", "0.75", "1.00"]
        # # Option 2. Relative to score max
        # interval = score_max / 4
        # legend_values = [round(i * interval,2) for i in range(0,4)] + [round(score_max, 2)]
        # legend_values = ["{:.2f}".format(v) for v in legend_values]

        # Tick labels: convert to path
        tick_labels = [text_to_path(text, size=font_size * 0.50, family=font_family) for text in legend_values]
        for label in tick_labels:
            legend_labels_hmax = max(legend_labels_hmax, label["h"])
            legend_labels_wmax = max(legend_labels_wmax, label["w"])
        title_label = text_to_path("Score", size=font_size * 0.75, family=font_family)
        legend_title_hmax = max(legend_title_hmax, label["h"])
        legend_title_wmax = max(legend_title_wmax, label["w"])
        

        # Width of legend + tick len + pad + label wmax + pad
        entry_w = heatmap_s + heatmap_pad + (heatmap_pad / 2) + legend_labels_wmax + heatmap_s
        if (legend_title_wmax + heatmap_s) > entry_w:
            entry_w = legend_title_wmax + heatmap_s

        for i_v, variant_type in enumerate(variant_types):
            box = {
                    "x": legend_x + (i_v * entry_w),
                    "y": legend_y,
                    "w": heatmap_s,
                    "h": (heatmap_s * 3) + (heatmap_pad * 2),
                    "fill": f"url(#{variant_type}_gradient)",
                }
            title = copy.deepcopy(title_label)
            title["x"] = box["x"]
            title["y"] = box["y"] - tip_pad
            for i_g,g in enumerate(title["glyphs"]):
                title["glyphs"][i_g]["x"] = title["x"] + g["x"]
                title["glyphs"][i_g]["y"] = title["y"]
            ticks = []
            for i_t,text in enumerate(reversed(legend_values)):
                x1 = box["x"] + box["w"]
                x2 = x1 + heatmap_pad
                y = legend_y + (i_t * (legend_h / (len(legend_values) - 1)))
                label = copy.deepcopy(list(reversed(tick_labels))[i_t])
                label["x"] = x2 + (heatmap_pad / 2)
                label["y"] = y

                # Adjust the label y position to center on first numeric char
                center_glyph_h = label["glyphs"][0]["h"]
                label["y"] = y + (center_glyph_h / 2)
                for i_g,g in enumerate(label["glyphs"]):
                    label["glyphs"][i_g]["x"] = label["x"] + g["x"]
                    label["glyphs"][i_g]["y"] = label["y"]

                tick = {"text":  text, "x1": x1, "x2":x2, "y": y, "label": copy.deepcopy(label) }
                ticks.append(tick)

            legend[variant_type] = {"title": title, "box": box, "ticks": ticks}
            legend_w += entry_w

    # -----------------------------------------------------------------------------
    # Draw Tree

    logging.info(f"Calculating final image dimensions.")
    # Final image dimensions
    if len(legend) > 0:
        logging.info("Final dimensions with legend.")
        width = math.ceil(legend_x + legend_w + margin)
        height = math.ceil(margin + variant_labels_wmax + tip_pad + tree_h + heatmap_s + margin)
    elif heatmap_w > 0 and heatmap_h > 0:
        logging.info("Final dimensions with heatmap.")
        width = math.ceil(heatmap_x + heatmap_w + margin)
        height = math.ceil(margin + variant_labels_wmax + tip_pad + tree_h + heatmap_s + margin)
    else:
        logging.info("Final dimensions without heatmap.")
        width = math.ceil(node_labels_x + node_labels_wmax + margin)
        height = math.ceil(tree_h + (2*margin))

    svg_path = os.path.join(output_dir, f"{prefix}plot.svg")
    logging.info(f"Rendering output svg ({width}x{height}): {svg_path}")
    with open(svg_path, 'w') as outfile:
        header = textwrap.dedent(
        f"""\
        <svg 
            version="1.1"
            xmlns="http://www.w3.org/2000/svg"
            xmlns:xlink="http://www.w3.org/1999/xlink"
            preserveAspectRatio="xMidYMid meet"
            width="{width}"
            height="{height}"
            viewbox="0 0 {width} {height}">
        """)
    
        outfile.write(header.strip() + "\n")

        branches = []
        tip_circles = []
        tip_dashes = []
        tip_labels = []
        variant_labels = OrderedDict()
        heatmap_boxes = OrderedDict()
        focal_boxes = OrderedDict()
        legend_entries = OrderedDict()

        # Variant Heatmap Labels
        for variant in variants:
            label = variants[variant]["label"]
            rx, ry = label["x"], label["y"]
            if variant not in variant_labels:
                variant_labels[variant] = []
                heatmap_boxes[variant]  = []
            for g in label["glyphs"]:
                line = f"<path transform='rotate(-90, {rx}, {ry}) translate({g['x']},{g['y']})' d='{g['d']}' />"
                variant_labels[variant].append(line)

        for node,data in plot_data.items():
            logging.debug(f"node: {node}, data: {data}")

            # Root branch, add custom little starter branch
            if not data["parent"]:
                cx, cy = data['x'], data['y']
                px, py = data['x'] - root_branch, cy
            # Non-root
            else:
                p_data = plot_data[data["parent"]]
                cx, cy, px, py = data['x'], data['y'], p_data['x'], p_data['y']

            # Draw tree lines
            if node in tree_nodes:
                line = f"<path d='M {cx} {cy} H {px} V {py}' style='stroke:black;stroke-width:2;fill:none'/>"
                branches.append(line)

            # Focal box
            if node in focal:
                # Draw the box to the midway point between it and the next label
                if len(plot_data) == 1:
                    rh = node_labels_hmax
                else:
                    # Figure out adjacent tip
                    tip_i = all_samples.index(node)
                    adjacent_tip = all_samples[tip_i+1] if tip_i == 0 else all_samples[tip_i-1]
                    rh = abs(plot_data[adjacent_tip]["y"] - data['y'])
                ry = data['y'] - (rh /2)
                rx = node_labels_x
                rw = node_labels_wmax + tip_pad + heatmap_w
                
                rect = f"<rect width='{rw}' height='{rh}' x='{rx}' y='{ry}' rx='1' ry='1' style='fill:grey;fill-opacity:.40'/>"
                focal_boxes[node] = [rect]

            # Dashed line: tree -> tip
            if node in tree_tips:
                lx = node_labels_x + (node_labels_wmax - data["label"]["w"])
                line = f"<path d='M {cx + 4} {cy} H {lx - 4}' style='stroke:grey;stroke-width:1;fill:none' stroke-dasharray='4 4'/>"
                tip_dashes.append(line)

            # Tip circles
            if node in focal:
                circle = f"<circle cx='{cx}' cy='{cy}' r='4' style='fill:black;stroke:black;stroke-width:1' />"
                tip_circles.append(circle)

            # Sample Labels
            if node in all_samples:
                for g in data["label"]["glyphs"]:
                    line = f"<path transform='translate({g['x']},{g['y']})' d='{g['d']}' />"
                    tip_labels.append(line)
            # Heatmap
            if node in all_samples:
                for variant,d in data["heatmap"].items():

                    hover_text = "\n".join([f"{k}: {v}" for k,v in d["hovertext"].items() if k not in all_samples and not k.endswith("_alt")])
                    rect = f"<rect width='{d['w']}' height='{d['h']}' x='{d['x']}' y='{d['y']}' rx='{d['r']}' ry='{d['r']}' style='fill:{d['fill']};fill-opacity:{d['fill_opacity']};stroke-width:1;stroke:{d['stroke']}'>"
                    title = f"<title>{hover_text}</title>"
                    heatmap_boxes[variant] += [rect, title, "</rect>"]

            # Legend
            for variant_type in legend:
                legend_entries[variant_type] = []
                # Gradient
                fill = palette[variant_type]
                gradient = f"""
                <linearGradient id="{variant_type}_gradient" x1="0" x2="0" y1="0" y2="1" gradientTransform="rotate(180 0.5 0.5)">
                    <stop offset="0%" stop-color="{fill}" stop-opacity="{score_min_alpha}" />
                    <stop offset="100%" stop-color="{fill}" stop-opacity="{score_max_alpha}" />
                </linearGradient>"""
                legend_entries[variant_type].append(gradient)
                # Legend Title
                title = legend[variant_type]["title"]
                for g in title["glyphs"]:
                    label = f"<path transform='translate({g['x']},{g['y']})' d='{g['d']}' />"
                    legend_entries[variant_type].append(label)
                # Legend box
                d = legend[variant_type]["box"]
                rect = f"<rect width='{d['w']}' height='{d['h']}' x='{d['x']}' y='{d['y']}' style='fill:{d['fill']};stroke-width:1;stroke:black'/>"
                legend_entries[variant_type].append(rect)
                # Legend Ticks
                for t in legend[variant_type]["ticks"]:
                    # Tick line
                    line = f"<path d='M {t['x1']} {t['y']} H {t['x2']} V {t['y']}' style='stroke:black;stroke-width:1;fill:none'/>"
                    legend_entries[variant_type].append(line)
                    # Tick Label
                    for g in t["label"]["glyphs"]:
                        label = f"<path transform='translate({g['x']},{g['y']})' d='{g['d']}' />"
                        legend_entries[variant_type].append(label)

        # Draw elements in groups
        indent = "    "
        outfile.write(f"{indent * 1}<g id='Plot'>" + "\n")

        # White canvas background
        background = f"<rect width='{width}' height='{height}' x='0' y='0' style='fill:white;stroke-width:1;stroke:white'/>"
        outfile.write(f"{indent * 2}<g id='Background'>\n{indent * 3}{background}\n{indent * 2}</g>" + "\n")

        # Focal boxes start
        outfile.write(f"{indent * 2}<g id='Focal'>" + "\n")
        # Focal Boxes
        outfile.write(f"{indent * 3}<g id='Focal Boxes'>\n")
        for sample,boxes in focal_boxes.items():
            outfile.write(f"{indent * 4}<g id='{sample}'>\n")
            outfile.write("\n".join([f"{indent * 5}{b}" for b in boxes]) + "\n")
            outfile.write(f"{indent * 4}</g>\n")
        outfile.write(f"{indent * 3}</g>\n")
        # Focal boxes end
        outfile.write(f"{indent * 2}</g>\n")
    
        # Tree Start
        outfile.write(f"{indent * 2}<g id='Tree'>" + "\n")
        # Branches
        outfile.write(f"{indent * 3}<g id='Branches'>\n")
        outfile.write("\n".join([f"{indent * 4}{b}" for b in branches]) + "\n")
        outfile.write(f"{indent * 3}</g>\n")
        # Tip Circles
        outfile.write(f"{indent * 3}<g id='Tip Circles'>\n")
        outfile.write("\n".join([f"{indent * 4}{c}" for c in tip_circles]) + "\n")
        outfile.write(f"{indent * 3}</g>\n")
        # Tip Dashes
        outfile.write(f"{indent * 3}<g id='Tip Dashes'>\n")
        outfile.write("\n".join([f"{indent * 4}{d}" for d in tip_dashes]) + "\n")
        outfile.write(f"{indent * 3}</g>\n")
        # Tip Labels
        outfile.write(f"{indent * 3}<g id='Tip Labels'>\n")
        outfile.write("\n".join([f"{indent * 4}{t}" for t in tip_labels]) + "\n")
        outfile.write(f"{indent * 3}</g>\n")     
        # Tree End
        outfile.write(f"{indent}</g>\n")

        # Heatmap Start
        outfile.write(f"{indent}<g id='Heatmap'>" + "\n")
        # Variant Labels
        outfile.write(f"{indent * 3}<g id='Variant Labels'>\n")
        for variant,paths in variant_labels.items():
            outfile.write(f"{indent * 4}<g id='{variant}'>\n")
            outfile.write("\n".join([f"{indent * 5}{p}" for p in paths]) + "\n")
            outfile.write(f"{indent * 4}</g>\n")
        outfile.write(f"{indent * 3}</g>\n")
        # Heatmap Boxes
        outfile.write(f"{indent * 3}<g id='Heatmap Boxes'>\n")
        for variant,boxes in heatmap_boxes.items():
            outfile.write(f"{indent * 4}<g id='{variant}'>\n")
            outfile.write("\n".join([f"{indent * 5}{b}" for b in boxes]) + "\n")
            outfile.write(f"{indent * 4}</g>\n")
        outfile.write(f"{indent * 3}</g>\n")
        # Heatmap End
        outfile.write(f"{indent * 2}</g>\n")

        # Legend
        if score_min != None and score_max != None:
            outfile.write(f"{indent * 2}<g id='Legend'>\n")
            for variant_type in legend_entries:
                outfile.write(f"{indent * 3}<g id='{variant_type}'>\n")
                for element in legend_entries[variant_type]:
                   outfile.write(f"{indent * 4}{element}\n")
                outfile.write(f"{indent * 3}</g>\n")
            outfile.write(f"{indent * 4}</g>\n")

                # outfile.write(f"{indent * 3}<g id='{variant_type}'>")
                # fill = palette[variant_type]
                # gradient = f"""
                # <linearGradient id="{variant_type}_gradient" x1="0" x2="0" y1="0" y2="1" gradientTransform="rotate(180 0.5 0.5)">
                #     <stop offset="0%" stop-color="{fill}" stop-opacity="{score_min_alpha}" />
                #     <stop offset="100%" stop-color="{fill}" stop-opacity="{score_max_alpha}" />
                # </linearGradient>"""
                # outfile.write(f"{indent * 4}{gradient}\n")
                # outfile.write(f"{indent * 3}</g>\n")
                # # Box
                # rect = f"<rect width='{legend_w}' height='{legend_h}' x='{legend_x}' y='{legend_y}' style='stroke:black;stroke-width:1;fill:url(#LegendGradient)'/>"
                # outfile.write(f"{indent * 3}{rect}\n")
            # # Ticks
            # for variant_type in legend_ticks:
            #     for t in legend_ticks[variant_type]:
            #         line = f"<path d='M {t['x1']} {t['y']} H {t['x2']} V {t['y']}' style='stroke:black;stroke-width:1;fill:none'/>"
            #         outfile.write(f"{indent * 3}{line}\n")
            #         # Label
            #         for g in t["label"]["glyphs"]:
            #             glyph = f"<path transform='translate({g['x']},{g['y']})' d='{g['d']}' />"
            #             outfile.write(f"{indent * 3}{glyph}\n")

        # Close out the plot group
        outfile.write(f"{indent * 1}</g>\n")
        outfile.write("</svg>")

    png_path = os.path.join(output_dir, f"{prefix}plot.png")
    logging.info(f"Rendering output png ({width}x{height}): {png_path}")
    cairosvg.svg2png(url=svg_path, write_to=png_path, output_width=width, output_height=height, scale=png_scale)

    return svg_path

def cli(args:str=None):

    # Parse input args from system or function input
    sys_argv_original = sys.argv
    if args != None:
        sys.argv = f"pangwas {args}".split(" ")
    command = " ".join(sys.argv)
    if len(sys.argv) == 1:
        sys.argv = ["pangwas", "--help"]

    # Handle sys exits by argparse more gracefully
    try:
        options = get_options()
    except SystemExit as exit_status:
        if f"{exit_status}" == "0":
            return 0
        elif f"{exit_status}" == "2":
            return exit_status
        else:
            msg = f"pangwas cli exited with status: {exit_status}"
            logging.error(msg)
            return exit_status

    # Restore original argv
    sys.argv = sys_argv_original

    if options.version:
        import importlib.metadata
        version = importlib.metadata.version("pangwas")
        print(f"pangwas v{version}")
        return 0

    logging.info("Begin")
    logging.info(f"Command: {command}")

    # Organize options as kwargs for functions
    kwargs = {k:v for k,v in vars(options).items() if k not in ["version", "subcommand"]}
    try:
        fn = globals()[options.subcommand]
        fn(**kwargs)
    except KeyError:
        logging.error(f"A pangwas subcommand is required (ex. pangwas extract).")
        return 1

    logging.info("Done")

if __name__ == "__main__":
    cli()
