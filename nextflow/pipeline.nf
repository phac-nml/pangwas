#!/usr/bin/env nextflow

// Specify the Domain-Specific Language (DSL) version
nextflow.enable.dsl = 2

// Some nf-core helper functions for start up and completion
include { PIPELINE_INITIALISATION } from "${projectDir}/subworkflows/local/utils_nfcore_pangwas_pipeline"
include { PIPELINE_COMPLETION     } from "${projectDir}/subworkflows/local/utils_nfcore_pangwas_pipeline"

// Custom processes
include { 
    ANNOTATE;
    EXTRACT; COLLECT; CLUSTER; DEFRAG; SUMMARIZE; 
    ALIGN;
    STRUCTURAL; SNPS; PRESENCE_ABSENCE; COMBINE;
    TREE; GWAS;
    MANHATTAN; HEATMAP;
} from "${projectDir}/modules/local/pangwas"

// Use Nextflow's built in hashing function
import nextflow.util.CacheHelper

// Main workflow
workflow {

    // Convert inputs to absolute paths
    def outdir = params.outdir
    if (outdir && !outdir.startsWith("/")){
       outdir = "${launchDir}/${outdir}"
    }
    def input = params.input
    if (input && !input.startsWith("/")){
       input = "${launchDir}/${input}"
    }
    def variants = params.variants
    if (variants && !variants.startsWith("/")){
       variants = "${launchDir}/${variants}"
    }
    def tree = params.tree
    if (tree && !tree.startsWith("/")){
       tree = "${launchDir}/${tree}"
    }
    def clusters = params.clusters
    if (clusters && !clusters.startsWith("/")){
       clusters = "${launchDir}/${clusters}"
    }
    def bed = params.bed
    if (bed && !bed.startsWith("/")){
       bed = "${launchDir}/${bed}"
    }

    // nf-core: Run initialisation tasks
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        outdir,
        input
    )
    
    // Set up the samplesheet
    ch_samplesheet = PIPELINE_INITIALISATION.out.samplesheet
    ch_samplesheet_file = Channel.of(input)
    ch_versions = Channel.empty()

    // Split samplesheet into separate channels for assembly and gff
    ch_assembly = ch_samplesheet.map   { meta -> [meta, meta.assembly] }
    ch_gff      = ch_samplesheet.filter{ meta -> meta.gff }.map{ meta -> [meta, meta.gff] }
    ch_empty    = Channel.empty()

    // Depending on the optional inputs, we can skip some steps
    skip_annotate         = (variants && tree) ? true : params.skip_annotate
    skip_extract          = ((variants && tree) || params.only_annotate) ? true : params.skip_extract
    skip_cluster          = (clusters) ? true : params.skip_cluster
    skip_presence_absence = (variants) ? true : params.skip_presence_absence
    skip_combine          = (variants) ? true : params.skip_combine
    skip_gwas             = (!params.trait) ? true: params.skip_gwas

    // ------------------------------------------------------------------------
    // Step 1. Annotate
    // ------------------------------------------------------------------------

    if (skip_annotate || !params.bakta_db) {
        ANNOTATE (ch_empty, ch_empty)
    }
    else {
        // Filter the assemblies for ones that don't have annotations
        ch_assembly_input = ch_samplesheet
            .filter { meta -> !meta.gff }
            .map{ meta -> [meta, meta.assembly] } 
        ANNOTATE (ch_assembly_input, file(params.bakta_db))
    }
    ch_versions = ANNOTATE.out.versions.mix(ch_versions)
    ch_gff      = ANNOTATE.out.gff.mix(ch_gff)

    // ------------------------------------------------------------------------
    // Stage 2: Cluster
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // Extract

    if (skip_extract) { EXTRACT ( ch_empty ) }
    else              { EXTRACT ( ch_gff )   }
    ch_versions = EXTRACT.out.versions.mix(ch_versions)

    // Collect file paths to table, hashing content so we know when to rerun
    ch_extract = EXTRACT.out.tsv
            .toSortedList( { a, b -> a[0].id <=> b[0].id } ).flatMap()
            .map { meta, tsv -> 
                def hash = CacheHelper.hasher([ tsv.text ]).hash().toString()
                "${tsv}\t${hash}\n"
            }
            .collectFile(name: "input.tsv", storeDir: "${outdir}/collect", sort: "index", cache: true)
            .map { tsv -> file(tsv) }   

    // ------------------------------------------------------------------------
    // Collect

    if (params.skip_collect) { COLLECT ( ch_empty ) }
    else                     { COLLECT ( ch_extract ) }
    ch_versions = COLLECT.out.versions.mix(ch_versions)
    ch_collect_fasta = COLLECT.out.fasta
    ch_collect_tsv = COLLECT.out.tsv

    // ------------------------------------------------------------------------
    // Cluster

    if (skip_cluster) { CLUSTER ( ch_empty ) }
    else              { CLUSTER ( ch_collect_fasta ) }
    ch_versions       = CLUSTER.out.versions.mix(ch_versions)    
    ch_clusters       = CLUSTER.out.clusters
    ch_representative = CLUSTER.out.representative

    // ------------------------------------------------------------------------
    // Defrag

    if (params.skip_defrag){
        DEFRAG (ch_empty, ch_empty )
        ch_clusters_final = ch_clusters
    } else { 
        DEFRAG ( ch_clusters, ch_representative )
        ch_clusters_final = DEFRAG.out.clusters
    }
    ch_versions  = DEFRAG.out.versions.mix(ch_versions)

    // ------------------------------------------------------------------------
    // Summarize

    if (params.skip_summarize) { SUMMARIZE (ch_emptyl, ch_empty) } 
    else                       { SUMMARIZE (ch_clusters_final, ch_collect_tsv) }
    ch_versions  = SUMMARIZE.out.versions.mix(ch_versions)
    ch_summarize = SUMMARIZE.out.clusters

    if (clusters) {
        ch_summarize = Channel.of(clusters)
    }

    // ------------------------------------------------------------------------
    // Stage 3: Align
    // ------------------------------------------------------------------------

    if (params.skip_align) { ALIGN ( ch_empty, ch_empty ) } 
    else                   { ALIGN ( ch_summarize, ch_collect_tsv ) }
    ch_versions   = ALIGN.out.versions.mix(ch_versions) 
    ch_alignments = ALIGN.out.cluster_alignments
    ch_alignment  = ALIGN.out.pangenome_alignment
    ch_consensus  = ALIGN.out.pangenome_consensus
    ch_bed        = ALIGN.out.pangenome_bed

    if (bed) {
        ch_bed = Channel.of(bed)
    }

    // ------------------------------------------------------------------------
    // Stage 4: Variants
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // Structural

    if (params.skip_structural) { STRUCTURAL ( ch_empty, ch_empty ) }
    else                        { STRUCTURAL ( ch_summarize, ch_alignments ) }
    ch_versions   = STRUCTURAL.out.versions.mix(ch_versions)
    ch_structural = STRUCTURAL.out.rtab.ifEmpty([])
    ch_variants   = ch_structural

    // ------------------------------------------------------------------------
    // SNPs

    if ( params.skip_snps ) { SNPS (ch_empty, ch_empty, ch_empty, ch_empty ) }
    else                    { SNPS (ch_alignment, ch_bed, ch_consensus, ch_structural) }
    ch_versions       = SNPS.out.versions.mix(ch_versions)
    ch_snps_rtab      = SNPS.out.rtab.ifEmpty([])
    ch_snps_alignment = SNPS.out.alignment_core
    ch_constant_sites = SNPS.out.constant_sites
    ch_variants       = ch_variants.concat(ch_snps_rtab)

    // ------------------------------------------------------------------------
    // Presence Absence

    if (skip_presence_absence) { PRESENCE_ABSENCE ( ch_empty ) }
    else                       { PRESENCE_ABSENCE ( ch_summarize ) }
    ch_versions         = PRESENCE_ABSENCE.out.versions.mix(ch_versions)
    ch_presence_absence = PRESENCE_ABSENCE.out.rtab.ifEmpty([])
    ch_variants         = ch_variants.concat(ch_presence_absence)

    // ------------------------------------------------------------------------
    // Combine

    if (skip_combine) { COMBINE (ch_empty) }
    else              { COMBINE (ch_variants.collect()) }
    ch_versions = COMBINE.out.versions.mix(ch_versions)
    ch_variants = COMBINE.out.rtab

    if (variants){
        ch_variants = Channel.of(variants)
    }

    // ------------------------------------------------------------------------
    // Stage 5: Tree

    if (tree){
        TREE ( ch_empty, ch_empty )
        ch_tree = Channel.of(tree)
    }
    else if (params.skip_tree) {
        TREE ( ch_empty, ch_empty )
        ch_tree = ch_empty
    }
    else {
        TREE ( ch_snps_alignment, ch_constant_sites )
        ch_tree = TREE.out.tree_rooted
    } 
    ch_versions = TREE.out.versions.mix(ch_versions)

    // ------------------------------------------------------------------------
    // Stage 6: GWAS

    if(skip_gwas) { GWAS ( ch_empty, ch_empty, ch_empty, ch_empty, ch_empty) }
    else          { GWAS ( params.trait, ch_samplesheet_file, ch_variants, ch_tree, ch_summarize.ifEmpty([]) ) }
    ch_gwas_all         = GWAS.out.locus_effects_all.flatMap()
    ch_gwas_significant = GWAS.out.locus_effects_significant.flatMap()
    ch_gwas_filtered    = GWAS.out.locus_effects_filtered.flatMap()
    ch_focal            = GWAS.out.focal.flatMap()
    ch_tree             = GWAS.out.tree.flatMap()
    ch_versions         = GWAS.out.versions.mix(ch_versions)


    // ------------------------------------------------------------------------
    // Stage 7: Plot

    if(params.skip_manhattan) { MANHATTAN (ch_empty, ch_empty ) }
    else                      { MANHATTAN (ch_gwas_all, ch_bed ) }
    ch_manhattan = MANHATTAN.out.svg
    ch_versions  = MANHATTAN.out.versions.mix(ch_versions)

    if      (params.variants_to_plot ==~ /.*filtered.*/)    { ch_gwas_to_plot = ch_gwas_filtered }
    else if (params.variants_to_plot ==~ /.*significant.*/) { ch_gwas_to_plot = ch_gwas_significant }
    else if (params.variants_to_plot ==~ /.*all.*/)         { ch_gwas_to_plot = ch_gwas_all }
    else                                                    { ch_gwas_to_plot = ch_empty }

    // TBD: Can we be certain that ch_gwas_to_plot and ch_focal will stay in the same order?
    if(params.skip_heatmap) { HEATMAP (ch_empty, ch_empty, ch_empty ) }
    else                    { HEATMAP (ch_gwas_to_plot, ch_tree.first(), ch_focal ) }
    ch_heatmap  = HEATMAP.out.svg
    ch_versions = HEATMAP.out.versions.mix(ch_versions)

    // ------------------------------------------------------------------------
    // COMPLETE
    // ------------------------------------------------------------------------

    // nf-core: Run completion tasks, remove all email/messaging parts
    PIPELINE_COMPLETION (
        null,
        null,
        null,
        outdir,
        params.monochrome_logs,
        null,
        Channel.empty()
    )
}
