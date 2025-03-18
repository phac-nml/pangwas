

process ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: "${params.outdir}/annotate/${prefix}", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"
    containerOptions "--writable-tmpfs"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*.gff3")        , emit: gff
    tuple val(meta), path("*.log")         , emit: log
    tuple val(meta), path("*.versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus == 1 ? 1 : task.cpus - 1
    def tmp_dir = workflow.containerEngine == 'singularity' ? '/bakta/tmp' : '/tmp'
    """
    pangwas annotate \\
      $args \\
      --fasta ${fasta} \\
      --db ${db} \\
      --sample ${prefix} \\
      --prefix ${prefix} \\
      2>&1 | tee ${prefix}.log

    major=\$(grep -m 1 major ${db}/version.json | cut -d " " -f 4 | sed -E 's/"|,//g')
    minor=\$(grep -m 1 minor ${db}/version.json | cut -d " " -f 4 | sed -E 's/"|,//g')
    type=\$(grep -m 1 type ${db}/version.json | cut -d " " -f 4 | sed -E 's/"|,//g')

    cat <<-END_VERSIONS > ${prefix}.versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
        bakta_db: \${major}.\${minor}-\${type}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    touch ${prefix}.log

    cat <<-END_VERSIONS > ${prefix}.versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}

process EXTRACT {
    label "process_single"

    publishDir path: "${params.outdir}/extract/${meta.id}", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.tsv"),    emit: tsv
    tuple val(meta), path("*.log"),    emit: log
    tuple val(meta), path("*.versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ""
    """
    pangwas extract $args --gff $gff --sample ${prefix} --prefix ${prefix} 2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > ${prefix}.versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ""
    """
    touch ${prefix}.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > ${prefix}.versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS    
    """
}

process COLLECT {
    label "process_single"

    publishDir path: "${params.outdir}/collect", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    val(tsv_paths) // TXT file of tsv paths from concatenated EXTRACT TSV

    output:
    path("*.tsv"),    emit: tsv
    path("*.fasta"),  emit: fasta
    path("*.log"),    emit: log
    path("*.versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    pangwas collect $args --tsv-paths $tsv_paths 2>&1 | tee collect.log

    cat <<-END_VERSIONS > collect.versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS    
    """

    stub:
    def args = task.ext.args ?: ""
    """
    touch collect.tsv
    touch collect.fasta
    touch collect.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS    
    """
}

process CLUSTER {
    label "process_high"

    publishDir path: "${params.outdir}/cluster", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(sequences)

    output:
    path("*clusters.tsv"),             emit: clusters
    path("*representative.fasta"),     emit: representative
    path("*versions.yml"),             emit: versions
    path("*.log"),                     emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    pangwas cluster \\
      $args \\
      --fasta $sequences \\
      --threads ${task.cpus} \\
      --memory '${task.memory}' \\
      2>&1 | tee cluster.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        mmseqs2: \$(mmseqs version)
    END_VERSIONS        
    """

    stub:
    """
    touch clusters.tsv
    touch representative.fasta
    touch cluster.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """
}

process DEFRAG {
    label "process_high"

    publishDir path: "${params.outdir}/defrag", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(clusters)       // Clusters in TSV format from CLUSTER
    path(representative) // Representative in FASTA format from CLUSTER

    output:
    path("*clusters.tsv"),         emit: clusters
    path("*align.tsv"),            emit: align
    path("*representative.fasta"), emit: representative
    path("*.log"),                 emit: log
    path("*versions.yml"),         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "defrag"
    def args = task.ext.args ?: ""
    """
    pangwas defrag \\
      $args \\
      --clusters $clusters \\
      --representative $representative \\
      --prefix ${prefix} \\
      --threads ${task.cpus} \\
      --memory '${task.memory}' \\
      2>&1 | tee defrag.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        mmseqs: \$(mmseqs version)
        ppanggolin: 2.2.0 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "defrag"
    """
    touch ${prefix}.clusters.tsv
    touch ${prefix}.representative.fasta
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        mmseqs: \$(mmseqs version)
        ppanggolin: 2.2.0
    END_VERSIONS
    """
}

process SUMMARIZE {
    label "process_single"

    publishDir path: "${params.outdir}/summarize", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(clusters)   // Clusters in TSV format from CLUSTER or DEFRAG
    path(regions)    // Sequence regions in TSV format from COLLECT

    output:
    path("*.clusters.tsv"),           emit: clusters
    path("*.phandango.csv"),          emit: phandango
    path("*.synteny.full.graphml"),   emit: graphml_full
    path("*.synteny.linear.graphml"), emit: graphml_linear
    path("*.synteny.full.gfa"),       emit: gfa_full
    path("*.synteny.linear.gfa"),     emit: gfa_linear
    path("*.log"),                    emit: log
    path("*versions.yml"),            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "summarize"
    def args = task.ext.args ?: ""
    """
    pangwas summarize \\
      $args \\
      --clusters $clusters \\
      --regions $regions \\
      --prefix ${prefix} \\
      2>&1 | tee summarize.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "summarize"
    """
    touch ${prefix}.clusters.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS    
    """
}

process ALIGN {
    label "process_high"

    publishDir path: "${params.outdir}/align", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(clusters)   // Clusters in TSV format from SUMMARIZE
    path(regions)    // Sequence regions in TSV format from COLLECT

    output:
    path("*sequences.tar.bz2")        , emit: cluster_sequences
    path("*alignments.tar.bz2")       , emit: cluster_alignments
    path("*consensus.tar.bz2")        , emit: cluster_consensus
    path("*representative.tar.bz2")   , emit: cluster_representative
    path("*pangenome.bed")            , emit: pangenome_bed
    path("*pangenome.aln")            , emit: pangenome_alignment
    path("*pangenome.consensus.fasta"), emit: pangenome_consensus
    path("*pangenome.gff3"),            emit: pangenome_gff
    path("*.log")                     , emit: log
    path("*versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    pangwas align \\
      $args \\
      --clusters $clusters \\
      --regions $regions \\
      --threads ${task.cpus} \\
      2>&1 | tee align.log

    for dir in alignments consensus sequences representative; do
      tar -cvjSf \${dir}.tar.bz2 \${dir}
      rm -rf \${dir}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        mafft: \$(mafft --version 2>&1 | cut -d ' ' -f 1 | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    """
    touch sequences.tar.xz
    touch alignments.tar.xz
    touch consensus.tar.xz
    touch representative.tar.xz
    touch pangenome.bed
    touch pangenome.representative.fasta
    touch pangenome.aln
    touch align.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        mafft: \$(mafft --version 2>&1 | cut -d ' ' -f 1 | sed 's/v//g')
    END_VERSIONS
    """
}

process STRUCTURAL {
    label "process_single"

    publishDir path: "${params.outdir}/structural", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(clusters)   // Clusters in TSV format from SUMMARIZE
    path(alignments) // Alignments directory in tar.bz2 format from ALIGN

    output:
    path("*.Rtab"), emit: rtab
    path("*.log") , emit: log
    path("*versions.yml") , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    mkdir -p alignments && tar -xvjSmf $alignments -C alignments --strip-components 1
    pangwas structural \\
      $args \\
      --clusters $clusters \\
      --alignments alignments \\
      2>&1 | tee structural.log
    rm -rf alignments

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    """
    touch structural.Rtab
    touch structural.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS    
    """
}

process SNPS {
    label 'process_single'

    publishDir path: "${params.outdir}/snps", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(alignment)  // Pangenome alignment in FASTA format from ALIGN
    path(bed)        // Pangenome bed from ALIGN
    path(consensus)  // Pangenome consensus in FASTA format from ALIGN    
    path(structural) // Optional strucutral Rtab from STRUCTURAL, set to [] to disable

    output: 
    path("*snps.all.fasta")         , emit: alignment_all
    path("*snps.core.fasta")        , emit: alignment_core
    path("*snps.all.tsv")           , emit: table_all
    path("*snps.core.tsv")          , emit: table_core
    path("*snps.all.vcf")           , emit: vcf_all
    path("*snps.core.vcf")          , emit: vcf_core
    path("*snps.Rtab")              , emit: rtab
    path("*snps.constant_sites.txt"), emit: constant_sites
    path("*.log")                   , emit: log
    path("*versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def structural_arg = structural ? "--structural $structural" : ""
    """
    pangwas snps \\
      $args $structural_arg \\
      --alignment $alignment \\
      --bed $bed \\
      --consensus $consensus \\
      2>&1 | tee snps.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    def prefix = "snps"
    """
    touch snps.all.fasta
    touch snps.core.fasta
    touch snps.all.tsv
    touch snps.core.tsv
    touch snps.Rtab
    touch snps.constant_sites.txt
    touch snps.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """
}

process PRESENCE_ABSENCE {
    label "process_single"

    publishDir path: "${params.outdir}/presence_absence", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(clusters)   // Clusters in TSV format from SUMMARIZE

    output:
    path("*.Rtab"),         emit: rtab
    path("*.log") ,         emit: log
    path("*versions.yml") , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "presence_absence"
    """
    pangwas presence_absence \\
      $args \\
      --clusters $clusters \\
      --output ${prefix}.Rtab \\
      2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    def prefix = "presence_absence"
    """
    touch ${prefix}.Rtab
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """
}

process COMBINE {
    label "process_single"

    publishDir path: "${params.outdir}/combine", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    val(rtab)

    output:
    path("*combine.Rtab"), emit: rtab
    path("*.log")        , emit: log
    path("*versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ""
    def rtab = rtab.join(" ")
    """
    pangwas combine $args --rtab $rtab 2>&1 | tee combine.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS    
    """

    stub:
    """
    touch combine.Rtab
    touch combine.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
    END_VERSIONS
    """
}

process TREE {
    label "process_high"

    publishDir path: "${params.outdir}/tree", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(alignment)
    path(constant_sites)

    output:
    path "*.treefile",                  emit: tree
    path "*.rooted.nwk",                emit: tree_rooted
    path "*.rooted.branch_support.tsv", emit: branch_support
    path "*.rooted.labelled_nodes.nwk", emit: tree_labelled_nodes
    path "*.rooted.plain.nwk",          emit: tree_plain
    path "*iqtree.log",                 emit: iqtree_log
    path "*.log",                       emit: log
    path "*versions.yml",               emit: versions
    path "*.model.gz",                  emit: model          , optional: true
    path "*.best_scheme",               emit: best_scheme    , optional: true
    path "*.best_scheme.nex",           emit: best_scheme_nex, optional: true
    path "*.bionj",                     emit: bionj          , optional: true
    path "*.mldist",                    emit: mldist         , optional: true

    script:
    def args   = task.ext.args ?: ""
    def constant_sites_arg = constant_sites ? "--constant-sites ${constant_sites}": ""
    """
    pangwas tree \\
      $args $constant_sites_arg \\
      --alignment ${alignment} \\
      --threads ${task.cpus} \\
      2>&1 | tee tree.pangwas.log

    mv tree.log iqtree.log
    mv tree.pangwas.log tree.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch tree.treefile
    touch tree.model.gz
    touch tree.best_scheme
    touch tree.best_scheme.nex
    touch tree.bionj
    touch tree.mldist
    touch tree.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
    END_VERSIONS
    """
}

process GWAS {
    label "process_high"

    publishDir path: "${params.outdir}/gwas/${trait}", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    val(trait)      // Variable to test with a GWAS, must match a column in the table
    path(table)     // Table of variables in TSV or CSV format (ex. samplesheet.csv)
    path(variants)  // Rtab of variants from (ex. from COMBINE)
    path(tree)      // Optional: Rooted tree in NEWICK format from TREE, set to [] to disable
    path(clusters)  // Optional: Clusters in TSV format from SUMMARIZE, set to [] to disable

    output:   
    // Locus Effects
    path("*.locus_effects.tsv"),                    emit: locus_effects_all
    path("*.locus_effects.patterns.txt"),           emit: locus_effects_patterns
    path("*.locus_effects.bonferroni.txt"),         emit: locus_effects_bonferonni
    path("*.locus_effects.significant.tsv"),        emit: locus_effects_significant
    path("*.locus_effects.significant.filter.tsv"), emit: locus_effects_filtered
    path("*.focal.txt"),                            emit: focal
    path("*.log"),                                  emit: log
    path("*versions.yml"),                          emit: versions
    // Distance matrices
    path("*kinship.tsv"),                           emit: kinship, optional: true
    path("*patristic.tsv"),                         emit: patristic, optional: true
    path("*.filter.nwk"),                           emit: tree, optional: true
    path("*midpoint.nwk"),                          emit: tree_midpoint, optional: true
    // Lineage Effects
    path("*.lineage_effects.tsv"),                  emit: lineage_effects, optional: true
    path("*.lineage_effects.significant.tsv"),      emit: lineage_effects_significant, optional: true
    // Other
    path("*.binarize.*"),                           emit: binarize, optional: true
    path("*.qq_plot.png"),                          emit: qq_plot , optional: true    

    script:
    def args = task.ext.args ?: ""
    def tree_arg = tree ? "--tree $tree" : ""
    def clusters_arg = clusters ? "--clusters $clusters" : ""
    """
    pangwas gwas \\
      $args $clusters_arg $tree_arg \\
      --column $trait \\
      --table $table \\
      --variants $variants \\
      --threads ${task.cpus} \\
      2>&1 | tee ${trait}.gwas.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        pyseer: \$( pyseer --version | cut -d " " -f 2 )
    END_VERSIONS
    """

    stub:
    """
    touch gwas.locus_effects.tsv
    touch gwas.locus_effects.patterns.txt
    touch gwas.locus_effects.bonferroni.txt
    touch gwas.locus_effects.significant.tsv
    touch gwas.locus_effects.significant.filter.tsv
    touch gwas.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        pyseer: \$( pyseer --version | cut -d " " -f 2 )
    END_VERSIONS
    """    
}

process HEATMAP {
    label "process_single"

    publishDir path: "${params.outdir}/heatmap/${params.trait}", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(gwas)  // GWAS results in TSV format from GWAS
    path(tree)  // Optional: Rooted phylogeny in NEWICK format from TREE, set to [] to disable
    path(focal) // Optional: TXT of focal sample IDs, set to [] to disable

    output:
    path("*.svg")        , emit: svg
    path("*.png")        , emit: png
    path("*.log")        , emit: log
    path("*versions.yml") , emit: versions

    script:
    def args = task.ext.args ?: ""
    def focal_arg = focal ? "--focal $focal" : ""
    def tree_arg = tree ? "--tree $tree" : ""
    """
    prefix=\$(basename $gwas | sed -E 's/\\..*//g')
    pangwas heatmap \\
        $args $focal_arg $tree_arg \\
        --gwas $gwas \\
        --prefix \${prefix} \\
        2>&1 | tee -a \${prefix}.plot.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
        svgpathstools: unknown
        cairo: \$(python -c "import cairo; print(cairo.version)")
        cairosvg: \$(python -c "import cairosvg; print(cairosvg.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch plot.svg
    touch plot.png
    touch plot.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
        svgpathstools: unknown
        cairo: \$(python -c "import cairo; print(cairo.version)")
        cairosvg: \$(python -c "import cairosvg; print(cairosvg.__version__)")
    END_VERSIONS
    """ 
}

process MANHATTAN {
    label "process_single"

    publishDir path: "${params.outdir}/manhattan/${params.trait}", mode: "copy", overwrite: true

    container "ghcr.io/phac-nml/pangwas:0.1.0"
    conda "pangwas=0.1.0"

    input:
    path(gwas)  // GWAS results in TSV format from GWAS
    path(bed)   // BED coordinates from ALIGN


    output:
    // plots will only be present if there are significant variants to plot
    path("*.svg")          , emit: svg,       optional: true
    path("*.png")          , emit: png,       optional: true
    path("*phandango.plot"), emit: phandango, optional: true
    path("*.log")          , emit: log
    path("*versions.yml")  , emit: versions

    script:
    def args = task.ext.args ?: ""
    """
    prefix=\$(basename $gwas | sed -E 's/\\..*//g')
    pangwas manhattan \\
        $args \\
        --gwas $gwas \\
        --bed $bed \\
        --prefix \${prefix} \\
        2>&1 | tee -a \${prefix}.plot.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        svgpathstools: unknown
        cairo: \$(python -c "import cairo; print(cairo.version)")
        cairosvg: \$(python -c "import cairosvg; print(cairosvg.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch plot.svg
    touch plot.png
    touch plot.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangwas: \$(pangwas --version | cut -d " " -f 2 | sed 's/v//g')
        svgpathstools: unknown
        cairo: \$(python -c "import cairo; print(cairo.version)")
        cairosvg: \$(python -c "import cairosvg; print(cairosvg.__version__)")
    END_VERSIONS
    """ 
}