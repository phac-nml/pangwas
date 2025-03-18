# Developers

## Testing

- **panGWAS** aims for >=95% code coverage in it's unit tests.
- Unit tests are defined in: `tests/config/test_pangwas.yaml`
- The following are some examples of typical unit tests.

```yaml
# -----------------------------------------------------------------------------
# extract

extract:

  # Basic unit test
  - name: extract
    params:
    function: pangwas.extract
    tags: ["default", "cov"]
    args:
        gff: "{data_dir}/{sample}.gff3"
        min_len: "{min_len}"
        outdir: "{outdir}"
    output:
        tsv: "{outdir}/{prefix}.tsv"
    variables:   
        data_dir: "data/test/gff"
        min_len: "10"
        outdir: "tests/observed/test/extract"
        prefix: "{sample}"
        sample:
        - "sample1"
        - "sample2"
        - "sample3"
        - "sample4"

  # Unit test that is expected to fail with a particular error
  - name: duplicate_contig_error
    params:
      tags: ["cov"]
      function: pangwas.extract
      args:
        gff: "{test_data}/{name}.gff3"
        outdir: "{outdir}"
        prefix: "{prefix}"        
      error_message: "Duplicate contig ID found: sample2_pesticin"
      variables:
        outdir: "tests/observed/test/extract/{name}"
        prefix: "{name}"
        test_data: "tests/data"

# -----------------------------------------------------------------------------
# cluster

cluster:

  # run the cluster function with various resources, and check that the
  # results are the same every single time.
  - name: reproducible
    params:
      tags: ["reproducible"]
      function: pangwas.cluster
      args:
        fasta: "{collect_dir}/sequences.fasta"
        outdir: "{outdir}"
        tmp: "{outdir}/tmp"
        memory: "{memory}"
        threads: "{threads}"
      output:
        clusters: "{outdir}/clusters.tsv"
        representative: "{outdir}/representative.fasta"
      variables:
        collect_dir: "tests/observed/test/collect"
        outdir: "tests/observed/test/cluster/t{threads}_mem{memory}"
        prefix: "t{threads}_mem{memory}"
        threads:
          - 1
          - 2
          - 4
        memory:
          - "1G"
          - "2G"
          - "4G"
```