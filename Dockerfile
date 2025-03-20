FROM mambaorg/micromamba:2.0.2-ubuntu24.10

# Transfer local files into container, runs as 'root', so we will
# put them in the opt dir for now, and then move them later
COPY conf /opt/pangwas/conf
COPY data /opt/pangwas/data
COPY Dockerfile /opt/pangwas/Dockerfile
COPY docs /opt/pangwas/docs
COPY environment.yml /opt/pangwas/environment.yml
COPY modules /opt/pangwas/modules
COPY nextflow.config /opt/pangwas/nextflow.config
COPY nextflow_schema.json /opt/pangwas/nextflow_schema.json
COPY pipeline.nf /opt/pangwas/pipeline.nf
COPY pyproject.toml /opt/pangwas/pyproject.toml
COPY README.md /opt/pangwas/README.md
COPY src /opt/pangwas/src
COPY subworkflows /opt/pangwas/subworkflows
COPY tests /opt/pangwas/tests

# Install the conda environment
# Transfer ownership to the non-root user (who is allowed to write to tmp)
# Install nextflow plugins
# Install python package and CLI
# Minimal test of CLI and import
# Change the default env to activate from 'base' to 'pangwas'
RUN micromamba create -y -f /opt/pangwas/environment.yml -n pangwas && \
    micromamba clean --all && \
    cp -r /opt/pangwas /tmp/pangwas && \
    cd /tmp/pangwas && micromamba run -n pangwas nextflow run . --help && \
    cd /tmp/pangwas && micromamba run -n pangwas pip install .[test] && \
    micromamba run -n pangwas pangwas --help && \
    micromamba run -n pangwas python -c "import pangwas; print(pangwas.__version__)"

# Change the default env to activate from 'base' to 'pangwas'
ENV ENV_NAME="pangwas"
WORKDIR /tmp/pangwas

# No default command, because we might use either CLI or nextflow