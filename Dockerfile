FROM mambaorg/micromamba:2.0.2-ubuntu24.10

# Install conda dependencies
COPY environment.yml /opt/pangwas/environment.yml
RUN micromamba create -y -f /opt/pangwas/environment.yml -n pangwas
RUN micromamba clean --all

# Transfer local files into container, runs as 'root', so we will
# put them in the opt dir for now, and then move them later
COPY data /opt/pangwas/data
COPY Dockerfile /opt/pangwas/Dockerfile
COPY docs /opt/pangwas/docs
COPY nextflow.config /opt/pangwas/nextflow.config
COPY nextflow /opt/pangwas/nextflow
COPY pyproject.toml /opt/pangwas/pyproject.toml
COPY README.md /opt/pangwas/README.md
COPY src /opt/pangwas/src
COPY tests /opt/pangwas/tests

# Transfer ownership to the non-root user (who is allowed to write to tmp)
RUN cp -r /opt/pangwas /tmp/pangwas

# Install nextflow plugins
RUN cd /tmp/pangwas && micromamba run -n pangwas nextflow run . --help

# Install python package and CLI
RUN cd /tmp/pangwas && micromamba run -n pangwas pip install .

# Minimal test of CLI and import
RUN micromamba run -n pangwas pangwas --help
RUN micromamba run -n pangwas python -c "import pangwas; print(pangwas.__version__)"

# Change the default env to activate from 'base' to 'pangwas'
ENV ENV_NAME="pangwas"
WORKDIR /tmp/pangwas

# No default command, because we might use either CLI or nextflow