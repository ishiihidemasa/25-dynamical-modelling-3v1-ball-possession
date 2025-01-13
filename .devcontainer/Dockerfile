# last update: 2025-01-10
FROM quay.io/jupyter/datascience-notebook:julia-1.11.1

WORKDIR /workdir
EXPOSE 8888

# install packages
COPY Project.toml Manifest.toml ./
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()'
# set default to the current project
ENV JULIA_PROJECT=@.
# enable multi-threading
ENV JULIA_NUM_THREADS=16
