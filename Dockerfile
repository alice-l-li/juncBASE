# Pre-download Miniconda installer in same directory

# Use a base image with ubuntu
FROM --platform=linux/amd64 ubuntu:24.04

LABEL authors="Alice Li" maintainer="ali5@fredhutch.org"


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory inside the container
WORKDIR /

# Update package lists, install the desired software (e.g., curl), and clean up
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    bzip2 \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy the pre-downloaded Miniconda installer into the image
# Make sure Miniconda3-latest-Linux-x86_64.sh is in the same folder as your Dockerfile
COPY Miniconda3-latest-Linux-x86_64.sh /tmp/miniconda.sh

RUN bash /tmp/miniconda.sh -b -p /opt/miniconda \
    && rm /tmp/miniconda.sh

# Initialize conda
ENV PATH=/opt/miniconda/bin:$PATH

# Make conda use system CA certificates
RUN conda config --system --set ssl_verify /etc/ssl/certs/ca-certificates.crt

# Create environment
COPY environment.yml /tmp/environment.yml
RUN conda config --system --remove channels defaults || true
RUN conda env create -f /tmp/environment.yml \
    && conda clean -afy

# Copy juncBASE scripts into a fixed location
WORKDIR /juncBASE
COPY . /juncBASE

# Auto-activate environment for interactive shells
RUN echo "source /opt/miniconda/etc/profile.d/conda.sh" >> /etc/bash.bashrc \
    && echo "conda activate juncBase-dependencies" >> /etc/bash.bashrc

ENTRYPOINT ["/bin/bash"]