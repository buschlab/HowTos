#!/bin/bash
export GOVERSION=1.19.1
export SINGULARITYVERSION=3.9.9

sudo apt-get update && sudo apt-get install -y \
build-essential \
uuid-dev \
libgpgme-dev \
squashfs-tools \
libseccomp-dev \
wget \
pkg-config \
git \
cryptsetup-bin

wget https://golang.org/dl/go$GOVERSION.linux-amd64.tar.gz
tar -C /usr/local -xzvf go$GOVERSION.linux-amd64.tar.gz && \
rm go$GOVERSION.linux-amd64.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
export GOPATH=${HOME}/go
export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

source ~/.bashrc && \
wget https://github.com/sylabs/singularity/archive/refs/tags/v${SINGULARITYVERSION}.tar.gz && \
tar -xzf v${SINGULARITYVERSION}.tar.gz && \
cd singularity-${SINGULARITYVERSION}
echo ${SINGULARITYVERSION} > VERSION
./mconfig && \
make -C ./builddir && \
make -C ./builddir install
cd .. && \
rm -rf singularity-${SINGULARITYVERSION} v$SINGULARITYVERSION.tar.gz
