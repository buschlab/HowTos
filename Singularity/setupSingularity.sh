#!/bin/bash
export GOVERSION=1.15.3
export SINGULARITYVERSION=3.6.4

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
wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITYVERSION}/singularity-${SINGULARITYVERSION}.tar.gz && \
tar -xzf singularity-${SINGULARITYVERSION}.tar.gz && \
cd singularity
./mconfig && \
make -C ./builddir && \
make -C ./builddir install
cd .. && \
rm -rf singularity singularity-$SINGULARITYVERSION.tar.gz
