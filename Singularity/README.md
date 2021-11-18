# Setup Singularity on your machine

For installing Singularity on Debian based systems like Debian, Ubuntu or Mint you can use this [setup file](https://github.com/buschlab/HowTos/blob/master/Singularity/setupSingularity.sh).

# Recipes for OMICS Singularity Containers

### List of available recipes

1. [ASERecipe](https://github.com/buschlab/HowTos/blob/master/Singularity/ASERecipe)
   Contains all Tools required for allele specific expression pipeline from *fastqc* to *phASER*.
2. [PipitsRecipe](https://github.com/buschlab/HowTos/blob/master/Singularity/PipitsRecipe)
   Container for Pipits pipeline. 
3. [CNVkit](https://github.com/buschlab/HowTos/blob/master/Singularity/cnvkit.def) Contains *CNVkit* for copy number analysis.
4. [BamSnap](https://github.com/buschlab/HowTos/blob/master/Singularity/bamsnap.def) Contains *BamSnap*, a visualization tool for BAM files.

### Pointers

###### Build Errors

```
sudo singularity build ASE_neu.simg ASERecipe 

Using container recipe deffile: ASERecipe
Sanitizing environment
Adding base Singularity environment to container
ERROR: 'Bootstrap' type not supported: docker
Cleaning up...
```

The cause of this error might be end of line artefacts that can be cleaned up with the *dos2unix* application.
