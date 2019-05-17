# Recipes for OMICS Singularity Containers

### List of available recipes

1. ASERecipe
   Contains all Tools required for allele specific expression pipeline from *fastqc* to *phASER*.
2. PipitsRecipe
   Container for Pipits pipeline. 

### Pointers

###### Build Errors

```sudo singularity build ASE_neu.simg ASERecipe 

Using container recipe deffile: ASERecipe
Sanitizing environment
Adding base Singularity environment to container
ERROR: 'Bootstrap' type not supported: docker
Cleaning up...
```

The cause of this error might be end of line artefacts that can be cleaned up with the *dos2unix* application.
