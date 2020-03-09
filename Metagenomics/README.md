# How to process and analyse Metagenomics

### Sections:
1. ### Kaiju
   Fast and sensitive taxonomic classification for metagenomics
2. ### Humann2 
   HUMAnN is a pipeline for efficient and accurate functional profiling of microbial 
   pathways in a community from metagenomic or metatranscriptomic sequencing data, aimed 
   to describe the metabolic potential of a microbial community and its members.
   

## Kaiju
0. [Prerequisites](#0.-prerequisites)
1. [Pre-Processing](#1.-pre-processing)
2. Kaiju databases
3. Run kaiju
4. Downstream - Krona
5. Downstream - R

#### 0. Prerequisites

This Workflow consists of 3 main steps.

1. ##### Raw data pre-processing
2. ##### running Kaiju (Protein-to-DNA classifier)
3. ##### Import into R for further processing

To that end you'll need to install [KneadData](https://bitbucket.org/biobakery/kneaddata/wiki/Home), [Kaiju](https://github.com/bioinformatics-centre/kaiju) and [R](https://www.r-project.org), as well as their respective requirements and additional packages.

#### 1. Pre-processing

Before running Kaiju the input-files have to be pre-processed and sanitised, i.e., perform quality control and remove all unwanted reads by mapping against the host organism. The tool for this job is 
[KneadData](https://bitbucket.org/biobakery/kneaddata/wiki/Home).

- multiple databases can be used by adding the `-db` flag and specifying additional paths
- at least on a computational cluster it was a good idea to specify the paths to both [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) by using the respective flags in the call to KneadData.
- The `--remove-intermediate-output` flag makes for a cleaner output directory and less storage overhead. 

In this case the hosts were B6-mice and for a stringent workflow the sequences were also mapped against a human reference to remove possible contamination.

```bash
# set variables to use
WORKDIR=$"/path/to/project/folder"
REFDB=$WORKDIR"/db"

cd $WORKDIR
# download DBs
# black 6 mouse reference-db
kneaddata_database --download mouse_C57BL bowtie2 $WORKDIR"/db/BT2_mouse/"
# human reference-db
kneaddata_database --download human_genome bowtie2 $SCRATCH"/databases/"

cd data/RAW_READS

# Kneaddata - perform trimming and removal of host sequences
ls *_1.fastq | sed 's/_/ /g' | awk -v workdir=$WORKDIR -v refdb=$REFDB '{print "kneaddata --input "$1"_1.fastq --input "$1"_2.fastq -db "refdb"/BT2_mouse -db "refdb"/BT2_human -t 10 --trimmomatic-options SLIDINGWINDOW:4:20:MINLEN:50 --bowtie2 /path/to/bowtie2/installation --trimmomatic /path/to/trimmomatic/installation --remove-intermediate-output --output " workdir "/data/results" }' > 01_kneaddata_PER.sh

# execution
bash 01_kneaddata_PER.sh

```

#### [2. Kaiju databases](#kaiju-two)

Several databases are available, incorporating data from bacterial to viral genomes. Refer to [Kaiju tutorial](https://github.com/bioinformatics-centre/kaiju) for an overview of DBs and content.vIn this case we chose the NCBI BLASTplus DB as is it is the most comprehensive. The DB can found and downloaded directly on the developers [website](http://kaiju.binf.ku.dk/server) or downloaded via the programs `kaiju-makedb` command.

```bash

mkdir kaijudb
cd kaijudb
kaiju-makedb -s nr_euk

```

#### [3. Run Kaiju](#kaiju-three)

Kaiju requires `nodes.dmp`and `kaiju_db.fmi` as well as the pre-processed input file.

```bash

kaiju -t nodes.dmp -f kaiju_db_*.fmi -i inputfile.fastq

```

Kaiju can work on paired-end files, as long as they are sorted in the same order. Also, the `z` flag can be used to execute in parallel.

```bash

# Storing the path to in variables makes the script easier to reaad and maintain.
NODES=$WORKDIR"/db/kaijudb/BLASTplus/nodes.dmp"
FMI=$WORKDIR"/db/kaijudb/BLASTplus/kaiju_db_nr_euk.fmi"

# List files and construct execution-script
ls *paired_1.sorted.fq.gz | sed 's/_/ /g' | awk -v threads=$THREADS -v nodes=$NODES -v fmi=$FMI '{print "kaiju -z "threads" -t "nodes" -f "fmi" -i "$1"_"$2"_"$3"_"$4" -j "$1"_"$2"_"$3"_2.sorted.fq.gz -o ../results/"$1"_BLASTplus_kaiju.out"}' > PtD_classifier_blast.sh

# Execute
bash PtD_classifier_blast.sh

```


#### [4. Downstream - Krona](#kaiju-four)

The Kaiju output can be converted to a format fitting the requirements of Krona, as described in 
[Kaiju github source](https://github.com/bioinformatics-centre/kaiju). And then imported into [Krona](https://github.com/marbl/Krona/wiki).

```bash

kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona

```


#### [5. Downstream - R](#kaiju-five)

First the output-files have to be summarised and combined with a complete taxonomic classification


```bash

# The following commands assume that the db-files are located in the file-root folder db - adjust accordingly.
# Input-file names look like - "ID_DBname_kaiju.out"
# The general command is:
kaiju2table -t nodes.dmp -n names.dmp -r genus -o kaiju_summary.tsv kaiju.out [kaiju2.out, ...]

## 1. List all BLASTplus-DB output-files, summarise and add classification
ls *.out | sed 's/_/ /g' | awk '{print "kaiju2table -t /db/kaijudb/BLASTplus/nodes.dmp -n /db/kaijudb/BLASTplus/names.dmp -l superkingdom,phylum,class,order,family,genus,species -r genus -o "$1"_"$2"_summary.tsv "$1"_"$2"_"$3}' > start_sum_kaiju.sh
## and execute
bash start_sum_kaiju.sh

```

Now import and join the sample files in R.

```r

### IMPORT & JOIN KAIJU Output
importKaiju <- function( filepath ) {
  ## list files in directory
  sample_id <- dir(filepath)
  input_files <- file.path(filepath, sample_id)
  ## reduce names to IDs - assuming that string before first underscore contains sample-IDs
  names(input_files) <- paste0(sapply(sample_id, function(x) unlist(strsplit(x, "_", fixed = FALSE))[1])) 
  ## start file import
  
  print("DEBUG: start importing")
  
  first <- TRUE
  for(idx in input_files) {
    if(first) {
      main_table <- data.table::fread(idx, select=c(3:5),fill=TRUE, data.table=FALSE) # ToDo: adjustable column-selection
      main_table$taxon_id <- as.character(main_table$taxon_id)
      
      colnames(main_table)[1] <- names(input_files[which(input_files %in% idx)])
      main_table <- main_table[,c(2,3,1)]
      first <- FALSE
      
    } else {
      tmp_table <- data.table::fread(idx, select=c(3:5),fill=TRUE, data.table=FALSE) # ToDo: adjustable column-selection
      tmp_table$taxon_id <- as.character(tmp_table$taxon_id)
      
      colnames(tmp_table)[1] <- names(input_files[which(input_files %in% idx)])
      print(paste("DEBUG: ID ", names(input_files[which(input_files %in% idx)]), " done"))
      
      # join tables
      main_table <- dplyr::full_join(main_table, tmp_table, by=c("taxon_id","taxon_name"))
      
    }
  }
  ## replace <NA> with zeros
  main_table[is.na(main_table)] <- 0
  
  return(main_table)
}

blastp_path <- paste(FILESOURCE, "results/kaiju/summaries/blastp", sep=""); sample_id <- dir(path)
## Import and join files
blastp_import <- importKaiju(blastp_path)
saveRDS(blastp_import, file=paste0(FILESOURCE,"rds/kaiju_blastp.rds"), compress = TRUE)

```

Format and filter to construct phyloseq-object.

```r

### CONSTRUCT phyloseq object
## check + drop duplicate ids - can be done with taxon_name as well
id_occur <- data.frame(table(blastp_import$taxon_id))
id_occur[id_occur$Freq > 1,]
## take a look
test_id <- blastp_import[blastp_import$taxon_id %in% id_occur$Var1[id_occur$Freq > 1],]
## drop
blastp_import <- blastp_import[blastp_import$taxon_id %in% id_occur$Var1[id_occur$Freq == 1],]
## split classification column
classification <- reshape2::colsplit(blastp_import$taxon_name,";",c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
classification$Species <- NA

# put humpdy-dumpdy together again
OTU = phyloseq::otu_table(data.matrix(blastp_import[,3:dim(blastp_import)[2]]), taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(as.matrix(classification))
physeq = phyloseq::phyloseq(OTU, TAX)
saveRDS(physeq, file=paste0(FILESOURCE,"rds/phyloseq_kaiju_blastp.rds"), compress = TRUE)

## ADD METADATA
### 1. load data
metadata <- readRDS(file=paste0(FILESOURCE,"rds/all_metadata.rds"), refhook = NULL)
physeq <- readRDS(file=paste0(FILESOURCE,"rds/phyloseq_kaiju_blastp.rds"), refhook = NULL)

### 2. combine data
# get all the sample names
(samples <- colnames(phyloseq::otu_table(physeq)))
rownames(metadata) <- metadata$ID
names(metadata)[1] <- "Group"

# load current meta-df 
# add metadata to phyloseq object
phyloseq::sample_data(physeq) <- metadata
saveRDS(physeq, file=paste0(FILESOURCE,"rds/phyloseq_kaiju_blastp.rds"), compress = TRUE)


```

## Humann2 
1. install which is tricky
2. run pipeline
3. handle the three given filetypes
4. magic
   
   

