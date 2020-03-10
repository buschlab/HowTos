# <a name="startofpage"></a> How to process and analyse Metagenomics

### Sections:
1. ### [Kaiju](#kaijumain)
   Fast and sensitive taxonomic classification for metagenomics
2. ### [Humann2](#humann2main) 
   HUMAnN is a pipeline for efficient and accurate functional profiling of microbial 
   pathways in a community from metagenomic or metatranscriptomic sequencing data, aimed 
   to describe the metabolic potential of a microbial community and its members.
   

## <a name="kaijumain"></a> Kaiju
0. [Prerequisites](#kaijuzero)
1. [Pre-Processing](#kaijuone)
2. [Kaiju databases](#kaijutwo)
3. [Run kaiju](#kaijuthree)
4. [Downstream - Krona](#kaijutfour)
5. [Downstream - R](#kaijufive)

#### <a name="kaijuzero"></a> 0. Prerequisites 

This Workflow consists of 3 main steps.

1. ##### Raw data pre-processing
2. ##### running Kaiju (Protein-to-DNA classifier)
3. ##### Import into R for further processing

To that end you'll need to install [KneadData](https://bitbucket.org/biobakery/kneaddata/wiki/Home), [Kaiju](https://github.com/bioinformatics-centre/kaiju) and [R](https://www.r-project.org), as well as their respective requirements and additional packages.

#### <a name="kaijuone"></a> 1. Pre-processing
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

#### <a name="kaijutwo"></a> 2. Kaiju databases

Several databases are available, incorporating data from bacterial to viral genomes. Refer to [Kaiju tutorial](https://github.com/bioinformatics-centre/kaiju) for an overview of DBs and content.vIn this case we chose the NCBI BLASTplus DB as is it is the most comprehensive. The DB can found and downloaded directly on the developers [website](http://kaiju.binf.ku.dk/server) or downloaded via the programs `kaiju-makedb` command.

```bash

mkdir kaijudb
cd kaijudb
kaiju-makedb -s nr_euk

```

#### <a name="kaijuthree"></a> 3. Run Kaiju

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

###### Supplementary

Sometimes the raw-data comes with weird header, wrong linebreak symbols or any number of formating issues that exacerbate the work.

An issue with Kaiju are the format of sequence identifiers in paired-end files and the oder of reads in both files. Luckily, the tool-kit [SeqKit](https://bioinf.shenwei.me/seqkit/) can help with all sorts of formatting and sorting issues. The following code checks seq-identifiers, sort and filters paired-end reads to prepare them for Kaiju.

```bash

### 1. sort & remove duplicates 
ls *.fastq | sed 's/\./ /g' | awk '{print "seqkit rmdup "$1"."$2"."$3" | seqkit replace --pattern \" .+\" --replacement \" 1\" | seqkit sample --proportion 0.9 --rand-seed 1 --out-file "$1".reads.fq.gz"}' > sortNremove.sh
bash sortNremove.sh

### 2. extract intersecting reads
ls *1.reads.fq.gz | sed 's/_/ /g' | awk '{print "seqkit --name --only-id "$1"_"$2"_"$3"_"$4"_"$5" "$1"_"$2"_"$3"_"$4"_2.reads.fq.gz | sort | uniq -d > "$1".txt"}' > uniqueReads.sh
bash uniqueReads.sh

### 3. extract matching reads to create new files
ls *.reads.fq.gz | sed 's/_/ /g' | awk '{print "seqkit grep --pattern-file "$1".txt "$1"_"$2"_"$3"_"$4"_"$5" -o "$1"_unique_"$4"_"$5}' > printReads.sh
bash printReads.sh

### 4. sort both files aaaand you're done
ls *_unique_* | sed 's/\./ /g' | awk '{print "gzip -d -c "$1"."$2".fq.gz | seqkit fx2tab | sort -k1,1 -T . | seqkit tab2fx | gzip -c > "$1".sorted.fq.gz"}' > finalSort.sh
bash finalSort.sh


```



#### <a name="kaijufour"></a> 4. Downstream - Krona

The Kaiju output can be converted to a format fitting the requirements of Krona, as described in 
[Kaiju github source](https://github.com/bioinformatics-centre/kaiju). And then imported into [Krona](https://github.com/marbl/Krona/wiki).

```bash

kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona

```


#### <a name="kaijufive"></a> 5. Downstream - R

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


## <a name="humann2main"></a> Humann2 

[HUMAnN2](http://huttenhower.sph.harvard.edu/humann) is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data.

0. [Prerequisites](#humann2zero)
1. [Prepare DBs](#humann2one)
2. [Run pipeline](#humann2two)
3. [Output](#humann2three)
4. [Downstream - R](#humann2four)


#### <a name="humann2zero"></a> 0. Prerequisites 
  
Since [HUMAnN2](http://huttenhower.sph.harvard.edu/humann) basically is a complete workflow for metagenomic-data, it incorporates a variety of tools and libraries. The probably easiest way to install is via pip, because it takes care of most dependencies ([Bowtie](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Diamond](https://diamond.readthedocs.io/en/latest/)) assuming that the correct flags have been set during installation. 

```bash

## HUMAnN2 new-install on an OMICS-Cluster with limited write permissions
#	`--user` to install in home directory
# 	`--build-diamond` compiles diamond from source
#	`--replace-dependencies-install`re-installs diamond and bowtie2 if already present

pip install humann2 --user --install-option='--build-diamond' --install-option='--replace-dependencies-install'

```

[MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) still has to be installed seperately. After Installation it is advisable to execute the provided test-scripts `humann2_test --run-functional-tests-tools` in order to ensure functionality of the pipeline. Please refer to the [HUMAnN2](http://huttenhower.sph.harvard.edu/humann) website to find detailed instructions.

Technically, no pre-processing of the input files is required. I recommend to do it anyhow, because:
1. If you intend to use other tools, such as Kaiju, you'll need to do it anyway
2. For my taste it's more dependable and format errors are easier to identify/fix than after a failed pipeline run.
3. Since [HUMAnN2](http://huttenhower.sph.harvard.edu/humann) makes no practical use of paired-end read information, it is recommended/required to concatenate both files into a single one. So you'll have to do some work in any case.



#### <a name="humann2one"></a> 1. Prepare DBs

[HUMAnN2](http://huttenhower.sph.harvard.edu/humann) requires both a nucleotide and and protein database to work properly. In this case we used ChocoPhlAn and UniRef90 DBs.

The provided function `humann2_databases` downloads the databases and adjusts the internal pathways of the humann2 framework. However, in case Cluster restrictions prevent you from using the work-directory for execution - or you want to download and check out, you can provide the DB-location on execution. The flags `--nucleotide-database` and `--protein-database` require the path to the respective DBs.

```bash

## Download nucleotide DB
humann2_databases --download chocophlan full $INSTALL_LOCATION

## Download protein DB - uniref90_diamond is the most comprehensive and thus largest DB
humann2_databases --download uniref uniref90_diamond $INSTALL_LOCATION

```


#### <a name="humann2two"></a> 2. Run pipeline

When all the requirements are satisfied and the appropriate DBs have been selected, it is time to run the pipeline. This is very easy, as you can see in the following comand. One note would be to supply the paths to `--metaphlan`, `--bowtie2` and `--diamond` manually. At least on an OMICS-cluster system the execution was more reliable that way.

```bash

## List files and prepare execution script.
ls *.fastq.gz | sed 's/_/ /g' | awk -v threads=$THREADS -v scratch=$SCRATCH '{print "humann2 --metaphlan /path/to/metaphlan2/ --diamond /path/to/diamond/ --bowtie2 /path/to/bowtie2/  --input "$1"_"$2"_"$3"_"$4" --output ../results/"$1" --threads "threads" --nucleotide-database /db/humanndb/chocophlan --protein-database /db/humanndb/uniref"}' > start_humann2.sh
## execute
bash start_humann2.sh

```


#### <a name="humann2three"></a> 3. Output

The pipeline creates 3 main output files.
1. ID_genefamilies.tsv - Gene families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
2. ID_pathabundance.tsv - The abundance of each pathway in the community as a function of the abundances of the pathway's component reactions, with each reaction's abundance computed as the sum over abundances of genes catalyzing the reaction.
3. ID_pathcoverage.tsv - Assigns a confidence score (I=[0,1]) to pathways and associated reactions. So it's more of a qualitytive measure for the detected pathways.

Those files can for example be used with the provided functions for downstream analyses. Next to those files, the pipeline also creates outputs for intermediate steps. Those can be used for downstream analyses as well as re-runs of parts of the pipeline. 

For information on integrated downstream analyses and intermediate output-files, please refer to the [Manual](https://github.com/biobakery/humann#initial-installation).


#### <a name="humann2four"></a> 4. Downstream - R

The following transformations are helpful for downstream analyses and facilitate the import to R. 

##### Preparation I

```bash

## Regroup genes to functional categories, i.e., all the reads that are associated to one gene/protein are grouped together - resulting in a total count for the gene/protein and the respective sub-classifications into single species (if possible)
# This step can be omitted, but it beats comming up with a personal strategy for aggregation of functional groupings.
ls *genefamilies.tsv | sed 's/_/ /g' | awk '{print "humann2_regroup_table --input "$1"_"$2"_"$3"_"$4"_"$5" --output "$1"_genefamilies_regrouped.tsv --groups uniref90_pfam"}' > H2_regroup.sh

## Normalize to copies per million 
# flag `--special n` excludes unmapped reads from normalization
# flag `--units relab`calculates relative abundance
ls *genefamilies_regrouped.tsv | sed 's/_/ /g' | awk '{print "humann2_renorm_table --special n --input "$1"_"$2"_"$3" --output "$1"_genefamilies_cpm.tsv --units cpm"}' > H2_cpm.sh

## Join the tables - input: folder containing normalized ouput files
humann2_join_tables --input pathways/ --output MG_genefamilies_cpm.tsv --file_name genefamilies_cpm

## Attach names to features, i.e., add description of function to the entities
humann2_rename_table --input MG_genefamilies_cpm.tsv --output MG_genefamilies_cpm_pfam.tsv --names pfam

## Fix ID-denominators
# flag `--update-snames` (not used here) can be used to keep track of the performed operations by adding them to the respective sample-IDs. In any case, the sID updated with the output filename when the pipeline finishes. 
# to fix sIDs you can simply use the following `sed` command or just employ `gsub` after import to R
sed 's/_unique_paired_concat_Coverage//g' MG_genefamilies_cpm_pfam.tsv > MG_genefamilies_cpm_pfam_fixedID.tsv

```

##### Preparation II
Similar preparations can be done for pathway-abundance or coverage files.

```bash
## Normalize pathway-abundance tables
# relative Abundance
ls *pathabundance.tsv | sed 's/_/ /g' | awk '{print "humann2_renorm_table --special n --input "$1"_"$2"_"$3"_"$4"_"$5" --output "$1"_pathway_relab.tsv --units relab --update-snames"}' > H2_relAbundance.sh
# copies per million
ls *pathabundance.tsv | sed 's/_/ /g' | awk '{print "humann2_renorm_table --special n --input "$1"_"$2"_"$3"_"$4"_"$5" --output "$1"_pathway_cpm.tsv --units cpm"}' > H2_CPM.sh

## Join pathway-abundance tables
# relative Abundance
humann2_join_tables --input pathways/ --output MG_pathabundance_relab.tsv --file_name relab
# copies per million
humann2_join_tables --input pathways/ --output MG_pathabundance_cpm.tsv --file_name cpm

### Fix ID-denominators
# relative Abundance
sed 's/_unique_paired_concat_Abundance//g' MG_pathabundance_relab.tsv > MG_Adina_pathabundance_relab_fixedID.tsv
# copies per million
sed 's/_unique_paired_concat_Abundance//g' MG_pathabundance_cpm.tsv > MG_pathabundance_cpm_fixedID.tsv

```bash

##### Import
And finally the import and re-format operations for processing in R.

```r

## Import any kind of textfile is most convenient with the `fread()` function.
genefam_cpm_pfam <- data.table::fread("MG_genefamilies_cpm_pfam.tsv", data.table=FALSE)
# Fix sample-IDs if not done already
colnames(genefam_cpm_pfam) <- gsub("_unique_paired_concat_Abundance-RPKs","",colnames(genefam_cpm_pfam))
# I like to save my imports in a separate file
saveRDS(genefam_cpm_pfam, file=paste0(FILESOURCE,"rds/humann2_genefam_cpm_nS_regrouped_pfam.rds"), compress = TRUE)

```

Some additional formatting is required, because - GeneID,Function,Genus,Species are concatenated in a single string. The subsequent function should accomplish the separation on basic output-files, regrouped genefamilies and annotated pathway-coverage/abundance files. But, no guarantees..

```r

## 'prepPWS' reformats humann2 pathabundance and pathcoverage tables
# @param pw_table: table to prepare
# @param dropZrows: set TRUE to remove rows that contain only zeros, i.e., no abundance or not confidently detected
# @return pw_table: the reformatted table
prepPWS <- function(i_table, dropZrows=TRUE) {
  if(dropZrows) {
    i_table <- i_table[apply(i_table[,-1], 1, function(x) !all(x==0)),]
  }
  if(any(colnames(i_table) %in% "# Pathway")) {
    # rename first column
    colnames(i_table)[1] <- "pathways"
    # separate the column
    tmp_split_1 <- data.frame(do.call('rbind', strsplit(as.character(i_table$pathways),'|',fixed=TRUE)), stringsAsFactors = FALSE)
    # replace denominator of unclassified counts with NAs
    tmp_split_1$X2 <- sapply(1:dim(tmp_split_1)[1], function(x) ifelse(tmp_split_1$X2[x] == tmp_split_1$X1[x], NA,tmp_split_1$X2[x]))
    # split first column into pathway denominator and task-desciption
    tmp_split_2 <- data.frame(do.call('rbind', strsplit(as.character(tmp_split_1$X1),': ',fixed=TRUE)), stringsAsFactors = FALSE)
    # replace with NAs
    tmp_split_2$X2 <- sapply(1:dim(tmp_split_2)[1], function(x) ifelse(tmp_split_2$X2[x] == tmp_split_2$X1[x], NA,tmp_split_2$X2[x]))
    # split taxonomy into genus and species
    tmp_split_3 <- data.frame(do.call('rbind', strsplit(as.character(tmp_split_1$X2),'.',fixed=TRUE)), stringsAsFactors = FALSE)
    tmp_split_3$X1 <- gsub("g__","",tmp_split_3$X1, fixed=TRUE) 
    tmp_split_3$X2 <- gsub("s__","",tmp_split_3$X2, fixed=TRUE)
    # put humpdy-dumpdy back together again
    tmp_frame <- data.frame("pathways"=tmp_split_2$X1, "task"=tmp_split_2$X2, "genus"=tmp_split_3$X1, 
      "species"=tmp_split_3$X2, stringsAsFactors = FALSE)
    tmp_frame <- cbind(tmp_frame, i_table[,-1])
  } else {
    colnames(i_table)[1] <- "gene_family"
    
    if(grepl(":",i_table[1,1])) {
      tmp_split_1 <- reshape2::colsplit(i_table[,1],": ",c("X1","X2"))
      tmp_split_2 <- reshape2::colsplit(tmp_split_1$X2,"\\|",c("X1","X2"))
      tmp_split_3 <- reshape2::colsplit(tmp_split_2$X2,"\\.",c("genus","species"))
      tmp_split_3$genus <- gsub("g__","",tmp_split_3$genus, fixed=TRUE) 
      tmp_split_3$species <- gsub("s__","",tmp_split_3$species, fixed=TRUE)
      
      tmp_frame <- data.frame("ID"=tmp_split_1$X1, "task"=tmp_split_2$X1, "genus"=tmp_split_3$genus, "species"=tmp_split_3$species, i_table[,-1], stringsAsFactors = FALSE)
      tmp_frame$genus[nchar(tmp_frame$genus) == 0] <- NA
      tmp_frame$species[nchar(tmp_frame$species) == 0] <- NA
    } else {
      # separate the column
      tmp_split_1 <- data.frame(do.call('rbind', strsplit(as.character(i_table$gene_family),'|',fixed=TRUE)), stringsAsFactors = FALSE)
      # replace denominator of unclassified counts with NAs
      tmp_split_1$X2 <- sapply(1:dim(tmp_split_1)[1], function(x) ifelse(tmp_split_1$X2[x] == tmp_split_1$X1[x], NA,tmp_split_1$X2[x]))
      # split first column into pathway denominator and task-desciption
      tmp_split_2 <- data.frame(do.call('rbind', strsplit(as.character(tmp_split_1$X1),': ',fixed=TRUE)), stringsAsFactors = FALSE)
      # replace with NAs
      tmp_split_2$X2 <- sapply(1:dim(tmp_split_2)[1], function(x) ifelse(tmp_split_2$X2[x] == tmp_split_2$X1[x], NA,tmp_split_2$X2[x]))
      # split taxonomy into genus and species
      tmp_split_3 <- data.frame(do.call('rbind', strsplit(as.character(tmp_split_1$X2),'.',fixed=TRUE)), stringsAsFactors = FALSE)
      tmp_split_3$X1 <- gsub("g__","",tmp_split_3$X1, fixed=TRUE)
      tmp_split_3$X2 <- gsub("s__","",tmp_split_3$X2, fixed=TRUE)
      # put humpdy-dumpdy back together again
      tmp_frame <- data.frame("ID"=tmp_split_2$X1, "task"=tmp_split_2$X2, "genus"=tmp_split_3$X1,
        "species"=tmp_split_3$X2, stringsAsFactors = FALSE)
      tmp_frame$ID <- gsub("^.*?_","",tmp_frame$ID)
      tmp_frame <- cbind(tmp_frame, i_table[,-1])
    }
  
    
    
  }
  return(tmp_frame)
}

```

The function is implemented in a way that allows you to filter for 'NA' in genus/taxonomy columns to get the total count of a pathway/gene/protein and drop the sub-classification.

```r

# keep only aggregated pw-counts
df_reduced <- df[is.na(df$taxonomy),]

# only keep aggregated protein-counts
df <- df[is.na(df$genus),]

```

From here on out you can basically apply every conceivable analysis.



###### Take me back to the [beginning](#startofpage)

