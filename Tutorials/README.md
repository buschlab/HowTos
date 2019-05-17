## CRAM to FASTQ conversion
Since the format only the difference between the reference and the aligned read, the reference genome is needed for the conversion to fastq. Usually, the information about the reference genome is requested from an **EBA** webservice during conversion. However, the webservice doesn't always respond and therefore produces empty fastq files. Hence, it is recommended to create a local reference database with **SAMtools** prio to conversion.

### Workflow
1. Extract name of reference genome from cram header and download **.fasta** file
2. Create reference database with ```/path/to/samtools/misc/seq_cache_populate.pl```
3. Set appropriate environment variables
4. Use **SAMtools sort** and **fastq** to create valid fastq files

### Example script

```#! /bin/bash
#SBATCH --cluster=omics
#SBATCH --partition=shortterm,longterm
#SBATCH --no-requeue
#SBATCH --mem=180GB
#SBATCH -c 48
#SBATCH --time 4320


path_data=/data/lied_egypt_genome/ega_data_EGAD00001001372


ls *cram | sed 's/\_/ /g' | uniq | sort -u | awk -v path_data=$path_data '{print "samtools sort -n $SCRATCH/cram/"$1"_"$2"_"$3" | samtools fastq -1 $SCRATCH/fastq/"$1".R1.fq -2 $SCRATCH/fastq/"$1".R2.fq -; gzip $SCRATCH/fastq/"$1".R*.fq; mv $SCRATCH/fastq/"$1".R*.fq.gz "path_data"/fastq"}' >cram2fastq.sh;


export XDG_CACHE_HOME=$SCRATCH
cp -r cram $SCRATCH
cp -r hs37d5.fa $SCRATCH
mkdir $SCRATCH/fastq
cd $SCRATCH


seq_cache_populate.pl -root $SCRATCH/cache $SCRATCH/hs37d5.fa
export REF_PATH=$SCRATCH/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
export REF_CACHE=$SCRATCH/cache/%2s/%2s/%s

~/workspace/toolbox/parallelize_cmds.pl 48 1 $path_data/cram2fastq.sh
```
