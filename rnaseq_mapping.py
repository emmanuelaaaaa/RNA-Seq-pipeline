#!/usr/bin/env python

## Mapping script

import fnmatch, os, sys, argparse, string, re
from subprocess import call
from sending_qsubjobs import submit_qsub_jobs

# at some point look at http://pymotw.com/2/optparse/index.html#module-optparse to organise the options
parser = argparse.ArgumentParser(description='This script will map the fastq files. The output directory and the build are necessary arguments. By default the directory with the fastq files is assumed to be a subfolder named "fastq_files" in the main directory. If it is different than that, it needs to be specified. It is recommended to have created a fastq directory (within the main directory) with symbolic links to the downloaded fastq files. Also it is easier if the names have been changed to a more informative name rather than the project number. If the files are paired-end the _1.fastq and _2.fastq suffixes (or _1/2.fastq.gz) are expected.')

parser.add_argument('-o', help='Required. Provide the root directory of the project where the files with the output will be created.', required=True)
parser.add_argument('-fastq', help='Provide the directory where the fastq files are.')
parser.add_argument('-build', choices=['hg38', 'hg19', 'hg18', 'dm3', 'mm9','mm10'], help='Required. Build to be used.', required=True)
parser.add_argument('-map',  default='star', choices=['star', 'tophat'], help='Aligner to be used. By default star.')
parser.add_argument('-sanger', action='store_const', default='phred', const='sanger', help='Flag for sanger score.')
parser.add_argument('-s', action='store_const', default='paired', const='single', help='Flag for single-end data.')
parser.add_argument('-strandedness', choices=['yes', 'no', 'reverse'], default='yes', help='Strandedness (yes/reverse/no). Only necessary for tophat. By default yes.')
parser.add_argument('-withcuff', action='store_true', help='Flag to be used if cufflinks will be run after STAR. Flag will set an extra argument in the mapping, which is necessary when running cufflinks.')
parser.add_argument('-fastqnozip', action='store_const', default='yes', const='no', help='Flag for not zipped fastq files.')
parser.add_argument('-user', default='erepapi', help='User to which the notifications will be sent (used when sending the jobs to the queue).')

args = parser.parse_args()

### now setting the rest of the initial values and directories

main_dir=args.o
if not main_dir.endswith('/'):
    main_dir=main_dir+'/'

if args.build=='hg19':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
    indeces_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STAR'
    fasta_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='hg18':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf'
    indeces_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/STAR'
    fasta_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa'    
elif args.build=='hg38':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
    indeces_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/STAR'
    fasta_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'    
elif args.build=='mm9':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf'
    indeces_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR'
    fasta_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='mm10':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
    indeces_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/STAR'
    fasta_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/mm10_sorted.fasta'
elif args.build=='dm3':
    gtf_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf'
    indeces_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/STAR'
    fasta_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa'

if args.fastq==None:
    fastq_dir= main_dir + 'fastq_files/'
else:
    fastq_dir=args.fastq
    if not fastq_dir.endswith('/'):
        fastq_dir=fastq_dir+'/'
  
log_dir = main_dir + 'logs/'

if args.map=='tophat':
    mapping_dir= main_dir + 'tophat_aligned/'
    aligned_dir= main_dir + 'tophat_links/' # only different for tophat
else:
    mapping_dir= main_dir + 'star_aligned/'
    aligned_dir= main_dir + 'star_aligned/' # only different for tophat

pairedorsingle=args.s
mapper=args.map 
encodedscores=args.sanger 
isitstranded=args.strandedness 
withcuff=args.withcuff
arefastqgz=args.fastqnozip
user=args.user

##########################
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if not os.path.exists(mapping_dir):
    os.makedirs(mapping_dir)

############################
### PART1: align fastq with tophat/star

# To check encoding do:
# fastq_scores.pl -i input.fastq

if encodedscores=='phred':
    solexa='solexa1.3-quals'
elif encodedscores=='sanger':
    solexa='solexa-quals'

if arefastqgz=='yes':
    fastq_suff= '.fastq.gz'
    star_option=' --readFilesCommand zcat'
else:
    fastq_suff= '.fastq'
    star_option=''
    
if pairedorsingle=='single':
    temp='*'
elif pairedorsingle=='paired':
    temp='*_1'
    
if isitstranded=='yes':
    tophat_strand=' --library-type fr-secondstrand'
elif isitstranded=='reverse':
    tophat_strand=' --library-type fr-firststrand'
else:
    tophat_strand=''

if withcuff:
    star_cuff=' --outSAMstrandField intronMotif'
else:
    star_cuff=''

os.chdir(fastq_dir)  
my_fastq = []

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, temp + fastq_suff):
        my_fastq.append(file)    

if not my_fastq:
    print 'No fastq files found in ' + fastq_dir
    print 'Check if the -s and -fastqnozip flags are set properly.'
    sys.exit()

for filename in my_fastq:
    if pairedorsingle=='single':
        sample=filename[0:len(filename)-len(fastq_suff)]
        #print sample
        if mapper=='tophat':
            command='/package/tophat/2.0.8b/bin/tophat --' + solexa + tophat_strand +' -o ' + mapping_dir + sample + ' --GTF ' + gtf_file + ' ' + indeces_file + ' ' + fastq_dir + filename 
        elif mapper=='star':
            command='module load rna-star \n' \
            + 'STAR --genomeDir ' + star_genome_dir + ' --outSAMtype BAM SortedByCoordinate --readFilesIn ' + fastq_dir + filename + star_option + star_cuff + ' --outFileNamePrefix ' + mapping_dir + sample + ' --runThreadN 4'
            
    elif pairedorsingle=='paired':
        sample=filename[0:len(filename)-len(fastq_suff)-2]
        fastq1=sample + '_1' + fastq_suff
        fastq2=sample + '_2' + fastq_suff        
        #print sample
        if mapper=='tophat':
            command='/package/tophat/2.0.8b/bin/tophat --' + solexa + tophat_strand +' -o ' + mapping_dir + sample + ' --GTF ' + gtf_file + ' ' + indeces_file + ' ' + fastq_dir +  fastq1 + ' ' + fastq_dir + fastq2 
        elif mapper=='star':
            command='module load rna-star \n' \
            + 'STAR --genomeDir ' + star_genome_dir + ' --outSAMtype BAM SortedByCoordinate --readFilesIn ' + fastq_dir + fastq1 + ' ' + fastq_dir + fastq2 + star_option + star_cuff + ' --outFileNamePrefix ' + mapping_dir + sample + ' --runThreadN 4'

    print command
    submit_qsub_jobs(command, nameqsub='qsubjob_' + mapper+ '_' + sample, my_dir=log_dir, namejob='mapping_' + mapper+ '_' + sample, logfile= sample + '.' + mapper + 'output.$JOB_ID', errfile= sample + '.' + mapper + 'error.$JOB_ID', user=user)        
    
