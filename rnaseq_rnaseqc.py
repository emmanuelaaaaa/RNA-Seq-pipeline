#!/usr/bin/env python

## Filtering script

import fnmatch, os, sys, argparse, string, re
from subprocess import call
from sending_qsubjobs import submit_qsub_jobs

# at some point look at http://pymotw.com/2/optparse/index.html#module-optparse to organise the options
parser = argparse.ArgumentParser(description='This script will do the rnaseqc on the mapped files (to follow rnaseq_mapping.py). The main directory is necessary argument (where the rnaseqc subfolder will be created). The directory with the mapped data can be specified. If not, it is expected to be called star_aligned or tophat_aligned (as the output of the rnaseq_mapping.py).')

parser.add_argument('-o', help='Required. Provide the root directory of the project where the files with the output will be created.', required=True)
parser.add_argument('-map', choices=['star', 'tophat'], help='Required. Aligner used for the mapping (tophat files are expected to be named "accepted_hits.bam" and for star to have the .bam suffix).', required=True)
parser.add_argument('-build', choices=['hg38','GRCh38','hg19', 'hg18', 'dm3', 'mm9','mm10','rn5','mm9_oldind'], help='Required. Build to be used.', required=True)
parser.add_argument('-mapdir', help='Provide the directory where the mapped files are in. If not specified, the directory is expected to be called star_aligned or tophat_aligned.')
parser.add_argument('-suffixbam', default= 'Aligned.out.bam', help='Suffix of mapped files (to be removed). Tophat files are expected to be named "accepted_hits.bam" and for star to have the "Aligned.out.bam" suffix). Otherwise it needs to be specified. (Default is "Aligned.out.bam")')
#parser.add_argument('-s', action='store_const', default='paired', const='single', help='Flag for single-end data.')
parser.add_argument('-group', default='Illumina', help='Group necessary for the rnaseqc. By default Illumina.')
parser.add_argument('-user', default='erepapi', help='User to whom the notifications will be sent (used when sending the jobs to the queue).')

args = parser.parse_args()

### now setting the rest of the initial values and directories


main_dir=args.o
if not main_dir.endswith('/'):
    main_dir=main_dir+'/'
  
if args.build=='hg19':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='GRCh38':
    gtf_file='/databank/igenomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Homo_sapiens.GRCh38.86.gtf'
    fasta_file='/databank/igenomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fasta'
elif args.build=='hg18':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='hg38':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='mm9':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='mm10':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/mm10_sorted.fasta'
elif args.build=='dm3':
    gtf_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='rn5':
    gtf_file='/databank/igenomes/Rattus_norvegicus/UCSC/rn5/Annotation/Genes/genes.gtf'
    fasta_file='/databank/igenomes/Rattus_norvegicus/UCSC/rn5/Sequence/WholeGenomeFasta/genome.fa'
elif args.build=='mm9_oldind':
    gtf_file='/databank/raw/gtf/genes_ucsc_mm9.gtf'
    fasta_file='/databank/raw/mm9_all/mm9_all.fa'

log_dir = main_dir + 'logs/'

if args.mapdir==None:
    if args.map=='tophat':
        mapping_dir= main_dir + 'tophat_aligned/'
        aligned_dir= main_dir + 'tophat_links/' # only different for tophat
    else:
        mapping_dir= main_dir + 'star_aligned/'
        aligned_dir= main_dir + 'star_aligned/' # only different for tophat
else:
    mapping_dir= args.mapdir
    if not mapping_dir.endswith('/'):
        mapping_dir=mapping_dir+'/'
    if args.map=='tophat':
        aligned_dir= main_dir + 'tophat_links/'
    else:
        aligned_dir= mapping_dir

rnaseqc_dir = main_dir + 'rnaseqc/'

#pairedorsingle=args.s
mapper=args.map 
aligneddatasuffix=args.suffixbam 
group_platform=args.group
user=args.user

##########################
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if not os.path.exists(mapping_dir) or not os.listdir(mapping_dir):
    print 'Cannot find the mapped files. Is the mapping directory there?'

if not os.path.exists(aligned_dir):
    os.makedirs(aligned_dir)

if not os.path.exists(rnaseqc_dir):
    os.makedirs(rnaseqc_dir)

############################
### Part2a: make symbolic links of each tophat-mapped bam/sam to a single dir

if mapper=='tophat' and not os.listdir(aligned_dir):
    os.chdir(mapping_dir)
    my_fastq = []
    for file in os.listdir('.'):
        my_fastq.append(file)       

    for filename in my_fastq:
        sample=filename
        command='ln -s ' + mapping_dir + filename + '/accepted_hits.bam ' + aligned_dir + sample +'.bam'
        call(command, shell=True) 

### Part2d: do the rna-seqc 

if mapper=='tophat':
    aligneddatasuffix='.bam'

suffixaligned='*' + aligneddatasuffix

os.chdir(aligned_dir)
my_aligned = []
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, suffixaligned):
        my_aligned.append(file)

#if pairedorsingle=='single':
#    rnaseqc_flag='-singleEnd '
#else:
#    rnaseqc_flag=''

main_command='module load picard-tools/1.105 \n' \
        + 'module load samtools/0.1.19 \n' \
        + 'module unload java \n' \
        + 'module load java/7u80 \n' \
        + '#RNA-SeQC v1.1.8.1 07/11/14 \n' \

samplefile=open(main_dir + 'samplefile_rnaseqc.txt','w')
samplefile.write('Sample ID\tBam File\tNotes\n')

for bam in my_aligned:
    bam=bam[0:len(bam)-len(aligneddatasuffix)]
    command_part2='AddOrReplaceReadGroups INPUT=' + aligned_dir + bam + aligneddatasuffix + ' OUTPUT=' + rnaseqc_dir + bam + '_h.bam SO=coordinate RGID=' + bam + ' RGLB=RNA_seq_lib RGPL=' + group_platform + ' RGSM=' + bam + ' RGPU=run1 \n' \
    + 'MarkDuplicates INPUT=' + rnaseqc_dir + bam + '_h.bam OUTPUT=' + rnaseqc_dir + bam + '_h_md.bam ASSUME_SORTED=true METRICS_FILE=' + rnaseqc_dir + bam + '_duplicates \n' \
    + 'ReorderSam I=' + rnaseqc_dir + bam + '_h_md.bam O=' + rnaseqc_dir + bam + '_h_md_s.bam R=' + fasta_file + ' \n' \
    + 'samtools index ' + rnaseqc_dir + bam + '_h_md_s.bam \n' \
    + 'rm ' + rnaseqc_dir + bam + '_h_md.bam ' + rnaseqc_dir + bam + '_h.bam \n'
    main_command=main_command + command_part2
    samplefile.write(bam + '\t' + rnaseqc_dir + bam + '_h_md_s.bam\tNA\n')

samplefile.close()

command_part3='rna-seqc -o ' + rnaseqc_dir + 'rnaseqc_output -r '+ fasta_file +' -t ' + gtf_file + ' -s '+ main_dir + 'samplefile_rnaseqc.txt'
main_command=main_command + command_part3
print main_command

submit_qsub_jobs(main_command, my_dir=log_dir, nameqsub='qsubjob_rnaseqc', namejob='rnaseqc', logfile= 'rnaseqc_output.$JOB_ID', errfile= 'rnaseqc_error.$JOB_ID', user=user)

