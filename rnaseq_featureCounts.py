#!/usr/bin/env python

## featureCounts script

import fnmatch, os, sys, argparse, string, re
from subprocess import call
from sending_qsubjobs import submit_qsub_jobs

# at some point look at http://pymotw.com/2/optparse/index.html#module-optparse to organise the options
parser = argparse.ArgumentParser(description='This script will provide the counts table to be used in the analysis. The main directory is necessary argument. The directory with the filtered data can be specified. If not, it is expected to be called filtered_bams (as the output of the rnaseq_filtering.py).')

parser.add_argument('-o', help='Required. Provide the root directory of the project where the files with the output will be created.', required=True)
parser.add_argument('-build', choices=['hg38', 'hg19', 'hg18', 'dm3', 'mm9','mm10'], help='Required. Build to be used.', required=True)
parser.add_argument('-filtereddir', help='Provide the directory where the (filtered) bams are in. If not specified, the files are expected to be in bamfiles/filtered_bams_bigwigs')
parser.add_argument('-s', action='store_const', default='paired', const='single', help='Flag for single-end data.')
parser.add_argument('-stranded', choices=['yes', 'no', 'reverse'], default='no', help='Is the sequencing stranded (yes/no/reverse). By default no.')
parser.add_argument('-t', default='5', help='Number of threads. 5 by default.')
parser.add_argument('-user', default='erepapi', help='User to whom the notifications will be sent (used when sending the jobs to the queue).')

args = parser.parse_args()

### now setting the rest of the initial values and directories

main_dir=args.o
if not main_dir.endswith('/'):
    main_dir=main_dir+'/'
  
if args.build=='hg19':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
elif args.build=='hg18':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf'
elif args.build=='hg38':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
elif args.build=='mm9':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf'
elif args.build=='mm10':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
elif args.build=='dm3':
    gtf_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf'

log_dir = main_dir + 'logs/'

if args.filtereddir==None:
    filtered_bams_dir = main_dir + 'bamfiles/filtered_bams_bigwigs/'
else:
    filtered_bams_dir= args.filtereddir
    if not filtered_bams_dir.endswith('/'):
        filtered_bams_dir=filtered_bams_dir+'/'

pairedorsingle=args.s
t=args.t
isitstranded=args.stranded 
user=args.user

##########################
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if not os.path.exists(filtered_bams_dir) or not os.listdir(filtered_bams_dir):
    print 'Cannot find the filtered files. Is the filtering directory there?'

############################

os.chdir(filtered_bams_dir)
my_bams_f = []

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.bam'):
        my_bams_f.append(file) 
      
if pairedorsingle=='single':
    paired_flag=''
else:
    paired_flag='-p '    

isitstranded=isitstranded.lower()

if isitstranded=='no':
    strandSpecific_flag=0
elif isitstranded=='reverse':
    strandSpecific_flag=2
elif isitstranded=='yes':
    strandSpecific_flag=1
else:
    print 'Are the files stranded? Specify the flag correctly (yes/no/reverse).'
    
f_bams_dir= ' ' + filtered_bams_dir

command='module load subread/1.5.0-p2 \n' \
        + 'featureCounts -T ' + t + ' -a ' + gtf_file + ' -t exon -g gene_id -o ../features_counted ' + paired_flag + ' -s ' + str(strandSpecific_flag) + f_bams_dir + f_bams_dir.join(my_bams_f)
submit_qsub_jobs(command, my_dir=log_dir, nameqsub='qsubjob_featureCounts', namejob='featureCounts', logfile= 'featureCountstimes.$JOB_ID', errfile= 'featureCountsoutput.$JOB_ID', user=user)
