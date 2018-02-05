#!/usr/bin/env python

## Filtering script

import fnmatch, os, sys, argparse, string, re
from subprocess import call
from sending_qsubjobs import submit_qsub_jobs

# at some point look at http://pymotw.com/2/optparse/index.html#module-optparse to organise the options
parser = argparse.ArgumentParser(description='This script will filter the mapped files (to follow rnaseq_mapping.py). The main directory is necessary argument (where the filtered_bams subfolder will be created). The directory with the mapped data can be specified. If not, it is expected to be called star_aligned or tophat_aligned (as the output of the rnaseq_mapping.py).')

parser.add_argument('-o', help='Required. Provide the root directory of the project where the files with the output will be created.', required=True)
parser.add_argument('-map', choices=['star', 'tophat'], help='Required. Aligner used for the mapping (tophat files are expected to be named "accepted_hits.bam" and for star to have the Aligned.sortedByCoord.out.bam suffix - otherwise it needs to be set).', required=True)
parser.add_argument('-mapdir', help='Provide the directory where the mapped files are in. If not specified, the directory is expected to be called star_aligned or tophat_aligned under the root directory specified above.')
parser.add_argument('-suffixbam', default= 'Aligned.out.bam', help='Suffix of mapped files (to be removed). Tophat files are expected to be named "accepted_hits.bam" and for star to have the "Aligned.out.bam" suffix). Otherwise it needs to be specified. (Default is "Aligned.out.bam")')
parser.add_argument('-s', action='store_const', default='paired', const='single', help='Flag for single-end data.')
parser.add_argument('-rmdup', action='store_const', default=False, const=True, help='Remove PCR duplicates (True/False). By default False.')
parser.add_argument('-rmmultimap', action='store_const', default=False, const=True, help='Remove multimapping (True/False). By default False.')
parser.add_argument('-user', default='erepapi', help='User to whom the notifications will be sent (used when sending the jobs to the queue).')

args = parser.parse_args()

### now setting the rest of the initial values and directories

main_dir=args.o
if not main_dir.endswith('/'):
    main_dir=main_dir+'/'
  
log_dir = main_dir + 'logs/'

if args.mapdir==None:
    if args.map=='tophat':
        mapping_dir= main_dir + 'tophat_aligned/'
        aligned_dir= main_dir + 'tophat_links/' # only different for tophat
        filtered_bams_dir = main_dir + 'bamfiles_bigwigs/filtered_bams_bigwigs/'
    else:
        mapping_dir= main_dir + 'star_aligned/'
        aligned_dir= main_dir + 'star_aligned/' # only different for tophat
        filtered_bams_dir = main_dir + 'bamfiles_bigwigs/filtered_bams_bigwigs/'
else:
    mapping_dir= args.mapdir
    if not mapping_dir.endswith('/'):
        mapping_dir=mapping_dir+'/'
    if args.map=='tophat':
        aligned_dir= main_dir + 'tophat_links/'
    else:
        aligned_dir= mapping_dir
    filtered_bams_dir = main_dir + 'bamfiles_bigwigs/filtered_bams_bigwigs/'

pairedorsingle=args.s
mapper=args.map 
aligneddatasuffix=args.suffixbam 
rmdup=args.rmdup
rmmultimap=args.rmmultimap
user=args.user

##########################
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if not os.path.exists(mapping_dir):
    print 'Cannot find the mapped files. Is the mapping directory there?'

if not os.path.exists(aligned_dir):
    os.makedirs(aligned_dir)

if not os.path.exists(main_dir + 'bamfiles_bigwigs/'):
    os.makedirs(main_dir + 'bamfiles_bigwigs/')

if not os.path.exists(filtered_bams_dir):
    os.makedirs(filtered_bams_dir)

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

### Part2b: filter the bams

if mapper=='tophat':
    aligneddatasuffix='.bam'

suffixaligned='*' + aligneddatasuffix

samtoolsflag_f='-b'

# not removing duplicates in single-end data is prefered because it biases-reduces the numbers of small and highly expressed genes
# doing the filtering here and not in featureCounts, since the samtools sort needs bams anyway

os.chdir(aligned_dir)
my_aligned = []
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, suffixaligned):
        my_aligned.append(file)

if pairedorsingle=='single':
    rmdupflag='-s '
else:
    rmdupflag=''
        
for bam in my_aligned:
    bam=bam[0:len(bam)-len(aligneddatasuffix)]
    command_part1='module load samtools/0.1.19 \n' 
    if rmmultimap==True:
        command_rmmultimap='samtools view -q 4 -@ 4 ' + samtoolsflag_f + ' -o '+ filtered_bams_dir + bam + '_multimap_filtered.bam ' + aligned_dir + bam + aligneddatasuffix + ' \n'
        multimapbam = '_multimap_filtered.bam '
        if rmdup==True:
            rmthemultimap = 'rm ' + filtered_bams_dir + bam + '_multimap_filtered.bam'
        else:
            rmthemultimap = ''
    else:
        command_rmmultimap=''
        multimapbam = aligneddatasuffix
        rmthemultimap = ''
        
    if rmdup==True:
        if rmmultimap==True:
            command_rmdup='samtools sort ' + filtered_bams_dir + bam + multimapbam  + ' ' + filtered_bams_dir + bam + '_sorted_filtered\n' \
            + 'samtools rmdup ' + rmdupflag + filtered_bams_dir + bam + '_sorted_filtered.bam'  + ' ' + filtered_bams_dir + bam + '_final_filtered.bam \n' \
            + 'rm ' + filtered_bams_dir + bam + '_sorted_filtered.bam \n'
        else:
            command_rmdup='samtools sort ' + aligned_dir + bam + multimapbam  + ' ' + filtered_bams_dir + bam + '_sorted_filtered\n' \
            + 'samtools rmdup ' + rmdupflag + filtered_bams_dir + bam + '_sorted_filtered.bam'  + ' ' + filtered_bams_dir + bam + '_final_filtered.bam \n' \
            + 'rm ' + filtered_bams_dir + bam + '_sorted_filtered.bam \n'
    else:
        command_rmdup=''
    command=command_part1 + command_rmmultimap + command_rmdup + rmthemultimap

    submit_qsub_jobs(command, my_dir=log_dir, nameqsub='qsubjob_filter_' + bam, namejob='filtering_' + bam, logfile= bam +'.filteringoutput.$JOB_ID', errfile= bam +'.filteringerror.$JOB_ID', user=user)
    print command
