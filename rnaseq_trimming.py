#!/usr/bin/env python

## Mapping script

import fnmatch, os, sys, argparse, string, re, math
from subprocess import call
from sending_qsubjobs import submit_qsub_jobs

# at some point look at http://pymotw.com/2/optparse/index.html#module-optparse to organise the options
parser = argparse.ArgumentParser(description='This script will map the fastq files. The output directory and the build are necessary arguments. By default the directory with the fastq files is assumed to be a subfolder named "fastq_files" in the main directory. If it is different than that, it needs to be specified. It is recommended to have created a fastq directory (within the main directory) with symbolic links to the downloaded fastq files. Also it is easier if the names have been changed to a more informative name rather than the project number. If the files are paired-end the _1.fastq and _2.fastq suffixes (or _1/2.fastq.gz) are expected.')

parser.add_argument('-o', help='Required. Provide the root directory of the project where the files with the output will be created.', required=True)
parser.add_argument('-fastq', help='Provide the directory where the fastq files are.')
parser.add_argument('-adapter', choices=['truseq','nextera', 'sequence'], default='truseq', help='Trim for the TruSeq adapter, the Nextera sequence (CTGTCTCTTATA) or provide your own sequence, in which case the sequence_3prime and/or sequence_5prime will need to be set.')
parser.add_argument('-q', help='Trim low-quality ends from reads before adapter removal.')
parser.add_argument('-sequence_3prime', help='Provide the sequence for the 3prime trimming.')
parser.add_argument('-sequence_5prime', help='Provide the sequence for the 5prime trimming')
parser.add_argument('-user', default='erepapi', help='User to which the notifications will be sent (used when sending the jobs to the queue).')
parser.add_argument('-s', action='store_true', help='Flag for single-end data.')
parser.add_argument('-NextSeq', action='store_true', help='Flag to specify files produced by NextSeq.')
parser.add_argument('-fastqnozip', action='store_const', default='yes', const='no', help='Flag for not zipped fastq files.')

args = parser.parse_args()

### now setting the rest of the initial values and directories

main_dir=args.o
singleend=args.s
user=args.user
NextSeq=args.NextSeq
adapter=args.adapter
qual=args.q
seq1=args.sequence_3prime
seq2=args.sequence_5prime
arefastqgz=args.fastqnozip

if arefastqgz=='yes':
    fastq_pattern= 'fastq.gz'
    star_option=' --readFilesCommand zcat'
else:
    fastq_pattern= 'fastq'
    star_option=''

##########################
### setting up the directories
if not main_dir.endswith('/'):
    main_dir=main_dir+'/'

if args.fastq==None:
    fastq_dir= main_dir + 'fastq_files/'
else:
    fastq_dir=args.fastq
    if not fastq_dir.endswith('/'):
        fastq_dir=fastq_dir+'/'
  
log_dir = main_dir + 'logs/'
cutadapt_dir= main_dir + 'cutadapt/'

##########################
### if the directories are missing, create them
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if not os.path.exists(cutadapt_dir):
    os.makedirs(cutadapt_dir)

##########################
## setting up the necessary parameters

if adapter=='truseq' and singleend:
    adapterflag= '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
elif adapter=='truseq':
    adapterflag= '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

if adapter=='nextera' and singleend:
    adapterflag= '-a CTGTCTCTTATA'
elif adapter=='nextera':
    adapterflag= '-b CTGTCTCTTATA '

if adapter=='sequence' and seq1==None and seq2==None:
    print 'No sequences are provided for the trimming. Please check your parameters.'
    sys.exit()
elif adapter=='sequence' and singleend and seq2==None:
    adapterflag= ' -a ' + seq1
elif adapter=='sequence' and singleend and seq1==None:
    adapterflag= ' -g ' + seq2
elif adapter=='sequence' and singleend and seq1!=None and seq2!=None:
    print 'Please check your flags, you have provided both the single end strand two sequences for trimming. This does not look correct.'
    sys.exit()
elif adapter=='sequence' and seq1!=None and seq2==None:
    adapterflag= ' -a ' + seq1 + ' -A ' + seq1
elif adapter=='sequence' and seq1!=None and seq2!=None:
    adapterflag= ' -a ' + seq1 + ' -G ' + seq2
    print 'Please check in the logs (and the cutadapt manual) that this is doing what you expect it to be doing!'
elif adapter=='sequence' and seq1==None and seq2!=None:
    adapterflag= ' -g ' + seq2 + ' -G ' + seq2
    print 'Please check in the logs (and the cutadapt manual) that this is doing what you expect it to be doing!'

if qual==None:
    qual_param=''
else:
    qual_param=' -q ' + qual

##########################
### get the fastq files

os.chdir(fastq_dir)  
my_fastq = []
samplefiles = []

for root, dirs, files in os.walk(fastq_dir, topdown=True):
    for name in files:
        if name.endswith(fastq_pattern):
            my_fastq.append(os.path.join(root, name))
            samplefiles.append(name)

num_files = len(my_fastq)

if not my_fastq:
    print 'No fastq files found in ' + fastq_dir
    print 'Check if the -s and -fastqnozip flags are set properly.'
    sys.exit()

listset = []
for i in samplefiles:
    if NextSeq:
        listset.append(i[0:len(i)-len('_L001_R1_001.' + fastq_pattern)])
    else:
        listset.append(i[0:len(i)-len('_1.' + fastq_pattern)])

uniqsamples = list(set(listset)) # to get the unique samples for paired-end or Next-Seq 
toprange = int(math.ceil(float(len(uniqsamples))/20)) # to get how many sets of 20 samples (the number of qsub to send)

#checking for problems in naming scheme
for i in uniqsamples:
    mycount=0
    wierdsamples = []
    for j in uniqsamples:
        if i in j:
            mycount += 1
            wierdsamples.append(j)
            if mycount>1:
                print 'Naming scheme problematic, sample "' + i + '" within samples ' 
                print wierdsamples
                print 'Check that the mapping and fastqc has been done properly (probably not). You may need to remap that sample manually.' 

############################
### check names - appropriate flags

nextseqpattern=re.search( r'L00[1-4]_R', samplefiles[0])
if nextseqpattern and not NextSeq:
    print('The files look like they may be coming from the NextSeq machine. You may want to add the -NextSeq flag.')
    sys.exit()

all_samplepattern1 = []
all_samplepattern2 = []
for mystring in samplefiles:
    samplepattern1=re.search( r'_2.f(ast)?q', mystring)
    samplepattern2=re.search( r'_R2_001.fastq', mystring)
    if (samplepattern1 or samplepattern2) and singleend:
        print('The files look like they are coming from paired-end data. You may want to remove the -s flag.')
        sys.exit()
    all_samplepattern1.append(samplepattern1)
    all_samplepattern2.append(samplepattern2)

if ((not any(all_samplepattern1)) and (not any(all_samplepattern2)) and (not singleend)):
    print('The files look like they are coming from single-end data. You may want to add the -s flag.')
    sys.exit()

############################
### cutadapt:

if NextSeq:
    temp_dir = main_dir + 'temp/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

main_command='module load cutadapt/1.8.3 \n' \
       
for i in range(1,toprange+1):
    for filename in uniqsamples[(i-1)*20:i*20]:
        foundfiles=[]
        for m in my_fastq:
            if filename in m:
                foundfiles.append(m)
        foundfiles.sort()
        if NextSeq:
            if singleend:
                mergedfiles = ' '.join(foundfiles)
                extra_command = 'cat ' + mergedfiles + ' > ' + temp_dir + filename + '.' + fastq_pattern + ' \n'
                extra_command2 = 'rm ' + temp_dir + filename + '.' + fastq_pattern + ' \n'
                mergedfiles =  temp_dir + filename + '.' + fastq_pattern
                outputfiles =  cutadapt_dir + filename + '.' + fastq_pattern
            else:
                mergedfiles1 = ' '.join(foundfiles[0:7:2])
                mergedfiles2 = ' '.join(foundfiles[1:8:2])
                extra_command = 'cat ' + mergedfiles1 + ' > ' + temp_dir + filename + '_1.' + fastq_pattern + ' \n' + 'cat ' + mergedfiles2 + ' > ' + temp_dir + filename + '_2.' + fastq_pattern + ' \n'
                extra_command2 = 'rm ' + temp_dir + filename + '_1.' + fastq_pattern + ' ' + temp_dir + filename + '_2.' + fastq_pattern + ' \n'
                mergedfiles1 = temp_dir + filename + '_1.' + fastq_pattern
                mergedfiles2 = temp_dir + filename + '_2.' + fastq_pattern
                outputfiles1 =  cutadapt_dir + filename + '_1.' + fastq_pattern
                outputfiles2 =  cutadapt_dir + filename + '_2.' + fastq_pattern
        else:
            if singleend:
                mergedfiles = foundfiles[0]
                extra_command = ''
                extra_command2 = ''
                outputfiles =  cutadapt_dir + filename + '.' + fastq_pattern
            else:
                mergedfiles1 = foundfiles[0]
                mergedfiles2 = foundfiles[1]
                extra_command = ''
                extra_command2 = ''
                outputfiles1 =  cutadapt_dir + filename + '_1.' + fastq_pattern
                outputfiles2 =  cutadapt_dir + filename + '_2.' + fastq_pattern
        if singleend:
            command_cutadapt=extra_command + 'cutadapt -m 10 ' + adapterflag + qual_param + ' -o ' +  outputfiles + ' ' + mergedfiles +' \n' + extra_command2
        else:
            command_cutadapt=extra_command + 'cutadapt -m 10 ' + adapterflag + qual_param + ' -o ' + outputfiles1 + ' -p ' + outputfiles2 + ' ' + mergedfiles1 + ' ' + mergedfiles2 + ' \n' + extra_command2
        main_command = main_command + command_cutadapt

    submit_qsub_jobs(main_command, nameqsub='qsubjob_trimming_' + str(i), my_dir=log_dir, namejob='cutadapt_' + str(i), logfile= 'cutadapt_batch' + str(i) + 'output.$JOB_ID', errfile= 'cutadapt_batch' + str(i) + 'error.$JOB_ID', user=user)        
    #print main_command
    main_command='module load cutadapt/1.8.3 \n' \
