#!/usr/bin/env python

## Mapping script

import fnmatch, os, sys, argparse, string, re, math
from subprocess import call
from sending_qsubjobs import submit_qsub_jobs

# at some point look at http://pymotw.com/2/optparse/index.html#module-optparse to organise the options
parser = argparse.ArgumentParser(description='This script will map the fastq files. The output directory and the build are necessary arguments. By default the directory with the fastq files is assumed to be a subfolder named "fastq_files" in the main directory. If it is different than that, it needs to be specified. It is recommended to have created a fastq directory (within the main directory) with symbolic links to the downloaded fastq files. Also it is easier if the names have been changed to a more informative name rather than the project number. If the files are paired-end the _1.fastq and _2.fastq suffixes (or _1/2.fastq.gz) are expected.')

parser.add_argument('-o', help='Required. Provide the root directory of the project where the files with the output will be created.', required=True)
parser.add_argument('-fastq', help='Provide the directory where the fastq files are. Note: symlinked directories do not work at the moment.')
parser.add_argument('-build', choices=['hg38', 'hg19', 'hg18', 'dm3', 'mm9','mm10'], help='Required. Build to be used.', required=True)
parser.add_argument('-map',  default='star', choices=['star', 'tophat'], help='Aligner to be used. By default star.')
parser.add_argument('-sanger', action='store_const', default='phred', const='sanger', help='Flag for sanger score.')
parser.add_argument('-strandedness', choices=['yes', 'no', 'reverse'], default='no', help='Strandedness (yes/reverse/no). Only necessary for tophat. By default no.')
parser.add_argument('-withcuff', action='store_true', help='Flag to be used if cufflinks will be run after STAR. Flag will set an extra argument in the mapping, which is necessary when running cufflinks.')
parser.add_argument('-fastqnozip', action='store_const', default='yes', const='no', help='Flag for not zipped fastq files.')
parser.add_argument('-user', default='erepapi', help='User to which the notifications will be sent (used when sending the jobs to the queue).')
parser.add_argument('-s', action='store_true', help='Flag for single-end data.')

#extra
parser.add_argument('-NextSeq', action='store_true', help='Flag to specify files produced by NextSeq.')
parser.add_argument('-usegtf', action='store_true', help='Flag for using the gtf with mapping.')
parser.add_argument('-useENCODEflags', action='store_true', help='Flag for using the options recommended by ENCODE for mapping long RNA-Seq (mRNAs (poly-A(+)), rRNA-depleted total RNA, or poly-A(-) RNA populations that are size-selected to be longer than approximately 200 bp).')
parser.add_argument('-ERCC',  action='store_true', help='Flag for specifying mapping with the ERCC, when the ERCC spike-ins were included in the samples.')
parser.add_argument('-adapterseq', help='Provide the sequence to be used for clipping with STAR  (using the clip3pAdapterSeq flag).')

parser.add_argument('-fastqc', action='store_true', help='Flag for running fastqc')
parser.add_argument('-doonlyfastqc', action='store_true', help='Flag for sending the scripts for fastqc only and exiting without doing the mapping.')

# have changed
parser.add_argument('-fastqscreen', action='store_true', help='Flag for running fastq_screen')
parser.add_argument('-doonlyfastqscreen', action='store_true', help='Flag for sending the scripts for fastq_screen only and exiting without doing the mapping or the fastqc.')

args = parser.parse_args()

### now setting the rest of the initial values and directories

main_dir=args.o
singleend=args.s
mapper=args.map 
encodedscores=args.sanger 
isitstranded=args.strandedness 
withcuff=args.withcuff
NextSeq=args.NextSeq
arefastqgz=args.fastqnozip
user=args.user

NextSeq=args.NextSeq
usegft=args.usegtf
useENCODEflags=args.useENCODEflags
ercc=args.ERCC
runfastqc=args.fastqc
runonlyfastqc=args.doonlyfastqc
runfastqscreen=args.fastqscreen
runonlyfastqscreen=args.doonlyfastqscreen

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
fastqc_dir= main_dir + 'fastqc/'
fastqscreen_dir= main_dir + 'fastq_screen/'

if args.map=='tophat':
    mapping_dir= main_dir + 'tophat_aligned/'
    aligned_dir= main_dir + 'tophat_links/' # only different for tophat
else:
    mapping_dir= main_dir + 'star_aligned/'
    aligned_dir= main_dir + 'star_aligned/' # only different for tophat

##########################
### if the directories are missing, create them
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if (runfastqc or runonlyfastqc) and not os.path.exists(fastqc_dir):
    os.makedirs(fastqc_dir)

if (runfastqscreen or runonlyfastqscreen) and not os.path.exists(fastqscreen_dir):
    os.makedirs(fastqscreen_dir)

if not os.path.exists(mapping_dir):
    os.makedirs(mapping_dir)

if args.map=='tophat' and not os.path.exists(aligned_dir):
    os.makedirs(aligned_dir)

##########################
## setting up the necessary parameters

if args.build=='hg19':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
    gtf_file_ERCC='/databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes_with_ERCC.gtf'
    indeces_file='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STAR'
    star_genome_dir_ERCC='/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STAR_with_ERCC'
elif args.build=='hg18':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Annotation/Genes/genes.gtf'
    gtf_file_ERCC=''
    indeces_file='/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/STAR'
    star_genome_dir_ERCC=''
elif args.build=='hg38':
    gtf_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
    gtf_file_ERCC='/databank/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes_with_ERCC.gtf'
    indeces_file='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/STAR'
    star_genome_dir_ERCC='/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/STAR_with_ERCC'
elif args.build=='mm9':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf'
    gtf_file_ERCC='/databank/igenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf'
    indeces_file='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR'
    star_genome_dir_ERCC='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR_with_ERCC'
elif args.build=='mm10':
    gtf_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
    gtf_file_ERCC='/databank/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes_with_ERCC.gtf'
    indeces_file='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/STAR'
    star_genome_dir_ERCC='/databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/STAR_with_ERCC'
elif args.build=='dm3':
    gtf_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf'
    gtf_file_ERCC=''
    indeces_file='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome' # for tophat
    star_genome_dir='/databank/igenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/STAR'
    star_genome_dir_ERCC=''

if ercc:
    gtf_file = gtf_file_ERCC
    star_genome_dir = star_genome_dir_ERCC
    if gtf_file=='' or star_genome_dir=='':
        print 'The genome files for with the ERCC are missing for this genome. Update the genomes and then update the pipeline.'
        sys.exit()

if args.adapterseq==None:
    adapterflag= ''
else:
    adapterflag= ' --clip3pAdapterSeq  ' + args.adapterseq

if arefastqgz=='yes':
    fastq_pattern= 'fastq.gz'
    star_option=' --readFilesCommand zcat'
else:
    fastq_pattern= 'fastq'
    star_option=''

if encodedscores=='phred':
    solexa='solexa1.3-quals'
elif encodedscores=='sanger':
    solexa='solexa-quals'
    
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

if usegft:
    gtfflag = ' --sjdbGTFfile ' + gtf_file
    genomeLoad = ' --genomeLoad NoSharedMemory'
else:
    gtfflag = ''
    genomeLoad = ' --genomeLoad LoadAndKeep'

if useENCODEflags:
    encodeflags=' --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 '
else:
    encodeflags=''

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
### PART1: run fastqc and fastq_screen

if NextSeq:
    temp_dir = main_dir + 'temp/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

if runfastqc or runonlyfastqc or runfastqscreen or runonlyfastqscreen:
    main_command_fastqcommands='module purge \n' + 'module load fastq_screen \n' + 'module load fastqc/0.11.4 \n' \
       
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
                else:
                    mergedfiles1 = ' '.join(foundfiles[0:7:2])
                    mergedfiles2 = ' '.join(foundfiles[1:8:2])
                    extra_command = 'cat ' + mergedfiles1 + ' > ' + temp_dir + filename + '_1.' + fastq_pattern + ' \n' + 'cat ' + mergedfiles2 + ' > ' + temp_dir + filename + '_2.' + fastq_pattern + ' \n'
                    extra_command2 = 'rm ' + temp_dir + filename + '_1.' + fastq_pattern + ' ' + temp_dir + filename + '_2.' + fastq_pattern + ' \n'
                    mergedfiles1 = temp_dir + filename + '_1.' + fastq_pattern
                    mergedfiles2 = temp_dir + filename + '_2.' + fastq_pattern
            else:
                if singleend:
                    mergedfiles = foundfiles[0]
                    extra_command = ''
                    extra_command2 = ''
                else:
                    mergedfiles1 = foundfiles[0]
                    mergedfiles2 = foundfiles[1]
                    extra_command = ''
                    extra_command2 = ''
            if (runfastqc and not runfastqscreen) or (runonlyfastqc and not runonlyfastqscreen):
                if singleend:
                    command_fastqc=extra_command + 'fastqc ' + mergedfiles + ' -o ' + fastqc_dir +' --extract \n' + extra_command2
                else:
                    command_fastqc=extra_command + 'fastqc ' + mergedfiles1 + ' -o ' + fastqc_dir +' --extract \n' + 'fastqc ' + mergedfiles2 + ' -o ' + fastqc_dir +' --extract \n' + extra_command2
                main_command_fastqcommands = main_command_fastqcommands + command_fastqc

            if (runfastqscreen and not runfastqc) or (runonlyfastqscreen and not runonlyfastqc):
                if singleend:
                    command_fastqscreen=extra_command + 'fastq_screen --aligner bowtie2 ' + mergedfiles + ' --outdir ' + fastqscreen_dir + ' \n' + extra_command2
                else:
                    command_fastqscreen=extra_command + 'fastq_screen --aligner bowtie2 ' + mergedfiles1 + ' --outdir ' + fastqscreen_dir +' \n' + extra_command2
                main_command_fastqcommands = main_command_fastqcommands + command_fastqscreen

            if (runfastqc and runfastqscreen) or (runonlyfastqc and runonlyfastqscreen):
                if singleend:
                    command_fastqc=extra_command + 'fastqc ' + mergedfiles + ' -o ' + fastqc_dir +' --extract \n'
                    command_fastqscreen='fastq_screen --aligner bowtie2 ' + mergedfiles + ' --outdir ' + fastqscreen_dir + ' \n' + extra_command2
                else:
                    command_fastqc=extra_command + 'fastqc ' + mergedfiles1 + ' -o ' + fastqc_dir +' --extract \n' + 'fastqc ' + mergedfiles2 + ' -o ' + fastqc_dir +' --extract \n' 
                    command_fastqscreen='fastq_screen --aligner bowtie2 ' + mergedfiles1 + ' --outdir ' + fastqscreen_dir + ' \n' + extra_command2
                main_command_fastqcommands = main_command_fastqcommands + command_fastqc +  command_fastqscreen

        submit_qsub_jobs(main_command_fastqcommands, nameqsub='qsubjob_fastq_' + str(i), my_dir=log_dir, namejob='fastq_' + str(i), logfile= 'fastq_batch' + str(i) + 'output.$JOB_ID', errfile= 'fastq_batch' + str(i) + 'error.$JOB_ID', user=user)        
        #print main_command_fastqcommands
        main_command_fastqcommands='module purge \n' + 'module load fastq_screen \n' + 'module load fastqc/0.11.4 \n' \

if runonlyfastqc or runonlyfastqscreen:
    print 'The fastqc or fastq_screen qsub file(s) been sent. It will now exit because you have added the -doonlyfastqc and/or -doonlyfastqscreen flag. '
    sys.exit()

##########################
### PART2: map fastq

# To check encoding do:
# fastq_scores.pl -i input.fastq


# to delete when done:
#for file in os.listdir('.'):
#    if fnmatch.fnmatch(file, temp + fastq_suff):
#        my_fastq.append(file)    

if len(uniqsamples)>30:
    pipeline='SC'
else:
    pipeline='bulk'

if pipeline=='bulk':
    #mapping is done in individual qsubs   
    for filename in uniqsamples:
        foundfiles=[]
        for m in my_fastq:
            if filename in m:
                foundfiles.append(m)
        foundfiles.sort()
        if singleend:
            if NextSeq:
                allfilenames = ','.join(foundfiles)
            else:
                allfilenames = foundfiles[0] # there is only one
        else:
            if filename.endswith('_2.' + fastq_pattern) or filename.endswith('_R2_001.' + fastq_pattern):
                continue
            if NextSeq:
                allfilenames = ','.join(foundfiles[0:7:2]) + ' ' + ','.join(foundfiles[1:8:2])
            else:
                allfilenames=' '.join(foundfiles)
        if mapper=='tophat':
            command='module load tophat/2.0.13 \n' \
                + 'tophat --' + solexa + tophat_strand +' -o ' + mapping_dir + filename + ' --GTF ' + gtf_file + ' ' + indeces_file + ' ' +  allfilenames
        elif mapper=='star':
            command='module load rna-star/2.4.2a \n' \
            + 'STAR --genomeDir ' + star_genome_dir + gtfflag + encodeflags + adapterflag + ' --outSAMtype BAM Unsorted' + genomeLoad + ' --readFilesIn ' + allfilenames + star_option + star_cuff + ' --outFileNamePrefix ' + mapping_dir + filename + ' --runThreadN 4'
        submit_qsub_jobs(command, nameqsub='qsubjob_' + mapper+ '_' + filename, my_dir=log_dir, namejob='mapping_' + mapper+ '_' + filename, logfile= filename + '.' + mapper + 'output.$JOB_ID', errfile= filename + '.' + mapper + 'error.$JOB_ID', user=user)        
else:
    #mapping will be done in sets of 20
    if mapper=='tophat':
        print 'Too many files to be mapping with tophat. This script is not configured for this.'
        sys.exit()
    main_command_star = 'module load rna-star/2.4.2a \n'
    for i in range(1,toprange+1):
        for filename in uniqsamples[(i-1)*20:i*20]:
            foundfiles=[]
            for m in my_fastq:
                if filename in m:
                    foundfiles.append(m)
            foundfiles.sort()
            if singleend:
                if NextSeq:
                    allfilenames = ','.join(foundfiles)
                else:
                    allfilenames = foundfiles[0] # there is only one
            else:
                if filename.endswith('_2.' + fastq_pattern) or filename.endswith('_R2_001.' + fastq_pattern):
                    continue
                if NextSeq:
                    allfilenames = ','.join(foundfiles[0:7:2]) + ' ' + ','.join(foundfiles[1:8:2])
                else:
                    allfilenames=' '.join(foundfiles)
            command_star='STAR --genomeDir ' + star_genome_dir + gtfflag + encodeflags + adapterflag  + ' --outSAMtype BAM Unsorted' + genomeLoad + ' --readFilesIn ' + allfilenames + star_option +  ' --outFileNamePrefix ' + mapping_dir + filename + ' --runThreadN 2 \n' 
            main_command_star = main_command_star + command_star
                
        removegenome = 'STAR --genomeDir ' + star_genome_dir + ' --genomeLoad Remove \n'
        main_command_star = main_command_star + removegenome
        submit_qsub_jobs(main_command_star, nameqsub='qsubjob_star_batch' + str(i), my_dir=log_dir, namejob='mapping_' + str(i), logfile= 'mapping_batch' + str(i) + 'output.$JOB_ID', errfile= 'mapping_batch' + str(i) + 'error.$JOB_ID', user=user)        
        #print main_command_star
        main_command_star='module load rna-star \n' \
