import os
from subprocess import call

TEMPLATE_SERIAL = """
#!/bin/sh
#####################################
#$ -cwd
#$ -q batchq 
#$ -M {user}
#$ -m eas
#$ -N {name}
#$ -e {errfile}
#$ -o {logfile}
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""
 
def submit_qsub_jobs(code, my_dir=os.getcwd() , nameqsub='qsubjob', namejob="job", logfile="output.$JOB_ID", errfile="error.$JOB_ID", cleanup=False, user='erepapi'):
    os.chdir(my_dir)
    open(nameqsub + '.qsub', 'wb').write(TEMPLATE_SERIAL.format(script=code, name=namejob, logfile=logfile, errfile=errfile, user=user))
    try:
        call('qsub ' + nameqsub + '.qsub', shell=True)
    finally:
        if cleanup:
            os.remove(base + '.qsub')
