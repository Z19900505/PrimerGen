#!/usr/bin/python
# Samuel Moijueh
# July 2016

import sys
import re
import time
import os
import datetime
from optparse import OptionParser
import errno
from array import array
import subprocess

sys.path.append('/usr/lib64/python2.7/')


###########################################################
'''
== SYNOPOSIS ==
This Program performs a series of steps for running the final script in parallel on a HPC:

1) initialize new directories for the final script of the pipeline
2) Split the bed file file from the previous script into smaller bed files no greater than 2500 lines
3) create Bash job scripts for each bed file.
    These Bash scripts will be run later as jobs on the Sun Grid Engine using PBS
4) create a Bash qsub submitter file

Running the jobs in parallel is a great way for automating a large number of (repetitve) tasks.
This kind of multiprocessing ensures that server resources are allocated efficiently.

USAGE: python 2_PrepareJobs.py

'''
#########################################################


# global variables
now = datetime.datetime.now()
base_directory = "/group/das-lab/Samuel/Jag_PrimerDesign_Project/" + str(now.month) + "-" + str(now.day) + "_InSilico_PrimerDesign"


def main():
    start_time = time.time()
    initialize_pipeline_directories()
    Alamut_ExonSequences_file = retrieve_alamut_sequences_file()
    splitAlamutFile(Alamut_ExonSequences_file)
    create_job_scripts()
    # missing_transcripts_fh.close();
    # create_qsub_submitter();
    end_time = time.time()
    print "Total time = " + str(end_time - start_time) + "\n"


# creates the directories required for the pipeline
def initialize_pipeline_directories():
    directories = [base_directory + "/PrimerGen/input/AlamutFiles/",
                   base_directory + "/Job_Scripts/",
                   base_directory + "/PrimerGen/output/JobLogs/",
                   base_directory + "/PrimerGen/output/Processed_Primers/",
                   base_directory + "/PrimerGen/input/AlamutFiles/",
                   base_directory + "/tool_configurations/BLAT/input/",
                   base_directory + "/tool_configurations/BLAT/output/",
                   base_directory + "/tool_configurations/Primer3/input/",
                   base_directory + "/tool_configurations/Primer3/output/"]

    for directory in directories:
        dir = os.path.dirname(directory.strip())
        if not os.path.exists(dir):
            os.makedirs(dir)


# opens the Alamut file (created in the previous script) that contains the gene transcripts and their exon sequences
def retrieve_alamut_sequences_file():
    usage = "python SplitAlamutFile.py -i AlamutSeqfile.bed"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input_File", dest="input", help="The path to the alamut file.")
    (options, args) = parser.parse_args()
    alamutfile = options.input

    # USAGE
    if(len(sys.argv) < 3):
        print("Usage: python revised_PrepareJobs.py -i 2_22_NM_transcript_primers.txt");
        sys.exit()
    if(not os.path.exists(options.input)):
        print("Error: %s file doesn't exist, please check and try again" % (options.input))
        sys.exit()

    return alamutfile

# splits the Alamut bed file into smaller bed files (<2500 line-increments)
def splitAlamutFile(AlamutFile):
    AlamutFile = open(AlamutFile, "r")
    lineCounter = 0
    FileCounter = 0
    for line in AlamutFile:
        lineCounter = lineCounter + 1
        if lineCounter == 1:
            FileCounter = FileCounter + 1
            NewFile = open(base_directory + '/PrimerGen/input/AlamutFiles/AlamutFile_' + str(FileCounter) + '.bed', 'w')
        NewFile.write(line)
        if lineCounter == 2500:  # this line determines how big each bed file should be
            lineCounter = 0


# creates the bash scripts that call revised_PrimerGen.py
def create_job_scripts():
    for filename in os.listdir(base_directory + "/PrimerGen/input/AlamutFiles/"):
        if filename.startswith("AlamutFile"):
            # print filename;
            SeqFileNumber = re.search(r'\d+', filename).group()

            mkdir_recursive(base_directory + "/PrimerGen/output/JobLogs/AlamutFile" + SeqFileNumber)

            job_header = "#!/bin/bash\n\n" + \
                "#PBS -N PrimerGen_AlamutFile_" + str(SeqFileNumber) + "\n" + \
                "#PBS -S /bin/bash\n" + \
                "#PBS -l walltime=250:00:00\n" + \
                "#PBS -l mem=16GB\n" + \
                "#PBS -o " + base_directory + "/PrimerGen/output/JobLogs/AlamutFile" + SeqFileNumber + "/STDOUT.txt\n" + \
                "#PBS -e " + base_directory + "/PrimerGen/output/JobLogs/AlamutFile" + SeqFileNumber + "/BLAT_err.txt\n" + \
                "#PBS -l nodes=1:ppn=1\n\n";

            #call_to_PrimerGen = "python /group/das-lab/Samuel/Jag_PrimerDesign_Project/old_scripts/4_PrimerGen_testing_copy.py ";
            call_to_PrimerGen = "/apps/compilers/python/2.7.6/bin/python2.7 /group/das-lab/Samuel/Jag_PrimerDesign_Project/3_PrimerGen.py "
            parameters = "-p %s -i %s -t %s" % \
                (base_directory, base_directory + "/PrimerGen/input/AlamutFiles/AlamutFile_" + str(SeqFileNumber) + ".bed", \
                 #"/group/das-lab/Samuel/Jag_PrimerDesign_Project/input/2.vcf");
                    "/group/das-lab/Samuel/Jag_PrimerDesign_Project/input/1000G_Phase_3_gpos.vcf")
            BatchFile = open(base_directory + '/Job_Scripts/Job_AlamutFile_' + str(SeqFileNumber) + '.sh', 'w')
            BatchFile.write(job_header + call_to_PrimerGen + parameters)  # + " > " + base_directory + "/errors.txt")
            BatchFile.close()


def create_qsub_submitter():
    mkdir_recursive(base_directory + "/Job_Scripts/qsub_submitter/")
    count = 1
    job_first_pass = True
    external_counter = 1
    qsub_submitter_fh = open(base_directory + "/Job_Scripts/qsub_submitter/JobSet" + str(count) + ".sh", "w")
    for filename in os.listdir(base_directory + "/PrimerGen/input/AlamutFiles/"):
        if (external_counter % 100):
            qsub_submitter_fh = open(base_directory + "/Job_Scripts/qsub_submitter/JobSet" + str(count) + ".sh", "w")
            count += 1
            job_first_pass = True
        if (job_first_pass):
            job_header = "#!/bin/bash\n\n" + \
                "#PBS -N PrimerGen_qsub_summitter_" + str(count) + "\n" + \
                "#PBS -S /bin/bash\n" + \
                "#PBS -l walltime=170:00:00\n" + \
                "#PBS -l mem=16GB\n" + \
                "#PBS -o " + base_directory + "/Job_Scripts/qsub_submitter/STDOUT.txt\n" + \
                "#PBS -e " + base_directory + "/Job_Scripts/qsub_submitter/BLAT_err.txt\n" + \
                "#PBS -l nodes=1:ppn=1\n\n";
            qsub_submitter_fh.write(job_header)
            job_first_pass = False
        qsub_submitter_fh.write("qsub " + filename + "\n")
        external_counter += 1


def mkdir_recursive(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def runProcess(exe):
    p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while(True):
        retcode = p.poll()  # returns None while subprocess is running
        line = p.stdout.readline()
        yield line
        if(retcode is not None):
            break



# Program entry point
if __name__ == "__main__":
    main()
