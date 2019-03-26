#!/usr/bin/python
# Samuel Moijueh
# July 2016

import sys
import re
import time
import random
from decimal import *
import sys
import os
import subprocess
import fileinput
from optparse import OptionParser
from random import randint
from subprocess import call
from Bio.Seq import Seq;
from Bio import Entrez, SeqIO;


###########################################################
'''
== SYNOPOSIS ==
This Program uses Primer3 and BLAT to design primers for exon sequences length is less than or equal to 400

We use the Thousannd Genomes VCF file as an additional check for unique binding

INPUT:  1) Alamut file with seq

OUTPUT: 1) Alamut file with Primer
        2) Alamut file where exon seq > 400


USAGE: python PrimerGen.py -i Alamut_withSeq.bed -o Alamut_withPrimer.bed > ../PrimerOutput/Alamutfile_seq400.txt

'''
###########################################################

######### GLOBAL PARAMETERS ######################
usage = "python revised_4_PrimerGen.py -i AlamutSeqfile.txt";
parser = OptionParser(usage = usage);
parser.add_option("-p", "--master_path", dest="master_path", help="The path to the working directory of this pipeline");
parser.add_option("-i", "--input_File", dest="input", help="The path to the target alamut BED file. These are the exon sequences we're creating Primers for.");
parser.add_option("-t", "--ThousandFile", dest = "Thousand", help = "The path to the ThousandG Variants VCF file.");
totalTimeStart = time.time();


(options, args) = parser.parse_args();
#print len(sys.argv);
if(len(sys.argv) < 7) :
	print("python revised_PrimerGen.py -p 2-18_InSilico_PrimerDesign -i AlamutFile_1.bed -t 1_1000G.vcf");
	sys.exit();
if(not os.path.exists(options.input)) :
	print("Error: %s file doesn't exist, please check and try again" %(options.input))
	sys.exit();

MASTER_PATH = options.master_path;
strandFile = options.input
AlamutFile = open(strandFile,'r')
ThousandGenomeFile = options.Thousand;

############### CONSTANTS #################
SEQUENCE_MAX_LENGTH = 400;
SEQUENCE_BUFFER		= 320;
ThousandGenomeDict = dict();
SeqFileNumber = re.findall('\d+', strandFile)
SeqFileNumber = str(SeqFileNumber[2])


############## FILE HANDLES ################
timeLog_fh = open(MASTER_PATH + "/PrimerGen/output/JobLogs/AlamutFile" + str(SeqFileNumber) + "/Job" + str(SeqFileNumber) + "_TimeLog.txt", "w");
log_fh = open(MASTER_PATH + "/PrimerGen/output/JobLogs/AlamutFile" + str(SeqFileNumber) + "/Job" + str(SeqFileNumber) + "_PrimerGenLog.txt", "w");
# remind myself to append to file
outputfile1 = open(MASTER_PATH + '/PrimerGen/output/Processed_Primers/JobSet' + str(SeqFileNumber) + '_primers_mySQL_short.txt', 'w');
notFoundOutputfile = open(MASTER_PATH + '/PrimerGen/output/JobLogs/AlamutFile' + str(SeqFileNumber) +'/Job' + str(SeqFileNumber) + '_notFound_output.txt', 'w');


############## GLOBAL VARIABLES ###########
primer_counter = -1;
PreviousGeneName = "";


################################################
def main():
	######### LOAD 1000 Genomes VCF File #######
	total_start = time.time();
	ThousandG_Start = time.time();
	ThousandGenomeDict = LoadThousandGenomeFile();
	ThousandG_End = time.time();
	timeLog_fh.write("LoadThousandGenomeFile(): " + str(ThousandG_End - ThousandG_Start) + " seconds\n");

	######### GENERATE PRIMER PAIRS and SELECT THE BEST ONE ########
	GeneratePrimers();
	total_end = time.time();
	log_fh.write("TOTAL TIME: " + str(total_end - total_start) + "\n");
	log_fh.close();
	outputfile1.close()
	notFoundOutputfile.close()
	totalTimeEnd = time.time();
	timeLog_fh.write("Completed 4_PrimerGen.py: " + str(totalTimeEnd - totalTimeStart) + " seconds\n");
	timeLog_fh.close();
	log_fh.close();


# reads the thousand genomes variant file into a dictionary where the
# KEY is a string of the 'Chromosome;Genomic Position' and the VALUE is a placeholder (1).
# return the dictionary; we're only concerned about where the variants are located.
# NOTE: splitting the VCF file by the bed file significantly speeds up this function
# the thousand genomes VCF file contains all the variants found in the 1000 genomes project
def LoadThousandGenomeFile():
	for line in fileinput.input([ThousandGenomeFile]):
		line = line.strip();
		if line.startswith("#"):	# skip comments
			continue;
		else:
			lineSplit = line.split('\t')
			ThousandGenomeDict[lineSplit[0]+';'+lineSplit[1].strip()] = 1
	return ThousandGenomeDict;


# generates primer pairs for sequences found in the alamut file
def GeneratePrimers():
	global primer_counter;
	global PreviousGeneName;
	PRIMER_LENGTH = SEQUENCE_MAX_LENGTH + SEQUENCE_BUFFER*2;  #length of the primer we would like to generate

	handle = open("/group/das-lab/Samuel/Jag_PrimerDesign_Project/input/hg19.fa", "rU");
	genomeSeq = dict();
	for record in SeqIO.parse(handle, "fasta"):
		genomeSeq[record.id]=record.seq;
	for line in AlamutFile:
		AllprimersToBlat = dict();
		columns = line.split('\t');
		Sequence = columns[len(columns)-1].replace("\n","");
		if len(Sequence) <= 1040:
			# store the variables of the alamut file
			ChromosomeNo, ExonStart, ExonEnd, Transcript, Exon, Gene, Strand, sequence = columns;

			if primer_counter == -1:
				PreviousGeneName = Gene;
				primer_counter = 0;

			start = time.time();
			PrimerSeq = callPrimer3(Sequence,200,-200,ChromosomeNo,Exon,Transcript,columns[1],columns[2],AllprimersToBlat,Strand,SeqFileNumber)
			end = time.time();
			print "callPrimer3(): " + str(end - start) + " seconds";
			start = time.time();
			SeqDetails = BLATCheck(AllprimersToBlat,ExonStart,ExonEnd,SeqFileNumber,ChromosomeNo, Transcript, Gene, Exon, Strand);
			log_fh.write("Strand = " + str(Strand));
			log_fh.write("\n== The final primer string is ==\n");

			end = time.time();
			timeLog_fh.write("BLATCheck(): " + str(end - start) + " seconds\n");
			if SeqDetails == 0:
				log_fh.write("Unable to find a primer pair for this exon.\n");
				log_fh.write("== The original exon sequence is ==\n");
				log_fh.write(line + "\n");
				notFoundOutputfile.write(line);
			else:
				list_of_primer_coordinates = [];
				primer_seq_mapper = dict();
				SeqDetails_test_splitSeq = SeqDetails.split("_");
				true_right_seq = SeqDetails_test_splitSeq[0];
				true_left_test_split = SeqDetails_test_splitSeq[1].split('>');
				true_left_seq = true_left_test_split[0];
				SeqDetails_test_splitAgain= SeqDetails.split('>');

				SeqDetails_test_splitCoord = SeqDetails_test_splitAgain[1].split('^');
				true_left_almost_coord = SeqDetails_test_splitCoord[1].split("-");

				true_right_almost_coord = SeqDetails_test_splitCoord[0].split("-")

				true_right_coord = SeqDetails_test_splitCoord[0];

				coordinate1 = true_left_almost_coord[0].split(";")[1]; coordinate2 = true_left_almost_coord[0].split("=")[0];
				coordinate3 = true_right_almost_coord[0].split(";")[1]; coordinate4 = true_right_almost_coord[0].split("=")[0];
				coordinate2 = coordinate2.split(";")[1]; coordinate4 = coordinate4.split(";")[1];


				log_fh.write("coordinate1 = " + coordinate1 + " | coordinate2 = " + coordinate2 + " | coordinate3 = " + coordinate3 + " | coordinate4 = " + coordinate4 + "\n");
				log_fh.write(SeqDetails + "\n");
				list_of_primer_coordinates = [int(coordinate1), int(coordinate2), int(coordinate3), int(coordinate4)];
				ProductSeq = genomeSeq[ChromosomeNo][min(list_of_primer_coordinates):max(list_of_primer_coordinates)]

				log_fh.write("NOTE: the genome coordinates are relative to the primer's position on the FORWARD strand\n\n");
				log_fh.write("== The original exon sequence is ==\n");
				log_fh.write(line + "\n");

				log_fh.write("==ProductSeq (relative to the forward strand)==\n" + str(ProductSeq) + "\n");
				GC_prodSeq = GCcontent(str(ProductSeq));
				GC_prodSeq = 0;
				end = time.time();
				timeLog_fh.write("GCcontent(): " + str(end - start) + " seconds\n");

				GCType = "";
				if GC_prodSeq >= 60:
					GCType = "GC_RICH"
				else:
					GCType = "HOTSTAR"

				if abs(int(min(list_of_primer_coordinates))-(int(ExonStart)-320)) > 1040:
					notFoundOutputfile.write(line)
				else:
					if Strand == "-":
						FPrimerSeq = "tgtaaaacgacggccagt"+true_left_seq+"_"+"caggaaacagctatgacc"+true_right_seq;
					else:
						SeqDetails_splitforSeq = SeqDetails.split("_");
						LeftPrimerSeq = "tgtaaaacgacggccagt"+SeqDetails_splitforSeq[0];
						RightPrimer_splitforSeq = SeqDetails_splitforSeq[1].split('>');
						RightPrimerSeq = "caggaaacagctatgacc"+RightPrimer_splitforSeq[0];
						log_fh.write("LeftPrimerSeq = " + LeftPrimerSeq + "\n");
						log_fh.write("RightPrimerSeq = " + RightPrimerSeq + "\n");
						FPrimerSeq = LeftPrimerSeq+"_"+RightPrimerSeq;
					if (Gene == PreviousGeneName):
						primer_counter = primer_counter + 1;
					else:
						primer_counter = 1;
						PreviousGeneName = Gene;
					PrimerName = Gene+"_IP"+str(primer_counter);
					PCR_settings = "";
					if GCType == "HOTSTAR":
						PCR_settings = "TD57";
					else:
						PCR_settings = "TD51_2";
					ProductCoordinates = str(ChromosomeNo) + ":" + str(min(list_of_primer_coordinates)) + "-" + str(max(list_of_primer_coordinates));

					mySQL_command_short = "INSERT INTO PrimerInfo VALUES(\'"+PrimerName+"\',NULL,\'"+"chr"+ProductCoordinates+"\', \'"+FPrimerSeq+"\', \'INSILICO\', NULL,\'"+GCType+"\', \'"+PCR_settings+"\', \'True\', \'SEQ M13\', \'SEQ M13\', \'PASS\', \'smoijueh\');";
					outputfile1.write(mySQL_command_short + "\n");
					log_fh.write("FPrimerSeq = " + FPrimerSeq + " | PrimerName = " + PrimerName + "\n");
					log_fh.write("ProductCoordinate = " + ProductCoordinates + "\n");
					log_fh.write("GCType = " + GCType + "\n");
		else:
			notFoundOutputfile.write(line);
			log_fh.write("Outside: Apparently the primer is problematic\n");

		log_fh.write("\n\t\t======== NEXT EXON ========\n");



# in this function we check each primer pair for homology
# we return the best scoring primer pair
def BLATCheck(AllPrimerstoBlat,ExonStart,ExonEnd,SeqFileNumber,ChromosomeNumber, transcriptID, gene, exon, strand):
	prepareBLATFileStart = time.time();
	FileName = PrepareBLATFile(AllPrimerstoBlat,gene,exon, SeqFileNumber)
	prepareBLATFileEnd = time.time();
	timeLog_fh.write("In BLATCheck()->prepareBLATFile(): " + str(prepareBLATFileEnd - prepareBLATFileStart) + " seconds\n");
	blatCallStart = time.time();

	# create directory for each Job
	my_dir = MASTER_PATH + '/tool_configurations/BLAT/input/Job'+str(SeqFileNumber);
	if not os.path.exists(my_dir):
		mkdir_recursive(my_dir);

	# create directory for each Job
	my_dir = MASTER_PATH + '/tool_configurations/BLAT/output/Job'+str(SeqFileNumber);
	if not os.path.exists(my_dir):
		mkdir_recursive(my_dir);

	subprocess.call(['/apps/software/blat/35x1/blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=100 -out=blast8 /group/das-lab/Samuel/Jag_PrimerDesign_Project/script4_output/output/PrimerInput/hg19.2bit ' + MASTER_PATH + '/tool_configurations/BLAT/input/Job'+str(SeqFileNumber)+'/'+'queryBatch_'+gene+'_'+exon+'.fa '+ MASTER_PATH + '/tool_configurations/BLAT/output/Job'+str(SeqFileNumber)+'/BLAToutput_'+gene+'_'+exon+'.psl'], shell=True);

	blatCallEnd = time.time();
	timeLog_fh.write("In BLATCheck()->calling blat: " + str(blatCallEnd - blatCallStart) + " seconds\n");

	BlatFile = open(MASTER_PATH+"/tool_configurations/BLAT/output/Job"+str(SeqFileNumber)+"/BLAToutput_"+gene+"_"+exon+".psl","r");

	InputPrimerDict = {};
	FinalPrimers = {};
	HitCountDict = {};
	# parsing the BLAT file. store the chromosome start and stop in a dictionary
	for line in BlatFile:
		lineSplit = line.split("\t");
		SeqID = lineSplit[0];
		Chr = lineSplit[1];
		# start-stop-percentMatch
		StartandStop = str(lineSplit[8])+"-"+str(lineSplit[9])+"-"+str(lineSplit[2]);
		StrtoStore = str(Chr)+":"+StartandStop;
		# if the SeqID is already in the dictionary then append the next primer StrtoStore
		if SeqID in InputPrimerDict.keys():
			InputPrimerDict[SeqID] = InputPrimerDict[SeqID]+";"+StrtoStore;
		else:
			InputPrimerDict[SeqID] = StrtoStore;

	log_fh.write("== In BLATCheck(): Processing Primer Pairs for Gene: " + gene + "; transcriptID: " + str(transcriptID) + "; EXON: " + str(exon) + " ==\n");
	log_fh.write("First level check: Looking for primer pairs with one and only one hit match for both primers\n");

	# loop every primer pair - 1st level
	for i in range(0,len(InputPrimerDict),2):
		if "SEQ"+str(i) in InputPrimerDict.keys() and "SEQ"+str(i+1) in InputPrimerDict.keys():
			# the primer coordinates are stored in a list
			LPrimerChrandCoordinates = InputPrimerDict["SEQ"+str(i)].split(';');
			RPrimerChrandCoordinates = InputPrimerDict["SEQ"+str(i+1)].split(';');
			LprimerHits = len(LPrimerChrandCoordinates);
			RPrimerHits = len(RPrimerChrandCoordinates);
			log_fh.write("SEQ"+str(i) + ">" + AllPrimerstoBlat["SEQ"+str(i)] + " & SEQ"+str(i+1) + ">" + AllPrimerstoBlat["SEQ"+str(i+1)] + "\n");
			log_fh.write('LPrimerChrandCoordinates = ' + ','.join(LPrimerChrandCoordinates) +  ' LPrimerHits = '  + str(LprimerHits) + '\n');
			log_fh.write('RPrimerChrandCoordinates = ' + ','.join(RPrimerChrandCoordinates) + ' RPrimerHits = ' + str(RPrimerHits) + '\n');
			# if there is a perfect match then
			if LprimerHits==1 and RPrimerHits==1:
				log_fh.write("Is there one and only one hit match for both primers? [Answer]: YES\n");
				# store the left primer positions
				Chrom_Position_Split = LPrimerChrandCoordinates[0].split(':');
				LPositionSplit = Chrom_Position_Split[1].split('-');
				# store the left primer chromosome
				FinalLChr = Chrom_Position_Split[0];
				# store the left primer start position
				FinalLStart = LPositionSplit[0];
				# store the left primer end position
				FinalLEnd = LPositionSplit[1];
				# store the left primer percent match
				LPrimerPercentMatch = LPositionSplit[2];
				# store the right primer positions
				Chrom_Position_Split = RPrimerChrandCoordinates[0].split(':');
				RPositionSplit = Chrom_Position_Split[1].split('-');
				# store the right primer chromosome
				FinalRChr = Chrom_Position_Split[0];
				# store the right primer start position
				FinalRStart = RPositionSplit[0];
				# store the right primer end position
				FinalREnd = RPositionSplit[1];
				# store the right primer percent match
				RPrimerPercentMatch = RPositionSplit[2];

				# check that the left and right primers are on the same chromosome
				if FinalLChr == FinalRChr:
					log_fh.write("Are the primers on the same chromosome? [Answer]: YES\n");
					FinalChr = FinalLChr;

					# check that the left and right primers have 100 percent match
					if float(LPrimerPercentMatch) == 100.00 and float(RPrimerPercentMatch) == 100.00:
						log_fh.write("Is there 100 percent match for BOTH primers? [Answer]: YES\n");

						# check that the primers are on the same chromosome as the exon
						if FinalChr == ChromosomeNumber:
							log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: YES\n");
							# store the positions
							Left_Coord_to_Check = str(FinalLStart)+';'+str(FinalLEnd);
							Right_Coord_to_Check = str(FinalREnd)+';'+str(FinalRStart);

							CheckThousandG_Start = time.time();
							ThousandGCheck = CheckThousandG(Left_Coord_to_Check,Right_Coord_to_Check,FinalChr)
							CheckThousandG_End = time.time();
							# if there are no Thousand Genome variants in our primer pairs then store it's information
							# we have found our match, no need to continue searching therefore we break
							if ThousandGCheck == 1:
								log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: NO. \n");");
								# check if the primer identified encompasses and spans the exon coordinates
								if (FinalLStart < ExonStart and FinalLEnd < ExonStart):
									if (FinalRStart > ExonEnd and FinalREnd > ExonEnd):
										log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: YES\n");
										log_fh.write("Store this Primer Pair in the FinalPrimers{} dictionary. We're DONE\n\n");
										HitCountDict["SEQ"+str(i)] = LprimerHits+RPrimerHits;
										FinalPrimers["SEQ"+str(i)] = str(FinalChr)+';'+str(FinalLStart)+'-'+str(FinalLEnd)+'='+str(LprimerHits);
										FinalPrimers["SEQ"+str(i+1)] = str(FinalChr)+';'+str(FinalRStart)+'-'+str(FinalREnd)+'='+str(RPrimerHits);
										break;
									else:
										log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
										log_fh.write("Go to the next primer pair\n");
								else:
									log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
									log_fh.write("Go to the next primer pair\n");
							else:
								log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: YES\n");
								log_fh.write("Go to the next primer pair\n");
					else:
						log_fh.write("LPrimerPercentMatch = " + str(LPrimerPercentMatch) + "% | RPrimerPercentMatch = " + str(RPrimerPercentMatch) + "%\n");
						log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: NO, go to the next primer pair!\n\n");
				else :
					log_fh.write("FinalLChr = " + FinalLChr + "; FinalRChr = " + FinalRChr + "\n");
					log_fh.write("Are the primers on the same chromosome? [Answer] NO, go to the next primer pair!\n\n");
			else:
				log_fh.write("\n");

	# if we haven't found a primer pair on the first level then proceed to the 2nd level
	if bool(FinalPrimers) == False:
		log_fh.write("We were unable to find a primer pair on the first level. Now we check the 2nd level\n");
		log_fh.write("== Second level: Looking for primer pairs with one and only one hit match for EITHER primer pair ==\n\n");
		# loop every primer pair - 2nd level
		for i in range(0,len(InputPrimerDict),2):

			if "SEQ"+str(i) in InputPrimerDict.keys() and "SEQ"+str(i+1) in InputPrimerDict.keys():

				# the primer coordinates are stored in a list
				LPrimerChrandCoordinates = InputPrimerDict["SEQ"+str(i)].split(';');
				RPrimerChrandCoordinates = InputPrimerDict["SEQ"+str(i+1)].split(';');

				LprimerHits = len(LPrimerChrandCoordinates);
				RPrimerHits = len(RPrimerChrandCoordinates);

				log_fh.write("SEQ"+str(i) + ">" + AllPrimerstoBlat["SEQ"+str(i)] + " & SEQ"+str(i+1) + ">" + AllPrimerstoBlat["SEQ"+str(i+1)] + "\n");
				log_fh.write('LPrimerChrandCoordinates = ' + ','.join(LPrimerChrandCoordinates) +  ' LPrimerHits = '  + str(LprimerHits) + '\n');
				log_fh.write('RPrimerChrandCoordinates = ' + ','.join(RPrimerChrandCoordinates) + ' RPrimerHits = ' + str(RPrimerHits) + '\n');

				# if either one of the primers have a perfect match
				if LprimerHits==1 or RPrimerHits==1:
					log_fh.write("Is there one and only one hit match for EITHER primer pair? [Answer]: YES\n");
					# store the left primer positions
					Chrom_Position_Split = LPrimerChrandCoordinates[0].split(':');
					LPositionSplit = Chrom_Position_Split[1].split('-');
					# store the left primer chromosome
					FinalLChr = Chrom_Position_Split[0];
					# store the left primer start
					FinalLStart = LPositionSplit[0];
					# store the left primer end
					FinalLEnd = LPositionSplit[1];
					# store the left primer percent match
					LPrimerPercentMatch = LPositionSplit[2];
					# store the right primer positions
					Chrom_Position_Split = RPrimerChrandCoordinates[0].split(':');
					RPositionSplit = Chrom_Position_Split[1].split('-');
					# store the right primer chromosome
					FinalRChr = Chrom_Position_Split[0];
					# store the right primer start
					FinalRStart = RPositionSplit[0];
					# store the right primer end
					FinalREnd = RPositionSplit[1];
					# store the right primer percent match
					RPrimerPercentMatch = RPositionSplit[2];

					# check that the left and right primers are on the same chromosome
					if FinalLChr == FinalRChr:
						log_fh.write("Are the primers on the same chromosome? [Answer]: YES\n");
						FinalChr = FinalLChr;

						if float(LPrimerPercentMatch) == 100.00 and float(RPrimerPercentMatch) == 100.00:
							log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: YES\n");

							if FinalLChr == ChromosomeNumber:
								log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: YES\n");
								# store the positions
								Left_Coord_to_Check = str(FinalLStart)+';'+str(FinalLEnd);
								Right_Coord_to_Check = str(FinalREnd)+';'+str(FinalRStart);

								CheckThousandG_Start = time.time();
								ThousandGCheck = CheckThousandG(Left_Coord_to_Check,Right_Coord_to_Check,FinalChr);
								CheckThousandG_End = time.time();
								# if there are no Thousand Genome variants in our primer pairs then store it's information
								if ThousandGCheck == 1:
									log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: NO.\n");

									# check if the primer identified encompasses and spans the exon coordinates
									if (FinalLStart < ExonStart and FinalLEnd < ExonStart):
										if (FinalRStart > ExonEnd and FinalREnd > ExonEnd):
											log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: YES\n");
											log_fh.write("Store this Primer Pair in the FinalPrimers{} dictionary. We're DONE\n\n");
											HitCountDict["SEQ"+str(i)] = LprimerHits+RPrimerHits;
											FinalPrimers["SEQ"+str(i)] = str(FinalChr)+';'+str(FinalLStart)+'-'+str(FinalLEnd)+'='+str(LprimerHits);
											FinalPrimers["SEQ"+str(i+1)] = str(FinalChr)+';'+str(FinalRStart)+'-'+str(FinalREnd)+'='+str(RPrimerHits);
											break;
										else:
											log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
											log_fh.write("Go to the next primer pair\n");
									else:
										log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
										log_fh.write("Go to the next primer pair\n");
								else:
									log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: YES.\n");
									log_fh.write("Go to the next primer pair\n\n");
							else:
								log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: NO.\n\n");
						else:
							log_fh.write("LPrimerPercentMatch = " + str(LPrimerPercentMatch) + "% | RPrimerPercentMatch = " + str(RPrimerPercentMatch) + "%\n");
							log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: NO, go to the next primer pair!\n\n");

					else:	# we need to loop through the primer pair hits and compare
						log_fh.write(">> The first primer pair were NOT on the same chromosome so now we're going to loop through the COMBOS:\n");
						# if there's one match in the left primer hits then we loop through the right primer hits
						if LprimerHits == 1:
							for primer_coordinates in RPrimerChrandCoordinates:
								coordinates = ''.join(primer_coordinates);
								log_fh.write("LPrimer = " + ''.join(LPrimerChrandCoordinates) + "; RPrimer = "); #+ coordinates + "\n");
								log_fh.write(primer_coordinates); log_fh.write("\n");
								RChrom_PositionSplit = primer_coordinates.split(':');
								FinalRChr = RChrom_PositionSplit[0];
								FinalRPositionSplit = RChrom_PositionSplit[1].split('-');
								FinalRStart = FinalRPositionSplit[0];
								FinalREnd = FinalRPositionSplit[1];
								RPrimerPercentMatch = FinalRPositionSplit[2];
								if FinalLChr == FinalRChr:
									FinalChr = FinalLChr;
									log_fh.write("Are the primers on the same chromosome? [Answer]: YES\n");

									if float(LPrimerPercentMatch) == 100.00 and float(RPrimerPercentMatch) == 100.00:
										log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: YES\n");

										if FinalChr == ChromosomeNumber:
											log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: YES\n");
											# store the positions
											Left_Coord_to_Check = str(FinalLStart)+';'+str(FinalLEnd);
											Right_Coord_to_Check = str(FinalREnd)+';'+str(FinalRStart);
											CheckThousandG_Start = time.time();
											ThousandGCheck = CheckThousandG(Left_Coord_to_Check,Right_Coord_to_Check,FinalChr);
											CheckThousandG_End = time.time();
											# if there are no Thusand Genome variants in our primer pairs then store it's information
											if ThousandGCheck == 1:
												log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: NO.\n");
												# check if the primer identified encompasses and spans the exon coordinates
												if (FinalLStart < ExonStart and FinalLEnd < ExonStart):
													if (FinalRStart > ExonEnd and FinalREnd > ExonEnd):
														log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: YES\n");
														log_fh.write("Store this Primer Pair in the FinalPrimers{} dictionary. We're DONE\n\n");
														HitCountDict["SEQ"+str(i)] = LprimerHits+RPrimerHits;
														FinalPrimers["SEQ"+str(i)] = str(FinalChr)+';'+str(FinalLStart)+'-'+str(FinalLEnd)+'='+str(LprimerHits);
														FinalPrimers["SEQ"+str(i+1)] = str(FinalChr)+';'+str(FinalRStart)+'-'+str(FinalREnd)+'='+str(RPrimerHits);
														break;
													else:
														log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
														log_fh.write("Go to the next primer pair\n");
												else:
													log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
													log_fh.write("Go to the next primer pair\n");
											else:
												log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: YES.\n");
												log_fh.write("Go the next primer pair\n");
										else:
											log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: NO.\n\n");
									else:
										log_fh.write("LPrimerPercentMatch = " + str(LPrimerPercentMatch) + "% | RPrimerPercentMatch = " + str(RPrimerPercentMatch) + "%\n");
										log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: NO, go to the next primer pair!\n\n");

							# ASSERT: the primers are on the same chromosome
							log_fh.write("FinalLChr = " + FinalLChr + "; FinalRChr = " + FinalRChr + "\n");
							log_fh.write("Are the primers on the same chromosome? [Answer]: NO, go to the next primer pair!\n\n");
						# if there's one match in the right primer hits then we loop through the left primer hits
						if RPrimerHits == 1:
							for primer_coordinates in LPrimerChrandCoordinates:
								coordinates = ''.join(primer_coordinates);
								log_fh.write("RPrimer = " + ''.join(RPrimerChrandCoordinates) + "; LPrimer = ");# + coordinates + "\n");
								log_fh.write(primer_coordinates); log_fh.write("\n");
								LChrom_PositionSplit = primer_coordinates.split(':');
								FinalLChr = LChrom_PositionSplit[0];
								FinalLPositionSplit = LChrom_PositionSplit[1].split('-');
								FinalLStart = FinalLPositionSplit[0];
								FinalLEnd = FinalLPositionSplit[1];
								LPrimerPercentMatch = FinalLPositionSplit[2];
								if FinalRChr == FinalLChr:
									FinalChr = FinalRChr;
									log_fh.write("Are the primers on the same chromosome? [Answer]: YES\n");

									if float(LPrimerPercentMatch) == 100.00 and float(RPrimerPercentMatch) == 100.00:
										log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: YES\n");

										if FinalChr == ChromosomeNumber:
											log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: YES\n");
											# store the positions
											Left_Coord_to_Check = str(FinalLStart)+';'+str(FinalLEnd);
											Right_Coord_to_Check = str(FinalREnd)+';'+str(FinalRStart);

											CheckThousandG_Start = time.time();
											ThousandGCheck = CheckThousandG(Left_Coord_to_Check,Right_Coord_to_Check,FinalChr);
											CheckThousandG_End = time.time();


											# if there are no Thousand Genome variants in our primer pairs then store it's information
											if ThousandGCheck == 1:
												log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: NO.\n");

												# check if the primer identified encompasses and spans the exon coordinates
												if (FinalLStart < ExonStart and FinalLEnd < ExonStart):
													if (FinalRStart > ExonEnd and FinalREnd > ExonEnd):
														log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: YES\n");
														log_fh.write("Store this Primer Pair in the FinalPrimers{} dictionary. We're DONE\n\n");
														HitCountDict["SEQ"+str(i)] = LprimerHits+RPrimerHits;
														FinalPrimers["SEQ"+str(i)] = str(FinalChr)+';'+str(FinalLStart)+'-'+str(FinalLEnd)+'='+str(LprimerHits);
														FinalPrimers["SEQ"+str(i+1)] = str(FinalChr)+';'+str(FinalRStart)+'-'+str(FinalREnd)+'='+str(RPrimerHits);
														break;
													else:
														log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
														log_fh.write("Go to the next primer pair\n");
												else:
													log_fh.write("Does the primer pair identified encompasses and spans the exon coordinates? [Answer]: NO\n");
													log_fh.write("Go to the next primer pair\n");
											else:
												log_fh.write("Is there a Thousand Genome variant overlapping the primer pairs? [Answer]: YES.\n");
												log_fh.write("Go to the next primer pair\n");
												#continue;
										else:
											log_fh.write("Are the primers on the same chromosome as the exon? [Answer]: NO.\n\n");
									else:
										log_fh.write("LPrimerPercentMatch = " + str(LPrimerPercentMatch) + "% | RPrimerPercentMatch = " + str(RPrimerPercentMatch) + "%\n");
										log_fh.write("Is it 100 percent match for BOTH primers? [Answer]: NO, go to the next primer pair!\n\n");


							# ASSERT: the primers are on the same chromosome
							log_fh.write("FinalLChr = " + FinalLChr + "; FinalRChr = " + FinalRChr + "\n");
							log_fh.write("Are the primers on the same chromosome? [Answer]: NO, go to the next primer pair!\n\n");
		# ASSERT: we have looped all the primer pairs for the 2nd level.
		if (len(FinalPrimers.keys()) == 2):
			log_fh.write("== keys in FinalPrimers dictionary ==\n");
			log_fh.write("strand = " + str(strand) + "\n");
			for key, value in FinalPrimers.iteritems():
				log_fh.write("key = " + key + " , value = " + value + "\n");

			if (int(re.search(r'\d+', FinalPrimers.keys()[0]).group()) % 2) == 0:
				return AllPrimerstoBlat[FinalPrimers.keys()[0]]+"_"+AllPrimerstoBlat[FinalPrimers.keys()[1]]+">"+FinalPrimers[FinalPrimers.keys()[0]]+"^"+FinalPrimers[FinalPrimers.keys()[1]];
			else:
				return AllPrimerstoBlat[FinalPrimers.keys()[1]]+"_"+AllPrimerstoBlat[FinalPrimers.keys()[0]]+">"+FinalPrimers[FinalPrimers.keys()[1]]+"^"+FinalPrimers[FinalPrimers.keys()[0]];
		else:
			return 0;

	else:
		log_fh.write("We HAVE identified a primer pair on the first level\n");
		log_fh.write("strand = " + str(strand) + "\n");
		for key, value in FinalPrimers.iteritems():
			log_fh.write(">"+key+"\t"+value+"\n");

		if strand == '-':
			return AllPrimerstoBlat[FinalPrimers.keys()[1].upper()]+"_"+AllPrimerstoBlat[FinalPrimers.keys()[0].upper()]+">"+FinalPrimers[FinalPrimers.keys()[1]]+"^"+FinalPrimers[FinalPrimers.keys()[0]];
		else:
			return AllPrimerstoBlat[FinalPrimers.keys()[0]]+"_"+AllPrimerstoBlat[FinalPrimers.keys()[1]]+">"+FinalPrimers[FinalPrimers.keys()[0]]+"^"+FinalPrimers[FinalPrimers.keys()[1]];


# parses the Primer3 output for each of the primer pairs
# here we're identifying primers pairs that we may want to BLAT
# and storing them in the dictionary AllprimersToBlat
def GetPrimersfromPrimer3Output(PrimerResults,InputSeq,ChrNo,Transcript,OriginalStart,
	OriginalChrStart,OriginalChrStop,Exon,OriginalSequence,AllprimersToBlat,Strand,SeqFileNumber):
	#log_fh.write("In GetPrimersfromPrimer3Output\n");
	CheckGCandProdSizeFlag = 0;
	# split the Primer3 results file into a list
	PrimerResults_LineSplit = PrimerResults.split('\n');
	PrimerResultsDict = {};
	# store the Primer3 results file into a dictionary
	for eachLine in PrimerResults_LineSplit:
		if '=' in eachLine and eachLine!='=':
			eachLine_split = eachLine.split("=");
			PrimerResultsDict[eachLine_split[0]] = eachLine_split[1];


	if ('PRIMER_PAIR_NUM_RETURNED' in PrimerResultsDict.keys()):
		NumberOfPrimers = PrimerResultsDict["PRIMER_PAIR_NUM_RETURNED"];
	else:
		return AllprimersToBlat;

	# looping each primer of this exon
	for i in range(0,int(NumberOfPrimers)):
		ThousandGCheck = 0;

		Len_ofOriginalSeq = len(OriginalSequence);
		Len_of_Exon = Len_ofOriginalSeq - 640;

		# the primers are 20 bps each
		# LeftPrimerstart is the the sequence position of the left primer start
		# RightPrimerStart is the sequence position of the right primer start

		LeftPrimerstart = OriginalSequence.find(PrimerResultsDict['PRIMER_LEFT_'+str(i)+'_SEQUENCE']);

		RightPrimerStart = OriginalSequence.find(revcom(PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE']));

		LeftprimerPositionforname = SEQUENCE_BUFFER-LeftPrimerstart;

		RightPrimerPositionforname = (RightPrimerStart-SEQUENCE_BUFFER)+len(PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE']);

		# check that the primer GC content is less than or equal to 60%
		PrimerCheckVar = CheckGCandProdSize(PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'],PrimerResultsDict['PRIMER_LEFT_'+str(i)+'_SEQUENCE']);

		# if the GC content of both the left and right primer pairs are less or equal to 60% then proceed
		if PrimerCheckVar == 1:
			CheckGCandProdSizeFlag = 1
			# append the primer pairs separated by an underscore
			Pseq = PrimerResultsDict['PRIMER_LEFT_'+str(i)+'_SEQUENCE']+'_'+PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE']

			# skip primers that are less than 20 bps; we want primers that are >= 20 bps
			if len(PrimerResultsDict['PRIMER_LEFT_'+str(i)+'_SEQUENCE'])<20 or len(PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'])<20:
				continue;
			# this is how we keep track of each primer that belongs to the exon
			NumberOftheprimerLeft = "SEQ"+str(len(AllprimersToBlat));
			NumberOftheprimerRight = "SEQ"+str(len(AllprimersToBlat)+1);

			if Strand == '-':
				# KEY is the unique name of the primer; VALUE is the primer sequence
				AllprimersToBlat[NumberOftheprimerLeft] = PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'].upper();
				AllprimersToBlat[NumberOftheprimerRight] = PrimerResultsDict['PRIMER_LEFT_'+str(i)+'_SEQUENCE'].upper();
			# if the transcript we're creating primers for is on the forward strand
			# then we BLAT the primer from left to right
			else:
				AllprimersToBlat[NumberOftheprimerLeft] = PrimerResultsDict['PRIMER_LEFT_'+str(i)+'_SEQUENCE'].upper();
				AllprimersToBlat[NumberOftheprimerRight] = PrimerResultsDict['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'].upper();

		# we're storing the primer pairs that pass the first three conditions in a dictionary
		# these primer pairs will later be processed by BLAT to search for homology

	if OriginalStart == 200:
		return callPrimer3(OriginalSequence,150,-150,ChrNo,Exon,Transcript,OriginalChrStart,OriginalChrStop,
			AllprimersToBlat,Strand,SeqFileNumber);
	if OriginalStart == 150:
		return callPrimer3(OriginalSequence,100,-100,ChrNo,Exon,Transcript,OriginalChrStart,OriginalChrStop,
			AllprimersToBlat,Strand,SeqFileNumber);
	if OriginalStart == 100:
		return callPrimer3(OriginalSequence,50,-50,ChrNo,Exon,Transcript,OriginalChrStart,OriginalChrStop,
			AllprimersToBlat,Strand,SeqFileNumber);
	if OriginalStart == 50:
		return callPrimer3(OriginalSequence,0,0,ChrNo,Exon,Transcript,OriginalChrStart,OriginalChrStop,
			AllprimersToBlat,Strand,SeqFileNumber);
	if OriginalStart == 0:
		return AllprimersToBlat;

def prepareInputFileforPrimer3(Sequence,ChromosomeNo,Exon,Transcript,Start,Stop,SeqFileNumber):
	Seqlen = len(Sequence)
	Lstart = '1'
	Llength = '70'
	RStart = str(Seqlen-70)
	Rlength = '70'
	SequenceID = "Chr"+ChromosomeNo+"-"+Transcript+"_ex"+Exon
	Text_toPrint = "PRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_INTERNAL_OLIGO=0\nPRIMER_PICK_RIGHT_PRIMER=1\nPRIMER_NUM_RETURN=10\nPRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=5\nPRIMER_PRODUCT_SIZE_RANGE=100-1000\nPRIMER_PRODUCT_OPT_SIZE=0\nPRIMER_PAIR_WT_PRODUCT_SIZE_LT=1.0\nPRIMER_PAIR_WT_PRODUCT_SIZE_GT=1.0\nPRIMER_MIN_SIZE=20\nPRIMER_INTERNAL_MIN_SIZE=20\nPRIMER_OPT_SIZE=20\nPRIMER_INTERNAL_OPT_SIZE=20\nPRIMER_MAX_SIZE=28\nPRIMER_INTERNAL_MAX_SIZE=28\nPRIMER_WT_SIZE_LT=1.0\nPRIMER_INTERNAL_WT_SIZE_LT=1.0\nPRIMER_WT_SIZE_GT=1.0\nPRIMER_INTERNAL_WT_SIZE_GT=1.0\nPRIMER_MIN_GC=20.0\nPRIMER_INTERNAL_MIN_GC=20.0\nPRIMER_OPT_GC_PERCENT=50.0\nPRIMER_MAX_GC=80.0\nPRIMER_INTERNAL_MAX_GC=80.0\nPRIMER_WT_GC_PERCENT_LT=0.0\nPRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0\nPRIMER_WT_GC_PERCENT_GT=0.0\nPRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0\nPRIMER_GC_CLAMP=0\nPRIMER_MAX_END_GC=5\nPRIMER_MIN_TM=58.0\nPRIMER_INTERNAL_MIN_TM=58.0\nPRIMER_OPT_TM=60.0\nPRIMER_INTERNAL_OPT_TM=60.0\nPRIMER_MAX_TM=63.0\nPRIMER_INTERNAL_MAX_TM=63.0\nPRIMER_PAIR_MAX_DIFF_TM=100.0\nPRIMER_WT_TM_LT=1.0\nPRIMER_INTERNAL_WT_TM_LT=1.0\nPRIMER_WT_TM_GT=1.0\nPRIMER_INTERNAL_WT_TM_GT=1.0\nPRIMER_PAIR_WT_DIFF_TM=0.0\nPRIMER_PRODUCT_MIN_TM=-1000000.0\nPRIMER_PRODUCT_OPT_TM=0.0\nPRIMER_PRODUCT_MAX_TM=1000000.0\nPRIMER_PAIR_WT_PRODUCT_TM_LT=0.0\nPRIMER_PAIR_WT_PRODUCT_TM_GT=0.0\nPRIMER_TM_FORMULA=0\nPRIMER_SALT_MONOVALENT=50.0\nPRIMER_INTERNAL_SALT_MONOVALENT=50.0\nPRIMER_SALT_DIVALENT=0.0\nPRIMER_INTERNAL_SALT_DIVALENT=0.0\nPRIMER_DNTP_CONC=0.0\nPRIMER_INTERNAL_DNTP_CONC=0.0\nPRIMER_SALT_CORRECTIONS=0\nPRIMER_DNA_CONC=50.0\nPRIMER_INTERNAL_DNA_CONC=50.0\nPRIMER_MAX_SELF_ANY=8.00\nPRIMER_INTERNAL_MAX_SELF_ANY=12.00\nPRIMER_PAIR_MAX_COMPL_ANY=8.00\nPRIMER_WT_SELF_ANY=0.0\nPRIMER_INTERNAL_WT_SELF_ANY=0.0\nPRIMER_PAIR_WT_COMPL_ANY=0.0\nPRIMER_MAX_SELF_END=3.00\nPRIMER_INTERNAL_MAX_SELF_END=12.00\nPRIMER_PAIR_MAX_COMPL_END=3.00\nPRIMER_WT_SELF_END=0.0\nPRIMER_INTERNAL_WT_SELF_END=0.0\nPRIMER_PAIR_WT_COMPL_END=0.0\nPRIMER_MAX_END_STABILITY=100.0\nPRIMER_WT_END_STABILITY=0.0\nPRIMER_MAX_NS_ACCEPTED=0\nPRIMER_INTERNAL_MAX_NS_ACCEPTED=0\nPRIMER_WT_NUM_NS=0.0\nPRIMER_INTERNAL_WT_NUM_NS=0.0\nPRIMER_MAX_POLY_X=5\nPRIMER_INTERNAL_MAX_POLY_X=5\nPRIMER_MIN_THREE_PRIME_DISTANCE=-1\nPRIMER_PICK_ANYWAY=0\nPRIMER_LOWERCASE_MASKING=0\nPRIMER_EXPLAIN_FLAG=0\nPRIMER_LIBERAL_BASE=0\nPRIMER_FIRST_BASE_INDEX=0\nPRIMER_MAX_TEMPLATE_MISPRIMING=-1.00\nPRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=-1.00\nPRIMER_WT_TEMPLATE_MISPRIMING=0.0\nPRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0\nPRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=1\nPRIMER_MAX_LIBRARY_MISPRIMING=12.00\nPRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00\nPRIMER_PAIR_MAX_LIBRARY_MISPRIMING=24.00\nPRIMER_WT_LIBRARY_MISPRIMING=0.0\nPRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0\nPRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0\nPRIMER_MIN_QUALITY=0\nPRIMER_INTERNAL_MIN_QUALITY=0\nPRIMER_MIN_END_QUALITY=0\nPRIMER_QUALITY_RANGE_MIN=0\nPRIMER_QUALITY_RANGE_MAX=100\nPRIMER_WT_SEQ_QUAL=0.0\nPRIMER_INTERNAL_WT_SEQ_QUAL=0.0\nPRIMER_PAIR_WT_PR_PENALTY=1.0\nPRIMER_PAIR_WT_IO_PENALTY=0.0\nPRIMER_INSIDE_PENALTY=-1.0\nPRIMER_OUTSIDE_PENALTY=0.0\nPRIMER_WT_POS_PENALTY=1.0\nPRIMER_SEQUENCING_LEAD=50\nPRIMER_SEQUENCING_SPACING=500\nPRIMER_SEQUENCING_INTERVAL=250\nPRIMER_SEQUENCING_ACCURACY=20\nPRIMER_WT_END_QUAL=0.0\nPRIMER_INTERNAL_WT_END_QUAL=0.0"
	Text_toPrint = Text_toPrint + "\n" + "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST="+Lstart+","+Llength+","+RStart+","+Rlength
	Text_toPrint = Text_toPrint + "\n" + "SEQUENCE_ID=" + SequenceID
	Text_toPrint = Text_toPrint + "\n" + "SEQUENCE_TEMPLATE=" + Sequence
	Text_toPrint = Text_toPrint + "\n" + "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/apps/software/primer3/2.3.6/src/primer3_config/"
	Text_toPrint = Text_toPrint + "\n" + "="
	print SequenceID + "_input"+"_AlamutFile_"+str(SeqFileNumber);

	Primer3Inputfile = open(MASTER_PATH+"/tool_configurations/Primer3/input/"+SequenceID+"_input"+"_AlamutFile_"+str(SeqFileNumber), 'w');
	Primer3Inputfile.write(Text_toPrint)
	Primer3Inputfile.close()
	return SequenceID+"_input"+"_AlamutFile_"+str(SeqFileNumber)

def callPrimer3(Sequence,Start,Stop,ChromosomeNo,Exon,Transcript,ChrStart,ChrStop,AllprimersToBlat,Strand,SeqFileNumber):
	if Stop == 0 :
		InitialSequence = Sequence
	else:
		InitialSequence = Sequence[Start:Stop]
	OriginalChrStart = ChrStart
	OriginalChrStop = ChrStop
	prepareInputFileforPrimer3_Start = time.time();
	FileName = MASTER_PATH + "/tool_configurations/Primer3/input/" + prepareInputFileforPrimer3(InitialSequence,ChromosomeNo,Exon,Transcript,Start,Stop,SeqFileNumber)
	prepareInputFileforPrimer3_End = time.time();
	timeLog_fh.write("In callPrimer3()->prepareInputFileforPrimer3(): " + str(prepareInputFileforPrimer3_End - prepareInputFileforPrimer3_Start) + " seconds\n");
	primer_results_fh = open(MASTER_PATH + '/tool_configurations/Primer3/output/' + str(SeqFileNumber) + '_Primer3_results.txt', 'a');
	primer_results_fh.write("File = " + FileName);
	print FileName;

	p = subprocess.Popen(['/apps/software/primer3/2.3.6/src/./primer3_core', FileName], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	PrimerResults, err = p.communicate()

	primer_results_fh.write("= Results = \nInitialSequence Length = " + str(len(InitialSequence)) + "\n" + PrimerResults + "=== POPO === \n");

	GetPrimersfromPrimer_Start = time.time();
	PrimerSeq = GetPrimersfromPrimer3Output(PrimerResults,InitialSequence,ChromosomeNo,Transcript,Start,OriginalChrStart,OriginalChrStop,Exon,Sequence,AllprimersToBlat,Strand,SeqFileNumber)
	GetPrimersfromPrimer_End = time.time();
	timeLog_fh.write("In callPrimer3()->GetPrimersfromPrimer3Output(): " + str(GetPrimersfromPrimer_End - GetPrimersfromPrimer_Start) + " seconds\n");

	return None;


#### HELPER FUNCTIONS ####

# returns the GC content percentage in the primer sequence
def GCcontent(ProdSeq):
        gcCount = 0
        totalBaseCount = 0
        gcCount += len(re.findall("[GC]", ProdSeq.upper()))
        totalBaseCount += len(re.findall("[GCTA]", ProdSeq.upper()))
        print ProdSeq,gcCount,totalBaseCount
        gcFraction = float(gcCount) / float(totalBaseCount)
        ProductGCpercent = gcFraction * 100
        return ProductGCpercent

# returns the complement of a genomic sequence
def complement(inputSeq):
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	letters = list(inputSeq)
	Seq = ""
	for eachbase in letters:
		base=eachbase.upper()
		if base=='N':
			Seq = Seq+'N'
		else:
			Seq = Seq + basecomplement[base]
	return Seq

# returns the reverse complement of a genomic sequence
def revcom(inputSeq):
	Seq = complement(inputSeq[::-1])
	return Seq

# check each genomic position of the primer pairs
# if the either primer has a SNP from the Thousand Genome then
# return 0 otherwise return 1
def CheckThousandG(LeftCoord,RightCoord,ChrNo):
	# split the start and end coordinates into a list
	LeftCoord_split = LeftCoord.split(';');
	RightCoord_split = RightCoord.split(';');
	# generates a list of numbers spanning the start and end coordinates of each primer
	LeftCoord_Range = range(int(LeftCoord_split[0]),int(LeftCoord_split[1])+1);
	RightCoord_Range = range(int(RightCoord_split[0]),int(RightCoord_split[1])+1);

	# check each coordinate within the left primer
	for eachCoord in LeftCoord_Range:
		CoordStr = ChrNo+';'+str(eachCoord);
		# if there's a SNP in that genomic position then return 0 otherwise continue searching
		if CoordStr in ThousandGenomeDict:
			return 0;
		else:
			continue;
	# check each coordinate within the right primer
	for eachCoord in RightCoord_Range:
		CoordStr = ChrNo+';'+str(eachCoord);
		# if there's a SNP in that genomic position then return 0 otherwise continue searching
		if CoordStr in ThousandGenomeDict:
			return 0;
		else:
			continue;
	# if we were unable to find a SNP within the primer pairs then return 1
	return 1;

# if the GC content of both the left AND right primers is
# less than or equal to 60 then return 1
# otherwise return 0
def CheckGCandProdSize(Leftprimer,RightPrimer):
	PrimerLeftGC = GCcontent(Leftprimer)
	PrimerRightGC = GCcontent(RightPrimer)
	PrimerGCFlag = 0
	if PrimerLeftGC<=60 and PrimerRightGC<=60:
		PrimerGCFlag = 1
	if PrimerGCFlag==1:
		return 1
	else:
		return 0

# creates the input file for BLAT. FASTA format
# >SEQ0
# atcgttga...atcg
def PrepareBLATFile(AllPrimerstoBlat,gene,exon, SeqFileNumber):
	# create directory for each Job
	my_dir = MASTER_PATH + '/tool_configurations/BLAT/input/Job'+str(SeqFileNumber);
	if not os.path.exists(my_dir):
		mkdir_recursive(my_dir);

	BLATInputfile = open(MASTER_PATH+"/tool_configurations/BLAT/input/Job"+str(SeqFileNumber)+"/queryBatch_"+gene+"_"+exon+".fa", 'w');
	for key, value in AllPrimerstoBlat.iteritems():
		BLATInputfile.write(">"+key+"\n"+value+"\n")
	BLATInputfile.close()
	return 1

def mkdir_recursive(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise



# Program entry point
if __name__ == "__main__":
    main()

