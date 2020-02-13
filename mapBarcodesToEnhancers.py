#!/home/unix/cgdeboer/bin/python3
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='Parses sorted SAM file and associates barcode sequences with the enhancers they mapped to.')
parser.add_argument('-is',dest='inSAM',	metavar='<inFile>',help='Input file of aligned reads to enhancers in sam format. must be sorted (unix sort)', required=True);
parser.add_argument('-ib',dest='inBarcodes',	metavar='<inFile>',help='Input file of barcodes in read\tbarcode format (must also be sorted)', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();

SAM_ID=0;
SAM_STARTPOS=3;
SAM_CIGAR=5;
SAM_SCORE=4;
SAM_TLEN=8;
SAM_REF=2;

inSAM=MYUTILS.smartGZOpen(args.inSAM,'r');
inBC=MYUTILS.smartGZOpen(args.inBarcodes,'r');


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

if (args.outFP is None):
	outFile= sys.stdout;
else:
	if args.verbose>0: warnings.warn("Outputting to file "+args.outFP);
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');


nextSAMLine = inSAM.readline().rstrip().split("\t");
numUnmatchedSAM = 0;
numUnmatchedBC = 0;
outFile.write("readID\tR1_REF\tR1_MAPQ\tR1_START\tR1_CIGAR\tR2_REF\tR2_MAPQ\tR2_START\tR2_CIGAR\tbarcode\n");
for line in inBC:
	if line is None or line == "": continue
	data=line.rstrip().split("\t");
	while data[0]>nextSAMLine[SAM_ID] and (nextSAMLine[0]!="" or len(nextSAMLine)!=1):
		nextSAMLine = inSAM.readline().rstrip().split("\t");
		numUnmatchedSAM+=1;
	if data[0]==nextSAMLine[SAM_ID]:
		curSAMLine = nextSAMLine;
		nextSAMLine = inSAM.readline().rstrip().split("\t");
		if data[0]!=nextSAMLine[SAM_ID]:
			sys.stderr..write("Next line in SAM file is not the same read ID as predecessor; skipping read: %s - %s\n" %(data[0], nextSAMLine[SAM_ID]));
		else: 
			if int(curSAMLine[SAM_TLEN])<0 and int(nextSAMLine[SAM_TLEN])>0: #swap if not F read first
				temp=curSAMLine;
				curSAMLine=nextSAMLine;
				nextSAMLine=temp;
			outFile.write(data[0]+"\t"+curSAMLine[SAM_REF]+"\t"+curSAMLine[SAM_SCORE]+"\t"+curSAMLine[SAM_CIGAR] + "\t" + curSAMLine[SAM_STARTPOS] + "\t"+nextSAMLine[SAM_REF]+"\t"+nextSAMLine[SAM_SCORE]+"\t"+nextSAMLine[SAM_CIGAR]+"\t" + nextSAMLine[SAM_STARTPOS]+"\t"+data[1]+"\n");
			nextSAMLine = inSAM.readline().rstrip().split("\t");
	else:
		numUnmatchedBC+=1;

sys.stderr.write("Completed.  Num unmatched barcodes = %i; Num unmatched SAM entries = %i\n"%(numUnmatchedBC, numUnmatchedSAM)); 
inBC.close();
inSAM.close();

outFile.close();
if (args.logFP is not None):
	logFile.close();
