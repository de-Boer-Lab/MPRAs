#!/home/unix/cgdeboer/bin/python3
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='DESCRIPTION.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of fastq', required=True);
parser.add_argument('-s',dest='startPos',	metavar='<startPos>',help='Where the barcode starts in fastq', required=True);
parser.add_argument('-n',dest='numBases',	metavar='<numBases>',help='Length of barcode', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();


inFile=MYUTILS.smartGZOpen(args.inFP,'r');

args.startPos = int(args.startPos);
args.numBases = int(args.numBases);

if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

if (args.outFP is None):
	outFile= sys.stdout;
else:
	if args.verbose>0: warnings.warn("Outputting to file "+args.outFP);
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');



#@M03102:139:000000000-AGUL4:1:1101:15230:1342 1:N:0:0
#TTGACCTA
#+
#CCCCCFFF

#CTGTTCCGCTATACTC

state = 0;
for line in inFile:
	if line is None or line == "": continue
	if state==0 and line[0]=="@":
		data = line[1:].split(" ");
		curID = data[0];
	elif state==1:
		outFile.write(curID+"\t"+line[args.startPos:(args.startPos+args.numBases)]+"\n");
	elif state==2 and line[0]=="+":
		pass
	elif state==3:
		state=-1;
	else:
		raise Exception("Reached bad state=%d for '%s' at line '%s'" %(state,id,line));
	state+=1;
		
inFile.close();
outFile.close();
if (args.logFP is not None):
	logFile.close();
