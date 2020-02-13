#!/home/unix/cgdeboer/bin/python3
import warnings
import MYUTILS
import sys
import argparse
parser = argparse.ArgumentParser(description='This program is intended to quantify MPRA RNA tags given sequencing data of tags and a reference map of tags - enhancers.')
parser.add_argument('-it',dest='inTagMap',	metavar='<inTagMap>',help='Input file of the mRNA tag -enhancer map', required=True);
parser.add_argument('-iq',dest='inFastq',	metavar='<inFastq>',help='Input file of fastq barcodes', required=True);
parser.add_argument('-c',dest='constT', metavar='<constT>',help='the constant region following the reads [default=None]', required=False);
parser.add_argument('-mc',dest='mmc', metavar='<mismatchesConst>',help='the number of mismatches to allow in constant region [default=1]', required=False,default = 1);
parser.add_argument('-mt',dest='mmt', metavar='<mismatchesTag>',help='the number of mismatches to allow in tag region [default=0]', required=False);
parser.add_argument('-o',dest='outFPre', metavar='<outFilePre>',help='Where to output results, prefix', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-nc',dest='noConstCheck', action='count',help='Ignore constant region matching', required=False, default=0);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();
args.mmc = int(args.mmc);
args.mmt = int(args.mmt);

if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

sys.stderr.write("Compiling possible mismatches to constant region...\n");
#adds mismatches to tag2tag hash
def addMMToTags(myHash, baseSeq, numToAdd, alphabet, ref):
	if numToAdd<=0:
		if baseSeq in myHash and myHash[baseSeq]!=ref: #collision
			if baseSeq!=myHash[baseSeq]:#other is not an exact match
				myHash[baseSeq]="NA";
		else:	
			myHash[baseSeq]=ref;
	else:
		for i in range(0,len(baseSeq)):
			for a in range(0,len(alphabet)):
				addMMToTags(myHash, baseSeq[0:i]+alphabet[a]+baseSeq[(i+1):len(baseSeq)], numToAdd-1, alphabet, ref);

#loads mismatches to constant region into hash with the provided alphabet
def loadMismatches(myHash, baseSeq, numToAdd, alphabet):
	if numToAdd<=0:
		myHash[baseSeq]=1;
	else:
		for i in range(0,len(baseSeq)):
			for a in range(0,len(alphabet)):
				loadMismatches(myHash, baseSeq[0:i]+alphabet[a]+baseSeq[(i+1):len(baseSeq)], numToAdd-1, alphabet);

possibleMismatchesConst = {}
alphabet = ["A","T","G","C","N"];
for i in range(0,args.mmc+1):
	loadMismatches(possibleMismatchesConst,args.constT,i,alphabet);



sys.stderr.write("Loading tag-enhancer map...\n");
tag2enh = {}
tag2tagCount = {}
tag2tag = {}
tagLength = -1;
inTagMap=MYUTILS.smartGZOpen(args.inTagMap,'r');
for line in inTagMap:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	tag2enh[data[1]] = data[0];
	tag2tagCount[data[1]] = 0;
	for i in range(0,args.mmt+1): #add mismatched barcodes to map
		addMMToTags(tag2tag,data[1],i,alphabet,data[1]);
	if tagLength==-1:
		tagLength=len(data[1]);
	elif tagLength!=len(data[1]):
		raise Exception("Error: not all tag lengths the same! previous = %i, now = %i for line %s" %(tagLength, len(data[1]), line));
inTagMap.close();

sys.stderr.write("The length of the tags is %i.\n"%(tagLength));
sys.stderr.write("Reading fastq and matching to tags...\n");
outFileUnmatched = MYUTILS.smartGZOpen(args.outFPre+"_unmatched.fastq.gz",'w');
outFileBadConst = MYUTILS.smartGZOpen(args.outFPre+"_badConstRegion.fastq.gz",'w');
inFastq=MYUTILS.smartGZOpen(args.inFastq,'r');
state=0;
mismatchedConstantRegion=0;
unmatchedTag=0;
tagCollision=0;
successfulTags=0;
tagMismatchOkay=0;
reads=0;
curOutLoc=0;
for line in inFastq:
	if line is None or line == "": continue
	line = line.rstrip();
	if state==0 and line[0]=="@":
		lastHead=line;
		pass;
	elif state==1:
		reads+=1;
		curTag = line[0:tagLength];
		constRegion = line[tagLength:len(line)];
		if len(constRegion)!=len(args.constT):
			raise Exception("Error: provided constant region not same length as non-tag component of fastq: %i vs %i for line %s" %(len(args.constT), len(constRegion), line));
		elif constRegion not in possibleMismatchesConst:
			mismatchedConstantRegion+=1;
			if args.noConstCheck==0:
				outFileBadConst.write(lastHead+"\n"+line+"\n");
				curOutLoc=outFileBadConst;
		elif curTag not in tag2tag:
			unmatchedTag+=1;
			outFileUnmatched.write(lastHead+"\n"+line+"\n");
			curOutLoc=outFileUnmatched;
		elif tag2tag[curTag]=="NA":
			tagCollision+=1;
		if curTag in tag2tag and tag2tag[curTag]!="NA" and not curOutLoc: # matched a tag successfully, not a collision tag,  and  not a mismatched constant region
			realTag = tag2tag[curTag];
			if realTag!=curTag:
				tagMismatchOkay+=1;
			tag2tagCount[realTag]+=1;
			successfulTags+=1;
	elif state==2 and line=="+":
		if curOutLoc:
			curOutLoc.write(line+"\n");
	elif state==3:
		state=-1;
		if curOutLoc:
			curOutLoc.write(line+"\n");
		curOutLoc=0;
	state+=1;
inFastq.close();
outFileUnmatched.close();
outFileBadConst.close();

sys.stderr.write("Outputting results...\n");
outFile = MYUTILS.smartGZOpen(args.outFPre+"_counts.txt.gz",'w');

tagsInLib = 0;
distinctNonNATagsObserved=0;
distinctNATagsObserved=0;
tagsMappingToNA=0;
outFile.write("tag\tenhancer\tcount\n");
for tag in sorted(tag2enh):
	outFile.write("%s\t%s\t%i\n"%(tag,tag2enh[tag],tag2tagCount[tag]));
	if tag2enh[tag][0:2]=="NA":
		tagsMappingToNA+=1;
		if tag2tagCount[tag]>0:
			distinctNATagsObserved+=1;
	else:
		if tag2tagCount[tag]>0:
			distinctNonNATagsObserved+=1;
	tagsInLib+=1;
outFile.close();

sys.stderr.write("Among %i reads:\n"%(reads));
sys.stderr.write("	%i failed to align due to mismatched constant regions (%i%%).\n"%(mismatchedConstantRegion,100*mismatchedConstantRegion/reads));
sys.stderr.write("	%i failed to align due to a tag collision (%i%%).\n"%(tagCollision,100*tagCollision/reads));
sys.stderr.write("	%i failed to align due to the tag matching nothing known (%i%%).\n"%(unmatchedTag,100*unmatchedTag/reads));
sys.stderr.write("	%i tags found a home (%i%%).\n"%(successfulTags,100*successfulTags/reads));
sys.stderr.write("	%i tags found a home because we allowed mismatches (%i%%).\n"%(tagMismatchOkay,100*tagMismatchOkay/reads));
sys.stderr.write("Among the library of %i tags :\n"%(tagsInLib));
sys.stderr.write("	%i map to nothing (are NA) (%i%%).\n"%(tagsMappingToNA,100*tagsMappingToNA/tagsInLib));
sys.stderr.write("Among the %i NA tags :\n"%(tagsMappingToNA));
sys.stderr.write("	%i were observed >=1 time in this data (%i%%).\n"%(distinctNATagsObserved,100*distinctNATagsObserved/tagsMappingToNA));
sys.stderr.write("Among the %i enhancer-mapping tags:\n"%(tagsInLib - tagsMappingToNA));
sys.stderr.write("	%i were observed >=1 time in this data (%i%%).\n"%(distinctNonNATagsObserved,100*distinctNonNATagsObserved/(tagsInLib-tagsMappingToNA)));
if (args.logFP is not None):
	logFile.close();

sys.stderr.write("Done!\n");
