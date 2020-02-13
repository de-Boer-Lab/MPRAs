#!/home/unix/cgdeboer/bin/python3
import warnings
import MYUTILS
import sys
import re
import argparse
parser = argparse.ArgumentParser(description='Filters out barcode-enhancer combinations and tallies associations.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of mapped enhancers-barcodes as made by mapBarcodesToEnhancers.py', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFile>',help='Where to output results prefix', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();


inFile=MYUTILS.smartGZOpen(args.inFP,'r');


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;


#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
F_RID=0;
F_R1_REF=1;
F_R1_MAPQ=2;
F_R1_CIGAR=3;
F_R1_START=4;
F_R2_REF=5;
F_R2_MAPQ=6;
F_R2_CIGAR=7;
F_R2_START=8;
F_BAR=9;


outFileUnmapped = MYUTILS.smartGZOpen(args.outFPre+".reads.unmapped.gz",'w');
outFileGapped = MYUTILS.smartGZOpen(args.outFPre+".reads.gapped.gz",'w');
outFileDifferentRef = MYUTILS.smartGZOpen(args.outFPre+".reads.diffRef.gz",'w');
outFileDifferentMAPQ = MYUTILS.smartGZOpen(args.outFPre+".reads.diffMAPQ.gz",'w');
outFileUngapped = MYUTILS.smartGZOpen(args.outFPre+".reads.ungapped.gz",'w');
outFileLowMAPQ = MYUTILS.smartGZOpen(args.outFPre+".reads.lowMAPQ.gz",'w');

def startAndCigar2GapPoss(cigar, start):
	if re.match("^[0-9]*M$",cigar):
		return "M"
	if cigar=="*":
		return "*";
	m = re.match("^([0-9]*)M(.*[^0-9])[0-9]*M$",cigar)
	if m:
		newSt = start + int(m.group(1));
		return "%iM%s"%(newSt,m.group(2));
	else:
		raise Exception("Bad cigar; could not find terminal matches %s"%cigar);
		

enhancerBC2CIGARs = {};
tag2enhancer = {};
gapsBC2Enh = {};
ambiguousBC2Enh = {};
unmappableBCs = {};
total=0;
diffRef=0;
notMapped=0;
diffMapq=0;
mapqLow=0;
noGaps=0;
hasGaps=0;
numNBarcodes = 0;
header = inFile.readline().rstrip().split("\t");
for line in inFile:
	total+=1;
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	if re.search("N", data[F_BAR]):
		numNBarcodes+=1;
		next;
	if data[F_R1_REF] != data[F_R2_REF]:
		diffRef+=1;
		unmappableBCs[data[F_BAR]]=1;
		outFileDifferentRef.write(line);
		continue;
	if data[F_R1_REF]=="*":
		notMapped+=1;
		unmappableBCs[data[F_BAR]]=1;
		outFileUnmapped.write(line);
		continue;
	enh_bar = data[F_R1_REF]+"\t"+data[F_BAR];
	r1GapDesc = startAndCigar2GapPoss(data[F_R1_CIGAR], int(data[F_R1_START]));
	r2GapDesc = startAndCigar2GapPoss(data[F_R2_CIGAR], int(data[F_R2_START]));
	
	if enh_bar not in enhancerBC2CIGARs:
		enhancerBC2CIGARs[enh_bar] = [{},{}]; # Read1 cigar, read2 cigar
	if r1GapDesc not in enhancerBC2CIGARs[enh_bar][0]:
		enhancerBC2CIGARs[enh_bar][0][r1GapDesc]=0;
	if r2GapDesc not in enhancerBC2CIGARs[enh_bar][1]:
		enhancerBC2CIGARs[enh_bar][1][r2GapDesc]=0;
	enhancerBC2CIGARs[enh_bar][0][r1GapDesc]+=1;
	enhancerBC2CIGARs[enh_bar][1][r2GapDesc]+=1;
	if data[F_R1_MAPQ]!=data[F_R2_MAPQ]:
		diffMapq+=1;
		outFileDifferentMAPQ.write(line);
	if int(data[F_R1_MAPQ])+int(data[F_R1_MAPQ])<12: #exclude because mapping is ambiguous
		mapqLow+=1;
		outFileLowMAPQ.write(line)
		if data[F_BAR] not in ambiguousBC2Enh:
			ambiguousBC2Enh[data[F_BAR]]={};
		if data[F_R1_REF] not in ambiguousBC2Enh[data[F_BAR]]:
			ambiguousBC2Enh[data[F_BAR]][data[F_R1_REF]]=1;
		else:
			ambiguousBC2Enh[data[F_BAR]][data[F_R1_REF]]+=1;
	elif r1GapDesc=="M" and r2GapDesc=="M": #mapping has no gaps
		noGaps+=1;
		outFileUngapped.write(line)
		if data[F_BAR] not in tag2enhancer:
			tag2enhancer[data[F_BAR]]={};
		if data[F_R1_REF] not in tag2enhancer[data[F_BAR]]:
			tag2enhancer[data[F_BAR]][data[F_R1_REF]]=1;
		else:
			tag2enhancer[data[F_BAR]][data[F_R1_REF]]+=1;
	else: # these have gaps in the alignment, excluded
		hasGaps+=1;
		outFileGapped.write(line);
		if data[F_BAR] not in gapsBC2Enh:
			gapsBC2Enh[data[F_BAR]]={};
		if data[F_R1_REF] not in gapsBC2Enh[data[F_BAR]]:
			gapsBC2Enh[data[F_BAR]][data[F_R1_REF]]=1;
		else:
			gapsBC2Enh[data[F_BAR]][data[F_R1_REF]]+=1;

inFile.close();

sys.stderr.write("Of all %i reads:\n"%(total));
sys.stderr.write("  barcodes containing Ns: %i (%i%%)\n"%(numNBarcodes, 100.0 * numNBarcodes/total));
sys.stderr.write("  unmapped reads: %i (%i%%)\n"%(notMapped, 100.0 * notMapped/total));
sys.stderr.write("  differing reference: %i (%i%%)\n"%(diffRef, 100.0 * diffRef/total));
sys.stderr.write("  differing MAPQ: %i (%i%%)\n"%(diffMapq, 100.0 * diffMapq/total));
sys.stderr.write("  low MAPQ: %i (%i%%)\n"%(mapqLow, 100.0 * mapqLow/total));
sys.stderr.write("  no gaps in alignment: %i (%i%%)\n"%(noGaps, 100.0 * noGaps/total));
sys.stderr.write("  gaps in alignment: %i (%i%%)\n"%(hasGaps, 100.0 * hasGaps/total));

total=0;
uniqueC1=0;
uniqueC2=0;
singleObs=0;
uniqueBoth=0;
dominantC1=0;
dominantC2=0;
dominantBoth=0;
perfectC1=0;
perfectC2=0;
perfectBoth=0;

outFileCIGARs = MYUTILS.smartGZOpen(args.outFPre+".cigars.gz",'w');
for enh_bar in enhancerBC2CIGARs:
	curLine = enh_bar;
	total+=1;
	lastUnique=False;
	lastPerfect=False;
	lastDom=False;
	for i in range(0,2): # for each read in pair
		curCigars = [];
		maxObs = 0;
		numObs = 0;
		for cigar in enhancerBC2CIGARs[enh_bar][i]:
			curCigars.append(cigar+":%i"%(enhancerBC2CIGARs[enh_bar][i][cigar]));
			maxObs = max([maxObs, enhancerBC2CIGARs[enh_bar][i][cigar]]);
			numObs+=enhancerBC2CIGARs[enh_bar][i][cigar];
			if cigar=="M": #perfect match
				if lastPerfect:
					perfectBoth+=1;
				if i==0:
					perfectC1+=1;
					lastPerfect=True;
				else:
					perfectC2+=1;
		if len(curCigars)==1: # only one cigar for this read
			if i==0:
				uniqueC1+=1;
				lastUnique=True;
			else:
				uniqueC2+=1;
				if lastUnique:
					uniqueBoth+=1;
		if numObs==1 and i==0: #only one cigar observed and the first read (so as not to count this twice since the second should also have only one obs (since they are part of a pair))
			singleObs+=1;
		if 1.0*maxObs/numObs>0.5: # the maximally observed cigar for the read is at least 50% of the cigars observed for the read
			if i==0:
				dominantC1+=1;
				lastDom=True;
			else:
				dominantC2+=1;
				if lastDom:
					dominantBoth+=1;
		curLine=curLine+"\t"+";".join(curCigars);
	outFileCIGARs.write(curLine+"\n");

outFileCIGARs.close();
		
sys.stderr.write("Of the %i enhancer-barcode pairs that mapped to anything (well or otherwise):\n"%(total));
sys.stderr.write("  have only one read supporting them: %i (%i%%)\n"%(singleObs, 100.0 * singleObs/total));
sys.stderr.write("  have only one CIGAR1: %i (%i%%)\n"%(uniqueC1, 100.0 * uniqueC1/total));
sys.stderr.write("  have only one CIGAR2: %i (%i%%)\n"%(uniqueC2, 100.0 * uniqueC2/total));
sys.stderr.write("  have both only one C1 and C2: %i (%i%%)\n"%(uniqueBoth, 100.0 * uniqueBoth/total));
sys.stderr.write("  have a dominant C1: %i (%i%%)\n"%(dominantC1, 100.0 * dominantC1/total));
sys.stderr.write("  have a dominant C2: %i (%i%%)\n"%(dominantC2, 100.0 * dominantC2/total));
sys.stderr.write("  have both a dominant C1 and C2: %i (%i%%)\n"%(dominantBoth, 100.0 * dominantBoth/total));
sys.stderr.write("  have an ungapped C1: %i (%i%%)\n"%(perfectC1, 100.0 * perfectC1/total));
sys.stderr.write("  have an ungapped C2: %i (%i%%)\n"%(perfectC2, 100.0 * perfectC2/total));
sys.stderr.write("  have both an ungapped C1 and C2: %i (%i%%)\n"%(perfectBoth, 100.0 * perfectBoth/total));

usabletag2enh={};
mmbc2bc = {};
for bc in unmappableBCs:
	usabletag2enh[bc]="NA_unmappable";
	for k in range(0,len(bc)):
		mmbc = bc[0:k]+"-"+bc[(k+1):len(bc)];
		if mmbc not in mmbc2bc:
			mmbc2bc[mmbc]=[];
		mmbc2bc[mmbc].append(bc);

#fill in 1 mismatch hash
for bc in tag2enhancer:
	for k in range(0,len(bc)):
		mmbc = bc[0:k]+"-"+bc[(k+1):len(bc)];
		if mmbc not in mmbc2bc:
			mmbc2bc[mmbc]=[];
		mmbc2bc[mmbc].append(bc);

	


totalTags = 0;
barcodeEnhPairs=0;
barcodeCollisions=0;
barcodeCollisions1MM=0;
barcodeErrors=0;
barcodeWasUnmappable=0;
barcodeWasAmbiguousWithOther=0;
barcodeMappedToSomethingElseWithGaps=0;
enh2usabletag = {};
for bc in tag2enhancer:
	totalTags+=1;
	#is this tag good?
	enh=list(tag2enhancer[bc].keys())[0];
	if len(tag2enhancer[bc])>1: # this tag mapped to more than one enhancer
		enh="NA_collision";
		barcodeCollisions+=1;
	elif bc in unmappableBCs:
		barcodeWasUnmappable+=1;
	elif bc in ambiguousBC2Enh and (len(ambiguousBC2Enh[bc])>1 or enh not in ambiguousBC2Enh[bc]): # bc was ambiguous with something (either more than one other thing, or not itself)
		enh="NA_ambiguous";
		barcodeWasAmbiguousWithOther+=1;
	elif bc in gapsBC2Enh and enh not in gapsBC2Enh[bc]: # barcode was in a read with gaps and never mapped to the correct enhancer with this gapped read
		enh="NA_indel"
		barcodeMappedToSomethingElseWithGaps+=1;
	usabletag2enh[bc]=enh;
	if enh not in enh2usabletag:
		enh2usabletag[enh]=[];
	enh2usabletag[enh].append(bc);

for bc in ambiguousBC2Enh:
	if bc not in usabletag2enh:
		mostHits=0;
		bestGuess="NA";
		for enh in ambiguousBC2Enh[bc]:
			if ambiguousBC2Enh[bc][enh]>mostHits:
				mostHits = ambiguousBC2Enh[bc][enh];
				bestGuess=enh;
		enh = "NA_ambiguous_"+enh;
		usabletag2enh[bc]=enh;
		for k in range(0,len(bc)):
			mmbc = bc[0:k]+"-"+bc[(k+1):len(bc)];
			if mmbc not in mmbc2bc:
				mmbc2bc[mmbc]=[];
			mmbc2bc[mmbc].append(bc);

for bc in gapsBC2Enh:
	if bc not in usabletag2enh:
		mostHits=0;
		bestGuess="NA";
		for enh in gapsBC2Enh[bc]:
			if gapsBC2Enh[bc][enh]>mostHits:
				mostHits = gapsBC2Enh[bc][enh];
				bestGuess=enh;
		enh = "NA_indel_"+enh;
		usabletag2enh[bc]=enh;
		for k in range(0,len(bc)):
			mmbc = bc[0:k]+"-"+bc[(k+1):len(bc)];
			if mmbc not in mmbc2bc:
				mmbc2bc[mmbc]=[];
			mmbc2bc[mmbc].append(bc);



for bc in usabletag2enh:
	enh=usabletag2enh[bc]
	for k in range(0,len(bc)):
		mmbc = bc[0:k]+"-"+bc[(k+1):len(bc)];
		if len(mmbc2bc[mmbc])!=1: # there are more than one barcodes that mapped to this mmmBC
			for i in range(0,len(mmbc2bc[mmbc])):
				if mmbc2bc[mmbc][i]!=bc and usabletag2enh[mmbc2bc[mmbc][i]]!=enh: #not the current barcode and barcode maps to a different enhancer
					barcodeCollisions1MM+=1;
				elif mmbc2bc[mmbc][i]!=bc and enh[0:2]!="NA" and usabletag2enh[mmbc2bc[mmbc][i]]==enh: #not mapping to the current tag, not NA, and mapping to the same enhancer
					barcodeErrors+=1;

sys.stderr.write("Of the %i barcodes that appeared in isolation to map to something uniquely:\n"%(totalTags));
sys.stderr.write("  mapped to >=2 different enhancers in different reads: %i (%i%%)\n"%(barcodeCollisions, 100.0 * barcodeCollisions/totalTags));
sys.stderr.write("  were identical to an ambiguously mapping read (and not the same enh): %i (%i%%)\n"%(barcodeWasAmbiguousWithOther, 100.0 * barcodeWasAmbiguousWithOther/totalTags));
sys.stderr.write("  were identical to a barcode whose read was unmappable: %i (%i%%)\n"%(barcodeWasUnmappable, 100.0 * barcodeWasUnmappable/totalTags));
sys.stderr.write("  were identical to a barcode whose read had gaps (and not same enh): %i (%i%%)\n"%(barcodeMappedToSomethingElseWithGaps, 100.0 * barcodeMappedToSomethingElseWithGaps/totalTags));
sys.stderr.write("  barcodes that were within 1 edit distance for the same enhancer: %i (%i%%)\n"%(barcodeErrors, 100.0 * barcodeErrors/totalTags));
sys.stderr.write("  barcodes that collided within 1 edit distance: %i (%i%%)\n"%(barcodeCollisions1MM, 100.0 * barcodeCollisions1MM/totalTags));

totalEnh = len(enh2usabletag)-1;#because of NA
singleton=0;
moreThanFive=0;
moreThan10=0;
moreThan20=0;
moreThan40=0;
moreThan80=0;
moreThan250=0;
moreThan500=0;
moreThan1000=0;
moreThan10000=0;
moreThan100000=0;
for enh in enh2usabletag:
	numTags = len(enh2usabletag[enh]);
	if numTags==1:
		singleton+=1;
	if numTags>=5:
		moreThanFive+=1;
	if numTags>=10:
		moreThan10+=1;
	if numTags>=20:
		moreThan20+=1;
	if numTags>=40:
		moreThan40+=1;
	if numTags>=80:
		moreThan80+=1;
	if numTags>=250:
		moreThan250+=1;
	if numTags>=500:
		moreThan500+=1;
	if numTags>=1000:
		moreThan1000+=1;
	if numTags>=10000:
		moreThan10000+=1;
	if numTags>=100000:
		moreThan100000+=1;

sys.stderr.write("Of the %i enhancers that passed filtering:\n"%(totalEnh));
sys.stderr.write("  have only one barcode: %i (%i%%)\n"%(singleton, 100.0 * singleton/totalEnh));
sys.stderr.write("  have more than five barcodes: %i (%i%%)\n"%(moreThanFive, 100.0 * moreThanFive/totalEnh));
sys.stderr.write("  have more than 10 barcodes: %i (%i%%)\n"%(moreThan10, 100.0 * moreThan10/totalEnh));
sys.stderr.write("  have more than 20 barcodes: %i (%i%%)\n"%(moreThan20, 100.0 * moreThan20/totalEnh));
sys.stderr.write("  have more than 40 barcodes: %i (%i%%)\n"%(moreThan40, 100.0 * moreThan40/totalEnh));
sys.stderr.write("  have more than 80 barcodes: %i (%i%%)\n"%(moreThan80, 100.0 * moreThan80/totalEnh));
sys.stderr.write("  have more than 250 barcodes: %i (%i%%)\n"%(moreThan250, 100.0 * moreThan250/totalEnh));
sys.stderr.write("  have more than 500 barcodes: %i (%i%%)\n"%(moreThan500, 100.0 * moreThan500/totalEnh));
sys.stderr.write("  have more than 1000 barcodes: %i (%i%%)\n"%(moreThan1000, 100.0 * moreThan1000/totalEnh));
sys.stderr.write("  have more than 10000 barcodes: %i (%i%%)\n"%(moreThan10000, 100.0 * moreThan10000/totalEnh));
sys.stderr.write("  have more than 100000 barcodes: %i (%i%%)\n"%(moreThan100000, 100.0 * moreThan100000/totalEnh));

outFileMap = MYUTILS.smartGZOpen(args.outFPre+".map.gz",'w');
for bc in usabletag2enh:
	enh = usabletag2enh[bc];
	if enh[0:2]!="NA":
		outFileMap.write("%s\t%s\t%i\n"%(enh,bc,tag2enhancer[bc][enh])); # last column indicates the number of times the barcode was observed in the data
	else:
		outFileMap.write("%s\t%s\tNA\n"%(enh,bc)); # last column indicates nothing but not to use this barcode
outFileMap.close();

if (args.logFP is not None):
	logFile.close();
