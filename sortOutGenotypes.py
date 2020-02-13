#!/usr/bin/python
import warnings
import MYUTILS
import re
import sys
import argparse
parser = argparse.ArgumentParser(description='Figure out the mapping of enhancer sequences to the SNP IDs and such contained within.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of excel', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFile>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();




if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;


#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));

I_VARID=0;
I_SNP=1;
I_CHR=2;
I_POS=3;
I_GC=4;
I_LEN=5;
I_J1=6;
I_REFA=7;
I_ALTA=8;
I_INCSNP=9;
I_SEQL=10;
I_SEQM=11;
I_SEQR=12;
I_SEQ=13;
I_PARTNERIDS=14;
I_PARTNERLOCS=15;
I_PARTNERAS=16;

snp2partner = {}; # dict containing primary SNP -> partnerSNPs
snp2all = {}; # dict containing SNP -> alleles
snp2ref = {}; #snp ID to ref allele 
snp2position = {};

inFile=MYUTILS.smartGZOpen(args.inFP,'r');
inFile.readline();
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	data[I_SEQ]=data[I_SEQ][7:(len(data[I_SEQ])-7)];
	curAltSNPs = data[I_PARTNERIDS].split(";")
	multiAllele = re.match("^(.*)_[0-9]+$",data[I_SNP])
	if multiAllele:
		data[I_SNP] = multiAllele.group(1);
	if data[I_SNP] not in snp2partner:
		snp2partner[data[I_SNP]]={};
		snp2all[data[I_SNP]] ={};
		snp2ref[data[I_SNP]] = data[I_REFA];
		snp2position[data[I_SNP]] = data[I_POS];
	snp2all[data[I_SNP]][data[I_SEQM]]=0;
	for i in range(0,len(curAltSNPs)):
		if curAltSNPs[i]!="" and curAltSNPs[i]!="-":
			multiAllele = re.match("^(.*)_[0-9]+$",curAltSNPs[i])
			if multiAllele:
				curAltSNPs[i] = multiAllele.group(1);
			snp2partner[data[I_SNP]][curAltSNPs[i]]=0;

snpseq2allele= {}
#go through the file again, recording any SNP_seqs that map to two different alleles (exclude these later);
inFile.seek(0,0);
inFile.readline();
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	data[I_SEQ]=data[I_SEQ][7:(len(data[I_SEQ])-7)];
	multiAllele = re.match("^(.*)_[0-9]+$",data[I_SNP])
	if multiAllele:
		data[I_SNP] = multiAllele.group(1);
	partnerSNPs = sorted(snp2partner[data[I_SNP]].keys());
	curAltSNPs = data[I_PARTNERIDS].split(";")
	curAltGenotypes = data[I_PARTNERAS].split(";");
	for i in range(0,len(partnerSNPs)):
		snp_seq = partnerSNPs[i]+"_"+data[I_SEQ];
		if partnerSNPs[i] not in snp2ref:
			allele="NA";
		else:
			allele = snp2ref[partnerSNPs[i]];
		for j in range(0,len(curAltSNPs)):
			multiAllele = re.match("^(.*)_[0-9]+$",curAltSNPs[i])
			if multiAllele:
				curAltSNPs[i] = multiAllele.group(1);
			if curAltSNPs[j]==partnerSNPs[i]: #override the previous allele if this instance is non-reference
				allele=curAltGenotypes[j];
		if snp_seq not in snpseq2allele:
			snpseq2allele[snp_seq]={};
		snpseq2allele[snp_seq][allele]=0;

outFilePrimary = MYUTILS.smartGZOpen(args.outFPre+"_byPrimary.txt.gz",'w');
outFileFull = MYUTILS.smartGZOpen(args.outFPre+"_all.txt.gz",'w');
outFilePrimary.write("VarID\tPrimarySNP\tchr\tpos\tstrand\tPrimaryAllele\tSequence\tPartnerSNPs\tPartnerAlleles\tPartnerPosition\n");
outFileFull.write("SNP\tchr\tpos\tstrand\tPosPrimary\tIsPrimary\tAllele\tSequence\n");
missingPartnerSNPs = {};
inFile.seek(0,0);
inFile.readline();
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	data[I_SEQ]=data[I_SEQ][7:(len(data[I_SEQ])-7)];
	strand = "+";
	if re.match("^.*_RC$",data[I_VARID]):
		strand="-";
	multiAllele = re.match("^(.*)_[0-9]+$",data[I_SNP])
	if multiAllele:
		data[I_SNP] = multiAllele.group(1);
	partnerSNPs = sorted(snp2partner[data[I_SNP]].keys());
	partnerAlleles = [];
	partnerPositions = [];
	curAltSNPs = data[I_PARTNERIDS].split(";")
	curAltGenotypes = data[I_PARTNERAS].split(";");
	outFileFull.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(data[I_SNP],data[I_CHR],data[I_POS], strand, data[I_POS],"T", data[I_SEQM], data[I_SEQ]));
	for i in reversed(range(0,len(partnerSNPs))): # throw out all the SNPs that aren't actually present in this sequence
		snp_seq = partnerSNPs[i]+"_"+data[I_SEQ];
		if len(snpseq2allele[snp_seq])>1:
			del snp2partner[data[I_SNP]][partnerSNPs[i]];
			del partnerSNPs[i];
	for i in range(0,len(partnerSNPs)):
		if partnerSNPs[i] not in snp2ref:
			missingPartnerSNPs[partnerSNPs[i]]=1;
			partnerAlleles.append("NA"); #default is the reference
			partnerPositions.append("NA"); #default is the reference
		else:
			partnerAlleles.append(snp2ref[partnerSNPs[i]]); #default is the reference
			partnerPositions.append(snp2position[partnerSNPs[i]]); #default is the reference
		for j in range(0,len(curAltSNPs)):
			multiAllele = re.match("^(.*)_[0-9]+$",curAltSNPs[i])
			if multiAllele:
				curAltSNPs[i] = multiAllele.group(1);
			if curAltSNPs[j]==partnerSNPs[i]:
				partnerAlleles[i]=curAltGenotypes[j];
		if partnerSNPs[i] not in snp2ref:
			outFileFull.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(partnerSNPs[i],data[I_CHR],"NA",strand, data[I_POS],"F", partnerAlleles[i], data[I_SEQ]));
		else:
			outFileFull.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(partnerSNPs[i],data[I_CHR],snp2position[partnerSNPs[i]],strand, data[I_POS],"F", partnerAlleles[i], data[I_SEQ]));
	outFilePrimary.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(data[I_VARID], data[I_SNP],data[I_CHR],data[I_POS],strand,data[I_SEQM],data[I_SEQ],";".join(partnerSNPs),";".join(partnerAlleles), ";".join(partnerPositions)));

sys.stderr.write("There were %i missing partner SNPs:\n"%(len(missingPartnerSNPs)));
sys.stderr.write("%s\n"%("\n".join(missingPartnerSNPs.keys())));
inFile.close();
outFileFull.close();
outFilePrimary.close();
if (args.logFP is not None):
	logFile.close();
