#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from Bio.SeqIO import parse
import os
import random
from collections import Counter

###########################################################################################################
## python3 ./fasta2bpfile.py chr_10.MSA.fasta ADR INT INV MAC 5000 0.2 4 0.00000001 0.0000000001 100000  ##
## argument 1: name of the fasta file (ex: chr_10.MSA.fasta)                                             ##
## argument 2: name for species A (ex: ADR)                                                              ##
## argument 3: name for species B (ex: INT)                                                              ##
## argument 4: name for species C (ex: INV)                                                              ##
## argument 5: name for species D (ex: MAC)                                                              ##
## argument 6: name of the outgroup (ex: truc)
## argument 7: size of the bin in bp (ex: 5000)                                                          ##
## argument 8: max proportion of missing data in the alignement (ex: 0.2)                                ##
## argument 9: minimum number of individuals in a given species at a given locus. (ex: 4)                ##
## argument 10: mutation rate mu per nucleotide and per generation (ex: 0.000000001)                      ##
## argument 11: recombinaration rate (ex: 0.0000000001 --> keep it small for big bins, i.e, 0.05mu       ##
## argument 12: size of the reference population (Nref) arbitrary fixed for ABC (ex: 100000)             ##
##                                                                                                       ##
## produces 3 output files:                                                                              ##
##         . general_infos.txt : informations about each bins (position, length, number of SNPs, etc...) ##
##         . observed_data.ms : the observed sequenced data in ms's format                               ##
##         . bpfile : used with observed_data.ms by mscalc to compute statistics                         ##
###########################################################################################################

bases = ['A', 'T', 'G', 'C']
def fasta2ms(aln, region, nameOut):
	nTot = len(aln['spA']) + len(aln['spB']) + len(aln['spC']) + len(aln['spD'])
	locus_length = len(aln['spA'][0])
	
	if region=='coding':
		liste_positions = range(2+3, locus_length-3, 3) # from the third position of the second codon, to the codon before the last one
	else:
		liste_positions = range(locus_length)

	nPos=len(liste_positions)
	nPos_valid = 0
	nPos_biallelic = 0
	positions = []	
	alignements = []
	# loop over all positions (or all 3rd coding positions)
	for pos in liste_positions:
		alleles_spA = [] # alleles from spA at position 'pos'
		alleles_spB = [] # alleles from spB at position 'pos'
		alleles_spC = [] # alleles from spC at position 'pos'
		alleles_spD = [] # alleles from spD at position 'pos'
		if nameOut!='NA':
			alleles_outgroup = []
		for ind in aln['spA']:
			nuc = ind[pos]
			if nuc in bases:
				alleles_spA.append(nuc)
		for ind in aln['spB']:
			nuc = ind[pos]
			if nuc in bases:
				alleles_spB.append(nuc)
		for ind in aln['spC']:
			nuc = ind[pos]
			if nuc in bases:
				alleles_spC.append(nuc)
		for ind in aln['spD']:
			nuc = ind[pos]
			if nuc in bases:
				alleles_spD.append(nuc)

		alleles_total = alleles_spA + alleles_spB + alleles_spC + alleles_spD

		# if there is an outgroup
		if nameOut!='NA':
			for ind in aln['nameOut']:
				nuc_out = ind[pos]
				if nuc_out in bases:
					alleles_outgroup.append(nuc_out)
			# randomly sample one concensus base from outgroup to orientate mutations
			if len(alleles_outgroup)>0:
				nuc_out=random.sample(population=alleles_outgroup, k=1)[0]
			else:
				nuc_out='NA'
		if len(alleles_total)==nTot:
			segregating_alleles = set(alleles_total)
			if nameOut!='NA': # if there is an outgroup
				if nuc_out!='NA': # if the concensus base is not 'NA'
					nPos_valid += 1
					segregating_alleles.add(nuc_out)
					if len(segregating_alleles) == 2: # if biallelic
						positions.append(pos/(1.0*locus_length))
						nPos_biallelic += 1
						alleles_total_msFormat = [ '0' if i==nuc_out else '1' for i in alleles_total ]
						alignements.append(alleles_total_msFormat)
			else: # if there is no outgroup
				nPos_valid += 1
				if len(segregating_alleles) == 2: # if biallelic
					positions.append(pos/(1.0*locus_length))
					nPos_biallelic += 1
					compteur = Counter(alleles_total)
					minor_allele = min(compteur, key=compteur.get)
					alleles_total_msFormat = [ '1' if i==minor_allele else '0' for i in alleles_total ] # minor allele==1 to count them more easily
					alignements.append(alleles_total_msFormat)

	res = {}
	res['nTot'] = nTot
	res['locus_length'] = locus_length
	res['nPos'] = nPos
	res['nPos_valid'] = nPos_valid
	res['nPos_biallelic'] = nPos_biallelic
	res['positions'] = positions
	res['nA']=len(aln['spA'])
	res['nB']=len(aln['spB'])
	res['nC']=len(aln['spC'])
	res['nD']=len(aln['spD'])
	
	if nPos_biallelic > 0:
		seq = '\n//\nsegsites: {0}\npositions: {1}\n'.format(nPos_biallelic, ' '.join([ str(i) for i in positions ]))
		for ind in range(nTot):
			for pos in range(nPos_biallelic):
				seq += '{0}'.format(alignements[pos][ind])
			seq += '\n'
	else:
		seq = '\n//\nsegsites: 0\n'
	res['ms'] = seq
    #print('ms alignment')
    #print(res['ms'])
	return(res)

nameOut='NA'
for tmp in sys.argv:
	arg=tmp.split('=')
	if arg[0]=='fastafile':
		fastafile=arg[1]
	if arg[0]=='region':
		region=arg[1] # coding ; noncoding
	if arg[0]=='spA':
		spA=arg[1]
	if arg[0]=='spB':
		spB=arg[1]
	if arg[0]=='spC':
		spC=arg[1]
	if arg[0]=='spD':
		spD=arg[1]
	if arg[0]=='nameOut':
		nameOut=arg[1]
	if arg[0]=='maxN': # max prop. of N tolerated in the alignement
		maxN=float(arg[1])
	if arg[0]=='nMin': # minimum number of individuals in a given species at a given locus
		nMin=int(arg[1])
	if arg[0]=='Lmin': # minimum number of valid positions to retain a locus
		Lmin=int(arg[1])
	if arg[0]=='mu':
		mu=float(arg[1]) # mutation rate per generation and per bp
	if arg[0]=='rec':
		rec=float(arg[1]) # recombination rate per generation and per bp
	if arg[0]=='Nref':
		Nref=int(arg[1]) # Nref : 100000
	if arg[0]=='datapath':
		datapath=arg[1] # Nref : 100000
	if arg[0]=='binpath':
		binpath=arg[1] # Nref : 100000
	if arg[0]=='simulationpath':
		simulationpath=arg[1] # Nref : 100000

#fastafile = sys.argv[1] # name of the fasta file
#spA = sys.argv[2] # name of species A
#spB = sys.argv[3]
#spC = sys.argv[4]
#spD = sys.argv[5]
#maxN = float(sys.argv[6]) # max prop. of N tolerated in the alignement
#nMin = int(sys.argv[7]) # minimum number of individuals in a given species at a given locus
#mu = float(sys.argv[8]) # mutation rate per generation and per bp
#rec = float(sys.argv[9]) # recombination rate per generation and per bp
#Nref = int(sys.argv[10]) # Nref : 100000


#### TEST
#fastafile='100_RNAseq.fasta'
#spA='E1'
#spB='W1'
#spC='W2'
#spD='W3'
#maxN=0.5
#nMin=4
#mu=0.000000001
#rec=0.00000001
#Nref=100000

###

data = {}
infile = parse(fastafile, 'fasta')
for i in infile:
	tmp = i.id.split('|')
	gene = tmp[0]
	species = tmp[1]
	individual = tmp[2]
	allele = tmp[3]
	
	if gene not in data:
		data[gene]={}
	
	if species not in data[gene]:
		data[gene][species]=[]
	
	if nameOut!='NA':
		if 'nameOut' not in data[gene]:
			data[gene]['nameOut']=[]
	if species==nameOut:
		data[gene]['nameOut'].append(i.seq)
	data[gene][species].append(i.seq)

outfile_infos = open('{datapath}/general_infos.txt'.format(datapath=datapath), 'w')
outfile_infos.write('gene\tdecision\tcomment\tlength\tn_treated_Sites\tn_valid_sites\tn_SNPs\n')

outfile_ms = open('{datapath}/observed_data.ms'.format(datapath=datapath), 'w')
outfile_ms.write('msnsam 4pop\n111 222 333\n\n')

bpfile_header = "#Nref:{0}\tmu:{1}\tr:{2}\t".format(Nref, mu, rec)
bpfile_L1=''# L
bpfile_L2=''# nA
bpfile_L3=''# nB
bpfile_L4=''# nC
bpfile_L5=''# nD
bpfile_L6=''# theta
bpfile_L7=''# rho
bpfile_L8=''# nSNPs 

nRetainedLoci = 0
for gene in data: # loop over genes
	bin_tmp = {}
	bin_tmp['spA'] = []
	bin_tmp['spB'] = []
	bin_tmp['spC'] = []
	bin_tmp['spD'] = []
	if nameOut!='NA':
		bin_tmp['nameOut'] = []
	
	nA = 0
	for seq in data[gene][spA]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nA += 1
			bin_tmp['spA'].append(seq_tmp)

	nB = 0
	for seq in data[gene][spB]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nB += 1
			bin_tmp['spB'].append(seq_tmp)

	nC = 0
	for seq in data[gene][spC]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nC += 1
			bin_tmp['spC'].append(seq_tmp)

	nD = 0
	for seq in data[gene][spD]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nD += 1
			bin_tmp['spD'].append(seq_tmp)

	if nameOut!='NA':
		for seq in data[gene]['nameOut']:
			seq_tmp = seq
			missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
			if missing_data < maxN:
				bin_tmp['nameOut'].append(seq_tmp)

	if nA>=nMin and nB>=nMin and nC>=nMin and nD>=nMin:
		res = fasta2ms(bin_tmp, region, nameOut)
		if res['nPos_valid']>=Lmin:
			nRetainedLoci += 1
			outfile_infos.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(gene, 'retained', 'ok', res['locus_length'], res['nPos'], res['nPos_valid'], res['nPos_biallelic']))
			outfile_ms.write(res['ms'])
			bpfile_L1 += '{0}\t'.format(res['nPos_valid'])
			bpfile_L2 += '{0}\t'.format(res['nA'])
			bpfile_L3 += '{0}\t'.format(res['nB'])
			bpfile_L4 += '{0}\t'.format(res['nC'])
			bpfile_L5 += '{0}\t'.format(res['nD'])
			bpfile_L6 += '{0}\t'.format(4*Nref*mu*res['nPos_valid'])
			bpfile_L7 += '{0}\t'.format(4*Nref*rec*res['nPos_valid'])
			bpfile_L8 += '{0}\t'.format(res['nPos_biallelic'])
		else:
			outfile_infos.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(gene, 'rejected', 'not enough valid sites', res['locus_length'], res['nPos'], res['nPos_valid'], res['nPos_biallelic']))
	else:
		missing_pop = []
		if nA<nMin:
			missing_pop.append(spA)
		if nB<nMin:
			missing_pop.append(spB)
		if nC<nMin:
			missing_pop.append(spC)
		if nD<nMin:
			missing_pop.append(spD)
		missing_pop = ','.join(missing_pop)
		outfile_infos.write('{0}\t{1}\t{2}\tNA\tNA\tNA\tNA\n'.format(gene, 'rejected', 'not enough sequences in ' + missing_pop))
outfile_ms.write('\n')
outfile_infos.close()
outfile_ms.close()

bpfile_header += 'nLoci:{nLoci}'.format(nLoci=nRetainedLoci)
bpfile = open('{datapath}/bpfile'.format(datapath=datapath), 'w')
bpfile.write(bpfile_header.strip() + '\n')
bpfile.write(bpfile_L1.strip() + '\n')
bpfile.write(bpfile_L2.strip() + '\n')
bpfile.write(bpfile_L3.strip() + '\n')
bpfile.write(bpfile_L4.strip() + '\n')
bpfile.write(bpfile_L5.strip() + '\n')
bpfile.write(bpfile_L6.strip() + '\n')
bpfile.write(bpfile_L7.strip() + '\n')
bpfile.write(bpfile_L8.strip() + '\n')

bpfile.close()

### get ABCstat.txt
if nameOut=='NA':
    folded='1'
else:
    folded='0'

commande = 'cat {datapath}/observed_data.ms | pypy {binpath}/mscalc_4pop.py datapath={datapath} simulationpath={simulationpath} folded={folded}'.format(binpath=binpath, datapath=datapath, simulationpath=simulationpath, folded=folded)
os.system(commande)


