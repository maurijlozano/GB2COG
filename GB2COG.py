#!/usr/bin/env python3
'''
GB2COG.py Extract Proteins from Genbank file.
M.J.Lozano 2020
'''
import argparse
import sys, os, re, subprocess,multiprocessing
import pandas as pd
from Bio import SeqIO

#functions
def CDD2COG(row):
	cdd = int(row['CDD'])
	cog = cdd2cog[cdd2cog['CDD'] == cdd].iloc[0,1]
	return cog

def COG2FUN(row):
	cog = row['COG']
	fun = cog2fun[cog2fun['COG'] == cog].iloc[0,1]
	fun=fun[0]
	return fun

def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='GB2COG extracts protein sequences and assigns COGs using rpsblast+..')
	parser.add_argument("-q", "--QueryFiles",help="Query genomes in genbank format...", dest="qfiles", action='append', nargs='+', required=False)
	parser.add_argument("-i", "--ID",help="Percent ID for rpsblast+...", dest="pid", action='store', required=False)
	parser.add_argument("-e", "--evalue",help="Evalue for rpsblast+...", dest="evalue", action='store', required=False)
	parser.add_argument("-c", "--queryCoverage",help="Query coverage for rpsblast+...", dest="qc", action='store', required=False)
	args = parser.parse_args()
	return args


args = parseArgs()
if not args.qfiles:
	print("Error: you must provide space separated Genebank files...")
	exit()

if not args.pid:
	pid = 0
else:
	pid = int(args.pid)

if not args.qc:
	qc = 0
else:
	qc = int(args.qc)
	
if not args.evalue:
	evalue = 0.001
else:	
	evalue = args.evalue

ncpu=multiprocessing.cpu_count()
cdd2cog = pd.read_table('cdd2cog.txt',names=['CDD','COG','Gene','Description','',])
cog2fun = pd.read_table('cog2fun.txt',names=['COG','func','name'])

for filename in args.qfiles[0]:
	print('Processsing '+str(filename)+'...')
	basename = os.path.splitext(filename)[0]
	i=1
	allCOGs=pd.DataFrame()
	for rec in SeqIO.parse(filename, "genbank"):
		n = 1
		description = rec.description
		if len(description) == 0:
			description = rec.name
		faaName = basename+'_'+str(i)+'.faa'
		faaCOG = basename+'_'+str(i)+'.cog'
		f = open(faaName, '+w')
		if rec.features:
			for feature in rec.features:
				if feature.type == "CDS":
					if not 'pseudo' in feature.qualifiers:
						feature_name = feature.qualifiers["locus_tag"][0]
						feature_product = re.sub('[^a-z^A-Z^0-9^-]','_',feature.qualifiers['product'][0])
						if 'protein_id' in feature.qualifiers:
							feature_protid = feature.qualifiers['protein_id'][0]
						else:
							feature_protid = ""
						feature_seq = feature.qualifiers['translation'][0]
						f.write('>'+feature_name+'|'+feature_product+'|'+feature_protid+'\n'+str(feature_seq)+'\n')
						n+=1
		f.close()
		subprocess.call(['rpsblast+','-query', faaName, '-db', 'Cog_LE/Cog', '-out', faaCOG, '-evalue', evalue,  '-num_threads', str(int(ncpu)), '-max_target_seqs' ,'1', '-outfmt', '6 qseqid sseqid pident qlen length evalue'])
		blastRes = pd.read_table(faaCOG,names=['qseqid','sseqid','pident','qlen','length','evalue'])
		blastRes = blastRes[blastRes['pident'] > pid]
		blastRes['qc'] = blastRes['qlen']/blastRes['length']*100
		blastRes = blastRes[blastRes['qc'] > qc]
		blastRes = blastRes[['qseqid','sseqid']]
		blastRes = blastRes.groupby('qseqid', as_index=False).first()
		#
		#search cdd on cdd2cog
		blastRes['CDD'] = blastRes.apply(lambda x: re.sub(r'gnl\|CDD\|','',x["sseqid"]), axis=1)
		blastRes['COG'] = blastRes.apply(lambda x: CDD2COG(x), axis=1)
		blastRes['FUN'] = blastRes.apply(lambda x: COG2FUN(x), axis=1)
		blastRes.to_csv(basename+'_'+str(i)+'_Full_Tab.csv')
		found = len(blastRes)
		COGdistrb = blastRes.groupby(['FUN'], as_index=False).aggregate({'COG':'count'})
		completeCOGS = pd.DataFrame(columns=['FUN','COG'])
		for k in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
			if any(COGdistrb['FUN']==k):
				completeCOGS = completeCOGS.append({'FUN':k,'COG':COGdistrb[COGdistrb['FUN']==k].iloc[0,1]},ignore_index=True)
			else:
				completeCOGS = completeCOGS.append({'FUN':k,'COG':0},ignore_index=True)
		#		
		completeCOGS = completeCOGS.append({'FUN':'Not Found','COG':n-found},ignore_index=True)
		completeCOGS = completeCOGS.set_index('FUN')
		completeCOGS = completeCOGS.rename(columns={'COG':description})
		completeCOGS.to_csv(faaCOG)
		#All in one table
		allCOGs = pd.concat([allCOGs,completeCOGS],axis=1)
		i+=1
	allCOGs.to_csv('COG_matrix.txt')

