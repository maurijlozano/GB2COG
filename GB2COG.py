#!/usr/bin/env python3
'''
GB2COG.py Extract Proteins from Genbank file.
M.J.Lozano 2020
'''
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



if len(sys.argv) < 2:
	print("Arguments required: space separated Genebank files...")
	exit()

ncpu=multiprocessing.cpu_count()
cdd2cog = pd.read_table('cdd2cog.txt',names=['CDD','COG','Gene','Description','',])
cog2fun = pd.read_table('cog2fun.txt',names=['COG','func','name'])

for filename in sys.argv[1:]:
	print('Processsing '+str(filename)+'...')
	basename = os.path.splitext(filename)[0]
	i=1
	for rec in SeqIO.parse(filename, "genbank"):
		n = 1
		faaName = basename+'_'+str(i)+'.faa'
		faaCOG = basename+'_'+str(i)+'.cog'
		f = open(faaName, '+w')
		if rec.features:
			for feature in rec.features:
				if feature.type == "CDS":
					if not 'pseudo' in feature.qualifiers:
						feature_name = feature.qualifiers["locus_tag"][0]
						feature_product = re.sub('[^a-z^A-Z^0-9^-]','_',feature.qualifiers['product'][0])
						feature_protid = feature.qualifiers['protein_id'][0]
						feature_seq = feature.qualifiers['translation'][0]
						f.write('>'+feature_name+'|'+feature_product+'|'+feature_protid+'\n'+str(feature_seq)+'\n')
						n+=1
		f.close()
		subprocess.call(['rpsblast+','-query', faaName, '-db', 'Cog_LE/Cog', '-out', faaCOG, '-evalue', '0.00001',  '-num_threads', str(int(ncpu)), '-max_target_seqs' ,'1', '-outfmt', '6'])
		blastRes = pd.read_table(faaCOG,names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
		blastRes = blastRes[['qseqid','sseqid']]
		blastRes = blastRes.groupby('qseqid', as_index=False).first()
		#
		#search cdd on cdd2cog
		blastRes['CDD'] = blastRes.apply(lambda x: re.sub(r'CDD:','',x["sseqid"]), axis=1)
		blastRes['COG'] = blastRes.apply(lambda x: CDD2COG(x), axis=1)
		blastRes['FUN'] = blastRes.apply(lambda x: COG2FUN(x), axis=1)
		blastRes.to_csv(basename+'_'+str(i)+'_Full_Tab.csv')
		found = len(blastRes)
		COGdistrb = blastRes.groupby(['FUN'], as_index=False).aggregate({'COG':'count'})
		COGdistrb = COGdistrb.append({'FUN':'Not Found','COG':n-found},ignore_index=True)
		COGdistrb.to_csv(faaCOG)
		i+=1
	

