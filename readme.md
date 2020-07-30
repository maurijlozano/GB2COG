#Genebank to COG function table
This python script reads Genebank files, extracts protein sequences, make rpsblast+ searches against COG_LE database and finally outputs a table of COG functional class vs number of genes.

#Requiremets
- Ncbi blast+ suite (rpsblast+)
- COG database, can be downloaded from (here)[ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/]
- Python modules
 - sys, os, re, subprocess,multiprocessing
- Python third party modules:
 - pandas
 - Biopython

The COG database must be on the same folder as GB2COG.py and with the path Cog_LE/Cog*
  
#Usage
```
Gb2COG.py <file1.gb> <file2.gb> ... <fileN.gb>
```

#Output
The script generates three files per genome entry: .faa protein fasta file, .csv COG table and .cog table with COG functional class vs counts.
