# Genebank to COG function table
This python script reads Genebank files, extracts protein sequences, make rpsblast+ searches against COG_LE database and finally outputs a table of COG functional class vs number of genes.

# Requiremets
- Ncbi blast+ suite (rpsblast+)
- COG database, can be downloaded from [here](ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/)
- Python modules
 - sys, os, re, subprocess,multiprocessing
- Python third party modules:
 - pandas
 - Biopython

The COG database must be on the same folder as GB2COG.py and with the path Cog_LE/Cog*
  
# Usage
```
usage: GB2COG.py [-h] [-q QFILES [QFILES ...]] [-i PID] [-e EVALUE] [-c QC]

GB2COG extracts protein sequences and assigns COGs using rpsblast+.

optional arguments:
  -h, --help            show this help message and exit
  -q QFILES [QFILES ...], --QueryFiles QFILES [QFILES ...]
                        Query genomes in genbank format...
  -i PID, --ID PID      Percent ID for rpsblast+
  -e EVALUE, --evalue EVALUE
                        Evalue for rpsblast+
  -c QC, --queryCoverage QC
                        Query coverage for rpsblast+
      
```

# Output
The script generates three files per genome entry: .faa protein fasta file, .csv COG table and .cog table with COG functional class vs counts.
