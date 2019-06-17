Various binaries optimized for certain platforms. Here's an example workflow:

# Create a BURST database out of the NCBI representative genomes
## Get the files using the standard way
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O- > assembly_summary.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O- >> assembly_summary.txt
awk -F "\t" '($5 == "representative genome" || $5 == "reference genome") && $14=="Full" && $11=="latest"{print $20}' assembly_summary.txt >> ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths | sed 's/ftp:/https:/' > ftpfilepaths
cat ftpfilepaths | xargs -n 1 -P 16 wget -q --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 99
```

## Rename the files to get rid of the trailing characters in the names
`for f in *.fna.gz; do FN=${f/_/-}; FN=${FN/_*/}; mv $f ${FN/-/_}.fna.gz; done`

## Lingenome without stripping anything (names by filenames and works on compressed files):
`lingenome . Rep94.fasta FILENAME`

## Make Burst db for queries up to 320bp in length
`burst_linux_DB15 -r Rep94.fasta -o Rep94.edx -a Rep94.acx -d DNA 320 -i 0.95 -t 48 -s 1500`

# Taxonomy 
## Ensure latest tid2gg.txt, and sorted version:
```
mkdir taxtmp && cd taxtmp
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
t2gg nodes.dmp names.dmp tid2gg.txt SUBONLY
sort -k1,1 tid2gg.txt > tid2gg.srt.txt
cd ..
```

## Get raw taxonomy by filename
`for f in *.fna.gz; do echo "${f/.fna.gz/}"; done | grep -F -f - assembly_summary.txt | cut -f 1,6,7,8 | sort -k2,2 > rawtax.tax`

## Get GG string for each taxon
`join -t $'\t' -12 -21 -e0 -o'1.1,2.2,1.4,0,1.3' rawtax.srt.tax taxtmp/tid2gg.srt.txt | sort -k2 > alltax.txt`

