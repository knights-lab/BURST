Readme_utils
This folder contians taxonomy, database, and tree-making utilities (built for 64-bit linux). 
All tools are written in C by Gabe A (send complaints and blame to algh0022@umn.edu)

This readme is presented as a walkthrough and doesn't cover all the functionality of these tools.
To learn more, see usage (run without parameters) and play around with them!

Database processing tools
-----------------
Example of getting rep84: 
# Viruses
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt -O- > assembly_summary.txt
awk -F "\t" '$11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths

# Archaea and Bacteria
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O- > assembly_summary.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O- >> assembly_summary.txt
awk -F "\t" '($5 == "representative genome" || $5 == "reference genome") && $14=="Full" && $11=="latest"{print $20}' assembly_summary.txt >> ftpdirpaths

# Consolidate and download (xargs used for parallel download/extraction and also to enable entire refseq download)
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
cat ftpfilepaths | xargs -n 1 -P 16 wget -q --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 99
printf '%s\0' *.gz | xargs -r0 -n 1 -P 16 gunzip

# Compile genome fragments into whole genomes, identify and separate plasmids, and rename headers to newick-compatible format
lingenome . AllGenomes.fasta AllPlasmids.fasta HEADFIX

# Note that, rarely, some plasmids and chromosomes are mixed up by NCBI (though not in rep to our knowledge).
# To avoid stripping plasmids, omit "AllPlasmids.fasta" from the lingenome command above

Taxonomy processing
-----------------
## Only do these first two steps once per database revision.
# Step 1: Get ncbi taxID (tid) to strain-level delimited string mapping
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
t2gg nodes.dmp names.dmp tid2gg.txt

# Step 2: Generate accession to tid mapping
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
gunzip *.gz
tail -n+2 nucl_gb.accession2taxid >> nucl_wgs.accession2taxid
a2gg_make tid2gg.txt nucl_wgs.accession2taxid a2gDB

# Step 3: Assign complete taxonomy string to any linearized fasta containing NCBI accession!
a2gg_parse /path/to/AllGenomes.fasta a2gDB outPrefix ">" "." FULL
# a2gg_parse is portable and relatively resource-light. It only requires a2gDB and a fasta file.


Tree tools
------------------
## Create a good genome-content Nei-Saitou long-branch-attraction-corrected Newick tree (for UniFrac, etc):
akmer94 AllGenomes.fasta AllGenomes.tre 12 ADJ DIRECT TREE
# The above settings take about a day on my personal machine. To speed it up at the cost of quality, you can do: 
# akmer94 AllGenomes.fasta AllGenomes.tre 11 HEUR4 ADJ DIRECT TREE

## To instead do de novo phylogeny from amplicon data, use combined_seqs.fna instead of AllGenomes.fasta:
akmer94 AllGenomes.fasta AllGenomes.tre HEUR ADJ DIRECT TREE
# (replace HEUR with HEUR4 or HEUR9 to speed up tree creation and reduce RAM use)
# To create clusters/rep set (arguably unnecessary to cluster if you have the full tree!):
# use your favorite tree-cutting tool to cut the tree at the desired height to create subtrees
# Choose the medpoint of each cut subtree as the representative sequence

## To instead measure distance between amplicon SAMPLES (rather than sequences) for direct denovo beta div:
# You'll need one fasta file per sample. QIIME's split_fasta can do this for a combined_seqs.fna
# Rarefy your samples to uniform depth (just call, eg., "head -1000" on each sample fasta)
lingenome mySampleFastas/ AllSamples.fasta HEADFIX
akmer94 AllSamples.fasta AllSamples.tre ADJ DIRECT TREE


Misc tools
----------------
"embalmulate" turns BURST output (b6 file) into two tables: an OTU table (or specific genome table) and interpolated taxonomy table.
"linfasta" is a much stripped-down version of lingenome. It operates on individual fasta files only to linearize them, nothing more.
