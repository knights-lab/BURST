echo "Searching with BURST using ALLPATHS..."
time /project/flatiron2/sop/burst12 -q seqs.fna -r /project/flatiron2/sop/gg97.edx -a /project/flatiron2/sop/gg97.acx -o gg97_allpaths.txt -i 0 -m ALLPATHS -b /project/flatiron2/sop/gg97.tax --taxacut 1000000

echo "Search with BURST using CAPITALIST..."
time /project/flatiron2/sop/burst12 -q seqs.fna -r /project/flatiron2/sop/gg97.edx -a /project/flatiron2/sop/gg97.acx -o gg97_capitalist.txt -i 0 -m CAPITALIST -b /project/flatiron2/sop/gg97.tax --taxacut 1000000

# build BT2 database and search
echo "Building BT2 database..."
mkdir bt2
/project/flatiron2/sop/bowtie2-build /project/flatiron/dan/gg_13_8/gg_13_8_otus/rep_set/97_otus.fasta bt2/97_otus --threads 96

echo "Running BT2..."
time /project/flatiron2/sop/bowtie2 -x bt2/97_otus -f seqs.fna --np 1 --rfg 0,1 --rdg 0,1 --mp 1,1 --score-min "L,0,-1" -S bt2.sam --no-head

echo "tallying hits..."
python tally_correct_hits_by_pct_id.py gg97_allpaths.txt gg97_capitalist.txt /project/flatiron/dan/gg_13_8/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt bt2.sam tally.txt
