# BURST
BURST (formerly known as embalmer) is an optimal, high-speed pairwise sequence aligner specialized in aligning many NGS short reads against large reference databases.  

## Why
As next-generation DNA sequencing data emerges faster than computational power can keep up, approximate heuristic solutions to the fundamental DNA alignment/mapping problem are increasingly used. Paradoxically, it seems, the more data we have, the less accurate the alignment algorithms used to analyze it. Algorithms with perfect sensitivity and specificity acheivable under mismatch constraints have been neglected in favor of techniques promising speedier alignment at the cost of absolute alignment quality (under some metrics of precision/accuracy/sensitivity/recall). BURST returns to the roots of provably optimal alignment algorithms, reinvigorating them with speedups as high as millions-fold without sacrificing any alignment quality whatsoever in the default operating modes. 

## What
BURST is a truly, mathematically optimal high-throughput end-to-end short-read DNA aligner. It supports:
- gapped, end-to-end alignment of variable-length short-reads (up to a few thousand bases) against arbitrary reference sequences
- guaranteed (first) best hit, guaranteed *all* tied hits, or guaranteed *all hits over specified identity cutoff*
- dual objective scoring function: 1) find lowest edit distance, 2) find highest BLAST ID among all possible least-cost paths
- optional optimal LCA taxonomy assignment (or other hierarchical variable if supplied) with customizable confidence cutoff
- optional optimal minimization of the number of unique references hit (tie-breaking by minimization)
- full IUPAC ambiguous base support in queries and references (with option to penalize N's in references)
- ease-of-use: no installation (just download and run) with concise commandline interface
- fast speed of operation. It can align 12 million 292-bp microbial amplicon sequences against the Greengenes 13.8 97% reference database in under 10 minutes on a quad E7-4850v2 server, or a couple of hours on a dual-core 2013 Macbook Air. It also aligns shotgun reads against very large databases (like the ~20GB IMG annotated bacterial genes database); 0.99.2 can align at a rate of *over 10,000 100bp reads per second* against a 31.5GB subset of RefSeq complete genomes at 98% alignment identity on a single 32-core Ivy Bridge server. 
- optional heuristic modes allowing faster alignments at lower similarity scores which, if used at your own risk, completely invalidates everything I just said about optimality

## What not
BURST does not currently implement the following, although all these are planned in future releases:
- clustering (as an output mode; it can use clustering in making its database)
- custom scoring matrices (although it supports any alphabet, including proteins, with option -x)
- local alignment (only end-to-end alignment is supported like in bowtie2/bwa/usearch default operation)
- finding very low identity matches with longer query sequences (it stops counting after accruing ~250 mismatches)
- Paired-end unstitched alignments (although this can be performed downstream by aligning both pairs in ALLPATHS mode and finding the reference mapped to by both pairs with acceptable orientation/distance). 

## What now
- Further speed improvements are in the works. Each speed improvement is guaranteed (mathematically) never to sacrifice alignment quality, even of a single alignment. That is, in the default optimal mode of operation. All bets are off using the heuristic modes, naturally.

## How
See [Releases](https://github.com/knights-lab/burst/releases) page for precompiled binaries for a variety of systems with no dependencies. Basically, just download one of the files on the releases page appropriate for your system (Windows, Linux, or Mac) and run it on the command line. If on Windows, pick an ".exe" version; if on macOS pick a ".mac" version; if on Linux pick a ".linux" version. If the default version (burst.exe, burst.mac, burst.linux) doesn't work, try the corresponding version with ".older" in the name, and if that still doesn't work, try the one with ".buzzard." Please let me know if you can't get the program to run on your system. 

### Easiest (Not for large reference databases or long reference sequences such as gigabase-length eukaryotic genomes):
`burst -r myRefs.fasta -q myQueries.fasta -o myAlignments.b6`

### Fastest (short version):
1. Create database
```
burst -r MyDB.fasta -a MyDB.acx -o MyDB.edx -d DNA -s
```

2. Search 

The default search mode, CAPITALIST, reports the smallest set of references necessary to explain all tied hits:
```
burst -q myQueries.fasta -a MyDB.acx -r MyDB.edx -o output.txt
```

Note that burst can also report LCA taxonomy for each query sequence if taxonomy is provided with `-b MyDB.tax` (a tab-delimited taxonomy file where the first column contains the entire sequence header for each sequence in the original fasta file, and the second column contains semi-colon-separated taxonomy). In this case the command above becomes:

```
burst -q myQueries.fasta -a MyDB.acx -r MyDB.edx -b MyDB.tax -o output.txt
```

**BEST** mode (report first best hit):
```
burst -q myQueries.fasta -a MyDB.acx -r MyDB.edx -m BEST -b MyDB.tax -o output.txt
```

**ALLPATHS** mode (larger output file; report all ties for best hit for every query sequence):
```
burst -q myQueries.fasta -a MyDB.acx -r MyDB.edx -m ALLPATHS -b MyDB.tax -o output.txt
```


### Fastest (longer, more detailed version):
1. Decide on the maximum lengths your queries will be, and the minimum identity you require of qualifying alignments. For example, for max query length of 320 bases and minimum identity of 0.97 (97%), you'd pass "-d DNA 320" and "-i 0.97" like below. *Note: databases assuming shorter maximum query length and higher minimum identities will run faster. If you only have HiSeq 125-bp data and you're only interested in alignments of 98% identity or better, you'd want to use something like "-d DNA 125" and "-i 0.98" instead.*
2. Run `burst -r MyDB.fasta -d DNA 320 -o MyDB.edx -a MyDB.acx -s 1 -i 0.97` to generate a database and accelerator. 
3. (optional, advanced) To refine the database, you can specify `-f` when building to enable fingerprint clustering. Conversely, if you have insufficient memory to make a database using `-d DNA`, consider using `-dp 2` or higher (partitions ease memory use) or use the non-compressive database mode `-d QUICK`.
4. Use the database for all future alignments: `burst -r MyDB.edx -a MyDB.acx -q MyQueries.fasta -o myAlignments.b6`

Other alignment modes, taxonomy parsing, tie-reporting, etc:
- Using "-m CAPITALIST" (the default) enables unique-reference minimization (reducing the number of unique references hit; useful for OTU picking or taxonomy assignment). 
- Using "-m BEST" produces the single highest BLAST-id alignment possible in the database, breaking ties by choosing the very first occurrence (in input order) within the original input fasta database. 
  - This has interesting implications if there is meaning to the order of the references (ordered by increasing taxonomic specificity or sequence abundance in another sample or a depth-first traversal of a clustogram). Otherwise it ensures consistency of best hit for the same input sequence.  
- Using "-m ALLPATHS" reports all tied best-hit alignments above the chosen identity threshold. 
  - It is up to the user to parse the resulting alignments in a meaningful way, such as running EM interpolation on the possibilities, informing a more intelligent, optimal "MAPQ"-like score, resolving coverage ambiguities in downstream programs, etc
- Using "-m FORAGE" reports all alignments, best or otherwise, that exist under the chosen identity threshold.
  - This is like ALLPATHS but doesn't just report the ties. It reports *everything* below the threshold. Use for primer identification, variable region sleuthing, exhaustive enumeration of all possible sequence-level relationships among organisms, running the input queries as their own references to use as a distance matrix, etc.
- Using "-b taxonomy.txt" along with "-m CAPITALIST" enables optimal LCA (lowest common ancestor) taxonomy assignment. 
  - The taxonomy file is a simple file with two columns separated by tab. The first column contains the name of a reference in the references, and the second contains its corresponding taxonomy (or other hierarchical assignment), with levels separated by semicolon. 
  - You can assign taxonomy more conservatively or speculatively by changing the --taxacut (-bc) parameter. Lower values are more speculative, higher values are more conservative. (You may also provide a confidence to --taxacut directly, e.g. --taxacut 0.85)
- Use "-y" during database creation and subsequent alignment to disable penalizing the base 'N' in queries and references. 


## Where
Output alignments are stored in the resulting .b6 file. This is a tab-delimited text file in [BLAST-6 column format](http://www.drive5.com/usearch/manual/blast6out.html). Columns 11 and 12 instead refer to total edit distance (number of differences between query and reference in total) and whether the query is an exact duplicate of the query above it (1 if so), respectively. If taxonomy is assigned (-m CAPITALIST -b taxonomy.txt), that particular read's (interpolated if CAPITALIST) taxonomy is reported in column 13. 

To find the latest version of BURST, see [How](#how) above.

## Who
Please contact Gabe Al-Ghalith or Dan Knights* (I'm sure you can find our contact info!)

## Woah (troubleshooting)
1. *I downloaded the program for my system but it won't run! Says "Permission denied" or "command not found":*
If on Linux or Mac, you may have to run the command "chmod +x" on the program first, and then run the program inside of the directory that contains it using a dot and slash before the name (for example, on Linux: "./burst.linux" if the file "burst.linux" is within the current working directory of the terminal). Another solution is to add the directory containing the program to the system PATH. This technique may vary by operating system and terminal type. 

2. *All queries align with 100% identity, many to the same strange reference sequence:*
Uh oh, looks like your database contains long series of "N"s (ambiguous bases). Because all ambiguities are resolved according to IUPAC standards, N actually matches perfectly to anything. For example, the nucleotide "K" matches "Y" but not "M," although "M" matches "Y". Although this opens up exciting new possibilities for leveraging ambiguity in aligning to SNP-aware databases, psuedo-clusters, and more, a stretch of 300 N's present in some poorly-curated databases will match any length-300 query perfectly. Disable this behavior by passing -n or --npenalize, which will force N's to be treated as mismatches against A, C, G, or T/U in the query. N will still be considered a match to any ambiguous nucleotides in the query. 

3. *I get "segmentation fault" (or other crash):*
This is likely a bug with BURST! Please contact me with no less than the following and I'll try to fix it:
  - The exact command-line used to run the program
  - The version of BURST used (run with -h to see help)
  - The operating system and amount of RAM (memory) in the computer running it
  - A minimalistic example of input and output to reproduce the problem. If it occurs using a DB (.edx), include the fasta file used to produce the DB. 

4. *I get no alignments with my amplicon reads, even though I know they're legit:*
Try reverse complementing (`-fr`). If that doesn't work, try removing sequencing platform adaptors and cleaning up and trimming the reads with [a QC pipeline](https://github.com/knights-lab/shi7). 

5. *Other program(s) give me more alignments; how can you say this is optimal?:*
First, more alignments doesn't mean correct alignments. Second, be careful when comparing technologies; BURST is a short-read aligner. It does not do local alignment like "BLAST" and hence does not do soft-trimming -- this is very much intentional and part of ensuring optimality of end-to-end alignments. An alignment of identity 97% spanning 97% of a query means that query is actually 97% x 97% = ~94% identical to its matched reference throughout. 

6. *It won't compile!* It is not recommended to compile this software yourself unless you have the Intel compiler and a lot of patience for profile-guided optimization using multi-pass compilation. Your binary probably won't be as fast as the one provided on the release page. If you nonetheless insist, you must use ICC or GNU GCC (NOT Apple/LLVM CLANG) and provide the additional compiler flag `-march=corei7` (or newer) and `-fwhole-program -O3` (or `-Ofast`)

## Cite
Al-Ghalith, Gabriel and Dan Knights. BURST enables optimal exhaustive DNA alignment for big data. DOI 2017:doi.org/10.5281/zenodo.806850

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1037553.svg)](https://doi.org/10.5281/zenodo.806850)
(please cite using DOI until manuscript is published)
