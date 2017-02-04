# embalmer
Where other optimal aligners go to die a slow death

## Why
As next-generation DNA sequencing data emerges faster than computational power can keep up, researchers are coming up with more and more approximate ("heuristic" or "guesswork") solutions to the fundamental DNA alignment problem. Paradoxically, the more data we have, the less accurate our alignments become as we sacrifice accuracy to the inherent biases of our approximation techniques. The dream of perfect sensitivity and specificity acheivable under mismatch constraints (sequencing/amplification noise, absence of perfectly-matching reference genes, etc) has long been forfeited in favor of techniques promising ever speedier alignment algorithms with misleading titles like "accurate" or "near optimal," when the reality shows each leap in speed comes with a loss of alignment quality (precision/accuracy/sensitivity/recall). EMBALMER returns to the roots of provably optimal alignment algorithms, reinvigorating them with speedups as high as millions-fold without sacrificing any alignment quality whatsoever. 

## What
EMBALMER is a truly, mathematically optimal high-throughput end-to-end short-read DNA aligner. It supports:
- gapped, end-to-end alignment of short-reads (up to a few thousand bases) against arbitrary reference sequences
- guaranteed (first) best hit, guaranteed *all* tied hits, or guaranteed *all hits over specified identity cutoff*
- dual objective scoring function: 1) find lowest edit distance, 2) find highest BLAST ID among all possible least-cost paths
- optional optimal LCA taxonomy assignment (or other hierarchical variable if supplied) with customizable confidence cutoff
- optional optimal minimization of the number of unique references hit (tie-breaking by minimization)
- full IUPAC ambiguous base support in queries and references (with option to penalize N's in references)
- ease-of-use: no installation (just download and run) with concise commandline interface
- fast speed of operation. It can align 12 million 292-bp microbial amplicon sequences against the Greengenes 13.8 97% reference database in under 10 minutes on a quad E7-4850v2 server, or a couple of hours on a dual-core 2013 Macbook Air. It also aligns shotgun reads against very large databases (like the ~20GB IMG annotated bacterial genes database), albeit at a slower speed: 45 minutes for 500,000 126-bp reads on the same server. 

## What not
EMBALMER does not currently implement the following, although all these are planned in future releases:
- clustering (as an output mode; it can use clustering in making its database)
- custom scoring matrices (although it supports any alphabet, including proteins, with option -x)
- local alignment (only end-to-end alignment is supported like in bowtie2/bwa/usearch default operation)
- finding very low identity matches with longer query sequences (it stops counting after accruing ~250 mismatches)

## What now
- Further speed improvements are in the works. Each speed improvement is guaranteed (mathematically) never to sacrifice alignment quality, even of a single alignment. 

## How
See [Releases](https://github.com/knights-lab/embalmer/releases) page for precompiled binaries for a variety of systems with no dependencies. Basically, just download one of the files on the releases page appropriate for your system (Windows, Linux, or Mac) and run it on the command line. If on Windows, pick an ".exe" version; if on macOS pick a ".mac" version; if on Linux pick a ".linux" version. If the default version (embalm.exe, embalm.mac, embalm.linux) doesn't work, try the corresponding version with ".older" in the name, and if that still doesn't work, try the one with ".buzzard." Please let me know if you can't get the program to run on your system. 

## Who
Please contact Gabe Al-Ghalith or Dan Knights* (I'm sure you can find our contact info!)

## Woah (troubleshooting)
If on Linux or Mac, you may have to run the command "chmod +x" on the program first, and then run the program inside of the directory that contains it with a dot (for example, on Linux: "./embalm.linux" if the file "embalm.linux" is within the current working directory of the terminal). Another solution is to add the directory containing the program to the system PATH. This technique may vary by operating system and terminal type. 
