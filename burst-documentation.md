# BURST (Better Ultrafast Reconciliation of Sequences Tool) Documentation

## Overview

BURST is an optimal, high-speed pairwise sequence aligner specialized in aligning many NGS short reads against large reference databases. It is designed to provide mathematically optimal alignments while maintaining exceptional speed, making it particularly useful for metagenomic studies and large-scale sequence analysis projects.

## Version

This documentation covers BURST version 1.0+.

## Features

1. Optimal end-to-end alignment of variable-length short reads (up to a few thousand bases) against arbitrary reference sequences
2. Gapped alignment support
3. Multiple alignment modes:
   - BEST: Report first best match by hybrid BLAST id
   - ALLPATHS: Report all ties with the same error profile
   - CAPITALIST: Minimize set of references AND interpolate taxonomy (default)
   - FORAGE: Report all matches above specified threshold
   - ANY: Report any valid hit above specified threshold
4. Optional optimal LCA taxonomy assignment with customizable confidence cutoff
5. Full IUPAC ambiguous base support in queries and references
6. Accelerator mode for faster alignment using k-mer hashing
7. Fingerprinting for additional filtering of potential matches
8. Support for reverse complement alignment
9. Database creation and management tools
10. Multithreading support for improved performance

## Usage

```
burst {options}
```

### Basic parameters:

- `--references (-r) <name>`: FASTA/edx DB of reference sequences [required]
- `--accelerator (-a) <name>`: Creates/uses a helper DB (acc/acx) [optional]
- `--queries (-q) <name>`: FASTA file of queries to search [required if aligning]
- `--output (-o) <name>`: Blast6/edb file for output alignments/database [required]

### Behavior parameters:

- `--forwardreverse (-fr)`: Also search the reverse complement of queries
- `--whitespace (-w)`: Write full query names in output (include whitespace)
- `--xalphabet (-x)`: Allow any alphabet and disable ambiguity matching
- `--nwildcard (-y)`: Allow N,X to match anything (in query and reference)
- `--taxonomy (-b) <name>`: Taxonomy map (to interpolate, use -m CAPITALIST)
- `--mode (-m) <name>`: Pick an alignment reporting mode (BEST, ALLPATHS, CAPITALIST, FORAGE, ANY)

### Performance parameters:

- `--dbpartition (-dp) <int>`: Split DB making into <int> chunks (lossy)
- `--taxacut (-bc) <num>`: Allow 1/<int> rank discord OR % conf
- `--taxa_ncbi (-bn)`: Assume NCBI header format '>xxx|accsn...' for taxonomy
- `--skipambig (-sa)`: Do not consider highly ambiguous queries (5+ ambigs)
- `--taxasuppress (-bs) [STRICT]`: Suppress taxonomic specificity by %ID
- `--id (-i) <decimal>`: Target minimum similarity (range 0-1)
- `--threads (-t) <int>`: How many logical processors to use
- `--shear (-s) [len]`: Shear references longer than [len] bases
- `--fingerprint (-f)`: Use sketch fingerprints to precheck matches (or cluster db)
- `--prepass (-p) [speed]`: Use ultra-heuristic pre-matching
- `--heuristic (-hr)`: Allow relaxed comparison of low-id matches
- `--noprogress`: Suppress progress indicator
- `--qbunch (-qb) <int>`: Pack QBUNCH with queries <int> divergent
- `--qbunch_max (-qm) <int>`: Max size of QBUNCH
- `--quickforage (-qf)`: Output FORAGE'd results inline
- `--cache (-c) <int>`: Performance tweaking parameter
- `--latency (-l) <int>`: Performance tweaking parameter

## Alignment Modes

1. BEST: Reports the first best match based on hybrid BLAST id.
2. ALLPATHS: Reports all ties with the same error profile.
3. CAPITALIST: Minimizes the set of references and interpolates taxonomy (default mode).
4. FORAGE: Reports all matches above the specified threshold.
5. ANY: Reports any valid hit above the specified threshold.

## Database Creation

BURST can create custom databases for faster alignment:

```
burst -r input_references.fasta -d [DNA|RNA|QUICK] [max_query_length] -o output_database.edx
```

Options:
- DNA/RNA: Creates a full database
- QUICK: Creates a faster, but potentially less sensitive database
- max_query_length: Optional parameter to specify the maximum expected query length

## Accelerator Creation

Note: BURST's accelerator formats are hard-coded for either prefixes of size 12 or 15. The version you're using is displayed in the BURST help string. The smaller size-12 prefix uses less memory but is slower (and is hence suitable for marker gene analysis). 

To create an accelerator file for even faster alignments:

```
burst -r input_references.fasta -d [options] -a output_accelerator.acx -o output_database.edx
```

## Input File Requirements

### Reference Sequences
- FASTA format
- Can be provided as raw FASTA or as a pre-built BURST database (.edx)

### Query Sequences
- FASTA or FASTQ format
- Gzipped input supported
- Maximum sequence length: 100MB (configurable)

### Taxonomy Mapping (optional)
- Tab-delimited text file
- Columns: sequence name, taxonomy string
- Taxonomy strings are semicolon-delimited

## Output Format

BURST outputs alignments in a modified BLAST-6 column format:

1. Query sequence name
2. Reference sequence name
3. Percent identity
4. Alignment length
5. Number of mismatches
6. Number of gap openings
7. Query start position
8. Query end position
9. Subject start position
10. Subject end position
11. E-value (set to -1 in BURST)
12. Bit score (used for other purposes in BURST)
13. Taxonomy (LCA-based, if provided and using CAPITALIST mode)

## Performance Considerations

1. Use the accelerator (-a) option for faster alignments on large databases
2. Increase the number of threads (-t) to utilize multiple CPU cores
3. Adjust the cache (-c) and latency (-l) parameters for fine-tuning performance
4. Use the fingerprint (-f) option for additional filtering of potential matches
5. Consider using the prepass (-p) option for ultra-fast, heuristic pre-matching

## Limitations

1. Local alignment is not supported (only end-to-end alignment)
2. Custom scoring matrices are not implemented
3. Paired-end unstitched alignments are not directly supported

## Error Handling

BURST includes basic error checking for input file formats and command-line arguments. It will print error messages and exit if it encounters issues like malformed input files or invalid options.

## Examples

1. Basic alignment:
   ```
   burst -r references.fasta -q queries.fasta -o alignments.b6 -i 0.97
   ```

2. Create a database and accelerator:
   ```
   burst -r references.fasta -d DNA 320 -a references.acx -o references.edx
   ```

3. Align using a pre-built database with taxonomy:
   ```
   burst -r references.edx -a references.acx -q queries.fasta -b taxonomy.txt -o alignments.b6 -m CAPITALIST
   ```

## Conclusion

BURST provides a powerful and flexible tool for optimal sequence alignment, particularly suited for metagenomic studies and large-scale sequence analysis projects. Its various output options and optimization features make it suitable for a wide range of applications, from simple best-hit reporting to complex taxonomic assignment tasks.
