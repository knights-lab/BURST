#define USAGE "Usage: LLsim input.lin.fna output.fna <numReads> <readLen> <numErrors>"
#define VERSION "0.92"
/* LLsim: High-performance short read simulator w/specified number of errors per read.
   INPUT: linearized FASTA file (alternating header and sequence)
   OUTPUT: FASTA of simulated reads with specified number of errors
   ARGS: numReads: how many reads to simulate from the references. 
         readLen: how long each read has to be (BEFORE potential additions/deletions)
		 numErrors: how many errors to introduce in the read relative to the reference
*/
// Note to Gabe: Add this to help string!

#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <time.h>
#include <sys/resource.h>

// Special random-number engine using primed xorshifts
static inline uint64_t LLRand64(uint64_t *x) {
	*x^=*x<<13; *x^=*x>>7; return *x^=*x<<17;}

// Convenience package to hold both a sequence ix and offset within it
typedef struct {uint64_t ix, off;} Ix_Off_pair;

// Set for modulo 32 ASCII values that aren't A-C-G-T-U
char bad[] = {1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1};

/* Function that indicates whether a sequence contains invalid (ambiguous)
   DNA characters. Y'know, for those of us who can't tolerate uncertainty. ;)
   INPUT: 
   	string containing DNA sequence (s), 
   	position within that sequence (off),
	length of the subsequence beginning at off (len)
   OUTPUT: boolean indicator of whether the span of s from off to len
    contains any invalid or ambiguous characters, post-masking.
*/
int isAmbig(char *s, uint64_t off, uint64_t len) {
	for (uint64_t i = off; i < off + len; ++i)
		if (bad[s[i] & 31]) return 1; //{puts("...BAD!"); return 1;}
	return 0;
}

/* Simple qcomp-compliant unsigned integer comparator. It's ugly at first
   sight, but the logic is simple -- if string less, return -1. Else return the
   result of a greater-than comparison (which itself can return only 0 or 1). 
   Thus a 0 result of this secondary test means 'equal' in this case since we 
   ruled out less-than on the LHS, and 1 means 'greater-than'. -1, 0, and 1
   coincide with the return codes required of a qcomp compliant comparator.
*/
static int u32cmp(const void *a, const void *b) {
	return *(uint32_t *)a < *(uint32_t *)b ? -1 : 
		*(uint32_t *)a > *(uint32_t *)b;
}

// Will hold buckets for the possible substitutions for each masked base
char *MUT[32] = {0};

/* Wide binary search function to return a specific entry with offset. 
   Intended to search through a list of numerical indices corresponding to
   starting locations of sequences. It finds the nearest preceding sequence
   and returns that along with how far into that sequence the key was.
   Because it's used internally over a list of known size, there is no need
   to check for cases 
   INPUT: 
    list of indices demarcating offsets (ixList)
    numeric key to look up in the list; can exist between valid offsets (key)
    length of the list of indices (range)
   OUTPUT: Ix_Off_pair (defined above) of: found base entry and offset
*/
static Ix_Off_pair uWBS(uint64_t *ixList, uint64_t key, uint64_t range) {
	uint64_t middle, low = 0, high = range;
	while (low <= high) {
		middle = low + ((high - low) >> 1);
		if (key > ixList[middle]) low = middle + 1;
		else if (key < ixList[middle]) high = middle - 1;
		else break; 
	}
	if (ixList[middle] > key) --middle;
	return (Ix_Off_pair){middle, key - ixList[middle]};
}

/* The brunt of the work happens in the main function. The required arguments are
   the commandline options specified in the usage information at the top of this file */
int main(int argc, char *argv[]) {
	puts("This is LLsim ["VERSION"] by Gabe and Lia");
	if (argc < 6) {puts(USAGE); exit(1);} // All arguments are required. See USAGE at top
	int numReads = atol(argv[3]), readLen = atol(argv[4]), numE = atol(argv[5]);
	if (numReads < 1 || readLen < 1 || numE < 0 || numE > readLen) {
		printf("Invalid read parameters: num %d, len %d, err %d\n",numReads,readLen,numE);
		exit(1);
	}
	
	// Prepare to directly read page table on disk. "struct stat" and "read" are 
	// POSIX opcalls. They work on all major operating systems.
	struct stat sb; int inf = open(argv[1],O_RDONLY); 
	// Unlike fopen, open produces the return code '-1' on failure condition.
	if (inf==-1) {printf("Cannot open input '%s'\n",argv[1]); exit(2);}
	fstat(inf,&sb); uint64_t sz_inf = sb.st_size; // opcall to get file size
	// Create a direct pagetable memory map to character representation of file
	char *mm_inf = mmap(0,sz_inf+4096,PROT_READ|PROT_WRITE,MAP_PRIVATE,inf,0);
	if (mm_inf[sz_inf-1] == 0) mm_inf[sz_inf-1] = '\n'; // support non-NL file end
	madvise(mm_inf,sz_inf+4096,MADV_SEQUENTIAL); // better caching for in-order reads
	if (!sz_inf) {puts("fasta file is empty. Exiting."); exit(1);}

	// Use standard I/O for output file with increased buffer length
	FILE *out = fopen(argv[2],"wb"); 
	if (!out) {printf("Cannot open output '%s'\n",argv[2]); exit(2);}
	setvbuf(out,0,_IOFBF,1<<20); // increase buffer

	// Store the byte offsets of all headers and sequences dynamically
	uint64_t cur_sz = 1000, ns = 0, total_bases = 0, msl = 0;
	char **Heads = malloc(cur_sz*sizeof(*Heads)), **Seqs = malloc(cur_sz*sizeof(*Seqs));
	if (!Heads || !Seqs) {puts("ERROR: Out of memory"); exit(3);}
	char *mm_ptr = mm_inf;
	while (*mm_ptr == '>') {  // continue until out of headers
		if (ns >= cur_sz) {   // dynamically resize data structures
			Heads = realloc(Heads,(cur_sz+=cur_sz)*sizeof(*Heads)),
			Seqs = realloc(Seqs,cur_sz*sizeof(*Seqs));
			if (!Heads || !Seqs) {puts("ERROR: Out of memory [2]"); exit(3);}
		}
		Heads[ns] = mm_ptr;  // current pointer must be header
		while (*mm_ptr && *mm_ptr != '\n') ++mm_ptr; // proceed to EOL
		Seqs[ns] = ++mm_ptr; // current pointer (post-EOL) must be sequence
		while (*mm_ptr && *mm_ptr != '\n') ++mm_ptr; // proceed to EOL
		uint64_t sl = mm_ptr - Seqs[ns]; // compute sequence length
		msl = sl > msl ? sl : msl;       // update max sequence length
		total_bases += sl;               // tally total number of bases
		++mm_ptr, ++ns;                  // proceed post-EOL, bump num seqs
	} 
	// Shrink data structures back down to exact (nl) size
	Heads = realloc(Heads,ns*sizeof(*Heads));
	Seqs = realloc(Seqs,ns*sizeof(*Seqs));
	printf("There were %llu sequences identified [max len %llu]\n",ns,msl);
	if (msl < readLen) {puts("ERROR: max ref length shorter than desired query length"); exit(1);}

	// Create list of sequence offsets for later wide binary search
	uint64_t *Offsets = malloc((ns+1)*sizeof(*Offsets));
	if (!Offsets) {puts("ERROR: Out of memory [3]"); exit(3);}
	*Offsets = 0; 
	for (uint64_t i = 1; i < ns; ++i) 
		Offsets[i] = Offsets[i-1] + (Heads[i] - Seqs[i-1] - 1);
	Offsets[ns] = total_bases;
	// Note to Gabe: place this in previous loop to avoid secondary LOC [TODO]

	// Preset mutation table (MUT is globally defined previously)
	MUT[0] = (uint8_t[]){'A','C','G','T'}; // insert (all bases accepted)
	MUT[1] = (uint8_t[]){'C','G','T'}; //A
	MUT[3] = (uint8_t[]){'A','G','T'}; //C
	MUT[7] = (uint8_t[]){'A','C','T'}; //G
	MUT[20]= (uint8_t[]){'A','C','G'}; //T
	MUT[21]= (uint8_t[]){'A','C','G'}; //U

	// Goal: simulate random number out of the valid bases, get ref and ix in it
	// If subsequence is valid (long enough, unambiguous), simulate from it, else reject
	// Note to Gabe: consider domain restriction instead of rejection sampling [TODO]
	uint64_t seed = 1; // TODO: Replace this with time(0) in production code?
	Ix_Off_pair loc;   // this will hold the found index and offset into the reference
	char *sr = calloc(2*readLen+1,1); // guarantee enough space for max (all-ins) expansion
	uint32_t *mIX = malloc(readLen*sizeof(*mIX)), *mShf = malloc(readLen*sizeof(*mShf));

	if (!sr || !mIX || !mShf) {puts("Memory error. Couldn't allocate data"); exit(3);}
	for (int i = 0; i < readLen; ++i) mShf[i] = i; // prepare linear permutation indices
	for (uint64_t i = 0; i < numReads; ++i) {
		uint64_t rand;
		do rand = LLRand64(&seed) % total_bases,
			loc = uWBS(Offsets, rand, ns-1);        // sample genetic region
		while (readLen + rand >= Offsets[loc.ix + 1] || 
			isAmbig(Seqs[loc.ix],loc.off,readLen)); // rejection conditions
		
		// Determine which indices (relative to reference) to target for errors
		for (int j = 0; j < numE; ++j) // sample numE indices from read
			mIX[j] = (LLRand64(&seed) % (readLen-j)) + j; 
		register uint32_t r,t;         // init temp variables (Knuth deck shuffle)
		for (int j = 0; j < numE; ++j) // sample without replacement (permute)
			t = mShf[j], r = mIX[j],
			mShf[j] = mShf[r], mShf[r] = t; // only need to sample numE cards/ixs
		// sort the indices so copying can proceed in order of bases in reference
		qsort(mShf,numE,sizeof(*mShf),u32cmp); // Get them in order
		// Note to Gabe: avoid expensive qsort syscall with sort network
		// ... or avoid sorting altogether in favor of select-N (small num E's)
		
		// Types of mutation: substitution, insertion, deletion. 
		// Randomly determine each mutation type by index 0-4. 
		// 0-2 = substitute. 3 = insert, 4 = delete
		uint32_t qix = 0, mix = 0, ni = 0, nd = 0; 
		char *qp = Seqs[loc.ix] + loc.off, typeS; // store ref pointer
		fwrite(Heads[loc.ix],1,Seqs[loc.ix]-Heads[loc.ix]-1,out); // write header
		//fprintf(out, ">%u [%u]", loc.ix, loc.off);
		fprintf(out, " @%u: ",loc.off+1); // embed original location into header
		fwrite(qp,1,readLen,out);         // also embed original sequence itself
		fprintf(out, " ");
		int li = 0, ld = 0;
		for (uint32_t j = 0; j < numE; ++j) {
			// Copy in all bases up to the current error index
			for (; qix + ni < mShf[j]; ++qix) 
				sr[mix++] = qp[qix];
			int type = LLRand64(&seed) % 5; // select type of mutation
			// Restrict successive adjacent indels (instead, select a substitution)
			if (j && mShf[j] == mShf[j-1] + 1 && (type == 3 && li) || (type == 4 && ld))
				type = LLRand64(&seed) % 3;
			if (type < 3) li = ld = 0, typeS = 'S', // substitution
				sr[mix++] = MUT[qp[qix++] & 31][type];
			else if (type == 3) ld = 1, li = 0, typeS = 'D', ++nd, ++qix; // deletion
			else li = 1, ld = 0, typeS = 'I', ++ni, // insertion
				sr[mix++] = MUT[0][LLRand64(&seed) % 4];
			fprintf(out, "%c%u", typeS, mShf[j]);   // add mutation type+loc to header
		}
		for (; qix < readLen; ++qix)    // copy in all bases after last edit index
			sr[mix++] = qp[qix];
		fprintf(out, "\n");             // move to sequence line
		fwrite(sr,1,readLen+ni-nd,out); // dump sequence
		fprintf(out,"\n");
	}
}

/* TODO: 
1. Add proper help menu (and real commandline parser) [easy]: 0%
2. Add ability to perform RC [trivial]: 0%
3. Detect and filter cases where errors nullify each other [hard]: 50%
*/