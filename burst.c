/* BURST aligner -- fast optimal aligner by Gabe. 
Copyright (C) 2015-2018 Knights Lab, Regents of the University of Minnesota.
This software is released under the GNU Affero General Public License (AGPL) v3.0.
*/
#define VER "v0.99.7LL"
#define _LARGEFILE_SOURCE_
#define FILE_OFFSET_BITS 64
#include <stdio.h>
#include <inttypes.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
// Ensure ability to read files > 4GB on Windows
#ifdef _WIN32
	#include <windows.h>
	#define fseek _fseeki64
	#define ftell _ftelli64
#endif
// OpenMP is the multi-threading protocol used by BURST
#ifdef _OPENMP
	#include <omp.h>
	#define HAS_OMP 1
	#define getNumberOfCores() omp_get_max_threads()
#else
	#define HAS_OMP 0
	#define omp_get_thread_num() (0)
	#define getNumberOfCores() (1)
	#define omp_get_wtime() ((double)clock()/CLOCKS_PER_SEC)
#endif
// This is the header that includes the assembly intrinsic functions
#include <immintrin.h>
#ifdef __AVX__
	#define ASMVER "AVX-128"
#elif __SSE4_1__
	#define ASMVER "SSE4.1"
#else
	// These redefinitions of SSE4+ intrinsics are necessary on all other (non-SSE4)
	// instruction sets. This allows functions that use SSE4 functions to still
	// run, albeit with slight performance penalty.
	#define _mm_popcnt_u64 __builtin_popcountll
	#define _mm_blendv_epi8(f,t,m) \
		_mm_xor_si128(_mm_and_si128(_mm_xor_si128(f,t),m),f) 
	#define _mm_blendv_ps(f,t,m) \
		_mm_xor_ps(_mm_and_ps(_mm_xor_ps(f,t),m),f)
	#define _mm_cvtepu8_epi32(x) \
		_mm_unpacklo_epi16(_mm_unpacklo_epi8(x,_mm_setzero_si128()),_mm_setzero_si128())
	#define _mm_extract_epi8(a,x) \
		(0xff & (_mm_extract_epi16(a,((x)>>1))>>(8*((x) & 0x01))))
	#define _mm_stream_load_si128 _mm_load_si128
	#ifdef __SSSE3__
		#define ASMVER "SSSE3"
	#else
		// With SSE2 (no longer officially supported), even more redefinitions
		// are necessary to perform assembly functions including blends and lt/gt
		#define ASMVER "SSE2"
		#define _mm_cmple_epu16(x,y) \
			_mm_cmpeq_epi16(_mm_subs_epu16(x, y), _mm_setzero_si128())
		#define _mm_blendv_si128 (x, y, mask) \
			_mm_or_si128(_mm_andnot_si128(mask, x), _mm_and_si128(mask, y))
		#define _mm_min_epu16(x,y) \
			_mm_blendv_si128(y, x, _mm_cmple_epu16(x, y))
		#define _mm_lddqu_si128 _mm_loadu_si128
	#endif
#endif
#include "math.h"

// This enum type defines all possible running modes in -m. These
// are tested later in functions that change their behavior based on
// which alignment mode is used (see do_alignments()).			
typedef enum {
	FORAGE, BEST, ALLPATHS, CAPITALIST, ANY,
	MATRIX, PAIRED, MAPPED, INLINE, COMMUNIST
} Mode;
// This enum type defines all of the database creation modes supported by BURST
typedef enum {DNA_16, DNA_8, DNA_4, PROT, DATA, QUICK} DBType;
Mode RUNMODE = CAPITALIST; // The default alignment run mode is now CAPITALIST; (-m)
uint32_t DO_PREPASS = 0;   // Prepass heuristic is not enabled by default; (-p)
uint32_t LATENCY = 16;     // The max length difference allowed between seqs within a DB cluster
int THREADS = 1;           // By default, BURST uses one thread, but this is overriden by OpenMP env vars; (-t)
int cacheSz = 150;         // How many rows of the alignment matrix to remember. More = more RAM bandwidth
int Xalpha = 0;            // Whether to align with arbitrary symbol exact matching; (-x)
int REBASE = 0;            // Size of shears to use when creating a database
int DO_FP = 0;             // Whether to use fingerprints for DB creation (or aligning); (-f)
int DO_ACCEL = 0;          // Whether to create or use an accelerator (DB creation/alignment); (-a)
int DO_HEUR = 0;           // Prevents demotion of sequences to non-accelerator alignment by identity; (-hr)

uint32_t TAXACUT = 10;     // Ratio (to 1) of taxonomic agreements required for CAPITALIST interpolation; (-bc)
float THRES = 0.97f;       // Identity/similarity threshold to restrict valid alignments; (-i)
long REBASE_AMT = 500, DB_QLEN = 500; // Stores shear length and max query length (DB creation); (-s); (-d [] qLen)
#define VECSZ 16
#ifndef SCOUR_N
#define SCOUR_N 15
#endif
#define SCOUR_R (32 - 2*SCOUR_N)
#define SCOUR_L 4
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define PRINT_USAGE() { \
	printf("\nBURST aligner (" VER "; DB%d version)\n", SCOUR_N); \
	printf("Compiled with " ASMVER " and%s multithreading\n", !HAS_OMP ? " WITHOUT" : "");\
	printf("\nBasic parameters:\n");\
	printf("--references (-r) <name>: FASTA/edx DB of reference sequences [required]\n");\
	printf("--accelerator (-a) <name>: Creates/uses a helper DB (acc/acx) [optional]\n"); \
	puts("--queries (-q) <name>: FASTA file of queries to search [required if aligning]");\
	puts("--output (-o) <name>: Blast6/edb file for output alignments/database [required]");\
	puts("\nBehavior parameters:"); \
	puts("--forwardreverse (-fr): also search the reverse complement of queries"); \
	puts("--whitespace (-w): write full query names in output (include whitespace)"); \
	/*puts("--npenalize (-n): Force A,C,G,T,U,N,X in query to mismatch X,N in reference");*/ \
	puts("--xalphabet (-x): Allow any alphabet and disable ambiguity matching");\
	puts("--nwildcard (-y): Allow N,X to match anything (in query and reference)"); \
	puts("--taxonomy (-b) <name>: taxonomy map (to interpolate, use -m CAPITALIST)");\
	puts("--mode (-m) <name>: Pick an alignment reporting mode by name. Available modes:");\
	puts("  BEST (report first best match by hybrid BLAST id)"); \
	puts("  ALLPATHS (report all ties with same error profile)"); \
	puts("  CAPITALIST (minimize set of references AND interpolate taxonomy) [default]"); \
	puts("  COMMUNIST (min. ref. set, taxa interpolation, and fair redistribution)"); \
	puts("  FORAGE (report all matches above specified threshold)"); \
	puts("  ANY (report any valid hit above specified threshold)"); \
	/*puts("  [not enabled]: PAIRED, MAPPED, INLINE");*/ \
	puts("--makedb (-d) [name qLen]: Create a database from input references");\
	printf("  [name]: Optional. Can be DNA, RNA, or QUICK [QUICK]\n");\
	printf("  [qLen]: Optional. Max query length to search in DB [%d]\n",DB_QLEN); \
	printf("\nPerformance parameters:\n"); \
	printf("--dbpartition (-dp) <int>: Split DB making into <int> chunks (lossy) [%u]\n",1); \
	printf("--taxacut (-bc) <num>: allow 1/<int> rank discord OR %% conf; 1/[%u]\n",TAXACUT); \
	printf("--taxa_ncbi (-bn): Assume NCBI header format '>xxx|accsn...' for taxonomy\n"); \
	printf("--skipambig (-sa): Do not consider highly ambiguous queries (5+ ambigs)\n"); \
	printf("--taxasuppress (-bs) [STRICT]: Suppress taxonomic specificity by %%ID\n"); \
	printf("--id (-i) <decimal>: target minimum similarity (range 0-1) [%.2f]\n",THRES);\
	printf("--threads (-t) <int>: How many logical processors to use [%u]\n",THREADS);\
	printf("--shear (-s) [len]: Shear references longer than [len] bases [%ld]\n",REBASE_AMT);\
	/*printf("--unique (-u): Dereplicate references (lossless preprocessing)\n");*/ \
	printf("--fingerprint (-f): Use sketch fingerprints to precheck matches (or cluster db)\n"); \
	printf("--prepass (-p) [speed]: use ultra-heuristic pre-matching\n"); \
	printf("  [speed]: Optional. Integer, maximum search effort [16]\n"); \
	printf("--heuristic (-hr): allow relaxed comparison of low-id matches\n"); \
	printf("--noprogress: suppress progress indicator\n"); \
	/*printf("--cache (-c) <int>: Performance tweaking parameter [%d]\n",cacheSz); */ \
	/* printf("--latency (-l) <int>: Performance tweaking parameter [%d]\n",LATENCY); */ \
	/*printf("--clustradius (-cr) <int>: Performance tweaking parameter [auto]\n");*/ \
	puts("\n--help (-h): Shows this help screen with version info"); \
	puts("\nExample: burst -r myRefs.fasta -q myQs.fasta -o outputs.txt -i 0.98"); \
	puts("\nLicensed under the GNU Affero General Public License v3.0"); \
	exit(1); \
}

// Gap penalty is worth 1 mismatch. This is not alterable currently, as any other
// value would violate later optimality edit distance assumptions.
#define GAP 1
// The width of a vector (in bytes) on the target machine. This is also immutable.
#define SCD 16
// Internal handler for database versions, to allow recognition. EDB is obsolete.
#define EDB_VERSION 2
#define EDX_VERSION 3
// There are two accelerator versions. One stores 16 million sheared references,
// and the other (big/LARGE-format) stores up to 256 million sheared references.
#define ACC_VERSION 0
#define ACC_VERSION_BIG 1
char Z = 1; // Whether to penalize N's (and also the penalty itself); (-n)
// BAK and RVC are obsolute debug-only reverse complement handlers
char *BAK = "\0ACGTNKMRYSWBVHD\0"; // 
char *RVC = "\0TGCANMKYRSWVBDH\0";
char RVT[] = {0,4,3,2,1,5,7,6,9,8,10,11,13,12,15,14}; // for reverse complementing
//			 0,T,G,C,A,N,M,K,Y,R,S, W, V, B, D, H
__m128i GAPV;
__m128i SCOREFAST[16] = {0};
char SCORENVedN[] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //.
				   //0  1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
				   //.  A C G T N K M R Y S W B V H D
					 -1,0,1,1,1,0,1,0,0,1,1,0,1,0,0,0, //A
					 -1,1,0,1,1,0,1,0,1,0,0,1,0,0,0,1, //C
					 -1,1,1,0,1,0,0,1,0,1,0,1,0,0,1,0, //G
					 -1,1,1,1,0,0,0,1,1,0,1,0,0,1,0,0, //T/U
					 -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //N/X
					 -1,1,1,0,0,0,0,1,1,1,1,1,0,1,1,0, //K
					 -1,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1, //M
					 -1,0,1,0,1,0,1,1,0,1,1,1,1,0,1,0, //R
					 -1,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1, //Y
					 -1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1, //S
					 -1,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0, //W
					 -1,1,0,0,0,0,0,1,1,0,0,1,0,1,1,1, //B
					 -1,0,0,0,1,0,1,0,0,1,0,1,1,0,1,1, //V
					 -1,0,0,1,0,0,1,0,1,0,1,0,1,1,0,1, //H
					 -1,0,1,0,0,0,0,1,0,1,1,0,1,1,1,0, //D
					};
char CHAR2NUM[128], CHAR2IX[32] = {0,1,12,2,15,5,5,3,14,5,5,6,5,7,5,5, 
								   5,5,8,10,4,4,13,11,5,9,5,0,0,0,0,0};

typedef union {
	__m128i v;
	__m128 vf;
	char c[16];
	uint8_t u8[16];
	uint16_t u16[8];
	uint32_t u32[4];
	uint64_t u64[2];
	int8_t i8[16];
	int16_t i16[8];
	int32_t i32[4];
	int64_t i64[2];
	float f[4];
	double d[2];
} DualCoil;

typedef struct {
	DualCoil sc, sh;
} SSMat;

typedef struct {
	float score[16];
} Scores16;

typedef struct {uint32_t v, i;} Split; 
typedef struct {char *s; uint32_t ix;} StrIxPair;
typedef struct {char *s; uint64_t ix;} StrIxPair64;

typedef struct {
	float score[16];
	uint32_t finalPos[16];
	uint8_t numGapR[16], numGapQ[16];
} MetaPack;

typedef union {
	uint8_t a[32];
	uint16_t h[16];
	uint32_t s[8];
	uint64_t w[4];
} Prince;

// Structure to store sequence fingerprints
typedef struct {
	Prince *P;   // Prince
	uint8_t *N;  // Pops
	void *initP; 
	uint32_t nf, *Ptrs;
} PackaPrince;

// One hash contains many references (accelerator lookup)
typedef struct {
	uint32_t len, cap;
	uint32_t *Refs;
} Accelerant;
typedef struct AccelNode AccelNode;
struct AccelNode {
	AccelNode *next;
	uint32_t ref;
};

// For sorting references by property
typedef struct {
	char *seq;
	uint32_t len, ix;
} Tuxedo;

// For intra-taxonomic cutoffs, use these similarity values
typedef struct {char *Head, *Tax;} TaxPair_t;
                           //K     P     C     O     F     G     S     SS+
float TAXLEVELS_STRICT[] = {.65f, .75f, .78f, .82f, .86f, .94f, .98f, .995f},
	 TAXLEVELS_LENIENT[] = {.55f, .70f, .75f, .80f, .84f, .93f, .97f, .985f},
	 *TAXLEVELS = TAXLEVELS_LENIENT;

// Primary storage unit for each unique query
typedef struct {
	char *Seq;    //8
	uint16_t div; //2
	uint32_t six; //4
	uint8_t rc;   //1
} UniBin; //15

// Primary storage unit for shared query (RC) attributes
typedef struct {
	uint32_t len; //4
	uint16_t ed;  //2
} ShrBin; //6

// Contains all info about references; for passing between functions
typedef struct {
	char **RefHead, **RefSeq;
	uint32_t *RefLen, *ClumpLen, *RefStart, *RefIxSrt, *TmpRIX, *RefDedupIx;
	uint32_t totR, origTotR, numRclumps, maxLenR;
	DualCoil **RefClump, **ProfClump;
	
	// Fingerprint tools
	Prince *Centroids;
	PackaPrince FingerprintsR;
	TaxPair_t *Taxonomy;  // To store taxon and index
	uint32_t taxa_parsed;
	uint32_t clustradius; // Cluster refinement passes
	uint32_t numRefHeads, *RefMap; // For easy reference name lookup

	uint32_t *BadList, badListSz;  // Stores ambiguous references in acx
	void **Accelerators;           // Used for constructing acx
	DBType dbType;                 // Corresponds to second arg of '-d'
	int cparts, skipAmbig;
} Reference_Data;

// Contains all info about queries; for passing between functions
typedef struct {
	char **QHead, *SeqDumpRC; 
	uint32_t totQ, numUniqQ, maxLenQ, *Offset, *QBins, 
		maxED, maxDiv, minLenQ;
	int incl_whitespace, taxasuppress, rc, skipAmbig, quiet;

	PackaPrince FingerprintsQ;
	UniBin *UniBins;
	ShrBin *ShrBins;
} Query_Data;

// Aligned memory allocation (for vectorized operations)
void * malloc_a(size_t algn, size_t size, void **oldPtr) {
    uintptr_t mask = ~(uintptr_t)(algn - 1);
	*oldPtr = malloc(size+algn-1);
    return (void *)(((uintptr_t)*oldPtr+algn-1) & mask);
}
// Aligned (zerod) memory allocation (for vectorized operations)
void * calloc_a(size_t algn, size_t size, void **oldPtr) {
    uintptr_t mask = ~(uintptr_t)(algn - 1);
	*oldPtr = calloc(size+algn-1,1);
    return (void *)(((uintptr_t)*oldPtr+algn-1) & mask);
}
// For sorting (comparing) structures by strings ('s') they contain
int cmpStrIx(const void *a, const void *b) {
	return strcmp(((StrIxPair*)a)->s,((StrIxPair*)b)->s); }

// A generic function that sorts by strings in parallel
// Linear-time aggregation of common prefices, then N*logN sort within each
// To use, pass in Pack, len, and define structure type, num letters per block, and nib function
#define PARALLEL_SORT_PROTOTYPE(STRUCT_TYPE, STR_MEMBER, NUMLET, NIBFUNC) { \
	static const uint32_t NLB = 1 << (NUMLET*4); \
	STRUCT_TYPE **Bins = calloc(NLB,sizeof(*Bins)); \
	uint32_t *BcountO = calloc(1+NLB,sizeof(*BcountO)), \
		*Bcount = BcountO + 1; \
	_Pragma("omp parallel for") \
	for (uint32_t i = 0; i < len; ++i) { \
		char *s = Pack[i].STR_MEMBER; \
		uint32_t nib = NIBFUNC; \
		_Pragma("omp atomic update") \
		++Bcount[nib]; \
	} \
	for (uint32_t i = 0; i < NLB; ++i) if (Bcount[i]) \
		Bins[i] = malloc(Bcount[i]*sizeof(*Bins[i])), \
		Bcount[i] = 0; \
	_Pragma("omp parallel for") \
	for (uint32_t i = 0; i < len; ++i) { \
		char *s = Pack[i].STR_MEMBER; \
		uint32_t nib = NIBFUNC; \
		uint32_t ix; \
		_Pragma("omp atomic capture") \
		ix = Bcount[nib]++; \
		Bins[nib][ix] = Pack[i]; \
	} \
	int cmpFunc(const void *a, const void *b) { \
		return strcmp(((STRUCT_TYPE *)a)->STR_MEMBER + NUMLET, \
		((STRUCT_TYPE *)b)->STR_MEMBER + NUMLET); \
	} \
	_Pragma("omp parallel for schedule(dynamic,1)") \
	for (uint32_t i = 0; i < NLB; ++i) if (Bins[i]) \
		qsort(Bins[i],Bcount[i],sizeof(*Bins[i]),cmpFunc); \
	--Bcount; \
	for (uint32_t i = 2; i <= NLB; ++i) Bcount[i] += Bcount[i-1]; \
	_Pragma("omp parallel for schedule(dynamic,1)") \
	for (uint32_t i = 0; i < NLB; ++i) { \
		uint32_t init = Bcount[i], ed = Bcount[i+1]; \
		for (uint32_t pt = init; pt < ed; ++pt) \
			Pack[pt] = Bins[i][pt - init]; \
		free(Bins[i]); \
	} \
	free(Bins), free(BcountO); \
}
// For initial sorting/splitting of strings into prefix buckets
#define NIB4 s[0] << 12 | s[1] << 8 | s[2] << 4 | s[3]
#define NIB5 s[0] << 16 | s[1] << 12 | s[2] << 8 | s[3] << 4 | s[4]

// Some instantiations of the generic sort for different structs
static inline void parallel_sort_strpack(StrIxPair *Pack, uint32_t len) 
	PARALLEL_SORT_PROTOTYPE(StrIxPair, s, 5, NIB5)
static inline void parallel_sort_unibin(UniBin *Pack, uint32_t len) 
	PARALLEL_SORT_PROTOTYPE(UniBin, Seq, 5, NIB5)
static inline void parallel_sort_tuxedo(Tuxedo *Pack, uint32_t len) {
	int tuxCmp(const void *a, const void *b) {
		Tuxedo *A = (Tuxedo *)a, *B = (Tuxedo *)b;
		char *s1 = A->seq, *s2 = B->seq;
		uint32_t len1 = A->len, len2 = B->len, ml = MIN(len1,len2);
		uint32_t i = 5; 
		for (; i < ml; ++i) if (s1[i] != s2[i]) break;
		if (i < ml) return s1[i] - s2[i];
		return len1 < len2 ? -1 : 1;
	}
	void tuxSort(Tuxedo *a, size_t n, size_t s) 
		{qsort(a,n,s,tuxCmp);}
	#define qsort(a,n,s,c) tuxSort(a,n,s)
	PARALLEL_SORT_PROTOTYPE(Tuxedo, seq, 5, NIB5)
	#undef qsort
}

// Taxonomy lookup by binary search
char NULLTAX[1] = {0};
static char * taxa_lookup_generic(char *key, uint32_t sz, TaxPair_t *Dict) {
	TaxPair_t *p = Dict;
	while (sz) {
		uint32_t w = sz >> 1; 
		char *ref_s = p[w+1].Head, *key_s = key;
		
		while (*ref_s == *key_s++) if (!*ref_s++) return p[w+1].Tax; 
		if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	char *ref_s = p->Head, *key_s = key;
	while (*ref_s == *key_s++) if (!*ref_s++) return p->Tax;
	return NULLTAX; //return p->Tax; 
}
// Specialized version of the above for NCBI assembly formats ('>xxx|accsn...')
static char * taxa_lookup_ncbi(char *key, uint32_t sz, TaxPair_t *Dict) {
	TaxPair_t *p = Dict;
	int in = 0;
	while (sz) {
		uint32_t w = sz >> 1; 
		char *ref_s = p[w+1].Head, *key_s = key + 4;
		while (*ref_s == *key_s++) if (!*ref_s++) return p[w+1].Tax; 
		if (*(key_s-1) == '.' && !*ref_s) return p[w+1].Tax;
		char keyLast = *(key_s-1) == '.' ? 0 : *(key_s-1);
		if (*ref_s < keyLast) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	char *ref_s = p->Head, *key_s = key + 4;
	while (*ref_s == *key_s++) if (!*ref_s++) return p->Tax;
	if (*(key_s-1) == '.' && !*ref_s) return p->Tax;
	return NULLTAX; //return p->Tax; 
}

// Dynamic function pointer for taxa lookup; based on whether '-bn' specified
char * (*taxa_lookup)(char *, uint32_t, TaxPair_t *) = taxa_lookup_generic;

// Parses tab-delimited taxonomy file line-by-line into Tax_Pair_t object 'Obj'
// (stores header and taxonomy string separately). Returns number parsed.
size_t parse_taxonomy(char *filename, TaxPair_t **Obj) {
	static const size_t linelen = 10000000; //1k entries, 10 meg line lines
	size_t cur_sz = 1000, ns = 0;
	FILE *file = fopen(filename,"rb");
	if (!file) 
		{fprintf(stderr,"Cannot open TAXONOMY file: %s.\n",filename); exit(2);}
	TaxPair_t *T = malloc(cur_sz * sizeof(*T));
	char *line = calloc(linelen,1), *lineO = line; 
	while (line = fgets(line, linelen, file)) {
		if (ns == cur_sz) { // double all data structures
			T = realloc(T, (cur_sz*=2)*sizeof(*T));
			if (!T) {fputs("OOM:T:parse_taxonomy\n",stderr); exit(3);}
		}
		uint32_t i, j; 
		for (i = 0; line[i] != '\t'; ++i) 
			if (!line[i]) {fprintf(stderr,"ERROR: invalid taxonomy [%u]\n",ns); exit(2);}
		char *temp = malloc(i+1);
		if (!temp) {fputs("OOM:temp:parse_taxonomy\n",stderr); exit(3);}
		memcpy(temp,line,i), temp[i] = 0;
		T[ns].Head = temp;
		for (j = ++i; line[j] && line[j] != '\n' && line[j] != '\r' && line[j] != '\t'; ++j); 
		temp = malloc(j - i + 1);
		if (!temp) {fputs("OOM:temp:parse_taxonomy\n",stderr); exit(3);}
		memcpy(temp,line + i, j-i), temp[j-i] = 0;
		T[ns].Tax = temp;
		//printf("[%u] Head = %s\nTax=%s\n",ns,T[ns].Head,T[ns].Tax);
		++ns;
	}
	free(lineO);
	*Obj = realloc(T, ns*sizeof(*T)); // shrink T
	fclose(file);
	return ns;
}

// Standard line-by-line fasta sequence parser. Grows a bunch of parallel
// arrays containing headers, sequences, and their lengths. Returns number parsed.
// Used in standard reference database creation. No longer used for query parsing
size_t parse_tl_fasta(char * filename, char ***HeadersP, char ***SeqsP, uint32_t **LengthsP) {
	static const size_t linelen = INT32_MAX; //1k entries, 2-gig lines
	FILE *file = fopen(filename,"rb");
	if (!file) { 
		fprintf(stderr,"Cannot open FASTA file: %s.\n",filename); exit(2); }
	size_t cur_sz = 1000, ns = -1, len16;
	char **Headers = malloc(cur_sz * sizeof(*Headers)),
		 **Seqs = malloc(cur_sz * sizeof(*Seqs));
	uint32_t *Lengths = malloc(cur_sz * sizeof(*Lengths));
	char *line = malloc(linelen), *lineO = line; 
	if (!Headers || !Seqs || !Lengths || !line) {
		fputs("OOM:parse_tl_fasta\n",stderr); exit(3); }
	int lastHd = 0;
	while (line = fgets(line, linelen, file)) {
		size_t len = strlen(line); // Kill newlines
		if (len && line[len-1] == '\n') --len;
		if (len && line[len-1] == '\r') --len;
		line[len] = 0;
		switch (*line) {
			case '>':  // We could be in the (a) header.
				if (lastHd) break;
				if (ns++ == cur_sz) { // double all data structures
					cur_sz += cur_sz;
					Headers = realloc(Headers, cur_sz*sizeof(*Headers));
					Seqs = realloc(Seqs, cur_sz*sizeof(*Seqs));
					Lengths = realloc(Lengths, cur_sz*sizeof(*Lengths));
					if (!Headers || !Seqs || !Lengths) {
						fputs("OOM:parse_tl_fasta\n",stderr); exit(3); }
				}
				lastHd = 1;
				Headers[ns] = memcpy(malloc(len), line+1, len);
				Lengths[ns] = 0; 
			case '\0': case ' ': break;
			default: // we're in sequence (hopefully!)
				lastHd = 0;
				uint32_t a = Lengths[ns] + len;
				len16 = a + (15 & (16 - (a & 15)));
				if (!Lengths[ns]) Seqs[ns] = malloc(len16+1);
				else Seqs[ns] = realloc(Seqs[ns],len16+1); 
				memcpy(Seqs[ns] + Lengths[ns],line,len);
				Lengths[ns] += len;
				memset(Seqs[ns]+Lengths[ns],'\0',len16-Lengths[ns]+1); // trailing nil
		}
	}
	if (lastHd) puts("WARNING: file ends on header. Skipping last sequence."), --ns;
	free(lineO);
	Headers = realloc(Headers,++ns*sizeof(*Headers)), Seqs = realloc(Seqs,ns*sizeof(*Seqs));
	Lengths = realloc(Lengths,ns*sizeof(*Lengths));
	*HeadersP = Headers; *SeqsP = Seqs; *LengthsP = Lengths;
	if (ns >= UINT32_MAX) puts("WARNING: >4 billion sequences processed.");
	return ns;
}

// Variation of the above that internally arranges sequences contiguously in
// memory. Instead of sequences themselves, an offset pointer is returned.
// This specialized function is intended for database creation which requires
// packed references for vectorization and deduplication purposes.
size_t parse_tl_fasta_db(char *ref_FN, char ***HeadersP, char **SeqsP, 
 char ***OffsP, uint64_t *numLet) {
	size_t linelen = INT32_MAX; 
	size_t cur_len = linelen, cur_sz = 1000;
	FILE *file = fopen(ref_FN,"rb");
	if (!file) 
		{ fprintf(stderr,"Cannot open FASTA file: %s.\n",ref_FN); exit(2); }
	char **Headers = malloc(cur_sz * sizeof(*Headers)),
		 *Seqs = malloc(cur_len * sizeof(*Seqs));
	char **Offsets = malloc(cur_sz * sizeof(*Offsets));
	uint64_t cur_off = -1, ns = -1;
	char *line = malloc(linelen), *lineO = line; 
	if (!line) {fputs("OOM:ptfdLn\n",stderr); exit(3);}
	int lastHd = 0;
	while (line = fgets(line, linelen, file)) {
		size_t len = strlen(line); // kill newlines
		if (line[len-1] == '\n') --len;
		if (line[len-1] == '\r') --len;
		line[len] = 0;
		switch (*line) {
			case '>': // We could be in the (a) header.
				if (lastHd) break;
				if (++ns == cur_sz) { // double all data structures
					cur_sz += cur_sz;
					Headers = realloc(Headers, cur_sz*sizeof(*Headers));
					Offsets = realloc(Offsets, cur_sz*sizeof(*Offsets));
					if (!Headers || !Offsets) {
						fputs("OOM:db_parse",stderr); exit(3); }
				}
				lastHd = 1;
				Headers[ns] = memcpy(malloc(len), line+1, len);
				Offsets[ns] = (char*)++cur_off; // insert null post-seq 
			case '\0': case ' ': break;
			default: // we're in sequence (hopefully!)
				lastHd = 0;
				if (cur_off + len + 1 >= cur_len) { //
					Seqs = realloc(Seqs, (cur_len*=2)*sizeof(*Seqs));
					if (!Seqs) {fputs("OOM:db_parse",stderr); exit(3);}
				} // guaranteed to be enough -- but what about a while loop
				// now dump the current sequence bits into the trunk
				//update cur_off
				memcpy(Seqs + cur_off,line,len+1); //+1 for null
				cur_off += len;
		}
	}
	// add tail offset at ++cur_off
	// reallocate data structures with one extra space (tail)
	if (lastHd) puts("WARNING: file ends on header. Skipping last sequence."), --ns;
	free(lineO);
	Headers = realloc(Headers, (++ns+1)*sizeof(*Headers));
	Offsets = realloc(Offsets, (ns+1)*sizeof(*Offsets));
	Seqs = realloc(Seqs,++cur_off+16); // ensure alignment of tail
	if (!Headers || !Offsets || !Seqs) {
		fputs("Error OOM terminus p1",stderr); exit(3); }
	memset(Seqs+cur_off,'\0',15); 
	for (uint64_t i = 0; i < ns; ++i)
		Offsets[i] = (uint64_t)Offsets[i] + Seqs;
	Offsets[ns] = Seqs + cur_off;
	Headers[ns] = 0; // just a tailcap
	
	// return the proper data now
	*HeadersP = Headers; *SeqsP = Seqs; *OffsP = Offsets;
	*numLet = cur_off;
	return ns;
} 
// Finds the next newline character. Only for when the sequence
// is guaranteed not to terminate before next newline.
// Accessory function for the fast fasta parser below.
static inline char * findNL_abs(char *a) {
	__m128i nl = _mm_set1_epi8('\n');
	while (1) {
		__m128i v = _mm_lddqu_si128((void*)a);
		__m128i e = _mm_cmpeq_epi8(v,nl);
		uint16_t x = _mm_movemask_epi8(e);
		if (x) return a + __builtin_ctz(x);
		a += 16;
	}
}
// General version of the above that stops also at nulls.
static inline char * findNL_or_eol(char *a) {
	__m128i nil = _mm_set1_epi8(0), nl = _mm_set1_epi8('\n');
	while (1) {
		__m128i v = _mm_lddqu_si128((void*)a);
		__m128i e = _mm_cmpeq_epi8(v,nl);
		uint16_t x = _mm_movemask_epi8(e);
		if (x) return a + __builtin_ctz(x);
		__m128i n = _mm_cmpeq_epi8(v,nil);
		uint16_t z = _mm_movemask_epi8(n);
		if (z) return a + __builtin_ctz(z);
		a += 16;
	}
}

// Faster fasta parser intended for large query files. Reads in entire
// queries file, then indexes into it directly to access queries.
size_t parse_tl_faster(char * filename, char ***HeadersP, char ***SeqsP, uint32_t **LengthsP) {
	FILE *file = fopen(filename,"rb");
	if (!file) { 
		fprintf(stderr,"Cannot open FASTA file: %s.\n",filename); exit(2); }
	uint64_t sz = 0; 
	fseeko(file,0,SEEK_END); sz = ftello(file); rewind(file);
	double wt = omp_get_wtime();
	char *dump = malloc(16+sz);
	uint64_t bytesRead = fread(dump, 1, sz, file);
	memset(dump+sz,0,16);
	if (*dump != '>') {fputs("ERROR: Malformatted FASTA file.\n",stderr); exit(1);}
	wt = omp_get_wtime();
	uint32_t numNL = 0, numLT = 0;
	#pragma omp parallel for simd reduction(+:numNL,numLT)
	for (uint64_t i = 0; i < sz; ++i) 
		numNL += dump[i]=='\n',
		numLT += dump[i]=='>'; 
	numNL += numNL & 1;
	if (numLT != numNL/2) {fputs("ERROR: line count != '>' * 2\n",stderr); exit(1);}

	char **Headers = malloc(numLT * sizeof(*Headers)),
		 **Seqs = malloc(numLT * sizeof(*Seqs));
	uint32_t *Lengths = malloc(numLT * sizeof(*Lengths));
	if (!Headers || !Seqs || !Lengths) {
		fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }\
	wt = omp_get_wtime();
	Headers[0] = dump+1;
	char *nl = findNL_abs(Headers[0]); 
	*nl = 0;
	Seqs[0] = nl + 1;
	nl = findNL_or_eol(Seqs[0]); 
	*nl = 0;
	Lengths[0] = nl - Seqs[0];
	uint64_t ix = 1;
	#pragma omp parallel for
	for (uint64_t i = (nl-dump); i < sz; ++i) if (dump[i]=='>') {
		uint32_t tix;
		#pragma omp atomic capture
		tix = ix++;
		char *s = dump + i + 1;
		Headers[tix] = s;
		char *nl = findNL_abs(s); 
		*nl = 0;
		if (*(nl-1) == '\r') *(nl-1) = 0;
		s = nl + 1;
		Seqs[tix] = s;
		nl = findNL_or_eol(s); 
		*nl = 0;
		if (*(nl - 1) == '\r') *--nl = 0;
		Lengths[tix] = nl - s;
	}
	if (ix > UINT32_MAX) puts("WARNING: greater than 4 billion queries");
	*HeadersP = Headers; *SeqsP = Seqs; *LengthsP = Lengths;
	return ix;
}

/// Prototypes for nucleotide-nucleotide scoring functions. 
// An alphabet-agnostic scorer (for "xalpha" mode) penalizing inequality
// Because of SSE vectorization constraints, must flip from "true" (-1)
// to zero-penalty (0) with overflow, and false (0) to penalty (1)
#define DIAGSC_XALPHA _mm_add_epi8(_mm_cmpeq_epi8(_mm_set1_epi8(qLet), \
	rChunk),_mm_set1_epi8(1))
// SSSE3 shuffle is used as a single-instruction score lookup per letter
#ifdef __SSSE3__
	#define DIAGSC_MAT16 _mm_shuffle_epi8(SCOREFAST[qLet], rChunk)
#else
	// If SSSE3 isn't available, operate on pre-computed score profiles
	#define DIAGSC_MAT16 _mm_load_si128((void*)(profile + (x-1)*SCD + qLet))
#endif

// Scoring function for hybrid edit-distance BLAST-id sequence comparison.
// Intended to "re-score" within an already-determined set of initial alignments,
// it is intended to break ties by BLAST-id as well as populate the alignment 
// result with more statistics about the alignment (number of gaps, etc).
// The setup of the alignment also differs from primary scorers (further below) 
// in that it is not designed to allow 'seeking' into the alignment matrix, so
// it optimizes calculation from the first row.
#define RESCOREM_PROTYPE(DIAG_FUNC) {\
	uint32_t y, x; \
	--query; --ref; ++qlen; /* ++rwidth; */ \
	/* __m128i maxEDv = _mm_set1_epi8(maxED+1 < 255 ? maxED+1 : 255); for checking if BAD via score >= maxED+1 */ \
	__m128i maxEDv = _mm_set1_epi8(maxED+1); \
	DualCoil *restrict prevSc = Matrix + width, *restrict prevSh = Shifts + width, \
		*restrict prevShR = ShiftR + width, *restrict curSc = Matrix, *restrict curSh = Shifts, \
		*restrict curShR = ShiftR; \
	uint32_t LB = 1, HB = rwidth, LBN = 1, HBN = rwidth; \
	{ /* Iteration 1 only */ \
		_mm_store_si128((void*)(curSh),_mm_setzero_si128()); \
		_mm_store_si128((void*)(curShR),_mm_set1_epi8(1)); \
		_mm_store_si128((void*)(curSc),_mm_set1_epi8(1)); /* column 0 */ \
		char qLet = query[1]; \
		for (x = 1; x < rwidth; ++x) {\
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); /* refRow += 16; */ \
			__m128i curRow_x_1 = _mm_load_si128((void*)(curSc+(x-1))); \
			__m128i score = DIAG_FUNC; \
			/* test: if I'm a 1 and the left is a zero, give me a shift of 1 else 0 */ \
			__m128i getshiftL = _mm_and_si128(_mm_cmpeq_epi8(score,_mm_set1_epi8(1)), \
				_mm_cmpeq_epi8(curRow_x_1,_mm_setzero_si128())); \
			__m128i shift = _mm_and_si128(getshiftL,_mm_set1_epi8(1)); /* its left shift will be 1 */ \
			_mm_store_si128((void*)(curSc+x),score); \
			_mm_store_si128((void*)(curSh+x),shift); \
			_mm_store_si128((void*)(curShR+x),_mm_setzero_si128()); \
		} \
	} \
	for (y=2; y < qlen; ++y) { \
		LB = LBN, HB = HBN; \
		LBN = 0; \
		char qLet = query[y]; \
		DualCoil *temp = curSc; curSc = prevSc; prevSc = temp; \
		temp = curSh; curSh = prevSh; prevSh = temp; \
		temp = curShR; curShR = prevShR; prevShR = temp; \
		_mm_store_si128((void*)(curSh),_mm_setzero_si128()); \
		__m128i newMin = _mm_set1_epi8(MIN(y,255)); \
		_mm_store_si128((void*)(curSc),newMin); /* col 0 */ \
		_mm_store_si128((void*)(curShR),newMin); \
		for (x = LB; x < HB; ++x) { \
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); /* refRow += 16; */ \
			__m128i prevRow_x = _mm_load_si128((void*)(prevSc+x)); \
			__m128i prevShf_x = _mm_load_si128((void*)(prevSh+x)); \
			__m128i prevShfR_x = _mm_load_si128((void*)(prevShR+x)); \
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prevSc+(x-1))); \
			__m128i prevShf_x_1 = _mm_load_si128((void*)(prevSh+(x-1))); \
			__m128i prevShfR_x_1 = _mm_load_si128((void*)(prevShR+(x-1))); \
			__m128i curRow_x_1 = _mm_load_si128((void*)(curSc+(x-1))); \
			__m128i curShf_x_1 = _mm_load_si128((void*)(curSh+(x-1))); \
			__m128i curShfR_x_1 = _mm_load_si128((void*)(curShR+(x-1))); \
			\
			__m128i diagSc = DIAG_FUNC; \
			__m128i scoreOld = prevRow_x_1; \
			__m128i shift = prevShf_x_1; \
			__m128i shiftR = prevShfR_x_1; \
			__m128i score = _mm_adds_epu8(scoreOld, diagSc); \
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); \
			__m128i shiftU = prevShf_x; \
			__m128i shiftRU = _mm_adds_epu8(prevShfR_x, _mm_set1_epi8(1)); \
			__m128i scoreM = _mm_min_epu8(scoreU,score); \
			__m128i shiftM = _mm_min_epu8(shiftU,shift); \
			__m128i shiftU_le_shift = _mm_cmpeq_epi8(shiftM,shiftU); \
			__m128i scoreU_eq_score = _mm_cmpeq_epi8(scoreU,score); \
			__m128i scoreM_eq_score = _mm_cmpeq_epi8(score,scoreM); /* O <= U */ \
			__m128i tiebreak = _mm_andnot_si128(shiftU_le_shift,scoreU_eq_score); \
			__m128i condition = _mm_andnot_si128(tiebreak,scoreM_eq_score); \
			shift = _mm_blendv_epi8(shiftU,shift,condition); \
			shiftR = _mm_blendv_epi8(shiftRU,shiftR,condition); \
			score = scoreM; \
			\
			/* consider L */ \
			__m128i scoreLold = curRow_x_1; \
			__m128i shiftLold = curShf_x_1; \
			__m128i shiftRLold = curShfR_x_1; \
			__m128i scoreL = _mm_adds_epu8(scoreLold,_mm_set1_epi8(GAP)); \
			__m128i shiftL = _mm_adds_epu8(shiftLold,_mm_set1_epi8(1)); \
			__m128i shiftRL = shiftRLold; \
			scoreM = _mm_min_epu8(scoreL,score); \
			shiftM = _mm_min_epu8(shiftL,shift); \
			__m128i shiftL_le_shift = _mm_cmpeq_epi8(shiftM,shiftL); \
			__m128i scoreL_eq_score = _mm_cmpeq_epi8(scoreL, score); \
			scoreM_eq_score = _mm_cmpeq_epi8(score,scoreM); /* O <= U */ \
			tiebreak = _mm_andnot_si128(shiftL_le_shift,scoreL_eq_score); \
			condition = _mm_andnot_si128(tiebreak,scoreM_eq_score); \
			\
			shift = _mm_blendv_epi8(shiftL,shift,condition); \
			shiftR = _mm_blendv_epi8(shiftRL,shiftR,condition); \
			score = scoreM; \
			\
			/* Do bounds handling */ \
			__m128i anyBad = _mm_cmpeq_epi8(maxEDv,_mm_min_epu8(maxEDv,score)); \
			score = _mm_or_si128(anyBad,score); \
			if (_mm_movemask_epi8(anyBad) != 0xFFFF) { \
				if (!LBN) LBN = x; \
				HBN = x; \
			} \
			_mm_store_si128((void*)(curSc+x),score); \
			_mm_store_si128((void*)(curSh+x),shift); \
			_mm_store_si128((void*)(curShR+x),shiftR); \
		} \
		if (!LBN) { \
			printf("\nCRITICAL ERROR: Truncation within known good path.\n"); \
			printf("--> maxED = %u, y = %u, qlen=%u\n", maxED, y, qlen); \
			exit(1); \
		} \
		LBN+=y>maxED, ++HBN; \
		_mm_store_si128((void*)(curSc+HBN),_mm_set1_epi8(-1)); /* max the right block */ \
		_mm_store_si128((void*)(prevSc+LBN - 1), _mm_set1_epi8(-1)); /* max left of new */ \
		HBN += HBN < rwidth; \
	} \
	\
	/* Do scores (hybrid) */ \
	__m128i curShfV = _mm_setzero_si128(), curShfRV = _mm_setzero_si128(), \
		minIntV = _mm_set1_epi8(-1); \
	for (uint32_t i = LB; i < HB; ++i) { \
		__m128i score = _mm_load_si128((void*)(curSc + i)), \
			shift = _mm_load_si128((void*)(curSh + i)), \
			shiftR = _mm_load_si128((void*)(curShR + i)); \
		__m128i scoreM = _mm_min_epu8(score, minIntV); \
		__m128i shiftM = _mm_min_epu8(shift, curShfV); \
		\
		__m128i shift_le_curShfV = _mm_cmpeq_epi8(shiftM,shift); \
		__m128i score_eq_minIntV = _mm_cmpeq_epi8(score, minIntV); \
		__m128i scoreM_eq_minIntV = _mm_cmpeq_epi8(scoreM,minIntV); \
		__m128i tiebreak = _mm_andnot_si128(shift_le_curShfV,score_eq_minIntV); \
		__m128i condition = _mm_andnot_si128(tiebreak,scoreM_eq_minIntV); \
		\
		curShfV = _mm_blendv_epi8(shift,curShfV,condition); \
		curShfRV = _mm_blendv_epi8(shiftR,curShfRV,condition); \
		minIntV = scoreM; \
	} \
	\
	__m128 QLm1 = _mm_set1_ps(qlen - 1); \
	__m128 sc = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(minIntV)); \
	__m128 sh = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(curShfV)); \
	sc = _mm_sub_ps(_mm_set1_ps(1), _mm_div_ps(sc,_mm_add_ps(QLm1,sh))); \
	_mm_stream_ps(M16->score,sc); \
	__m128 sc2 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_srli_si128(minIntV,4))); \
	__m128 sh2 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_srli_si128(curShfV,4))); \
	sc2 = _mm_sub_ps(_mm_set1_ps(1), _mm_div_ps(sc2,_mm_add_ps(QLm1,sh2))); \
	_mm_stream_ps(M16->score+4,sc2); \
	__m128 sc3 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_srli_si128(minIntV,8))); \
	__m128 sh3 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_srli_si128(curShfV,8))); \
	sc3 = _mm_sub_ps(_mm_set1_ps(1), _mm_div_ps(sc3,_mm_add_ps(QLm1,sh3))); \
	_mm_stream_ps(M16->score+8,sc3); \
	__m128 sc4 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_srli_si128(minIntV,12))); \
	__m128 sh4 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_srli_si128(curShfV,12))); \
	sc4 = _mm_sub_ps(_mm_set1_ps(1), _mm_div_ps(sc4,_mm_add_ps(QLm1,sh4))); \
	_mm_stream_ps(M16->score+12,sc4); \
	\
	/* calculate alignment index */ \
	__m128i t = _mm_set1_epi32(255); \
	__m128i I1 = _mm_set1_epi32(-1), I2 = _mm_set1_epi32(-1), \
		I3 = _mm_set1_epi32(-1), I4 = _mm_set1_epi32(-1); \
	for (uint32_t i = LB; i < HB; ++i) { \
		__m128i score = _mm_load_si128((void*)(curSc + i)), \
			shift = _mm_load_si128((void*)(curSh + i)); \
		__m128i isGood = _mm_and_si128(_mm_cmpeq_epi8(score,minIntV), \
			_mm_cmpeq_epi8(shift,curShfV)); \
		__m128i set1 = _mm_cmpeq_epi32(_mm_cvtepu8_epi32(isGood),t); \
		__m128i set2 = _mm_cmpeq_epi32(_mm_cvtepu8_epi32(_mm_srli_si128(isGood,4)),t); \
		__m128i set3 = _mm_cmpeq_epi32(_mm_cvtepu8_epi32(_mm_srli_si128(isGood,8)),t); \
		__m128i set4 = _mm_cmpeq_epi32(_mm_cvtepu8_epi32(_mm_srli_si128(isGood,12)),t); \
		I1 = _mm_blendv_epi8(I1,_mm_set1_epi32(i),set1); \
		I2 = _mm_blendv_epi8(I2,_mm_set1_epi32(i),set2); \
		I3 = _mm_blendv_epi8(I3,_mm_set1_epi32(i),set3); \
		I4 = _mm_blendv_epi8(I4,_mm_set1_epi32(i),set4); \
	} \
	_mm_stream_si128((void*)(M16->finalPos),I1); \
	_mm_stream_si128((void*)(M16->finalPos+4),I2); \
	_mm_stream_si128((void*)(M16->finalPos+8),I3); \
	_mm_stream_si128((void*)(M16->finalPos+12),I4); \
	_mm_stream_si128((void*)M16->numGapR,curShfRV); \
	_mm_stream_si128((void*)M16->numGapQ,curShfV); \
}

// Two implementations of primary re-scoring aligner for the alphabet-sensitive
// and alphabet agnostic cases, respectively
inline void reScoreM_mat16(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, 
 DualCoil *Matrix, DualCoil *Shifts, DualCoil *ShiftR, uint32_t maxED, DualCoil *profile, MetaPack *M16) 
	RESCOREM_PROTYPE(DIAGSC_MAT16)

inline void reScoreM_xalpha(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, 
 DualCoil *Matrix, DualCoil *Shifts, DualCoil *ShiftR, uint32_t maxED, DualCoil *profile, MetaPack *M16) 
	RESCOREM_PROTYPE(DIAGSC_XALPHA)

// Aligner specialized for the heuristic mode ('-p'). Differs from primary
// aligners (futher below) in that it doesn't enable seeking into the alignment
// matrix. It is hence a hybrid of 'rescoreM' (above) and 'aded' (below) aligners 
inline uint32_t prune_ed_mat16(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, 
 uint32_t width, DualCoil *Matrix, DualCoil *profile, uint32_t maxED, uint8_t *MinA) {
	uint32_t y, x; 
	__m128i maxEDv = _mm_set1_epi8(maxED+1); // ensure: maxED <= 254! /* BAD if score >= maxED+1 */ 
	--query, --ref, ++qlen, ++rwidth; 
	//maxED = maxED < qlen ? maxED : qlen; 
	uint32_t LB = 1, HB = rwidth + maxED - qlen + 2, LBN, HBN;
	HB = MIN(rwidth, HB);
	
	DualCoil *restrict prev = Matrix + width, *restrict cur = Matrix; 
	_mm_store_si128((void*)(cur),_mm_set1_epi8(1)); 
	char qLet = query[1]; 
	for (x = 1; x < HB; ++x) { 
		__m128i rChunk = _mm_load_si128((void*)(ref+x)); 
		__m128i score = DIAGSC_MAT16; 
		_mm_store_si128((void*)(cur+x),score); 
	}
	_mm_store_si128((void*)(cur+HB),_mm_set1_epi8(-1));
	HB += HB < rwidth;
	for (y = 2; y <= maxED; ++y) { 
		char qLet = query[y]; 
		DualCoil *temp = cur; cur = prev, prev = temp; 
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); /* column 0 */ 
		for (x = 1; x < HB; ++x) { 
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); 
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); 
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); 
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); 
			__m128i diagSc = DIAGSC_MAT16; 
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); 
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreU,score); 
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreL,score); 
			_mm_store_si128((void*)(cur+x),score); 
		} 
		_mm_store_si128((void*)(cur+HB),_mm_set1_epi8(-1)); /* kill the right block */ 
		HB += HB < rwidth;
	} 
	HBN = HB; LBN = 1;
	for (; y < qlen; ++y) { 
		LB = LBN, HB = HBN; 
		LBN = 0; 
		char qLet = query[y]; 
		DualCoil *temp = cur; cur = prev, prev = temp; 
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); 
		for (x = LB; x < HB; ++x) { 
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); 
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); 
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); 
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); 
			
			__m128i diagSc = DIAGSC_MAT16; 
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); 
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreU,score); 
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreL,score); 
			
			/* Do bounds handling */ 
			__m128i anyBad = _mm_cmpeq_epi8(maxEDv,_mm_min_epu8(maxEDv,score)); 
			score = _mm_or_si128(anyBad,score); 
			if (_mm_movemask_epi8(anyBad) != 0xFFFF) { 
				if (!LBN) LBN = x; 
				HBN = x; 
			} 
			_mm_store_si128((void*)(cur+x),score); 
		} 
		
		if (!LBN) return -1; 
		++LBN; /* we guarantee the bounds will close in */ 
		++HBN; /* since it'll bound the 'for (...x < rightBound; ++x)' */ 
		_mm_store_si128((void*)(cur+HBN),_mm_set1_epi8(-1)); /* kill the right block */ 
		/*if (LBN > 1)*/ 
		_mm_store_si128((void*)(prev+LBN - 1), _mm_set1_epi8(-1)); /* kill left of new */ 
		HBN += HBN < rwidth; 
	} 
	/* Vectorized min reduction */ 
	__m128i minIntV = _mm_set1_epi8(-1); 
	for (uint32_t i = LB; i < HB; ++i) 
		minIntV = _mm_min_epu8(minIntV,_mm_load_si128((void*)(cur+i))); 
	//__m128i mc = _mm_min_epu8(minIntV,_mm_set1_epi8(err_ceil));
	//*minA = _mm_movemask_epi8(_mm_cmpeq_epi8(minIntV,mc));
	
	_mm_store_si128((void*)MinA,minIntV); 
	__m128i c = _mm_srli_si128(minIntV,8); 
	minIntV = _mm_min_epu8(minIntV,c); 
	c = _mm_srli_si128(minIntV,4); 
	minIntV = _mm_min_epu8(minIntV,c); 
	c = _mm_srli_si128(minIntV,2); 
	minIntV = _mm_min_epu8(minIntV,c); 
	c = _mm_srli_si128(minIntV,1); 
	minIntV = _mm_min_epu8(minIntV,c); 
	return _mm_extract_epi8(minIntV,0); /* save min of mins as new bound */ 
}

// The primary high-speed alignment core in BURST. This is the stripped, speed-
// oriented aligner using the synergistic pruning and caching construct in 
// the manuscript. It is designed to allow arbitrary seeking into the alignment
// matrix to bypass search space, in addition to narrowing the search space as
// alignment progresses. It does not record any alignment statistics, necessitating
// the rescoreM aligner further above when best candidate(s) are identified here.
#define ADED_PROTOTYPE(DIAG_FUNC) {\
if (startQ > *LoBound) return -1; /* truncation signal */ \
	uint32_t y, x; \
	__m128i maxEDv = _mm_set1_epi8(maxED+1); /* < 255 ? maxED+1 : 255); /* BAD if score >= maxED+1 */ \
	--query; --ref; ++qlen; \
	maxED = maxED < qlen ? maxED : qlen; \
	DualCoil *restrict cur = Matrix + startQ * width, *restrict prev = cur - width; \
	for (y = startQ; y <= maxED; ++y) { \
		LoBound[y+1] = 1, HiBound[y+1] = rwidth; \
		char qLet = query[y]; \
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); /* column 0 */ \
		for (x = 1; x < rwidth; ++x) { \
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); /* refRow += 16; */ \
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); \
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); \
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); \
			\
			__m128i diagSc = DIAG_FUNC; \
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); \
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); \
			score = _mm_min_epu8(scoreU,score); \
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); \
			score = _mm_min_epu8(scoreL,score); \
			\
			_mm_store_si128((void*)(cur+x),score); \
		} \
		DualCoil *temp = cur; \
		cur = (y <= cacheSz) ? cur + width : prev; \
		prev = (y <= cacheSz) ? prev + width : temp; \
	} \
	for (; y < qlen; ++y) { \
		LoBound[y+1] = 0; \
		char qLet = query[y]; \
		/*printf("---Entering LIM loop: y = %u [let %u] / maxED = %u, LoBound[y]=%u, LoBound[y+1]=%u,HiBound[y]=%u,HiBound[y+1]=%u...", \
			y,qLet,maxED, LoBound[y], LoBound[y+1], HiBound[y], HiBound[y+1]);*/ \
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); \
		for (x = LoBound[y]; x < HiBound[y]; ++x) { \
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); \
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); \
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); \
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); \
			\
			__m128i diagSc = DIAG_FUNC; \
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); \
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); \
			score = _mm_min_epu8(scoreU,score); \
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); \
			score = _mm_min_epu8(scoreL,score); \
			\
			/* Do bounds handling */ \
			__m128i anyBad = _mm_cmpeq_epi8(maxEDv,_mm_min_epu8(maxEDv,score)); \
			score = _mm_or_si128(anyBad,score); \
			if (_mm_movemask_epi8(anyBad) != 0xFFFF) { \
				if (!LoBound[y+1]) LoBound[y+1] = x; \
				HiBound[y+1] = x; \
			} \
			_mm_store_si128((void*)(cur+x),score); \
		} \
		\
		if (!LoBound[y+1]) { \
			*LoBound = y; /* truncation location */ \
			return -1; \
		} \
		++LoBound[y+1]; /* we guarantee the bounds will close in */ \
		++HiBound[y+1]; /* since it'll bound the 'for (...x < rightBound; ++x)' */ \
		_mm_store_si128((void*)(cur+HiBound[y+1]),_mm_set1_epi8(-1)); /* kill the right block */ \
		\
		DualCoil *temp = cur; \
		cur = y <= cacheSz ? cur + width : prev; \
		prev = y <= cacheSz ? prev + width : temp; \
		/*if (LoBound[y+1] > 1) */\
		_mm_store_si128((void*)(cur+LoBound[y+1] - 1), _mm_set1_epi8(-1)); /* kill left of new */ \
		HiBound[y+1] += HiBound[y+1] < rwidth; \
	} \
	\
	/* Vectorized min reduction */ \
	DualCoil *last = prev; \
	__m128i minIntV = _mm_set1_epi8(-1); \
	for (uint32_t i = LoBound[qlen-1]; i < HiBound[qlen-1]; ++i) \
		minIntV = _mm_min_epu8(minIntV,_mm_load_si128((void*)(last+i))); \
	_mm_stream_si128((void*)MinA,minIntV); \
	__m128i c = _mm_srli_si128(minIntV,8); \
	__m128i b = _mm_min_epu8(minIntV,c); \
	c = _mm_srli_si128(b,4); \
	b = _mm_min_epu8(b,c); \
	c = _mm_srli_si128(b,2); \
	b = _mm_min_epu8(b,c); \
	c = _mm_srli_si128(b,1); \
	b = _mm_min_epu8(b,c); \
	\
	*LoBound = -1; \
	return _mm_extract_epi8(b,0); /* save min of mins as new bound */ \
}
/// Implementation of the prototype "aded' aligner above. 
inline uint32_t aded_mat16(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, DualCoil *Matrix,
 DualCoil *profile, uint32_t maxED, uint32_t startQ, uint32_t *LoBound, uint32_t *HiBound, DualCoil *MinA) ADED_PROTOTYPE(DIAGSC_MAT16)
inline uint32_t aded_xalpha(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, DualCoil *Matrix,
 DualCoil *profile, uint32_t maxED, uint32_t startQ, uint32_t *LoBound, uint32_t *HiBound, DualCoil *MinA) 
 ADED_PROTOTYPE(DIAGSC_XALPHA)

// Variant of standard 'aded' above, which further optimizes search by restricting
// the reference length analyzed based on maximum allowed errors introduced at the
// ends of queries by sliding off the reference early in the alignment
inline uint32_t aded_mat16L(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, uint32_t minlen, 
 DualCoil *Matrix, DualCoil *profile, uint32_t maxED, uint32_t startQ, uint32_t *LoBound, uint32_t *HiBound, DualCoil *MinA) {
	if (startQ > *LoBound || minlen > rwidth + maxED) return -1; /* truncation signal */ 
	uint32_t y, x; 
	__m128i maxEDv = _mm_set1_epi8(maxED+1); /* < 255 ? maxED+1 : 255); /* BAD if score >= maxED+1 */ 
	--query; --ref; ++qlen; 
	maxED = maxED < qlen ? maxED : qlen; 
	if (startQ == 1) LoBound[1] = 1, HiBound[1] = rwidth + maxED - minlen + 1;
	HiBound[1] = MIN(HiBound[1],rwidth);
	
	DualCoil *restrict cur = Matrix + startQ * width, *restrict prev = cur - width; 
	for (y = startQ; y <= maxED; ++y) { 
		char qLet = query[y]; 
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); /* column 0 */ 
		for (x = 1; x < HiBound[y]; ++x) { 
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); /* refRow += 16; */ 
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); 
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); 
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); 
			
			__m128i diagSc = DIAGSC_MAT16; 
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); 
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreU,score); 
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreL,score); 
			
			_mm_store_si128((void*)(cur+x),score); 
		} 
		LoBound[y+1] = 1; 
		HiBound[y+1] = HiBound[y] + (HiBound[y] < rwidth); 
		_mm_store_si128((void*)(cur+HiBound[y]),_mm_set1_epi8(-1)); /* kill the right block */ 
		DualCoil *temp = cur; 
		cur = (y <= cacheSz) ? cur + width : prev; 
		prev = (y <= cacheSz) ? prev + width : temp; 
	} 
	for (; y < qlen; ++y) { 
		LoBound[y+1] = 0; 
		char qLet = query[y]; 
		/*printf("---Entering LIM loop: y = %u [let %u] / maxED = %u, LoBound[y]=%u, LoBound[y+1]=%u,HiBound[y]=%u,HiBound[y+1]=%u...", 
			y,qLet,maxED, LoBound[y], LoBound[y+1], HiBound[y], HiBound[y+1]);*/ 
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); 
		for (x = LoBound[y]; x < HiBound[y]; ++x) { 
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); 
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); 
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); 
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); 
			
			__m128i diagSc = DIAGSC_MAT16; 
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); 
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreU,score); 
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); 
			score = _mm_min_epu8(scoreL,score); 
			
			/* Do bounds handling */ 
			__m128i anyBad = _mm_cmpeq_epi8(maxEDv,_mm_min_epu8(maxEDv,score)); 
			score = _mm_or_si128(anyBad,score); 
			if (_mm_movemask_epi8(anyBad) != 0xFFFF) { 
				if (!LoBound[y+1]) LoBound[y+1] = x; 
				HiBound[y+1] = x; 
			} 
			_mm_store_si128((void*)(cur+x),score); 
		} 
		
		if (!LoBound[y+1]) { 
			*LoBound = y; /* truncation location */ 
			return -1; 
		} 
		++LoBound[y+1]; /* we guarantee the bounds will close in */ 
		++HiBound[y+1]; /* since it'll bound the 'for (...x < rightBound; ++x)' */ 
		_mm_store_si128((void*)(cur+HiBound[y+1]),_mm_set1_epi8(-1)); /* kill the right block */ 
		
		DualCoil *temp = cur; 
		cur = y <= cacheSz ? cur + width : prev; 
		prev = y <= cacheSz ? prev + width : temp; 
		/*if (LoBound[y+1] > 1) */
		_mm_store_si128((void*)(cur+LoBound[y+1] - 1), _mm_set1_epi8(-1)); /* kill left of new */ 
		HiBound[y+1] += HiBound[y+1] < rwidth; 
	} 
	
	/* Vectorized min reduction */ 
	DualCoil *last = prev; 
	__m128i minIntV = _mm_set1_epi8(-1); 
	for (uint32_t i = LoBound[qlen-1]; i < HiBound[qlen-1]; ++i) 
		minIntV = _mm_min_epu8(minIntV,_mm_load_si128((void*)(last+i))); 
	_mm_stream_si128((void*)MinA,minIntV); 
	__m128i c = _mm_srli_si128(minIntV,8); 
	__m128i b = _mm_min_epu8(minIntV,c); 
	c = _mm_srli_si128(b,4); 
	b = _mm_min_epu8(b,c); 
	c = _mm_srli_si128(b,2); 
	b = _mm_min_epu8(b,c); 
	c = _mm_srli_si128(b,1); 
	b = _mm_min_epu8(b,c); 
	
	*LoBound = -1; 
	return _mm_extract_epi8(b,0); /* save min of mins as new bound */ 
}

// Non-vectorized function that looks up small integers to replace DNA bases
inline void translateNV(char* string, size_t len) { 
	for (size_t i = 0; i < len; ++i) string[i] = CHAR2NUM[string[i]]; }

// Vectorized version of the above. Requires dual-shuffle and mask to get
// around the 16-value lookup limit (extend to masked 2x16 and blend)
//#ifdef __SSSE3__
static inline __m128i TWOSHUFAL(char *query) {
	__m128i d2 = _mm_set_epi8( 0, 0, 0, 0, 0, 5, 9, 5,11,13, 4, 4,10, 8, 5, 5); // P,Q... rev
	__m128i d1 = _mm_set_epi8( 5, 5, 7, 5, 6, 5, 5,14, 3, 5, 5,15, 2,12, 1, 0); // null,A...rev
	__m128i q = _mm_load_si128((void*)query); 
	__m128i f = _mm_set1_epi8(15);
	__m128i qf = _mm_and_si128(q,f);
	__m128i s = _mm_set1_epi8(16);
	__m128i L1f = _mm_shuffle_epi8(d1,qf);
	__m128i L2f = _mm_shuffle_epi8(d2,qf);
	__m128i qs = _mm_and_si128(q,s);
	__m128i eqs = _mm_cmpeq_epi8(qs,s);
	return _mm_blendv_epi8(L1f,L2f,eqs);
}
static inline void translate16aln(char* string, size_t len) { // makes string into nums
	size_t i=0; for (; i < len; i+= 16) 
		_mm_store_si128((void*)(string+i),TWOSHUFAL(string+i));
}
//#else
	#define translate16aln translateNV
//#endif

// Called before alignment so the scoring structures can be populated with
// The correct values. Necessary because can't define vector constants and
// wanted to change ambiguity handling based on '-n' or '-y' flags. 
void setScore() {
/*{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //. [0]
	//. A C G T N K M R Y S W B V H D
	 -1,0,1,1,1,Z,1,0,0,1,1,0,1,0,0,0, //A [1]
	 -1,1,0,1,1,Z,1,0,1,0,0,1,0,0,0,1, //C [2]
	 -1,1,1,0,1,Z,0,1,0,1,0,1,0,0,1,0, //G [3]
	 -1,1,1,1,0,Z,0,1,1,0,1,0,0,1,0,0, //T/U [4]
	 -1,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z, //N/X [5]
	 -1,1,1,0,0,Z,0,1,1,1,1,1,0,1,1,0, //K [6]
	 -1,0,0,1,1,Z,1,0,1,1,1,1,1,0,0,1, //M [7]
	 -1,0,1,0,1,Z,1,1,0,1,1,1,1,0,1,0, //R [8]
	 -1,1,0,1,0,Z,1,1,1,0,1,1,0,1,0,1, //Y [9]
	 -1,1,0,0,1,Z,1,1,1,1,0,1,0,0,1,1, //S [10]
	 -1,0,1,1,0,Z,1,1,1,1,1,0,1,1,0,0, //W [11]
	 -1,1,0,0,0,Z,0,1,1,0,0,1,0,1,1,1, //B [12]
	 -1,0,0,0,1,Z,1,0,0,1,0,1,1,0,1,1, //V [13]
	 -1,0,0,1,0,Z,1,0,1,0,1,0,1,1,0,1, //H [14]
	 -1,0,1,0,0,Z,0,1,0,1,1,0,1,1,1,0, //D [15]
	}; */
	if (Z) SCORENVedN[1*16 + 5] = Z, 
		   SCORENVedN[2*16 + 5] = Z, 
		   SCORENVedN[3*16 + 5] = Z, 
		   SCORENVedN[4*16 + 5] = Z,
		   SCORENVedN[5*16 + 5] = Z,
		   SCORENVedN[6*16 + 5] = Z,
		   SCORENVedN[7*16 + 5] = Z,
		   SCORENVedN[8*16 + 5] = Z,
		   SCORENVedN[9*16 + 5] = Z,
		   SCORENVedN[10*16 + 5] = Z,
		   SCORENVedN[11*16 + 5] = Z,
		   SCORENVedN[12*16 + 5] = Z,
		   SCORENVedN[13*16 + 5] = Z,
		   SCORENVedN[14*16 + 5] = Z,
		   SCORENVedN[15*16 + 5] = Z;
		   // Also penalize N in queries (new!)
		   SCORENVedN[5*16 + 1] = Z;
		   SCORENVedN[5*16 + 2] = Z;
		   SCORENVedN[5*16 + 3] = Z;
		   SCORENVedN[5*16 + 4] = Z;
		   SCORENVedN[5*16 + 6] = Z;
		   SCORENVedN[5*16 + 7] = Z;
		   SCORENVedN[5*16 + 8] = Z;
		   SCORENVedN[5*16 + 9] = Z;
		   SCORENVedN[5*16 + 10] = Z;
		   SCORENVedN[5*16 + 11] = Z;
		   SCORENVedN[5*16 + 12] = Z;
		   SCORENVedN[5*16 + 13] = Z;
		   SCORENVedN[5*16 + 14] = Z;
		   SCORENVedN[5*16 + 15] = Z;
		   
	#define BAD_IX 0
	for (int i = 0; i < 65; ++i) CHAR2NUM[i] = BAD_IX;
	for (int i = 65; i < 91; ++i) CHAR2NUM[i] = 5; //N
	for (int i = 91; i < 97; ++i) CHAR2NUM[i] = BAD_IX;
	for (int i = 97; i < 122; ++i) CHAR2NUM[i] = 5; //N
	for (int i = 122; i < 128; ++i) CHAR2NUM[i] = BAD_IX;
	CHAR2NUM['a'] = CHAR2NUM['A'] = 1;
	CHAR2NUM['c'] = CHAR2NUM['C'] = 2;
	CHAR2NUM['g'] = CHAR2NUM['G'] = 3;
	CHAR2NUM['t'] = CHAR2NUM['T'] = 4;
	CHAR2NUM['u'] = CHAR2NUM['U'] = 4; // 5=N, set above
	CHAR2NUM['k'] = CHAR2NUM['K'] = 6;
	CHAR2NUM['m'] = CHAR2NUM['M'] = 7;
	CHAR2NUM['r'] = CHAR2NUM['R'] = 8;
	CHAR2NUM['y'] = CHAR2NUM['Y'] = 9;
	CHAR2NUM['s'] = CHAR2NUM['S'] = 10;
	CHAR2NUM['w'] = CHAR2NUM['W'] = 11;
	CHAR2NUM['b'] = CHAR2NUM['B'] = 12;
	CHAR2NUM['v'] = CHAR2NUM['V'] = 13;
	CHAR2NUM['h'] = CHAR2NUM['H'] = 14;
	CHAR2NUM['d'] = CHAR2NUM['D'] = 15;
	
	GAPV = _mm_set1_epi8(GAP);
	SCOREFAST[0] = _mm_set1_epi8(-1); //.
	//							                        X U
	//							 D H V B W S Y R M K N T G C A .
	SCOREFAST[1]  = _mm_set_epi8(0,0,0,1,0,1,1,0,0,1,Z,1,1,1,0,-1); //A
	SCOREFAST[2]  = _mm_set_epi8(1,0,0,0,1,0,0,1,0,1,Z,1,1,0,1,-1); //C
	SCOREFAST[3]  = _mm_set_epi8(0,1,0,0,1,0,1,0,1,0,Z,1,0,1,1,-1); //G
	SCOREFAST[4]  = _mm_set_epi8(0,0,1,0,0,1,0,1,1,0,Z,0,1,1,1,-1); //T or U
	//SCOREFAST[5]  = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,Z,0,0,0,0,-1); //N or X
	SCOREFAST[5]  = _mm_set_epi8(Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,-1); //N or X
	SCOREFAST[6]  = _mm_set_epi8(0,1,1,0,1,1,1,1,1,0,Z,0,0,1,1,-1); //K
	SCOREFAST[7]  = _mm_set_epi8(1,0,0,1,1,1,1,1,0,1,Z,1,1,0,0,-1); //M
	SCOREFAST[8]  = _mm_set_epi8(0,1,0,1,1,1,1,0,1,1,Z,1,0,1,0,-1); //R
	SCOREFAST[9]  = _mm_set_epi8(1,0,1,0,1,1,0,1,1,1,Z,0,1,0,1,-1); //Y
	SCOREFAST[10] = _mm_set_epi8(1,1,0,0,1,0,1,1,1,1,Z,1,0,0,1,-1); //S
	SCOREFAST[11] = _mm_set_epi8(0,0,1,1,0,1,1,1,1,1,Z,0,1,1,0,-1); //W
	SCOREFAST[12] = _mm_set_epi8(1,1,1,0,1,0,0,1,1,0,Z,0,0,0,1,-1); //B
	SCOREFAST[13] = _mm_set_epi8(1,1,0,1,1,0,1,0,0,1,Z,1,0,0,0,-1); //V
	SCOREFAST[14] = _mm_set_epi8(1,0,1,1,0,1,0,1,0,1,Z,0,1,0,0,-1); //H
	SCOREFAST[15] = _mm_set_epi8(0,1,1,1,0,1,1,0,1,0,Z,0,0,1,0,-1); //D
}

// Compares string pointer arrays by containing strings
static int cpcmp(const void *a, const void *b) 
	{ return strcmp(**(char ***)a, **(char ***)b); }
// Compares Tuxedo objects (sequence + length) by length
static int cmpPackLen(const void *first, const void *second) {
	Tuxedo *a = (Tuxedo *)first, *b = (Tuxedo *)second;
	return a->len < b->len ? -1 : a->len > b->len;

// Compares Tuxedo objects (sequence + length) by sequence
static int cmpPackSeq(const void *first, const void *second) {
	Tuxedo *a = (Tuxedo *)first, *b = (Tuxedo *)second;
	return strcmp(a->seq, b->seq);
}

// Wide binary search. Returns lower index of range containing key
static inline uint64_t uWBS(uint64_t *ixList, uint64_t key, uint64_t maxIX) {
	uint64_t middle, low = 0, high = maxIX;
	while (low <= high) {
		middle = low + ((high - low) >> 1);
		if (key > ixList[middle]) low = middle + 1;
		else if (key < ixList[middle]) high = middle - 1;
		else break; 
	}
	return ixList[middle] > key ? middle - 1: middle;
}

// Fingerprint: patterned Jaccard hash.
// Add 4-mers by 2-bit hash that occur after a pattern (here, 'A'). 
// Effectively samples every 1/4^5 bases into a 256-bit bucket by presence,
// which can be efficiently compared to other 256-bit buckets by Jaccard.
// Allows rough similarity clustering as well as exclusion by dissimilarity.
// Effectively an optimal lower-bound error profiler when the references 
// receive patterned hashing and queries receive exhaustive hashing
#define NL 4
                    //0 A C G T N K M R Y S W B V H D
uint8_t A_COMPAT[] = {0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,1}; 
uint8_t A_COMPNN[] = {0,1,0,0,0,0,0,1,1,0,0,1,0,1,1,1}; 
uint8_t AMBIGMAP[] = {4,0,1,2,3,5,4,4,4,4,4,4,4,4,4,4};
uint8_t AMBIGS[][16] = {{-1},{0,-1},{1,-1},{2,-1},{3,-1},     //0ACGT
					 {0,1,2,3,-1},{2,3,-1},{0,1,-1},{0,2,-1}, //NKMR
					 {1,3,-1},{1,2,-1},{0,3,-1},{1,2,3,-1},   //YSWB
					 {0,1,2,-1},{0,1,3,-1},{0,2,3,-1}};       //VHD
void setAmbigPrince(Prince *P, char *S, uint8_t w, uint8_t ix) {
	if (ix == NL) P->a[w >> 3] |= 1 << (w & 7);
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		setAmbigPrince(P, S, w << 2 | AMBIGS[S[ix]][i], ix + 1);
}
// Creates fingerprints according to ambiguity handling and 
static inline PackaPrince create_fingerprints(char **Seqs, uint32_t N, uint32_t *Lens, uint32_t *Mask, uint8_t isRef, uint8_t dualAmbig) {//, char **Names) {
	uint8_t doDR = isRef && dualAmbig;
	if (doDR && (uint64_t)N * 2ull >= UINT32_MAX) {puts("CRITICAL ERROR FP21_b"); exit(4);}
	uint8_t *Pops = 0, bad = isRef ? UINT8_MAX : 0;
	void *F_init; //Prince *F = calloc_a(64,(N + doDR*N)*sizeof(*F),&F_init);
	Prince *F = calloc(N + doDR*N, sizeof(*F)); F_init = F;
	if (!F_init) {fputs("OOM:fingerprint",stderr); exit(3);}
	uint8_t UA[32] = {0};
	UA[1] = UA[2] = UA[3] = UA[4] = 1; // only A-T ix's are unambig
	Prince AllMatch; 
	for (int i = 0; i < sizeof(AllMatch); ++i) AllMatch.a[i] = bad;
	
	uint32_t aix = N, *A = malloc(doDR? N*sizeof(*A):0); // extra ambig space (refs+dualAmbig)
	if (doDR) 
		#pragma omp parallel for 
		for (uint32_t i = 0; i < N; ++i) {
		uint8_t *S = (uint8_t *)Seqs[Mask ? Mask[i] : i];
		uint32_t L = Lens[Mask ? Mask[i] : i], N_found = 0;

		for (uint32_t j = 0; j + NL < L; ++j) {
			if (A_COMPAT[S[j]]) setAmbigPrince(F+i,S+j+1,0,0);
			if (S[j]==5) N_found = 1;
		}
		
		if (N_found) {
			#pragma omp atomic capture
			A[i] = aix++;
			for (uint32_t j = 0; j + NL < L; ++j) {
				if (A_COMPNN[S[j]]) {
					if (S[j+1] == 5) {j += 1; continue;}
					else if (S[j+2] == 5) {j += 2; continue;}
					else if (S[j+3] == 5) {j += 3; continue;}
					else if (S[j+4] == 5) {j += 4; continue;}
					setAmbigPrince(F+A[i],S+j+1,0,0);
				}
			}
		} else A[i] = i;
	}
	else if (isRef) 
		#pragma omp parallel for
		for (uint32_t i = 0; i < N; ++i) {
			uint8_t *S = (uint8_t *)Seqs[Mask ? Mask[i] : i];
			uint32_t L = Lens[Mask ? Mask[i] : i];
			for (uint32_t j = 0; j + NL < L; ++j) 
				if (A_COMPAT[S[j]]) setAmbigPrince(F+i,S+j+1,0,0);
		}
	else if (dualAmbig) 
		#pragma omp parallel for
		for (uint32_t i = 0; i < N; ++i) {
			uint8_t *S = (uint8_t *)Seqs[Mask ? Mask[i] : i];
			uint32_t L = Lens[Mask ? Mask[i] : i];
			for (uint32_t j = 0; j + NL < L; ++j) if (S[j]==1) {
				if (S[j+1] > 4) {j += 1; continue;} // do UA instead of > 4?
				else if (S[j+2] > 4) {j += 2; continue;}
				else if (S[j+3] > 4) {j += 3; continue;}
				else if (S[j+4] > 4) {j += 4; continue;}
				uint8_t w = (S[j+1]-1) << 6 | (S[j+2]-1) << 4 | 
					(S[j+3]-1) << 2 | (S[j+4]-1);
				F[i].a[w >> 3] |= 1 << (w & 7);
				j += 4; 
			}
		}
	else 
		#pragma omp parallel for
		for (uint32_t i = 0; i < N; ++i) {
			uint8_t *S = (uint8_t *)Seqs[Mask ? Mask[i] : i];
			uint32_t L = Lens[Mask ? Mask[i] : i];
			if (L <= NL) continue;
			for (uint32_t j = 0; j < L; ++j) 
				if (!UA[S[j]]) {F[i] = AllMatch; goto NEXT_FP_ITER;}
			
			uint8_t w = S[0] - 1;
			for (int j = 1; j < NL; ++j) w <<= 2, w |= S[j] - 1; // prologue
			for (uint32_t j = NL, lix = -1; j < L; ++j) {
				w <<= 2, w |= S[j] - 1;
				if (S[j-NL] == 1 && j > NL + lix) //printf("%u,",w),
					F[i].a[w >> 3] |= 1 << (w & 7), lix = j;
			}
			NEXT_FP_ITER:NULL;
		}
	
	if (!isRef) {
		Pops = malloc(N*sizeof(*Pops));
		if (!Pops) {fputs("OOM:fingerprint",stderr); exit(3);}
		for (uint32_t i = 0; i < N; ++i) {	 
			// count totals
			int tot = 0;
			for (int j = 0; j < sizeof(Prince)/sizeof(uint64_t); ++j)
				tot += _mm_popcnt_u64(F[i].w[j]);
			Pops[i] = tot;
		}
	}
	//if (doDR) realloc(F_init,)
	if (doDR) F = realloc(F, aix*sizeof(*F)), F_init = F;
	return (PackaPrince){F, Pops, F_init, aix, A};
}
static inline PackaPrince create_fingerprintsQ(ShrBin *ShrBins, UniBin *UniBins, uint32_t N) {
	uint8_t *Pops = malloc(N*sizeof(*Pops));
	void *F_init; 
	Prince *F = calloc(N, sizeof(*F)); F_init = F;
	if (!F_init || !Pops) {fputs("OOM:fingerprintQ",stderr); exit(3);}
	#pragma omp parallel for
	for (uint32_t i = 0; i < N; ++i) {
		uint8_t *S = (uint8_t *)UniBins[i].Seq;
		uint32_t L = ShrBins[UniBins[i].six].len;
		for (uint32_t j = 0; j + NL < L; ++j) if (S[j]==1) {
			if (S[j+1] > 4) {j += 1; continue;} 
			else if (S[j+2] > 4) {j += 2; continue;}
			else if (S[j+3] > 4) {j += 3; continue;}
			else if (S[j+4] > 4) {j += 4; continue;}
			uint8_t w = (S[j+1]-1) << 6 | (S[j+2]-1) << 4 | 
				(S[j+3]-1) << 2 | (S[j+4]-1);
			F[i].a[w >> 3] |= 1 << (w & 7);
			j += 4; 
		}
	}
	#pragma omp parallel for
	for (uint32_t i = 0; i < N; ++i) {	 
		int tot = 0;
		for (int j = 0; j < sizeof(Prince)/sizeof(uint64_t); ++j)
			tot += _mm_popcnt_u64(F[i].w[j]);
		Pops[i] = tot > 255 ? 255 : tot;
	}
	return (PackaPrince){F, Pops, F_init, N, 0};
}
#undef NL


static inline uint64_t VPOP(__m256i v) {
	/* #ifdef __AVX2__
		__m256i BIT4 = _mm256_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4, 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
		__m256i MSK4 = _mm256_set1_epi8(0x0F);
		__m256i h = _mm256_shuffle_epi8(BIT4,_mm256_and_si256(MSK4,v));
		__m256i l = _mm256_shuffle_epi8(BIT4,_mm256_and_si256(MSK4,_mm256_srli_epi16(v,4)));
		v = _mm256_add_epi8(h,l);
		v = _mm256_sad_epu8(v,_mm256_setzero_si256());
		return _mm256_extract_epi64(v,0) + 
			_mm256_extract_epi64(v,1) + 
			_mm256_extract_epi64(v,2) + 
			_mm256_extract_epi64(v,3);
	#elif __AVX__ */
		return _mm_popcnt_u64(_mm256_extract_epi64(v, 0)) +
			_mm_popcnt_u64(_mm256_extract_epi64(v, 1)) + 
			_mm_popcnt_u64(_mm256_extract_epi64(v, 2)) +
			_mm_popcnt_u64(_mm256_extract_epi64(v, 3));
	//#endif
	// try a version that adds and shifts, followed by a final extraction?

}
static inline uint64_t FP_pop(Prince *a) { // TODO: try non-pointer instead?
	// #ifdef __AVX__
		// return VPOP(_mm256_load_si256((void*)a));
	// #else
		return _mm_popcnt_u64(a->w[0]) +
			_mm_popcnt_u64(a->w[1]) +
			_mm_popcnt_u64(a->w[2]) +
			_mm_popcnt_u64(a->w[3]);
	// #endif
}

static inline uint64_t FP_intersect(Prince a, Prince b) {
	// #ifdef __AVX2__
		// return VPOP(_mm256_or_si256(a.v,b.v));
	// #elif __AVX__
		// return VPOP(_mm256_castpd_si256(_mm256_or_pd(a.v,b.v)));
	// #else
		return _mm_popcnt_u64(a.w[0] & b.w[0]) +
			_mm_popcnt_u64(a.w[1] & b.w[1]) +
			_mm_popcnt_u64(a.w[2] & b.w[2]) +
			_mm_popcnt_u64(a.w[3] & b.w[3]);
	// #endif
}

static inline uint64_t FP_union(Prince a, Prince b) {
	// #ifdef __AVX2__
		// return VPOP(_mm256_or_si256(a.v,b.v));
	// #elif __AVX__
		// return VPOP(_mm256_castpd_si256(_mm256_or_pd(a.v,b.v)));
	// #else
		return _mm_popcnt_u64(a.w[0] | b.w[0]) +
			_mm_popcnt_u64(a.w[1] | b.w[1]) +
			_mm_popcnt_u64(a.w[2] | b.w[2]) +
			_mm_popcnt_u64(a.w[3] | b.w[3]);
	// #endif
}
static inline uint64_t FP_dist(Prince a, Prince b) {
	// #ifdef __AVX2__
		// return VPOP(_mm256_or_si256(a.v,b.v));
	// #elif __AVX__
		// return VPOP(_mm256_castpd_si256(_mm256_or_pd(a.v,b.v)));
	// #else
		return _mm_popcnt_u64(a.w[0] ^ b.w[0]) +
			_mm_popcnt_u64(a.w[1] ^ b.w[1]) +
			_mm_popcnt_u64(a.w[2] ^ b.w[2]) +
			_mm_popcnt_u64(a.w[3] ^ b.w[3]);
	// #endif
}

static inline Prince FP_combine(Prince a, Prince b) {
	Prince x;
	// #ifdef __AVX2__
		// x.v = _mm256_or_si256(a.v,b.v);
	// #elif __AVX__
		// x.v = _mm256_or_pd(a.v,b.v);
	// #else
		x.w[0] = a.w[0] | b.w[0];
		x.w[1] = a.w[1] | b.w[1];
		x.w[2] = a.w[2] | b.w[2];
		x.w[3] = a.w[3] | b.w[3];
	// #endif
	return x;
}

static inline Split FP_union16x2_NV(Prince *P1, Prince *P2) {
	Prince *P = P1, p[8], s;
	p[0] = FP_combine(P[0],P[1]);
	p[1] = FP_combine(P[2],P[3]);
	p[2] = FP_combine(P[4],P[5]);
	p[3] = FP_combine(P[6],P[7]);
	p[4] = FP_combine(P[8],P[9]);
	p[5] = FP_combine(P[10],P[11]);
	p[6] = FP_combine(P[12],P[13]);
	p[7] = FP_combine(P[14],P[15]);
	
	p[0] = FP_combine(p[0],p[1]);
	p[1] = FP_combine(p[2],p[3]);
	p[2] = FP_combine(p[4],p[5]);
	p[3] = FP_combine(p[6],p[7]);
	
	p[0] = FP_combine(p[0],p[1]);
	p[1] = FP_combine(p[2],p[3]);
	s = FP_combine(p[0],p[1]);
	
	P = P2;
	p[0] = FP_combine(P[0],P[1]);
	p[1] = FP_combine(P[2],P[3]);
	p[2] = FP_combine(P[4],P[5]);
	p[3] = FP_combine(P[6],P[7]);
	p[4] = FP_combine(P[8],P[9]);
	p[5] = FP_combine(P[10],P[11]);
	p[6] = FP_combine(P[12],P[13]);
	p[7] = FP_combine(P[14],P[15]);
	
	p[0] = FP_combine(p[0],p[1]);
	p[1] = FP_combine(p[2],p[3]);
	p[2] = FP_combine(p[4],p[5]);
	p[3] = FP_combine(p[6],p[7]);
	
	p[0] = FP_combine(p[0],p[1]);
	p[1] = FP_combine(p[2],p[3]);
	
	p[0] = FP_combine(p[0],p[1]);
	
	return (Split){FP_pop(&s),FP_pop(p)};
}
static inline Split FP_union16x2_ZA(Prince *P1, Prince *P2) {
	Prince p1 = *P1, p2 = *P2;
	for (int i = 1; i < 16; ++i)
		p1 = FP_combine(p1,P1[i]),
		p2 = FP_combine(p2,P2[i]);
	return (Split){FP_pop(&p1),FP_pop(&p2)};
}
#ifdef __AVX__
static inline Split FP_union16x2_AVX(Prince *P1, Prince *P2) {
	#ifdef __AVX2__
		#define REG __m256i
		#define LOAD(x) _mm256_load_si256((void*)(x))
		#define STORE(m,r) _mm256_store_si256((void*)(m),r)
		#define OR(x,y) _mm256_or_si256(x,y)
	#else
		#define REG __m256
		#define LOAD(x) _mm256_load_ps((void*)(x))
		#define STORE(m,r) _mm256_store_ps((void*)(m),r)
		//#define STORE(m,r) _mm256_store_si256((void*)(m),_mm256_castps_si256(r))
		#define OR(x,y) _mm256_or_ps(x,y)
	#endif
	Prince *P = P1;
	REG p1 = OR(LOAD(P+0),LOAD(P+1));
	REG p2 = OR(LOAD(P+2),LOAD(P+3));
	REG p3 = OR(LOAD(P+4),LOAD(P+5));
	REG p4 = OR(LOAD(P+6),LOAD(P+7));
	REG p5 = OR(LOAD(P+8),LOAD(P+9));
	REG p6 = OR(LOAD(P+10),LOAD(P+11));
	REG p7 = OR(LOAD(P+12),LOAD(P+13));
	REG p8 = OR(LOAD(P+14),LOAD(P+15));
	p1 = OR(p1,p2), p2 = OR(p3,p4), p3 = OR(p5,p6), p4 = OR(p7,p8);
	p1 = OR(p1,p2), p2 = OR(p3,p4);
	Prince s;
	STORE(&s,OR(p1,p2));
	P = P2;
	p1 = OR(LOAD(P+0),LOAD(P+1));
	p2 = OR(LOAD(P+2),LOAD(P+3));
	p3 = OR(LOAD(P+4),LOAD(P+5));
	p4 = OR(LOAD(P+6),LOAD(P+7));
	p5 = OR(LOAD(P+8),LOAD(P+9));
	p6 = OR(LOAD(P+10),LOAD(P+11));
	p7 = OR(LOAD(P+12),LOAD(P+13));
	p8 = OR(LOAD(P+14),LOAD(P+15));
	p1 = OR(p1,p2), p2 = OR(p3,p4), p3 = OR(p5,p6), p4 = OR(p7,p8);
	p1 = OR(p1,p2), p2 = OR(p3,p4);
	Prince t;
	STORE(&t,OR(p1,p2));
	return (Split){FP_pop(&s),FP_pop(&t)};
}
	#define FP_union16x2 FP_union16x2_AVX
#else
	#define FP_union16x2 FP_union16x2_NV
#endif
inline uint64_t QRand64(uint64_t *x) {
	*x^=*x<<13; *x^=*x>>7; return *x^=*x<<17;}
	//(*x)^=(*x)>>11; (*x)^=(*x)<<37; (*x)^=(*x)>>4; return *x;}
			
static inline void create_sse2_profiles(Reference_Data *Rd) {
	if (Xalpha) return;
	puts("WARNING: program was built without SSSE3. This may be slow.");
	uint32_t *ClumpLen = Rd->ClumpLen, numRclumps = Rd->numRclumps;
	DualCoil **RefClump = Rd->RefClump;
	double wtime = omp_get_wtime();
	printf("Preparing SSE2 profiles...\n"); //ProfClump
	DualCoil **ProfClump = malloc(numRclumps*sizeof(*ProfClump));
	if (!ProfClump) {fputs("OOM:ProfClump\n",stderr); exit(3);}
	for (uint32_t i = 0; i < numRclumps; ++i)  // for each ref
		ProfClump[i] = malloc(SCD*ClumpLen[i]*sizeof(*ProfClump[i]));
	#pragma omp parallel for
	for (uint32_t i = 0; i < numRclumps; ++i) { // for each ref
		for (uint32_t j = 0; j < ClumpLen[i]; ++j) // for each letter in the ref
			for (int k = 0; k < SCD; ++k) // for each possibility of query letter
				for (int z = 0; z < VECSZ; ++z) // against each parallel ref letter
					ProfClump[i][j*SCD+k].u8[z] = 
						SCORENVedN[RefClump[i][j].u8[z]*SCD + k];
	}
	printf("SSE2 Prep time: %f\n", omp_get_wtime() - wtime);
	Rd->ProfClump = ProfClump;
}

// New -- Cucasort mk2 inplace
inline int str_d_gt_0(char *A, char *B, uint32_t maxD) {
	//for (uint32_t i = 0; i < maxD; ++i)
	//	if (A[i] != B[i]) return A[i] > B[i];
	for (uint32_t i = 0; i < maxD; i+=16) {
		__m128i x = _mm_lddqu_si128((void*)A + i),
			y = _mm_lddqu_si128((void*)B + i);
		__m128i e = _mm_cmpeq_epi8(x,y);
		uint16_t t = _mm_movemask_epi8(e);
		if (t != 0xFFFF) {
			uint32_t p = __builtin_ctz(t^-1);
			if (i+p >= maxD) return 0;
			return A[i+p] > B[i+p]; //~t
		}
	}
	return 0;
}
inline void iSort(char **s, uint32_t N, uint32_t d, uint32_t maxD) { 
	for (uint32_t i = 1, j; i < N; ++i) {
		char *key = s[i];
		for (j = i; j && str_d_gt_0(s[j-1]+d,key+d,maxD-d); --j) s[j] = s[j-1]; 
		//if (j != i)
			//memmove(s + j + 1, s + j, sizeof(*s) * (i - j)),
			s[j] = key;
	}
}
void Cucasort(char **S, uint64_t N, uint32_t d, uint32_t maxD) {
	// Recursive termination
	if (N < 1 || d >= maxD) return;
	if (N < 32) iSort(S,N,d,maxD);
	uint64_t Bounds[16];
	{ // define scope compartment for temp vars
		uint64_t CacheO[17] = {0}, *Cache = CacheO + 1; //, Bounds[16];
		for (uint64_t i = 0; i < N; ++i) 
			++Cache[S[i][d]];
		*Bounds = *Cache;
		for (int i = 1; i < 16; ++i)
			Cache[i] += Cache[i-1],
			Bounds[i] = Cache[i]; // sep?
		--Cache;
		for (int i = 0; i < 16;) {
			if (Cache[i] >= Bounds[i]) {++i; continue;}
			char *cur = S[Cache[i]], *t;
			while (cur[d] != i)  // sync needed
				t = S[Cache[cur[d]]],
				S[Cache[cur[d]]++] = cur,
				cur = t;
			S[Cache[i]++] = cur;
		}
	} // scope variables destroyed (some)
	// second pass: scan and recurse on identical non-zero ranges
	for (int i = 1; i < 16; ++i) {
		char **NS = S+Bounds[i-1];
		uint64_t len = Bounds[i] - Bounds[i-1];
		if (len > 1) Cucasort(NS,len,d+1,maxD);
	}
}

//////////////////QS3
#define CUTOFF 32
#define MEDCUT 100
//#define SWAP(s,i,j) { char *t = s[i]; s[i] = s[j]; s[j] = t; }
inline void swap(char **s, int i, int j) 
	{ char *t = s[i]; s[i] = s[j]; s[j] = t; }
inline int m3(char **s, int ia, int ib, int ic, int d) {
	int va, vb, vc;
	if ((va=s[ia][d]) == (vb=s[ib][d])) return ia;
	if ((vc=s[ic][d]) == va || vc == vb) return ic;
	return va < vb ?
		(vb < vc ? ib : (va < vc ? ic : ia ) ) : 
		(vb > vc ? ib : (va < vc ? ia : ic ) ); 
} 
void Qs3w(char **s, unsigned n, int d, int maxD) {
	if (d >= maxD) return;
	if (n < CUTOFF) { iSort(s,n,d,maxD); return; } 
	unsigned pl = 0, pm = n >> 1u, z;
	int le, lt, gt, ge, v, pn = n-1;
	if (n > MEDCUT) {
		z = n >> 3u;
		pl = m3(s, pl, pl+z, pl + (z << 1u), d);
		pm = m3(s, pm-z, pm, pm+z, d);
		pn = m3(s, pn - (z << 1u), pn-z, pn, d);
	}
	pm = m3(s, pl, pm, pn, d);
	swap(s, 0, pm);
	v = s[0][d]; 
	for (le = 1; le < n && s[le][d] == v; le++);  
	if (le == n) { if (v) Qs3w(s, n, d+1, maxD); return; }
	lt = le; gt = ge = n-1;
	for (;;) {
		for ( ; lt <= gt && s[lt][d] <= v; lt++)
			if (s[lt][d] == v) swap(s, le++, lt);
		for ( ; lt <= gt && s[gt][d] >= v; gt--)
			if (s[gt][d] == v) swap(s, gt, ge--);
		if (lt > gt) break;
		swap(s, lt++, gt--);
	}
	int r = MIN(le, lt-le); 
	for (int i = 0; i < r; ++i) swap(s,i,lt-r+i);
	r = MIN(ge-gt, n-ge-1);
	for (int i = 0; i < r; ++i) swap(s,lt+i,n-r+i);
	Qs3w(s, lt-le, d, maxD);
	if (v) Qs3w(s + lt-le, le + n-ge-1, d+1, maxD);
	Qs3w(s + n-(ge-gt), ge-gt, d, maxD); 
}
/////////////////QS3
inline uint32_t whereDiff(char *A, char *B, uint32_t len) {
	for (uint32_t i = 0; i < len; i+=16) {
		__m128i x = _mm_lddqu_si128((void*)A + i),
			y = _mm_lddqu_si128((void*)B + i);
		__m128i e = _mm_cmpeq_epi8(x,y);
		uint16_t t = _mm_movemask_epi8(e);
		if (t != 0xFFFF) {
			uint32_t p = __builtin_ctz(t^-1);
			return i+p; //A[i+p] - B[i+p]; //~t
		}
	}
	return len;
}
inline uint32_t whereDiffMsk(char *A, char *B, uint32_t len) {
	for (uint32_t i = 0; i < len; i+=16) {
		__m128i x = _mm_lddqu_si128((void*)A + i),
			y = _mm_lddqu_si128((void*)B + i);
		__m128i e = _mm_cmpeq_epi8(_mm_and_si128(x,_mm_set1_epi8(15)),
			_mm_and_si128(y,_mm_set1_epi8(15)));
		uint16_t t = _mm_movemask_epi8(e);
		if (t != 0xFFFF) {
			uint32_t p = __builtin_ctz(t^-1);
			return i+p; //A[i+p] - B[i+p]; //~t
		}
	}
	return len;
}

// curate is 0 (no dedupe), 1 (dedupe), 2 (make db so skip sse2 profiles)
static inline void process_references(char *ref_FN, Reference_Data *Rd, uint32_t maxLenQ, int curate) {
	// variables for the db code
	char *SeqDump = 0; uint64_t numLet = 0; // *Offs = 0;
	if (Rd->dbType == QUICK) Rd->totR = parse_tl_fasta(ref_FN, &Rd->RefHead, &Rd->RefSeq, &Rd->RefLen);
	else Rd->totR = parse_tl_fasta_db(ref_FN, &Rd->RefHead, &SeqDump, &Rd->RefSeq, &numLet);
	printf("Parsed %u references.\n",Rd->totR);

	// Translate nucleotides by parallel register lookups
	if (!Xalpha && Rd->dbType == QUICK) for (uint32_t i = 0; i < Rd->totR; ++i) 
		translate16aln(Rd->RefSeq[i],Rd->RefLen[i]);
	else if (!Xalpha && Rd->dbType == DNA_16) translate16aln(SeqDump, numLet);

	// Shear the references to desired length + tail overlap
	uint32_t origR = Rd->totR, numRrebase = 0, *origRefLen = Rd->RefLen, *ReRefIx;
	char **origRefSeq = Rd->RefSeq, **origRefHead = Rd->RefHead;
	if (REBASE) {
		uint32_t minShear = maxLenQ / THRES,
			shear = minShear > REBASE_AMT ? minShear : REBASE_AMT,
			ov = minShear; 
		shear = shear < minShear ? minShear : shear; // at least minShear
		printf("\nInitiating database shearing procedure [shear %u, ov %u].\n",
			shear, ov);
		if (Rd->dbType == DNA_16) {
			#define NL 13
			static const uint64_t NLB = (uint64_t)1 << (NL*2);
			#if NL==16
			#define NIB (s[0]-1)<<30|(s[1]-1)<<28|(s[2]-1)<<26|(s[3]-1)<<24|(s[4]-1)<<22|(s[5]-1)<<20|\
				(s[6]-1)<<18|(s[7]-1)<<16|(s[8]-1)<<14|(s[9]-1)<<12|(s[10]-1)<<10|(s[11]-1)<<8|(s[12]-1)<<6|\
				(s[13]-1)<<4|(s[14]-1)<<2|(s[15]-1)
			#elif NL==15
			#define NIB (s[0]-1)<<28|(s[1]-1)<<26|(s[2]-1)<<24|(s[3]-1)<<22|(s[4]-1)<<20|\
				(s[5]-1)<<18|(s[6]-1)<<16|(s[7]-1)<<14|(s[8]-1)<<12|(s[9]-1)<<10|(s[10]-1)<<8|(s[11]-1)<<6|\
				(s[12]-1)<<4|(s[13]-1)<<2|(s[14]-1)
			#elif NL==14
			#define NIB (s[0]-1)<<26|(s[1]-1)<<24|(s[2]-1)<<22|(s[3]-1)<<20|\
				(s[4]-1)<<18|(s[5]-1)<<16|(s[6]-1)<<14|(s[7]-1)<<12|(s[8]-1)<<10|(s[9]-1)<<8|(s[10]-1)<<6|\
				(s[11]-1)<<4|(s[12]-1)<<2|(s[13]-1)
			#elif NL==13
			#define NIB (s[0]-1)<<24|(s[1]-1)<<22|(s[2]-1)<<20|\
				(s[3]-1)<<18|(s[4]-1)<<16|(s[5]-1)<<14|(s[6]-1)<<12|(s[7]-1)<<10|(s[8]-1)<<8|(s[9]-1)<<6|\
				(s[10]-1)<<4|(s[11]-1)<<2|(s[12]-1)
			#elif NL==12
			#define NIB (s[0]-1)<<22|(s[1]-1)<<20|\
				(s[2]-1)<<18|(s[3]-1)<<16|(s[4]-1)<<14|(s[5]-1)<<12|(s[6]-1)<<10|(s[7]-1)<<8|(s[8]-1)<<6|\
				(s[9]-1)<<4|(s[10]-1)<<2|(s[11]-1)
			#endif
			uint32_t shear16p5 = shear + ov; //((shear+ov) / 16u) * 16u + 5u;
			Rd->cparts += !Rd->cparts;
			uint32_t cp_range = origR / Rd->cparts + (origR % Rd->cparts != 0);
			printf("Using compressive optimization (%u partitions; ~%u refs).\n",Rd->cparts,cp_range);
			//printf("-->CompDB: shearOv conv = %u --> %u [NLB=%llu]\n",shear+ov,shear16p5,NLB);

			uint64_t *BcountO = calloc(NLB+3,sizeof(*BcountO));
			if (!BcountO) {fputs("OOM:BcountO\n",stderr); exit(3);}
			
			uint64_t maxChain = 0, maxSh = 0;
			double wtime;
			for (uint32_t rix = 0, pass = 0; rix < origR; rix += cp_range, ++pass) {
				wtime = omp_get_wtime();
				memset(BcountO,0,(NLB+3)*sizeof(*BcountO));
				uint64_t *Bcount = BcountO + 2;
				uint32_t red = MIN(origR,rix+cp_range);
				#pragma omp parallel for schedule(dynamic)
				for (uint32_t i = rix; i < red; ++i) { 
					char *so = origRefSeq[i];
					uint32_t len = origRefSeq[i+1] - origRefSeq[i] - 1;
					//printf("Seq i=%u, len=%u\n",i,len);
					if (len < shear16p5) continue;
					len -= shear16p5;
					for (uint32_t j = 0; j < len; ++j) {
						uint8_t *s = so + j;
						for (int k=0; k<NL; ++k) if (!s[k] || s[k] > 4) goto L4_P6;
						uint32_t nib = NIB;
						//if (nib >= NLB) printf("WARNING: nib %u >= max %u!\n",nib,NLB);
						#pragma omp atomic update
						++Bcount[nib]; 
						L4_P6:NULL;
					}
				}
				printf("[%u] First pass: populated bin counters [%f]\n",
					pass, omp_get_wtime()-wtime);

				wtime = omp_get_wtime();
				for (uint64_t i = 1; i <= NLB; ++i) Bcount[i] += Bcount[i-1];
				uint64_t totCnt = Bcount[NLB];
				printf("--> Out of the %llu original places, %llu are eligible.\n",numLet, totCnt);
				--Bcount;
				char **Ptrs = malloc(totCnt*sizeof(*Ptrs));
				if (!Ptrs) {fputs("OOM:Ptrs_X\n",stderr); exit(3);}
				printf("[%u] Second pass: shift count vector [%f]\n",
					pass, omp_get_wtime()-wtime);

				wtime = omp_get_wtime();
				#pragma omp parallel for schedule(dynamic)
				for (uint32_t i = rix; i < red; ++i) { 
					char *so = origRefSeq[i];
					uint32_t len = origRefSeq[i+1] - origRefSeq[i] - 1; // NEW
					if (len < shear16p5) continue;
					len -= shear16p5;
					for (uint32_t j = 0; j < len; ++j) {
						uint8_t *s = so + j;
						for (int k=0; k<NL; ++k) if (!s[k] || s[k] > 4) goto L4_P7;
						uint32_t nib = NIB;
						uint64_t ix;
						#pragma omp atomic capture
						ix = Bcount[nib]++;
						Ptrs[ix] = s; 
						L4_P7:NULL;
					}
				}
				printf("[%u] Third pass: add pointers to Ptr vector [%f]\n",
					pass, omp_get_wtime() - wtime);

				wtime = omp_get_wtime();
				--Bcount; 
				#pragma omp parallel for schedule(dynamic,1)
				for (uint64_t i = 0; i < NLB; ++i) {
					uint64_t len = Bcount[i+1]-Bcount[i];
					char **section = Ptrs + Bcount[i];
					Qs3w(section,len,NL,shear16p5);
					//Cucasort(section, len, 5, shear16p5);
				}
				printf("[%u] Fourth pass: sort segments [%f]\n", 
					pass, omp_get_wtime() - wtime);
				
				uint32_t eqlen = shear16p5 - NL, nibLen = 24 - NL;
				if (maxChain == 0 && maxSh == 0) {
					// Phase 2: duplicate checks
					uint64_t dupes = 0, tdupes = 0, useful = 3;
					wtime = omp_get_wtime();
					#pragma omp parallel for schedule(dynamic) reduction(+:dupes,tdupes) reduction(max:maxChain,maxSh)
					for (uint64_t i = 0; i < NLB; ++i) {
						uint64_t chain = 0, sh = 0;
						for (uint64_t j = Bcount[i]+1; j < Bcount[i+1]; ++j) {
							uint32_t where = whereDiff(Ptrs[j-1]+NL,Ptrs[j]+NL,eqlen);
							if (where >= nibLen) ++sh;
							else if (sh > maxSh) maxSh = sh;
						
							if (where >= eqlen) ++chain;
							else {
								dupes += chain >= useful;
								tdupes += chain > 0;
								if (chain > maxChain) maxChain = chain;
								chain = 0;
							}
						}
					}
					printf("[%u] Fifth pass: tally %llu=%f (%llu=%f) [%f]\n", pass, dupes,
						(double)dupes/totCnt, tdupes, (double)tdupes/totCnt, omp_get_wtime()-wtime);
					printf("--> Max chain = %llu, max sh = %llu\n",maxChain,maxSh);
				}
				
				// dynamic range compression and marking 
				uint64_t sh1 = sqrt(maxSh)/2, sh2 = sh1*4/3, sh3 = sh1*3;
				printf("--> Bounds: max %llu -> sh3 %llu -> sh2 %llu -> sh1 %llu\n",
					maxSh,sh3,sh2,sh1);
				
				wtime = omp_get_wtime();
				#pragma omp parallel for schedule(dynamic)
				for (uint64_t i = 0; i < NLB; ++i) {
					uint64_t chain = 0, sh = 0;
					for (uint64_t j = Bcount[i]+1; j < Bcount[i+1]; ++j) {
						uint32_t where = whereDiffMsk(Ptrs[j-1]+NL,Ptrs[j]+NL,eqlen);
						
						// First consider local matches
						if (where >= nibLen) ++sh;
						else {
							if (sh > sh1) {
								uint8_t conv = sh>=sh3 ? 3 : sh>=sh2 ? 2 : 1;
								for (uint64_t k = j; k >= j - sh; --k)
									*Ptrs[k-1] |= conv << 4;
							}
							sh = 0;
						}
						// Consider full dupes (takes priority)
						if (where >= eqlen) ++chain;
						else {
							if (chain) {
								uint64_t t = chain*2048/maxChain; 
								t = MIN(2048,t);
								uint8_t conv = 31 - __builtin_clz((uint32_t)t) + 4;
								for (uint64_t k = j; k >= j - chain; --k)
									*Ptrs[k-1] |= conv << 4;
							}
							chain = 0;
						}
					}
				}
				free(Ptrs);
				printf("[%u] Sixth pass: populate duplicates [%f]\n\n",
					pass, omp_get_wtime()-wtime);
			}
			printf("-->CompDB: All sequences tallied.\n");
			free(BcountO);  
			
			// Finally, use the markings to guide shearing
			wtime = omp_get_wtime();
			// goal: minimal shear size such that we start a shear at the max-ranked dupe
			// need to create reflens out of the offsets between files. 

			uint32_t *RefLen = Rd->RefLen = malloc(origR*sizeof(*RefLen));
			for (uint32_t i = 0; i < origR; ++i) {
				RefLen[i] = origRefSeq[i+1]-origRefSeq[i] - 1;
				long unit = (long)RefLen[i] - (long)ov;
				unit = unit <= 0 ? 1 : unit;
				numRrebase += unit / ov + (unit % ov != 0);
			}
			numRrebase *= 2; // worst-case with leeway
			//ReRefIx = malloc(numRrebase*sizeof(*ReRefIx));
			uint32_t *RefStart = malloc(numRrebase*sizeof(*RefStart));
			// Define (and check) new RefLen, RefSeq, and RefHead too for inline addition
			uint32_t *newRefLen = malloc(numRrebase*sizeof(*newRefLen));
			char **newRefSeq = malloc(numRrebase*sizeof(*newRefSeq)),
				**newRefHead = malloc(numRrebase*sizeof(*newRefHead));
			if (!RefStart || !newRefHead) {fputs("OOM:RefStart\n",stderr); exit(3);}
			printf("-->CompDB: Beginning rebase...\n");
			int x = 0; // "real" rebase ix
			for (uint32_t i = 0; i < origR; ++i) {
				//if (RefLen[i] > 300000000) {printf("\nWarning: Reference %u is long [%u]\n",i,RefLen[i]);}
				uint8_t *ref = (uint8_t *)origRefSeq[i];
				uint32_t bstFlgPos = 0, end = 0, bstFlg = *ref >> 4;
				while (end < RefLen[i]) {
					uint32_t ix;
					//#pragma omp atomic capture
					ix = x++;
					if (ix >= numRrebase) {fputs("ERROR: Rebase overflow.\n",stderr); exit(4);}
					RefStart[ix] = bstFlgPos;
					newRefSeq[ix] = ref + bstFlgPos;
					newRefHead[ix] = origRefHead[i];

					//uint32_t bstFlgO = bstFlg;

					// scan for best flag in this range
					uint32_t bf=0;
					uint32_t bi, maxIX = MIN(RefLen[i],bstFlgPos+shear);
					for (uint32_t j = bstFlgPos+1; j < maxIX; ++j)
						if (ref[j] >> 4 >= bf) bf = ref[j]>>4, bi = j;
					if (bf > bstFlg /* && (bi > bstFlgPos + LATENCY || bf > 3) */)  // terminate here + ov
						bstFlgPos = bi;
					else bstFlgPos += shear;
					end = bstFlg > 3 ? MIN(maxIX+ov,RefLen[i]) : MIN(bstFlgPos + ov,RefLen[i]);
					//if (bstFlgPos >= RefLen[i]) {printf("\rWARNING: Poor flag position at %u [%u / %u]",i,bstFlgPos,RefLen[i]);}
					if (bstFlgPos < RefLen[i]) bstFlg = ref[bstFlgPos] >> 4;
					//end = MIN(bstFlgPos + ov,RefLen[i]); // end of this shear
					newRefLen[ix] = end - RefStart[ix];
					//if (!i) printf("Shear %u: start %u [%u] -> bstpos %u end %u [%u]\n", ix, RefStart[ix],bstFlgO,bstFlgPos,end,bstFlg);
				}
			}
			printf("\nRebased [%u orig] --> [%u rebase] --> [%u adj]\n",origR,numRrebase/2,x);
			
			//ReRefIx = realloc(ReRefIx,x*sizeof(*ReRefIx));
			numRrebase = x;
			newRefLen = realloc(newRefLen,x*sizeof(*newRefLen));
			newRefSeq = realloc(newRefSeq,x*sizeof(*newRefSeq));
			newRefHead = realloc(newRefHead,x*sizeof(*newRefHead));
			RefStart = realloc(RefStart,x*sizeof(*RefStart));
			Rd->totR = numRrebase;
			Rd->RefLen = newRefLen;
			Rd->RefSeq = newRefSeq;
			Rd->RefHead = newRefHead;
			Rd->RefStart = RefStart;
			Rd->maxLenR = shear + ov;
			
			// reporting (DEBUG ONLY)
			/*printf("DEBUG: reporting each reference [%u, total %llu] and the number of 'high bases' in it\n",origR,numLet);
			uint64_t total_length_debug = 0, total_high_debug = 0, total_afterimage_debug = 0;
			for (uint32_t i = 0; i < origR; ++i) {
				uint32_t high = 0;
				for (uint32_t j = 0; j < RefLen[i]; ++j) high += (uint8_t)origRefSeq[i][j] >= 16;
				//printf("Reference %u [len %u]: %u high\n",i, RefLen[i],high);
				for (uint32_t j = 0; j < RefLen[i]; ++j) if ((uint8_t)origRefSeq[i][j] >= 16) 
					origRefSeq[i][j] = (uint8_t)origRefSeq[i][j] & 15;
				total_high_debug += high;
				high = 0;
				for (uint32_t j = 0; j < RefLen[i]; ++j) high += (uint8_t)origRefSeq[i][j] >= 16;
				//printf("-->Afterward: %u high\n",high);
				total_afterimage_debug += high;
				total_length_debug += RefLen[i];
			}

			printf("Total length = %llu, num high = %llu, afterimage = %llu\n",
				total_length_debug,total_high_debug,total_afterimage_debug);
			
			total_high_debug = 0;
			printf("\n\n\nDEBUG: reporting each sheared reference and the number of 'high bases' in them\n");
			for (uint32_t i = 0; i < numRrebase; ++i) {
				uint32_t high = 0;
				for (uint32_t j = 0; j < newRefLen[i]; ++j) high += (uint8_t)newRefSeq[i][j] >= 16;
				//printf("SHEARED Reference %u [len %u]: %u high\n",i, newRefLen[i],high);
				total_high_debug += high;
			}
			printf("High afterward = %llu\n",total_high_debug);
			*/
			
			#pragma omp parallel for
			for (uint64_t i = 0; i < numLet; ++i) SeqDump[i] &= 15;

			free(RefLen);
			printf("Duplicate-guided shearing [%f]\n",omp_get_wtime()-wtime);
			
			//exit(190);
			
			// Final check: defined numRrebase, all RefStart populated; origRefLen+origRefSeq and Rd pointers to them
		}
		else {
			for (uint32_t i = 0; i < origR; ++i) {
				long unit = (long)Rd->RefLen[i] - (long)ov;
				unit = unit < 0 ? 1 : unit;
				numRrebase += unit / shear + (unit % shear != 0);
			}
			ReRefIx = malloc(numRrebase*sizeof(*ReRefIx)); 
			Rd->RefStart = malloc(numRrebase*sizeof(*Rd->RefStart));
			if (!ReRefIx || !Rd->RefStart) {fputs("OOM:ReRefIx\n",stderr); exit(3);}
			for (uint32_t i = 0, x = 0; i < origR; ++i) {
				long unit = (long)Rd->RefLen[i] - (long)ov;
				unit = unit < 0 ? 1 : unit;
				for (uint32_t j = 0; j < unit; j+= shear) 
					ReRefIx[x] = i, Rd->RefStart[x++] = j ? j : 0;
			}
			
			printf("Shorn refs: %u, rebased clumps: %u\n",numRrebase, numRrebase/16 + (numRrebase%16 != 0));
		
			// Redefine existing variables to sheared versions
			Rd->totR = numRrebase;
			Rd->RefLen = malloc(Rd->totR*sizeof(*Rd->RefLen));
			Rd->RefSeq = malloc(Rd->totR*sizeof(*Rd->RefSeq));
			Rd->RefHead = malloc(Rd->totR*sizeof(*Rd->RefHead));
			if (!Rd->RefLen || !Rd->RefSeq || !Rd->RefHead) 
				{fputs("OOM:RefLen\n",stderr); exit(3);}
			Rd->maxLenR = shear + ov; // may actually be less...
			for (uint32_t i = 0; i < Rd->totR; ++i) 
				Rd->RefHead[i] = origRefHead[ReRefIx[i]],
				Rd->RefSeq[i] = origRefSeq[ReRefIx[i]] + Rd->RefStart[i],
				Rd->RefLen[i] = origRefLen[ReRefIx[i]] - Rd->RefStart[i],
				Rd->RefLen[i] = Rd->RefLen[i] > Rd->maxLenR ? Rd->maxLenR : Rd->RefLen[i];
			free(ReRefIx);
		}
		

	}

	Rd->RefIxSrt = malloc((1+Rd->totR) * sizeof(*Rd->RefIxSrt));
	if (!Rd->RefIxSrt) {fputs("OOM:RefIxSrt\n",stderr); exit(3);}
	// Sort (shorn) references within length-adjusted pods
	if (LATENCY) {
		// sort refs by length (TODO: Better sort)
		Tuxedo *SeqLenPair = malloc(Rd->totR * sizeof(*SeqLenPair));
		if (!SeqLenPair) {fputs("OOM:SeqLenPair\n",stderr); exit(3);}
		for (uint32_t i = 0; i < Rd->totR; ++i) {
			SeqLenPair[i].seq = Rd->RefSeq[i]; 
			SeqLenPair[i].len = Rd->RefLen[i];
			SeqLenPair[i].ix = i;
		}
		qsort(SeqLenPair, Rd->totR, sizeof(*SeqLenPair),cmpPackLen);

		// Figure out the max reference length while at it
		Rd->maxLenR = SeqLenPair[Rd->totR-1].len;

		//Now sort refs lexicographically within tolerance range. 
		uint32_t curTol = SeqLenPair[0].len, prev_ix = 0, lat = DO_FP ? 0 : LATENCY;
		for (uint32_t i = 1; i < Rd->totR; ++i) {
			if (SeqLenPair[i].len > curTol + lat) {
				curTol = SeqLenPair[i].len;
				//qsort(SeqLenPair + prev_ix, i - prev_ix, sizeof(*SeqLenPair),cmpPackSeq);
				if (i - prev_ix > 1) {
					Tuxedo *thisTux = SeqLenPair + prev_ix;
					uint32_t range = i - prev_ix;
					//printf("Sorting on i=%u, prev = %u (range=%u)\n",i,prev_ix,range);
					if (range > 256) parallel_sort_tuxedo(thisTux, range);
					else qsort(thisTux, range, sizeof(*SeqLenPair),cmpPackSeq);
				}
				prev_ix = i;
			}
		}
		if (prev_ix < Rd->totR-1)
			//qsort(SeqLenPair + prev_ix, Rd->totR - prev_ix, sizeof(*SeqLenPair),cmpPackSeq);
			if (Rd->totR - prev_ix > 1) //printf("Sorting final (range=%u)\n",Rd->totR - prev_ix);
				parallel_sort_tuxedo(SeqLenPair + prev_ix, Rd->totR - prev_ix);
		
		for (uint32_t i = 0; i < Rd->totR; ++i) Rd->RefIxSrt[i] = SeqLenPair[i].ix;
		free(SeqLenPair);
		
	} else for (uint32_t i = 0; i < Rd->totR; ++i) 
		Rd->maxLenR = Rd->RefLen[i] > Rd->maxLenR ? Rd->RefLen[i] : Rd->maxLenR,
		Rd->RefIxSrt[i] = i;
	Rd->origTotR = Rd->totR, Rd->TmpRIX = Rd->RefIxSrt;

	if (curate) {   // dedupe references (redefine RefIxSrt using RefDedupIx map)
		uint64_t tot_divR = 0, rdupes = 0;
		uint32_t uix = 0; // dupe-map, map[uniq_ix] -> orig_ix
		printf("Producing edb scaffold...\n");
		uint32_t totR = Rd->totR, *RefIxSrt = Rd->RefIxSrt, *RefLen = Rd->RefLen,
			*RefDedupIx = calloc(Rd->totR+1,sizeof(*Rd->RefDedupIx));
		char **RefSeq = Rd->RefSeq;
		if (!RefDedupIx) {fputs("OOM:RefDedupIx\n",stderr); exit(3);}
		for (uint32_t i = 1, j; i < totR; ++i) {
			uint32_t curIx = RefIxSrt[i], prevIx = RefIxSrt[i-1],
				maxL = MIN(RefLen[curIx],RefLen[prevIx]);
			for (j = 0; j < maxL; ++j) 
				if (RefSeq[curIx][j] != RefSeq[prevIx][j]) break;
			tot_divR += j;
			if (j == RefLen[curIx] && j == RefLen[prevIx]) ++rdupes; 
			else RefDedupIx[++uix] = i;
		}
		RefDedupIx[++uix] = totR; // add an end-cap for closure
		Rd->RefDedupIx = RefDedupIx;
		printf("Avg. R Divergence: %f (%u dupes, %u uniq)\n",(double)tot_divR/totR,rdupes,uix);

		// Ensure the first db hit is maintained as the unique representative
		for (uint32_t i = 0; i < uix; ++i) { // not nec. if stable sort?
			uint32_t bix = Rd->RefIxSrt[Rd->RefDedupIx[i]];
			for (uint32_t mi = Rd->RefDedupIx[i] + 1; mi < Rd->RefDedupIx[i+1]; ++mi) 
				if (Rd->RefIxSrt[mi] < bix) bix = Rd->RefIxSrt[mi],
					Rd->RefIxSrt[mi] = Rd->RefIxSrt[Rd->RefDedupIx[i]],
					Rd->RefIxSrt[Rd->RefDedupIx[i]] = bix;
		}

		// redefine totR, RefIxSrt (the mask into RefSeq and RefLen)
		Rd->TmpRIX = malloc((1+uix)*sizeof(*Rd->TmpRIX));
		for (uint32_t i = 0; i < uix; ++i) Rd->TmpRIX[i] = Rd->RefIxSrt[Rd->RefDedupIx[i]];
			// can theoretically not make a TmpRIX, just use old on both sides ^
		uint32_t *tmp = Rd->RefIxSrt;
		Rd->RefIxSrt = Rd->TmpRIX, Rd->TmpRIX = tmp; // now RefIxSrt has latest, old is backup
		//Rd->TmpRIX = Rd->RefIxSrt, Rd->RefIxSrt = 0; // can recalc during alignment
		Rd->totR = uix;
	}
	
	// Create vector clumps of references
	uint32_t numRclumps = Rd->totR/VECSZ, totRC = numRclumps + (numRclumps*VECSZ < Rd->totR);
	printf("There are %u references and hence %u clumps (+%u)\n",
		Rd->totR, numRclumps, numRclumps*VECSZ < Rd->totR);
	
	if (DO_FP) {
		// The goal: Ensure each successive 16 (VECSZ) refs are reasonably similar to one another
		PackaPrince Fp = create_fingerprints(Rd->RefSeq, Rd->totR, Rd->RefLen, Rd->RefIxSrt, 1 ,1);
		//Prince *restrict p = __builtin_assume_aligned(Fp.P,64);
		Prince *p = Fp.P;
		printf("How many fingerprints made = %u\n",Fp.nf);
		if (Z) { // Cluster by N-free profile if N-penalize
			puts("Note: penalizing N's in fingerprints may slow unpenalized alignments.");
			Prince t; // New: swap to error-free fingerprints anyway (higher quality clusters for N-case and without)
			for (uint32_t i = 0; i < Rd->totR; ++i) t = p[i], p[i] = p[Fp.Ptrs[i]], p[Fp.Ptrs[i]] = t; // swap FPs
		}
		
		// prelim stats: individual coverage.
		uint32_t *RefIxSrt = Rd->RefIxSrt, *RefLen = Rd->RefLen, totR = Rd->totR;
		
		uint64_t tot_cl = 0;
		double similarity = 0.06;
		#pragma omp parallel for reduction(+:tot_cl)
		for (uint32_t i = 0; i < totR; ++i) tot_cl += FP_pop(p+i);
		uint32_t sig_thres = Rd->clustradius; // + similarity * (double)tot_cl/totR;
		//if (!Rd->clustradius && sig_thres < 5) sig_thres = 5;
		//else if (sig_thres > 999) similarity = (double)(sig_thres-1000)/1000, sig_thres = 0;
		#define GMIN 250
		//sig_thres = MIN(GMIN,sig_thres); // make this MING
		printf("Average coverage (atomic) = %llu [%f], cluster refinements: %u\n",tot_cl,(double)tot_cl/totR, sig_thres);
		//#define PTHRES2 240
		uint32_t IXBIN[THREADS];
		
		
		// The direct fingerprint way
		double wtime = omp_get_wtime();
		#define BANDED_FP 1048576
		// Heuristic: sort FPs by first few bits (Prince.h), then iterate over bounded range
		if (BANDED_FP) {
			uint32_t *WordRange = malloc(totR*sizeof(*WordRange));
			uint32_t *RefLen = Rd->RefLen;
			uint32_t b = 0, z = totR;
			//for (uint32_t z = 1, b = 0; z <= totR; ++z) {
				//if (z == totR || RefLen[RefIxSrt[z]] > RefLen[RefIxSrt[b]] + LATENCY) {
					size_t pow = 24, binCnt = 1 << pow, rem = 32-pow;
					uint32_t *Pops = calloc(binCnt+2,sizeof(*Pops)), *Cnt_p = Pops + 1;
					//#pragma omp parallel for
					for (uint32_t i = b; i < z; ++i)
						//#pragma omp atomic
						++Cnt_p[*p[i].s>>rem];
					--Cnt_p;
					for (uint32_t i = 1; i < binCnt+2; ++i) 
						Cnt_p[i] += Cnt_p[i-1];
					for (uint32_t i = b; i < z; ++i)
						WordRange[Cnt_p[*p[i].s>>rem]+++b] = i;
					b = z;
					free(Pops);
					puts("Prepared cluster bands");
				//}
			//}
			// Reorganize fingerprints as well as the original references (RIX, TmpRIX, dedup)
			if (curate) { // generate new RefIxSrt, RefDedupIx, Original (TmpRIX)
				uint32_t *TmpRIX = Rd->TmpRIX, *NewOrig = malloc(Rd->origTotR*sizeof(*NewOrig)),
					*NewDedup = Rd->RefIxSrt, *RefDedupIx = Rd->RefDedupIx;
				if (!NewOrig) {fputs("OOM:NewOrig:BFP\n",stderr); exit(3);}
				uint32_t j = 0; for (uint32_t i = 0; i < totR; ++i) {
					NewDedup[i] = j;
					for (uint32_t z = RefDedupIx[WordRange[i]]; z < RefDedupIx[WordRange[i]+1]; ++z)
						NewOrig[j++] = TmpRIX[z];
				}
				if (j < Rd->origTotR) puts("Coalescing origin..."), NewOrig[j] = TmpRIX[Rd->origTotR-1];
				NewDedup[totR] = Rd->origTotR;
				free(Rd->TmpRIX);
				Rd->RefDedupIx = NewDedup, Rd->TmpRIX = NewOrig;
				for (uint32_t i = 0; i < totR; ++i) RefDedupIx[i] = NewOrig[NewDedup[i]];
				Rd->RefIxSrt = RefDedupIx;
			}
			else { // just generate a new RefIxSrt
				uint32_t *TmpRIX = malloc(totR*sizeof(*TmpRIX));
				if (!TmpRIX) {fputs("OOM:TmpRIX:BFP\n",stderr); exit(3);}
				for (uint32_t i = 0; i < totR; ++i) 
					TmpRIX[i] = RefIxSrt[WordRange[i]];
				free(RefIxSrt);	
				Rd->TmpRIX = TmpRIX, Rd->RefIxSrt = TmpRIX;
				RefIxSrt = TmpRIX;
			}

			// Deal with scattering the FP's and their pointers
			uint32_t totFPs = Fp.nf;
			uint32_t *NewPix = 0;
			Prince *NewP = malloc(totFPs*sizeof(*NewP));
			if (!NewP) {fputs("OOM:NewP:BFP\n",stderr); exit(3);}

			//if (Z) totFPs = totR; //, Fp.nf = 0;
			//else {
				NewPix = malloc(totR*sizeof(*NewPix));
				if (!NewPix) {fputs("OOM:NewPix:BFP\n",stderr); exit(3);}
				for (uint32_t i = 0; i < totR; ++i)
					NewPix[i] = Fp.Ptrs[WordRange[i]] >= totR ? Fp.Ptrs[WordRange[i]] : i;
				for (uint32_t i = totR; i < totFPs; ++i) NewP[i] = p[i];
			//}
			for (uint32_t i = 0; i < totR; ++i) NewP[i] = p[WordRange[i]];
			free(Fp.initP); 
			free(Fp.Ptrs); 
			p = NewP; 
			Fp.initP = Fp.P = p;
			Fp.Ptrs = NewPix;
			Rd->FingerprintsR = Fp;

			printf("Time to preprocess: %f\n",omp_get_wtime()-wtime);
			wtime = omp_get_wtime();
			// TODO: use just this grouping as "fast clustering"
		}

		// Even newer clusterer! ME-based in pseudorandom local neighborhood
		// Idea: Build cluster using ME/convergence 2^n members at a time over n iters.

		uint32_t tot16 = totR + (15 & (16 - (totR & 15))); //mod16
		uint32_t *IxArray = malloc(tot16*sizeof(*IxArray));
		uint32_t *ShfIx = malloc(tot16*sizeof(*ShfIx));
		uint32_t *Cache = malloc(tot16*sizeof(*Cache));
		if (!IxArray || !ShfIx || !Cache) {fputs("OOM:ShfIx\n",stderr); exit(3);}
		printf("TotR = %u, tot16 = %u\n",totR, tot16);
		wtime = omp_get_wtime();
		// Setup: populate IxArray and ShfIx in order
		for (uint32_t i = 0; i < tot16; ++i) IxArray[i] = i;
		printf("Preparing L(A(t)) L(8)...%f\n",omp_get_wtime() - wtime);

		uint64_t totalPop = 0;

		#define doRand 1
		//uint32_t *Rand_Ix = malloc((tot16/2+UINT16_MAX+2)*sizeof(*Rand_Ix));
		uint64_t mseed = rand();
		
		// prepare new FP cache for greedy clusterer
		void *P_init = 0; Prince *P = calloc_a(64,(tot16+1)*sizeof(*P),&P_init);
		if (!P_init) {fputs("OOM:P_new\n",stderr); exit(3);}
		for (uint32_t i = 0; i < totR; ++i) 
			P[i] = p[i];
		
		uint32_t RI[THREADS]; for (int z = 0; z < THREADS; ++z) RI[z] = -1;
		//void *PC_init = 0; Prince *PC = calloc_a(64,(tot16 >> 4)*sizeof(*PC),&PC_init);
		Prince *PC = calloc(tot16 >> 4,sizeof(*PC));
		if (!PC) {fputs("OOM:PC_init\n",stderr); exit(3);}
		// naive greedy clusterer
		{
		double wtime = omp_get_wtime();
		#ifdef CUBICLUST
		// cubic zirconiu(m;^3)
		uint32_t *T_Pack[THREADS], *T_Ppk[THREADS], 
			hshSz = (MIN(totR,BANDED_FP) >> 3) + 7;
		uint8_t *T_Hash[THREADS];
		#pragma omp parallel
		T_Hash[omp_get_thread_num()] = calloc(hshSz,1),
		T_Pack[omp_get_thread_num()] = calloc(16,sizeof(uint32_t)),
		T_Ppk[omp_get_thread_num()] = calloc(16,sizeof(uint32_t));
		
		for (uint32_t i = 0; i < totR; i+=16) {
			uint32_t gmin = -1, gdst = -1, 
				range1 = MIN(totR,i+THREADS), range2 = MIN(totR, i+BANDED_FP);
			#pragma omp parallel
			{
				int tid = omp_get_thread_num();
				uint8_t *Hash = T_Hash[tid];
				uint32_t *Pack = T_Pack[tid], *Ppk = T_Ppk[tid];
				uint32_t pmin = -1, pdst = -1, *t;
				
				#pragma omp for reduction(min:gmin)
				for (uint32_t j = i; j < range1; ++j) {
					Pack[0] = j, Hash[(j-i)>>3] = 1 << ((j-i) & 7);
					Prince centroid = P[j];
					uint32_t min = -1, dst = -1, totD = 0; 
					for (int z = 1; z < 16; ++z) { 
						min = -1, dst = -1; 
						for (uint32_t k = i; k < range2; ++k) {
							uint32_t m = k-i;
							if (Hash[m>>3] && Hash[m>>3] & (1 << (m&7))) continue;
							uint32_t t = FP_union(centroid,P[k]), td;
							if (t < min) min = t, Pack[z] = k, dst = FP_dist(centroid,P[k]);
							else if (t == min && (td=FP_dist(centroid,P[k])) < dst)
								Pack[z] = k, dst = td;
						}
						if (min != -1) centroid = FP_combine(centroid,P[Pack[z]]), totD += dst,
							Hash[(Pack[z]-i)>>3] |= 1 << ((Pack[z]-i) & 7);
						else Pack[z] = i;
					}
					for (int z = 0; z < 16; ++z) Hash[(Pack[z]-i)>>3] = 0;
					
					if (min < gmin) 
						gmin = pmin = min, pdst = totD, 
						t = Ppk, Ppk = Pack, Pack = t;
					else if (min == gmin && totD < pdst) 
						pdst = totD, t = Ppk, Ppk = Pack, Pack = t;
				}
				T_Pack[tid] = Pack, T_Ppk[tid] = Ppk;
				RI[tid] = pmin == gmin ? pdst : -1;
			}
			uint32_t whichZ = -1;
			for (int z = 0; z < THREADS; ++z) if (RI[z] != -1) { //if (RI[z]) { 
				//if (whichZ == -1 || *T_Ppk[z] < *T_Ppk[whichZ]) whichZ = z;
				if (RI[z] < gdst) gdst = RI[z], whichZ = z; 
				else if (RI[z] == gdst && (*T_Ppk[z] < *T_Ppk[whichZ])) whichZ = z;
			}
			uint32_t *Score = T_Ppk[whichZ];
			for (int z = 0; z < 16; ++z) {
				if (i+z >= totR) break;
				Prince tp = P[i+z]; P[i+z] = P[Score[z]]; P[Score[z]] = tp;
				uint32_t x = IxArray[i+z]; IxArray[i+z] = IxArray[Score[z]]; 
				IxArray[Score[z]] = x;
			}
			if (i & 255) printf("\rClustered %u...\r",i);
		}
		for (uint32_t i = 0; i < tot16; i+=16) {
			Prince charming = P[i];
			for (uint32_t j = 1; j < 16; ++j) 
				charming = FP_combine(charming,P[i+j]);
			PC[i >> 4] = charming;
		}
		#else
		Prince centroid = *P;
		for (uint32_t j = 1; j < totR; ++j) {
			uint32_t min = -1, mix = -1, dst = -1; 
			#pragma omp parallel
			{
				uint32_t tix = -1, tmin = -1, tdst = -1;
				#pragma omp for reduction(min:min) 
				for (uint32_t k = j; k < MIN(totR,j+BANDED_FP); ++k) {
					uint32_t t = FP_union(centroid,P[k]), td;
					if (t < min) min = tmin = t, tix = k, tdst = FP_dist(centroid,P[k]);
					else if (t == min && (td=FP_dist(centroid,P[k])) < tdst)
						tix = k, tdst = td;
				}
				RI[omp_get_thread_num()] = tmin == min ? tix : -1;
			}
			uint32_t t;
			for (uint32_t k = 0; k < THREADS; ++k)
				if (RI[k] != -1 && ((t=FP_dist(centroid,P[RI[k]])) < dst || 
					(t==dst && RI[k] < mix))) mix = RI[k];
			centroid = FP_combine(centroid,P[mix]); // amend centroid
			Prince tp = P[j]; P[j] = P[mix]; P[mix] = tp;
			uint32_t x = IxArray[j]; IxArray[j] = IxArray[mix]; IxArray[mix] = x;
			if (!((j+1) & 15)) PC[j >> 4] = centroid, centroid = P[j+1];
			if (j & 255) printf("\rClustered %u...\r",j);
			if (totR < tot16) PC[totR >> 4] = centroid;
		}
		#endif
		
		printf("\nDone with greedy (%f)\n", omp_get_wtime() - wtime);
		totalPop = 0;
		for (uint32_t i = 0; i < tot16 >> 4; ++i) totalPop += FP_pop(PC+i);
		printf("Ended at %llu [%f] density.\n",totalPop,(double)(totalPop << 4)/totR);
		#pragma omp parallel for
		for (uint32_t i = 0; i < tot16; ++i) 
			P[i] = IxArray[i] >= totR ? (Prince){0} : p[IxArray[i]];
		uint32_t *ClusFPs = malloc((tot16 >> 4)*sizeof(*ClusFPs));
		#pragma omp parallel for
		for (uint32_t i = 0; i < tot16; i+=16) {
			Prince charming = P[i];
			for (uint32_t j = 1; j < 16; ++j) 
				charming = FP_combine(charming,P[i+j]);
			ClusFPs[i >> 4] = FP_pop(&charming);
		}
		totalPop = 0;
		for (uint32_t i = 0; i < tot16 >> 4; ++i) totalPop += ClusFPs[i];
		free(ClusFPs);
		printf("Raw check: %llu [%f].\n",totalPop,(double)(totalPop << 4)/totR);
		
		totalPop = 0;
		#pragma omp parallel for reduction(+:totalPop)
		for (uint32_t i = 0; i < tot16; ++i) totalPop += FP_pop(P+i);
		printf("Individual {DEBUG} score = %llu [%f]\n",totalPop,(double)totalPop/totR);
		
		}
		// end of greedy clusterer
		
		if (sig_thres) { // Cluster refinement loop
			uint32_t *ClusFPs = malloc((tot16 >> 4)*sizeof(*ClusFPs));
			#pragma omp parallel for
			for (uint32_t i = 0; i < tot16; ++i) 
				P[i] = IxArray[i] >= totR ? (Prince){0} : p[IxArray[i]];
			#pragma omp parallel for
			for (uint32_t i = 0; i < tot16; i+=16) {
				Prince charming = P[i];
				for (uint32_t j = 1; j < 16; ++j) 
					charming = FP_combine(charming,P[i+j]);
				ClusFPs[i >> 4] = FP_pop(&charming);
			}
			totalPop = 0;
			for (uint32_t i = 0; i < tot16 >> 4; ++i) totalPop += ClusFPs[i];
			printf("Raw check: %llu [%f].\n",totalPop,(double)(totalPop << 4)/totR);
		
			#pragma omp parallel for 
			for (uint32_t i = 0; i < tot16 >> 4; ++i) ClusFPs[i] = FP_pop(PC+i);
			
			/* for (uint32_t i = 0; i < tot16; i+=16) {
				Prince charming = P[i];
				for (uint32_t j = 1; j < 16; ++j) 
					charming = FP_combine(charming,P[i+j]);
				ClusFPs[i >> 4] = FP_pop(&charming);
			} */
			totalPop = 0;
			for (uint32_t i = 0; i < tot16 >> 4; ++i) totalPop += ClusFPs[i];
			printf("Starting at %llu [%f] for final round.\n",totalPop,(double)(totalPop << 4)/totR);
			
			#pragma omp parallel for
			for (uint32_t i = 0; i < tot16; ++i) ShfIx[i] = i;
			
			uint32_t tot2 = (tot16 >> 4) - ((tot16 >> 4) & 1);
			wtime = omp_get_wtime();
			for (uint32_t i = 0; i < sig_thres; ++i) {
				#pragma omp parallel
				{
					uint64_t seed = omp_get_thread_num() + 1 + mseed;
					#pragma omp for
					for (uint32_t z = 0; z < (tot16 >> 4)-1; ++z) {
						uint32_t r = QRand64(&seed);
						Cache[z] = r % ((tot16 >> 4)-z) + z; 
					}
					#pragma omp single
					{
						mseed = seed;
						for (uint32_t z = 0; z < (tot16 >> 4); ++z) {
							register uint32_t t = ShfIx[z], r = Cache[z]; 
							ShfIx[z] = ShfIx[r], ShfIx[r] = t;
						}
					}
					#pragma omp for schedule(dynamic,1)
					for (uint32_t j = 0; j < tot2; j+=2) {
						uint32_t c1 = ShfIx[j], c2 = ShfIx[j+1], 
							c1o = c1 << 4, c2o = c2 << 4;
						// Exhaustive 256-way swap
						uint32_t r1 = MIN(totR,c1o+16),
							r2 = MIN(totR,c2o+16);
						for (uint32_t k = c1o; k < r1; ++k) {
							for (uint32_t m = c2o; m < r2; ++m) {
								Prince tp = P[k]; P[k] = P[m]; P[m] = tp;
								Split re = FP_union16x2(P+c1o, P+c2o);
								if (re.v + re.i < ClusFPs[c1] + ClusFPs[c2]) { 
									ClusFPs[c1] = re.v, ClusFPs[c2] = re.i;
									uint32_t t = IxArray[k]; IxArray[k] = IxArray[m]; IxArray[m] = t;
								} else P[m] = P[k], P[k] = tp; // swap back!
							}
						}
					}
				}
				if (!(i & 511)) {
					uint64_t sv = 0;
					#pragma omp parallel for simd reduction(+:sv)
					for (uint32_t i = 0; i < tot16>>4; ++i) sv += ClusFPs[i];
					printf("[Round %u] savings = %llu, cur = %llu\n",i,totalPop-sv,sv);
					totalPop = sv;
				}
			}
			totalPop = 0;
			#pragma omp parallel for reduction(+:totalPop)
			for (uint32_t i = 0; i < tot16 >> 4; ++i) totalPop += ClusFPs[i];
			printf("Finished at %llu [%f] for final round [%f].\n",totalPop,(double)(totalPop << 4)/totR,omp_get_wtime()-wtime);
			
			free(ClusFPs);
		}
		uint32_t *ClusIX = IxArray; // reassign!
		Prince *UnionPrince = PC;   // reassign!
		free(P_init); 
		free(ShfIx); free(Cache);
		

		// New heirarchical clusterer: faster, more ram-friendly 
		// ... removed 
		
		//#pragma omp sections
		{
			//#pragma omp section
			if (curate) { // generate new RefIxSrt, RefDedupIx, Original (TmpRIX)
				uint32_t *TmpRIX = Rd->TmpRIX, *NewOrig = malloc(Rd->origTotR*sizeof(*NewOrig)),
					*NewDedup = Rd->RefIxSrt, *RefDedupIx = Rd->RefDedupIx;
				if (!NewOrig) {fputs("OOM:NewOrig\n",stderr); exit(3);}
				uint32_t j = 0; for (uint32_t i = 0; i < totR; ++i) {
					NewDedup[i] = j;
					for (uint32_t z = RefDedupIx[ClusIX[i]]; z < RefDedupIx[ClusIX[i]+1]; ++z)
						NewOrig[j++] = TmpRIX[z];
				}
				/* printf("totR = %u, Final value of j = %u; NewDedup[totR-1] = %u\n", totR, j, NewDedup[totR-1]);
				printf("NewOrig[j-1] = %u, NewOrig[j] = %u\n",NewOrig[j-1], NewOrig[j]);
				uint8_t *BinUp = calloc(Rd->origTotR,1);
				for (uint32_t i = 0; i < Rd->origTotR-1; ++i) {
					++BinUp[NewOrig[i]];
				}
				for (uint32_t i = 0; i < Rd->origTotR; ++i) {
					if (!BinUp[i]) printf("Lacking representation at %u\n",i);
				}
				printf("last TmpRIX = %u\n",TmpRIX[Rd->origTotR-1]); */
				if (j < Rd->origTotR) puts("Coalescing..."), NewOrig[j] = TmpRIX[Rd->origTotR-1];
				NewDedup[totR] = Rd->origTotR;
				free(Rd->TmpRIX);
				Rd->RefDedupIx = NewDedup, Rd->TmpRIX = NewOrig;
				for (uint32_t i = 0; i < totR; ++i) RefDedupIx[i] = NewOrig[NewDedup[i]];
				Rd->RefIxSrt = RefDedupIx;
			}
			else { // just generate a new RefIxSrt
				uint32_t *TmpRIX = malloc(totR*sizeof(*TmpRIX));
				if (!TmpRIX) {fputs("OOM:TmpRIX2\n",stderr); exit(3);}
				for (uint32_t i = 0; i < totR; ++i) 
					TmpRIX[i] = RefIxSrt[ClusIX[i]];
				free(RefIxSrt);	
				Rd->TmpRIX = TmpRIX, Rd->RefIxSrt = TmpRIX;
			}
			//#pragma omp section
			//{
				// Deal with scattering the FP's and their pointers
				uint32_t totFPs = Fp.nf;
				uint32_t *NewPix = 0;
				Prince *NewP = calloc(totFPs + VECSZ-1,sizeof(*NewP));
				if (!NewP) {fputs("OOM:NewP\n",stderr); exit(3);}
				//if (Z) totFPs = totR, Fp.nf = 0;
				//else {
					NewPix = malloc(totR*sizeof(*NewPix));
					if (!NewPix) {fputs("OOM:NewPix\n",stderr); exit(3);}
					for (uint32_t i = 0; i < totR; ++i)
						NewPix[i] = Fp.Ptrs[ClusIX[i]] >= totR ? Fp.Ptrs[ClusIX[i]] : i;
					for (uint32_t i = totR; i < totFPs; ++i) NewP[i] = p[i];
				//}
				for (uint32_t i = 0; i < totR; ++i) NewP[i] = p[ClusIX[i]];
				free(Fp.initP); 
				free(Fp.Ptrs); 
				p = NewP; 
				Fp.initP = Fp.P = p;
				Fp.Ptrs = NewPix;
				if (Z) {
					Prince t; // New: swap back error-free fingerprints; re-calc centroids
					for (uint32_t i = 0; i < totR; ++i) t = p[i], p[i] = p[Fp.Ptrs[i]], p[Fp.Ptrs[i]] = t;
					for (uint32_t i = 0; i < totRC; ++i) {
						Prince NC = p[i*VECSZ];
						for (uint32_t j = i*VECSZ+1; j < (i+1)*VECSZ; ++j) 
							NC = FP_combine(NC,p[j]);
						UnionPrince[i] = NC;
					}
				}
				
				Rd->FingerprintsR = Fp;
				Rd->Centroids = UnionPrince;
			//}
		}
		
	}
	// get local maxima within clumps
	uint32_t *ClumpLen = malloc((numRclumps + 1) * sizeof(*ClumpLen)); 
	if (!ClumpLen) {fputs("OOM:ClumpLen\n",stderr); exit(3);}
	for (uint32_t i = 0, clen; i < numRclumps; ++i) {
		clen = 0;
		for (int j = 0; j < VECSZ; ++j) 
			if (Rd->RefLen[Rd->RefIxSrt[i*VECSZ+j]] > clen) 
				clen = Rd->RefLen[Rd->RefIxSrt[i*VECSZ+j]];
		ClumpLen[i] = clen;
	}
	ClumpLen[numRclumps] = 0; // endpiece
	for (uint32_t i = numRclumps*VECSZ; i < Rd->totR; ++i) {
		if (Rd->RefLen[Rd->RefIxSrt[i]] > ClumpLen[numRclumps]) 
			ClumpLen[numRclumps] = Rd->RefLen[Rd->RefIxSrt[i]];
	}
	// report the average clump length. 
	uint64_t tot_cl = 0;
	for (uint32_t i = 0; i < totRC; ++i) tot_cl += ClumpLen[i];
	printf("Average R pack length = %f\n", (double)tot_cl/totRC);

	// add references to clumps (TODO: NUMA-aware, first-touch/threadlocal?)
	DualCoil **RefClump = malloc(totRC*sizeof(*RefClump));
	if (!RefClump) {fputs("OOM:RefClump\n",stderr); exit(3);}
	
	//#pragma omp parallel for schedule(dynamic)
	for (uint32_t i = 0; i < numRclumps; ++i) { // put this in main loop (thread-local)
		RefClump[i] = malloc(ClumpLen[i]*sizeof(*RefClump[i])); // align?
		if (!RefClump[i]) {fputs("OOM:RefClump[i]\n",stderr); exit(3);}
		for (uint32_t j = 0; j < ClumpLen[i]; ++j) {
			for (uint32_t k = 0; k < VECSZ; ++k) { // parse letter j across refs
				if (Rd->RefLen[Rd->RefIxSrt[i*VECSZ+k]] >= j) 
					RefClump[i][j].u8[k] = Rd->RefSeq[Rd->RefIxSrt[i*VECSZ+k]][j];
				else RefClump[i][j].u8[k] = BAD_IX;
			}
		}
	}
	// clump endpiece (any stragglers not fitting into a clump)
	if (totRC > numRclumps) {
		RefClump[numRclumps] = malloc(ClumpLen[numRclumps]*sizeof(*RefClump[numRclumps]));
		if (!RefClump[numRclumps]) {fputs("OOM:RefClump[numRclumps]\n",stderr); exit(3);}
		for (uint32_t j = 0; j < ClumpLen[numRclumps]; ++j) {
			for (uint32_t k = 0; k < Rd->totR-numRclumps*VECSZ; ++k) { // parse letter j across refs
				if (Rd->RefLen[Rd->RefIxSrt[numRclumps*VECSZ+k]] >= j) 
					RefClump[numRclumps][j].u8[k] = Rd->RefSeq[Rd->RefIxSrt[numRclumps*VECSZ+k]][j];
				else RefClump[numRclumps][j].u8[k] = BAD_IX;
			}
			for (uint32_t k = Rd->totR-numRclumps*VECSZ; k < VECSZ; ++k) {
				RefClump[numRclumps][j].u8[k] = BAD_IX;
			}
		}
	}

	Rd->numRclumps = numRclumps = totRC;
	Rd->ClumpLen = ClumpLen;
	Rd->RefClump = RefClump;

	//If in SSE2 mode, prepare profiles (TODO: move one profile generation per loop per thread)
	#ifndef __SSSE3__
	if (curate < 2) create_sse2_profiles(Rd);
	#endif
	if (DO_ACCEL) return;
	if (!DO_FP) {
		if (Rd->dbType == QUICK) 
			for (uint32_t i = 0; i < origR; ++i) free(origRefSeq[i]);
		else free(SeqDump);
	}
	if (REBASE) free(origRefLen), free(origRefSeq), free(origRefHead);
}

static inline void dump_edb(FILE *output, Reference_Data *Rd) {
	fputc(0, output); // will be padded later
	uint64_t totRefHeadLen = 0; uint32_t shear = DB_QLEN / THRES;
	fwrite(&totRefHeadLen,sizeof(totRefHeadLen),1,output);
	fwrite(&shear,sizeof(uint32_t),1,output); // V3
	fwrite(&Rd->totR,sizeof(Rd->totR),1,output); 
	fwrite(&Rd->origTotR,sizeof(Rd->origTotR),1,output);
	fwrite(&Rd->numRclumps,sizeof(Rd->numRclumps),1,output);
	fwrite(&Rd->maxLenR,sizeof(Rd->maxLenR),1,output);

	// New: sort and dedupe
	StrIxPair *Str_Ix = malloc(Rd->origTotR*sizeof(*Str_Ix));
	for (uint32_t i = 0; i < Rd->origTotR; ++i) 
		Str_Ix[i] = (StrIxPair){Rd->RefHead[i],i};
	qsort(Str_Ix,Rd->origTotR,sizeof(*Str_Ix),cmpStrIx);
	uint32_t *RefMap = malloc(Rd->origTotR*sizeof(*RefMap));
	uint32_t nix = 0;
	char *cur = Str_Ix[nix].s;
	totRefHeadLen += 1 + fprintf(output, "%s", cur), fputc(0,output);
	for (uint32_t i = 0; i < Rd->origTotR; ++i) {
		if (strcmp(cur,Str_Ix[i].s)) {
			cur = Str_Ix[i].s;
			totRefHeadLen += 1 + fprintf(output, "%s", cur), fputc(0,output);
			++nix;
		}
		RefMap[Str_Ix[i].ix] = nix;	
	}
	nix++;
	fwrite(&nix, sizeof(nix), 1, output); // read in as numRefHeads
	fwrite(RefMap, sizeof(*RefMap),Rd->origTotR, output);

	if (REBASE) 
		fwrite(Rd->RefStart, sizeof(*Rd->RefStart), Rd->origTotR, output);
	if (Rd->totR != Rd->origTotR) //puts("Dupewriting"), 
		fwrite(Rd->RefDedupIx, sizeof(*Rd->RefDedupIx), Rd->totR+1, output);

	// Write (sheared/deduped) clumps and mask into originals
	fwrite(Rd->TmpRIX, sizeof(*Rd->TmpRIX), Rd->origTotR, output); // can infer RefIxSrt
	fwrite(Rd->ClumpLen, sizeof(*Rd->ClumpLen), Rd->numRclumps, output);
	
	// New
	uint64_t oldTotRPackLen = 0, totRPackLen = 0;
	uint32_t *ClumpLen = Rd->ClumpLen; //, maxLenR = 0;
	#pragma omp simd reduction(+:oldTotRPackLen,totRPackLen)
	for (uint32_t i = 0; i < Rd->numRclumps; ++i) {
		oldTotRPackLen += ClumpLen[i];
		//if (ClumpLen[i] > maxLenR) maxLenR = ClumpLen[i];
		totRPackLen += ClumpLen[i]/2u + (ClumpLen[i] & 1);
	}
	printf(" --> neoEDB: Original pack table: %llu [now %llu]\n",
			oldTotRPackLen, totRPackLen);
	
	uint64_t pix = 0, totWritten = 0;
	for (uint32_t i = 0; i < Rd->numRclumps; ++i) {
		DualCoil *rc = Rd->RefClump[i];
		for (uint32_t j = 0; j < ClumpLen[i]; j+=2) {
			__m128i i1 = _mm_load_si128((void*)(rc+j)); //(DC_Dump + pix++));
			__m128i i2 = _mm_setzero_si128();
			++pix;
			if (j + 1 < ClumpLen[i]) i2 = _mm_load_si128((void*)(rc+j+1)), ++pix; //(DC_Dump + pix++));
			i2 = _mm_and_si128(_mm_slli_epi16(i2,4),_mm_set1_epi8(0xF0));
			DualCoil x;
			_mm_store_si128((void*)&x, _mm_or_si128(i1, i2));
			fwrite(&x, sizeof(x), 1, output);
			++totWritten;
		}
	}
	printf(" --> neoEDB: Number read: %llu, written: %llu\n", pix, totWritten);
	//exit(101);
	if (DO_FP) { // write centroid and nf prints
		uint32_t nf = Rd->FingerprintsR.nf;
		fwrite(Rd->Centroids, sizeof(*Rd->Centroids), Rd->numRclumps, output);
		fwrite(&nf, sizeof(nf), 1, output);
		if (nf) fwrite(Rd->FingerprintsR.Ptrs, sizeof(*Rd->FingerprintsR.Ptrs), Rd->totR, output);
		else nf = Rd->totR; // if nf == 0, it means the FPs are unambiguous (N-penalized)
		fwrite(Rd->FingerprintsR.P, sizeof(*Rd->FingerprintsR.P), nf, output);
	}
	rewind(output);
	fputc(1 << 7 | REBASE << 6 | DO_FP << 5 | Xalpha << 4 | EDX_VERSION, output);
		//(EDB_VERSION - Xalpha), output); // 5: xalpha, 6: ?
	fwrite(&totRefHeadLen, sizeof(totRefHeadLen), 1, output);
}

static inline uint32_t read_edb(char *ref_FN, Reference_Data *Rd) {
	FILE *in = fopen(ref_FN, "rb"), *output = 0;
	if (!in) {fputs("ERROR: cannot parse EDB",stderr); exit(1);}
	uint8_t cb = fgetc(in), dbVer = cb & 0xF, doConv2_3 = 0, newVer = 0;
	if (dbVer != EDX_VERSION && dbVer != EDB_VERSION) {
		fprintf(stderr,"ERROR: invalid database version %u\n", dbVer);
		exit(1);
	}
	else if (dbVer == EDB_VERSION) 
		fprintf(stderr,"ERROR: Old DB version. Converting...\n"), 
		doConv2_3 = 1;
	REBASE = (uint8_t)(cb << 1) >> 7;
	DO_FP = (uint8_t)(cb << 2) >> 7 && DO_FP;
	if (!DO_FP) printf(" --> EDB: Fingerprints are DISABLED\n");
	int dbX = (uint8_t)(cb << 3) >> 7;
	if (dbX && !Xalpha || !dbX && Xalpha) {
		fprintf(stderr,"ERROR: DB made with%s Xalpha; queries %s use Xalpha.\n",
			dbX? "" : "out", Xalpha? "can't" : "must");
		exit(1);
	}
	uint64_t totRefHeadLen = 0;
	uint32_t totR = 0, origTotR = 0, numRclumps = 0, maxLenR = 0, shear = 0;
	fread(&totRefHeadLen, sizeof(totRefHeadLen), 1, in);
	if (dbVer > 1) fread(&shear, sizeof(shear), 1, in); // V2+
	fread(&totR, sizeof(totR), 1, in); 
	fread(&origTotR, sizeof(origTotR), 1, in);
	fread(&numRclumps, sizeof(numRclumps), 1, in);
	fread(&maxLenR, sizeof(maxLenR), 1, in);

	char *RH_dump = malloc(totRefHeadLen), // *RH_orig = RH_dump,
		**RefHead = malloc(origTotR * sizeof(*RefHead));
	uint32_t *RefMap = 0, *RefStart = 0, *RefDedupIx = 0, 
		*TmpRIX = malloc(origTotR * sizeof(*TmpRIX)), 
		*ClumpLen = malloc(numRclumps * sizeof(*ClumpLen));
	DualCoil **RefClump = malloc(numRclumps * sizeof(*RefClump));
	if (!RH_dump || !RefHead || !TmpRIX || !ClumpLen || !RefClump)
		{fputs("OOM:read_edb\n",stderr); exit(3);}
	fread(RH_dump, 1, totRefHeadLen, in);
	*RefHead = RH_dump;
	uint32_t numRefHeads = origTotR; 
	if (dbVer == EDX_VERSION) fread(&numRefHeads,sizeof(numRefHeads),1,in);

	for (uint32_t i = 1; i <numRefHeads; ++i) {
		while (*RH_dump++);
		RefHead[i] = RH_dump;
	}
	if (dbVer == EDX_VERSION) { // new unique header scheme
		puts(" --> EDB: Parsing compressed headers");
		RefMap = malloc(origTotR*sizeof(*RefMap));
		if (!RefMap) {fputs("OOM:read_edb\n",stderr); exit(3);}
		fread(RefMap, origTotR, sizeof(*RefMap), in);
		char **NewRefHead = malloc(origTotR*sizeof(*NewRefHead)), *t;
		for (uint32_t i = 0; i < origTotR; ++i) 
			NewRefHead[i] = RefHead[RefMap[i]];
		free(RefHead);
		//Rd->DedupRefHead = RefHead;
		Rd->RefMap = RefMap; // for de-duping alignments later
		RefHead = NewRefHead;
	}
	if (REBASE) {
		printf(" --> EDB: Sheared database ");
		if (dbVer > 1) printf("(shear size = %u)\n",shear);
		else printf("(WARNING: old DB with unknown shear size)\n");
		RefStart = malloc(origTotR * sizeof(*RefStart));
		if (!RefStart) {fputs("OOM:read_edb\n",stderr); exit(3);}
		fread(RefStart, sizeof(*RefStart), origTotR, in);
	}
	if (totR != origTotR) {
		puts(" --> EDB: Unique-reference database");
		RefDedupIx = malloc((totR+1)*sizeof(*RefDedupIx));
		if (!RefDedupIx) {fputs("OOM:read_edb\n",stderr); exit(3);}
		fread(RefDedupIx, sizeof(*RefDedupIx), totR + 1, in);
	}
	fread(TmpRIX, sizeof(*TmpRIX), origTotR, in);
	fread(ClumpLen, sizeof(*ClumpLen), numRclumps, in);

	uint64_t totRPackLen = 0;
	#pragma omp simd reduction(max:maxLenR) reduction(+:totRPackLen)
	for (uint32_t i = 0; i < numRclumps; ++i) {
		if (ClumpLen[i] > maxLenR) maxLenR = ClumpLen[i];
		totRPackLen += ClumpLen[i]/2u + (ClumpLen[i] & 1);
	}
	if (doConv2_3) {
		int slen = strlen(ref_FN);
		if (slen >= 4096) {puts("Filename error."); exit(1);}
		char newName[4096] = {0};
		strcpy(newName,ref_FN);
		if (!strcmp(newName+slen-4,".edb")) newName[slen-1] = 'x';
		else strcpy(newName+slen,".edx");
		output = fopen(newName,"wb");
		if (!output) {fprintf(stderr,"ERROR: Cannot write '%s'\n",newName); exit(1);}
		fprintf(stderr," --> ConvEDB: Re-writing database as '%s'\n",newName);
		newVer = 1 << 7 | REBASE << 6 | ((uint8_t)(cb << 2) >> 7) << 5 | Xalpha << 4 | EDX_VERSION;
		fputc(newVer, output);
		totRefHeadLen = 0;
		fwrite(&totRefHeadLen, sizeof(totRefHeadLen), 1, output);
		fwrite(&shear, sizeof(shear), 1, output);
		fwrite(&totR, sizeof(totR), 1, output);
		fwrite(&origTotR, sizeof(origTotR), 1, output);
		fwrite(&numRclumps, sizeof(numRclumps), 1, output);
		fwrite(&maxLenR, sizeof(maxLenR), 1, output);
		//fwrite(*RefHead, 1, totRefHeadLen, output);
		// Sort and dedupe.
		StrIxPair *Str_Ix = malloc(origTotR*sizeof(*Str_Ix));
		for (uint32_t i = 0; i < origTotR; ++i) 
			Str_Ix[i] = (StrIxPair){RefHead[i],i};
		qsort(Str_Ix,origTotR,sizeof(*Str_Ix),cmpStrIx);
		uint32_t *RefMap = malloc(origTotR*sizeof(*RefMap));
		uint32_t nix = 0;
		char *cur = Str_Ix[nix].s;
		totRefHeadLen += 1 + fprintf(output, "%s", cur), fputc(0,output);
		for (uint32_t i = 0; i < origTotR; ++i) {
			if (strcmp(cur,Str_Ix[i].s)) {
				cur = Str_Ix[i].s;
				totRefHeadLen += 1 + fprintf(output, "%s", cur), fputc(0,output);
				++nix;
			}
			RefMap[Str_Ix[i].ix] = nix;	
		}
		nix++;
		fwrite(&nix, sizeof(nix), 1, output); // read in as numRefHeads
		fwrite(RefMap, sizeof(*RefMap),origTotR, output);

		if (REBASE) fwrite(RefStart, sizeof(*RefStart), origTotR, output);
		if (totR != origTotR) fwrite(RefDedupIx, sizeof(*RefDedupIx), totR + 1, output);
		fwrite(TmpRIX, sizeof(*TmpRIX), origTotR, output);
		fwrite(ClumpLen, sizeof(*ClumpLen), numRclumps, output);

		uint64_t oldTotRPackLen = 0;
		#pragma omp simd reduction(+:oldTotRPackLen)
		for (uint32_t i = 0; i < numRclumps; ++i) 
			oldTotRPackLen += ClumpLen[i];
		printf(" --> ConvEDB: Original pack table: %llu [now %llu]\n",
			oldTotRPackLen, totRPackLen);
		DualCoil *DC_Dump = malloc((oldTotRPackLen+1)*sizeof(*DC_Dump));
		if (!DC_Dump) {fputs("OOM:read_edb\n",stderr); exit(3);}
		fread(DC_Dump, sizeof(*DC_Dump), oldTotRPackLen, in);
		uint64_t pix = 0, totWritten = 0;
		for (uint32_t i = 0; i < numRclumps; ++i) {
			for (uint32_t j = 0; j < ClumpLen[i]; j+=2) {
				__m128i i1 = _mm_load_si128((void*)(DC_Dump + pix++));
				__m128i i2 = _mm_setzero_si128();
				if (j + 1 < ClumpLen[i]) i2 = _mm_load_si128((void*)(DC_Dump + pix++));
				i2 = _mm_and_si128(_mm_slli_epi16(i2,4),_mm_set1_epi8(0xF0));
				DualCoil x;
				_mm_store_si128((void*)&x, _mm_or_si128(i1, i2));
				fwrite(&x, sizeof(x), 1, output);
				++totWritten;
			}
		}
		printf(" --> ConvEDB: Number read: %llu, written: %llu\n", pix, totWritten);
	}
	else { // read new format edb (TODO FIX)
		DualCoil *DC_Dump = malloc((totRPackLen+1)*sizeof(*DC_Dump));
		if (!DC_Dump) {fputs("OOM:read_edb\n",stderr); exit(3);}
		fread(DC_Dump, sizeof(*DC_Dump), totRPackLen, in);
		size_t accum = 0;
		//#pragma omp simd reduction(+:accum)
		*RefClump = DC_Dump;
		for (uint32_t i = 1; i < numRclumps; ++i) 
			RefClump[i] = RefClump[i-1] + ClumpLen[i-1]/2u + (ClumpLen[i-1] & 1);
	}
	Prince *Centroids = 0; 
	void *init_FP = 0; Prince *FullPrince = 0;
	uint32_t R=0, numFPs=0, unAmbig = 0, *FP_ptrs = 0;
	if (DO_FP || (doConv2_3 && (uint8_t)(cb << 2) >> 7)) {
		printf(" --> %sEDB: Fingerprinting is enabled ", doConv2_3 ? "Conv" : "");
		Centroids = malloc(numRclumps*sizeof(*Centroids));
		if (!Centroids) {fputs("OOM:read_edb\n",stderr); exit(3);}
		fread(Centroids, sizeof(*Centroids), numRclumps, in);
		if (dbVer) { // DB version > 0 [check individual bit?]
			fread(&numFPs, sizeof(numFPs), 1, in);
			if (!numFPs) numFPs = totR, R = 1, printf("[RESTRICTED DB] ");
			else {
				FP_ptrs = malloc(totR*sizeof(*FP_ptrs));
				if (!FP_ptrs) {fputs("OOM:read_edb\n",stderr); exit(3);}
				fread(FP_ptrs, sizeof(*FP_ptrs), totR, in);
			}
			FullPrince = calloc_a(64,(numFPs+VECSZ-1)*sizeof(*FullPrince),&init_FP);
			if (!init_FP) {fputs("OOM:read_edb\n",stderr); exit(3);}
			printf("(%u FPs)\n",numFPs);
			fread(FullPrince, sizeof(*FullPrince), numFPs, in);
		}
		else printf("(WARNING: old DB without individual FPs)\n");
		if (doConv2_3) {
			fwrite(Centroids, sizeof(*Centroids), numRclumps, output);
			fwrite(&numFPs, sizeof(numFPs), 1, output);
			if (numFPs) fwrite(FP_ptrs, sizeof(*FP_ptrs), totR, output);
			fwrite(FullPrince, sizeof(*FullPrince), numFPs, output);
			fprintf(stderr," --> ConvEDB: complete. Re-run with the new .edx\n");
		}
	}
	if (doConv2_3) {
		rewind(output);
		fputc(newVer,output);
		fwrite(&totRefHeadLen, sizeof(totRefHeadLen), 1, output);
		exit(0);
	}
	fclose(in);

	Rd->RefHead = RefHead, Rd->ClumpLen = ClumpLen, Rd->TmpRIX = TmpRIX,
	Rd->RefDedupIx = RefDedupIx;
	Rd->totR = totR, Rd->origTotR = origTotR, Rd->numRclumps = numRclumps;
	Rd->RefClump = RefClump, Rd->maxLenR = maxLenR;
	Rd->RefStart = RefStart;
	
	Rd->numRefHeads = numRefHeads;
	Rd->Centroids = Centroids;
	Rd->FingerprintsR.P = FullPrince;
	Rd->FingerprintsR.Ptrs = FP_ptrs;
	Rd->FingerprintsR.nf = R ? 0 : numFPs; 
	Rd->FingerprintsR.initP = init_FP;
	printf(" --> EDB: %u refs [%u orig], %u clumps, %u maxR\n",
		totR, origTotR, numRclumps, maxLenR);
	return shear;
}

static inline void process_queries(char *query_FN, Query_Data *Qd) {
	char **QHead, **QSeq; //size_t
	uint32_t *QLen, *Divergence;
	uint32_t totQ, maxLenQ = 0, minLenQ = UINT32_MAX;
	double wt = omp_get_wtime();
	totQ = parse_tl_faster(query_FN, &QHead, &QSeq, &QLen);
	if (!totQ) {fputs("ERROR: No queries found.",stderr); exit(1);}
	if (!Qd->incl_whitespace) {
		#pragma omp parallel for
		for (uint32_t i = 0; i < totQ; ++i) {
			char *q = QHead[i];
			while (*q && *q != ' ' && *q != '\t') ++q; *q = 0;
		}
	}
	printf("Parsed %u queries (%f). Calculating minMax... \n",totQ,omp_get_wtime()-wt);
	wt = omp_get_wtime();
	#pragma omp parallel for reduction(max:maxLenQ) reduction(min:minLenQ)
	for (uint32_t i=0; i < totQ; ++i) {
		if (QLen[i] > maxLenQ) maxLenQ = QLen[i];
		if (QLen[i] < minLenQ) minLenQ = QLen[i]; 
	}
	if (maxLenQ > (1 << 16)) fputs("WARNING: Max query length is very long\n",stderr);
	if (minLenQ < 5) fputs("WARNING: Min query length is less than 5 bases\n",stderr);
	printf("Found min %u, max %u (%f).\n",minLenQ, maxLenQ, omp_get_wtime()-wt);

	if (!Xalpha) {
		printf("Converting queries... ");
		wt = omp_get_wtime();
		#pragma omp parallel for
		for (uint32_t i = 0; i < totQ; ++i) translate16aln(QSeq[i],QLen[i]);
		printf("Converted (%f)\n", omp_get_wtime()-wt);
	}
	printf("Copying queries... ");
	wt = omp_get_wtime();
	StrIxPair *Query_Ix = malloc(totQ * sizeof(*Query_Ix));
	if (!Query_Ix) {fputs("OOM:Query_Ix\n",stderr); exit(3);}
	#pragma omp parallel for
	for (uint32_t i = 0; i < totQ; ++i) Query_Ix[i]=(StrIxPair){QSeq[i],i};
	printf("Copied (%f)\n",omp_get_wtime() - wt);
	printf("Sorting queries... ");
	wt = omp_get_wtime();
	parallel_sort_strpack(Query_Ix, totQ);
	printf("Sorted (%f)\n",omp_get_wtime()-wt);

	printf("Copying indices... ");
	wt = omp_get_wtime();
	uint32_t *NewIX = malloc(totQ*sizeof(*NewIX));
	if (!NewIX) {fputs("OOM:NewIX\n",stderr); exit(3);}
	#pragma omp parallel for 
	for (uint32_t i = 0; i < totQ; ++i) NewIX[i] = Query_Ix[i].ix;
	free(Query_Ix), Query_Ix = 0;
	printf("Copied (%f)\n",omp_get_wtime() - wt);

	// Now determine uniqueness in parallel
	printf("Determining uniqueness... ");
	wt = omp_get_wtime();
	uint8_t *IsUniq = calloc(totQ,1);
	if (!IsUniq) {fputs("OOM:IsUniq\n",stderr); exit(3);}
	IsUniq[0] = 1;
	uint32_t numUniqQ = 1, numUniqQ_rc = 0;
	#pragma omp parallel for reduction(+:numUniqQ)
	for (uint32_t i = 1; i < totQ; ++i) 
		if (strcmp(QSeq[NewIX[i-1]],QSeq[NewIX[i]]))
			IsUniq[i] = 1, ++numUniqQ;
	printf("Done (%f). Number unique: %u\n",omp_get_wtime() - wt, numUniqQ);

	printf("Collecting unique sequences... ");
	wt = omp_get_wtime();
	uint32_t *Offset = malloc((numUniqQ+1)*sizeof(*Offset));
	if (!Offset) {fputs("OOM:Offset\n",stderr); exit(3);}
	uint32_t uix = 0;
	for (uint32_t i = 0; i < totQ; ++i) 
		if (IsUniq[i]) Offset[uix++] = i;
	Offset[numUniqQ] = totQ;

	char **SrtQHead = malloc(totQ*sizeof(*SrtQHead));
	if (!SrtQHead) {fputs("OOM:SrtQHead\n",stderr); exit(3);}
	#pragma omp parallel for
	for (uint32_t i = 0; i < totQ; ++i) 
		SrtQHead[i] = QHead[NewIX[i]];
	free(QHead); QHead = SrtQHead;
	printf("Done (%f)\n",omp_get_wtime() - wt);

	printf("Creating data structures... ");
	wt = omp_get_wtime();
	numUniqQ_rc = numUniqQ + ((Qd->rc && !DO_PREPASS) ? numUniqQ : 0);
	UniBin *UniBins = calloc(numUniqQ_rc,sizeof(*UniBins));
	ShrBin *ShrBins = calloc(numUniqQ,sizeof(*ShrBins));
	if (!UniBins || !ShrBins) {fputs("OOM:Uni/ShrBins\n",stderr); exit(3);}
	float reqID = 1/THRES - 1;
	uint64_t uniqTotLen = 0; uint32_t maxED = 0;
	#pragma omp parallel for reduction(+:uniqTotLen) reduction(max:maxED)
	for (uint32_t i = 0; i < numUniqQ; ++i) {
		uint32_t rix = Offset[i];
		uint32_t len = QLen[NewIX[rix]], ed = reqID * len;
		ShrBins[i].len = len;
		ShrBins[i].ed = MIN(254,ed);
		UniBins[i].Seq = QSeq[NewIX[rix]];
		UniBins[i].six = i;
		uniqTotLen += len;
		maxED = ShrBins[i].ed > maxED ? ShrBins[i].ed : maxED;
	}
	free(IsUniq); free(NewIX); free(QSeq); free(QLen);
	IsUniq = 0, NewIX = 0, QSeq = 0, QLen = 0;
	printf("Done (%f) [maxED: %u]\n",omp_get_wtime() - wt, maxED);
	
	char *RCDump = 0;
	if (Qd->rc && !DO_PREPASS) {
		if (numUniqQ >= UINT32_MAX/2) {fprintf(stderr,"ERROR: too many queries for RC\n"); exit(1);}
		printf("Processing reverse complements... ");
		wt = omp_get_wtime();
		RCDump = malloc(uniqTotLen + numUniqQ + 1);
		if (!RCDump) {fputs("OOM:RCDump\n",stderr); exit(3);}
		char *ptr = RCDump;
		#pragma omp parallel for
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			uint32_t len = ShrBins[i].len, nix = numUniqQ + i;
			char *rcs, *org = UniBins[i].Seq;
			#pragma omp atomic capture
			{rcs = ptr; ptr += len + 1;}
			for (uint32_t j = 0; j < len; ++j)
				rcs[j] = RVT[org[len-j-1]];
			rcs[len] = 0;
			UniBins[nix].Seq = rcs;
			UniBins[nix].rc = 1;
			UniBins[nix].six = i;
		}
		printf("Processed (%f) [offset %llu / %llu]\n",omp_get_wtime() - wt,ptr-RCDump,uniqTotLen+numUniqQ);
	}

	uint32_t *QBins = calloc(5,sizeof(*QBins));
	uint32_t  nq_clear, nq_ambig, nq_bad;
	if (!DO_PREPASS && DO_ACCEL) {
		printf("Determining query ambiguity... ");
		wt = omp_get_wtime();
		uint32_t newUniqQ = numUniqQ_rc;
		// re-sort queries: ambig (0), unambig (1), bad (2) [tooManyErrors = 3?]
		uint8_t *QStat = malloc((1+newUniqQ)*sizeof(*QStat));
		//Umap = malloc(newUniqQ*sizeof(*Umap));
		if (!QStat /*|| !Umap*/) {fputs("OOM:QStat\n",stderr); exit(3);}
		memset(QStat,1,newUniqQ*sizeof(*QStat));
		QStat[newUniqQ] = 0; // endcap JIC
		#pragma omp parallel for
		for (uint32_t i = 0; i < newUniqQ; ++i) {
			//uint32_t ai = i >= numUniqQ ? RCMap[i-numUniqQ] : i;
			UniBin Ub = UniBins[i]; ShrBin Sb = ShrBins[Ub.six];

			uint32_t rowN = 0, totN = 0, len = Sb.len;
			char *s = Ub.Seq; 
			if (len < SCOUR_N || (!DO_HEUR && Sb.ed >= len/SCOUR_N)) 
				QStat[i] = 2;
		
			else for (uint32_t j = 0; j < len; ++j) {
				if ((totN += s[j] > 4 + Z) > 5) {QStat[i] = 2; break;} // New! +Z 'cuz N not ambig'
				else if (s[j] > 4) QStat[i] = 0; // but they still demote to ambig bin
				//if (s[j] > 4) {
				//		QStat[i] = 0; // set as ambig-containing
				//	if (++totN > 5 || ++rowN > 3) {QStat[i] = 2; break;}
				//} else rowN = 0;
			}
		}
		printf("Determined (%f)\n",omp_get_wtime() - wt);
		printf("Creating bins... ");
		wt = omp_get_wtime();
		{
			uint32_t *P = QBins + 1;
			#pragma omp parallel for
			for (uint32_t i = 0; i < newUniqQ; ++i) 
				#pragma omp atomic update
				++P[QStat[i]];
			nq_ambig = P[0], nq_clear = P[1], nq_bad = P[2];
			--P;
			for (uint32_t i = 1; i < 5; ++i) P[i] += P[i-1];
			
			uint32_t QIX[4] = {0}; // = {P[1], P[2], P[3], P[4], P[4]};
			for (uint32_t i = 0; i < 3; ++i) QIX[i] = P[i+1]; 
			QIX[3] = -1; // the last [faux] bin is never incrementable!
			//printf("%u -> %u -> %u -> %u\n",P[0],P[1],P[2],P[3]);
			//printf("QIX: %u -> %u -> %u -> %u\n",QIX[0],QIX[1],QIX[2],QIX[3]);
			uint8_t cur_bin = 0; //uint64_t totScan = 0;
			for (uint32_t i = 0; i < newUniqQ; ++i) {
				UniBin tb; uint8_t ts;
				//printf("[%u]: %u",i,QStat[i]);
				while (i >= QIX[cur_bin]) i = P[++cur_bin]; //, printf("->%u(b%u)",i,cur_bin); 
				while (QStat[i] > cur_bin) { // swap local and remote
					tb = UniBins[i], ts = QStat[i];
					//printf(" <%u->[%u]:%u>", ts, P[ts], QStat[P[ts]]);
					UniBins[i] = UniBins[P[ts]], UniBins[P[ts]] = tb;
					QStat[i] = QStat[P[ts]], QStat[P[ts]] = ts;
					++P[ts]; // remote increment
					//++totScan;
				}
				++P[cur_bin];
				//++totScan;
				//printf("\n");
			}
			// DEBUG: double-check that we are indeed in ascending order (at least in the QStat bins)
			//#pragma omp parallel for
			//for (uint32_t i = 1; i < newUniqQ; ++i) 
			//	if (QStat[i] < QStat[i-1]) {puts("ERROR IN QSTAT"); exit(10);}
			free(QStat);
			//printf("Tot: %llu, P[0,1,2,3] = %u,%u,%u,%u\n",totScan, P[0],P[1],P[2],P[3]);
			
		}
		/* Split *Ix_Code = malloc(newUniqQ*sizeof(*Ix_Code));
		for (uint32_t i = 0; i < newUniqQ; ++i) Ix_Code[i] = (Split){QStat[i],i};
		int SpltCmp(const void *a, const void *b) {
			Split *A = (Split *)a, *B = (Split *)b;
			return A->i < B->i ? -1 : B->i < A->i;
		}
		qsort(Ix_Code, newUniqQ, sizeof(*Ix_Code), SpltCmp);
		UniBin *NUniBins = malloc(newUniqQ*sizeof(*NUniBins));
		for (uint32_t i = 0; i < newUniqQ; ++i) NUniBins[i] = UniBins[Ix_Code[i].i];
		free(UniBins); free(Ix_Code);
		UniBins = NUniBins; */
		//++QBins; // TEST ONLY
		
		printf("Created (%f); Unambig: %u, ambig: %u, super-ambig: %u [%u,%u,%u]\n",omp_get_wtime() - wt,
			nq_clear, nq_ambig, nq_bad, QBins[0], QBins[1], QBins[2]);
	}
	if (!DO_PREPASS && (Qd->rc || (DO_ACCEL && nq_clear != numUniqQ_rc))) {
		printf("Re-sorting... ");
		wt = omp_get_wtime();
		if (!DO_ACCEL) parallel_sort_unibin(UniBins,numUniqQ_rc);
		else parallel_sort_unibin(UniBins,QBins[0]),
			parallel_sort_unibin(UniBins+QBins[0],QBins[1]-QBins[0]),
			parallel_sort_unibin(UniBins+QBins[1],QBins[2]-QBins[1]);
		printf("Re-sorted (%f)\n",omp_get_wtime() - wt);
	}
	printf("Calculating divergence... ");
	wt = omp_get_wtime();
	uint64_t tot_div = 1;
	uint32_t maxDiv = 0;
	UniBins[0].div = 1;
	if (THRES <= 0.75f) cacheSz /= 2;  // limit thrashing
	cacheSz = MIN(maxLenQ+2, cacheSz); // limit space
	#pragma omp parallel for reduction(+:tot_div) reduction(max:maxDiv)
	for (uint32_t i = 1; i < numUniqQ_rc; ++i) {
		uint32_t div = 1;
		char *prv = UniBins[i-1].Seq, *cur = UniBins[i].Seq;
		while (*prv && *prv++ == *cur++) ++div;
		div = MIN(cacheSz,div);
		div = MIN(ShrBins[UniBins[i-1].six].len,div);
		UniBins[i].div = div;
		tot_div += div;
		maxDiv = div > maxDiv ? div : maxDiv;
	}
	if (QBins[1] > QBins[0]) UniBins[QBins[0]].div = 1;
	if (QBins[2] > QBins[1]) UniBins[QBins[1]].div = 1;
	printf("Calculated (%f) [%f avg div; %u max]\n", omp_get_wtime() - wt,
		(double)tot_div/numUniqQ_rc, maxDiv);
	
	// DEBUG ONLY
	/* for (uint32_t i = 0; i < 30; ++i) {
		printf("%u: %s:",i,QHead[Offset[UniBins[i].six]]);
		for (uint32_t j = 0; j < ShrBins[UniBins[i].six].len; ++j)
			printf("%c",BAK[UniBins[i].Seq[j]]);
		printf("\n");
	} */
	
	Qd->QHead = QHead; Qd->SeqDumpRC = RCDump; Qd->Offset = Offset;
	Qd->totQ = totQ, Qd->numUniqQ = numUniqQ, Qd->maxDiv = maxDiv; 
	Qd->UniBins = UniBins; Qd->ShrBins = ShrBins; Qd->maxED = maxED;
	Qd->QBins = QBins; Qd->maxLenQ = maxLenQ; Qd->minLenQ = minLenQ;
	//exit(1);
}

//#define SCOUR_NB 3
//uint32_t getAmbigScour(char *S, uint32_t w, uint8_t ix) {
static void setAmbigScour20(uint16_t *Hash, Split *Refs, uint32_t *nref, uint32_t mmatch,
 void **Forest, char *S, uint32_t qdim, uint32_t t, uint8_t ix, uint32_t *Cache, uint32_t *cix) {
	if (ix == SCOUR_N) {
		for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=5) {
			uint64_t PX = *(uint64_t *)PXp;
			uint32_t px1 = PX & 0xFFFFF, px2 = (PX >> 20) & 0xFFFFF;
			if (!Hash[px1]) Cache[(*cix)++] = px1;
			if (Hash[px1] == mmatch) {
				Hash[px1] = qdim + mmatch + 1;
				Refs[(*nref)++].v = px1;
			} 
			else if (Hash[px1] < UINT16_MAX) ++Hash[px1];
			if (PXp + 3 >= Forest[t+1]) break;
			if (!Hash[px2]) Cache[(*cix)++] = px2;
			if (Hash[px2] == mmatch) {
				Hash[px2] = qdim + mmatch + 1;
				Refs[(*nref)++].v = px2;
			} 
			else if (Hash[px2] < UINT16_MAX) ++Hash[px2];
		}
	}
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		setAmbigScour20(Hash, Refs, nref, mmatch, Forest, S, qdim,
			t << 2 | AMBIGS[S[ix]][i], ix + 1, Cache, cix);
}

static void setAmbigScour24(uint16_t *Hash, Split *Refs, uint32_t *nref, uint32_t mmatch,
 void **Forest, char *S, uint32_t qdim, uint32_t w, uint8_t ix, uint32_t *Cache, uint32_t *cix) {
	if (ix == SCOUR_N) {
		for (void *PXp = Forest[w]; PXp < Forest[w+1]; PXp+=3) {
			uint32_t PX = *(uint32_t *)PXp & 0xFFFFFF;
			if (!Hash[PX]) Cache[(*cix)++] = PX;
			if (Hash[PX] == mmatch) 
				Hash[PX] = qdim + mmatch,
				Refs[(*nref)++].v = PX;
			else if (Hash[PX] < UINT16_MAX) ++Hash[PX];
		}
	}
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		setAmbigScour24(Hash, Refs, nref, mmatch, Forest, S, qdim,
			w << 2 | AMBIGS[S[ix]][i], ix + 1, Cache, cix);
}

static void setUnambigScour20(uint16_t *Hash, Split *Refs, uint32_t *nref, uint32_t mmatch,
 void **Forest, char *s, uint32_t qdim, uint32_t len, uint32_t *Cache, uint32_t *cix) {
  	uint32_t w = 0; // cache k-indices
	for (uint32_t k = 0; k + 1 < SCOUR_N; ++k)
		w <<= 2, w |= s[k]-1;
	for (uint32_t k = SCOUR_N-1; k < len; ++k) {
		w <<= 2, w |= s[k]-1;
		uint32_t t = (w << SCOUR_R) >> SCOUR_R;
		for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=5) {
			uint64_t PX = *(uint64_t *)PXp;// & 0xFFFFFF;
			uint32_t px1 = PX & 0xFFFFF, px2 = (PX >> 20) & 0xFFFFF;
			if (!Hash[px1]) Cache[(*cix)++] = px1;
			if (Hash[px1] == mmatch) {
				Hash[px1] = qdim + mmatch + 1;
				Refs[(*nref)++].v = px1;
			} 
			else if (Hash[px1] < UINT16_MAX) ++Hash[px1];
			if (PXp + 3 >= Forest[t+1]) break;
			if (!Hash[px2]) Cache[(*cix)++] = px2;
			if (Hash[px2] == mmatch) {
				Hash[px2] = qdim + mmatch + 1;
				Refs[(*nref)++].v = px2;
			} 
			else if (Hash[px2] < UINT16_MAX) ++Hash[px2];
		}
	}
}
static void setUnambigScour24(uint16_t *Hash, Split *Refs, uint32_t *nref, uint32_t mmatch,
 void **Forest, char *s, uint32_t qdim, uint32_t len, uint32_t *Cache, uint32_t *cix) {
	  uint32_t w = 0; // cache k-indices
	for (uint32_t k = 0; k + 1 < SCOUR_N; ++k)
		w <<= 2, w |= s[k]-1;
	for (uint32_t k = SCOUR_N-1; k < len; ++k) {
		w <<= 2, w |= s[k]-1;
		uint32_t t = (w << SCOUR_R) >> SCOUR_R;
		for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=3) {
			uint32_t PX = *(uint32_t *)PXp & 0xFFFFFF;
			if (!Hash[PX]) Cache[(*cix)++] = PX;
			if (Hash[PX] == mmatch) {
				Hash[PX] = qdim + mmatch + 1;
				Refs[(*nref)++].v = PX;
			} 
			else if (Hash[PX] < UINT16_MAX) ++Hash[PX];
		}
	}
}
static void (*setAmbigScour)(uint16_t *, Split *, uint32_t *, uint32_t, void **, 
	char *, uint32_t, uint32_t, uint8_t, uint32_t *, uint32_t *) = setAmbigScour20;
static void (*setUnambigScour)(uint16_t *, Split *, uint32_t *, uint32_t,
	void **, char *, uint32_t, uint32_t, uint32_t *, uint32_t *) = setUnambigScour20;

static inline void storeUnambigWords(Split *Word_Ix, char *s, uint32_t len, uint32_t *wix, uint32_t j) {
	uint32_t w = *s - 1;
	for (uint32_t k = 1; k + 1 < SCOUR_N; ++k)
		w <<= 2, w |= s[k]-1;
	for (uint32_t k = SCOUR_N-1; k < len; ++k) {
		w <<= 2, w |= s[k]-1;
		uint32_t t = (w << SCOUR_R) >> SCOUR_R;
		Word_Ix[(*wix)++] = (Split){t,j};
	}
}

static void storeAmbigWords(Split *Word_Ix, char *S, uint32_t *wix, uint32_t j, uint32_t w, uint32_t ix) {
	if (ix == SCOUR_N) Word_Ix[(*wix)++] = (Split){w, j};
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		storeAmbigWords(Word_Ix, S, wix, j, w << 2 | AMBIGS[S[ix]][i], ix + 1);
}

static void postScour24(uint16_t *Hash, void **Forest, Split *Word_Ix,
 uint32_t wix, uint32_t *Cache, uint32_t *cix) {
	uint32_t dLast = 0, max = 0;
	for (uint32_t i = 1; i <= wix; ++i) {
		if (i == wix || Word_Ix[i].v != Word_Ix[dLast].v) {
			max = (i-dLast) > max ? (i-dLast) : max;
			uint32_t t = Word_Ix[i-1].v;
			for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=3) {
				uint32_t PX = *(uint32_t *)PXp & 0xFFFFFF;
				if (!Hash[PX]) Cache[(*cix)++] = PX;
				Hash[PX] = MIN(max+Hash[PX],UINT16_MAX);
			}
			dLast = i, max = 0;
		}
		else if (Word_Ix[i].i != Word_Ix[dLast].i) {
			max = (i-dLast) > max ? (i-dLast) : max;
			dLast = i; 
		}
	}
}
static void postScour20(uint16_t *Hash, void **Forest, Split *Word_Ix,
 uint32_t wix, uint32_t *Cache, uint32_t *cix) {
	uint32_t dLast = 0, max = 0;
	for (uint32_t i = 1; i <= wix; ++i) {
		if (i == wix || Word_Ix[i].v != Word_Ix[dLast].v) {
			max = (i-dLast) > max ? (i-dLast) : max;
			uint32_t t = Word_Ix[i-1].v;
			for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=5) {
				uint64_t PX = *(uint64_t *)PXp;// & 0xFFFFFF;
				uint32_t px1 = PX & 0xFFFFF, px2 = (PX >> 20) & 0xFFFFF;
				//printf("-----ref1: %u, ref2: %u\n",px1,px2);
				if (!Hash[px1]) Cache[(*cix)++] = px1;
				Hash[px1] = MIN(max+Hash[px1],UINT16_MAX);
				if (PXp + 3 >= Forest[t+1]) break;
				if (!Hash[px2]) Cache[(*cix)++] = px2;
				Hash[px2] = MIN(max+Hash[px2],UINT16_MAX);
			}
			dLast = i, max = 0;
		}
		else if (Word_Ix[i].i != Word_Ix[dLast].i) {
			max = (i-dLast) > max ? (i-dLast) : max;
			dLast = i; 
		}
	}
}

static void (*postScour)(uint16_t *, void **, Split *,
 uint32_t , uint32_t *, uint32_t *) = postScour20;
	
void addAmbigScour(Accelerant *F, uint32_t r, char *S, uint32_t w, uint8_t ix) {
	if (ix == SCOUR_N) {
		if (F[w].len == F[w].cap) // check on success of realloc here?
			F[w].Refs = realloc(F[w].Refs,(F[w].cap+=(F[w].cap>>3)+1)*sizeof(*F[w].Refs));
		if (!F[w].len || F[w].Refs[F[w].len-1] != r) F[w].Refs[F[w].len++] = r;
	}
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		addAmbigScour(F, r, S, w << 2 | AMBIGS[S[ix]][i], ix + 1);
}
void countAmbigScour(uint8_t *C, uint32_t *W, uint32_t *I, char *S, uint32_t w, uint8_t ix) {
	if (ix == SCOUR_N) { if (!(C[w >> 3] & (1 << (w & 7)))) 
		W[(*I)++] = w, C[w >> 3] |= (1 << (w & 7)); }
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		countAmbigScour(C, W, I, S, w << 2 | AMBIGS[S[ix]][i], ix + 1);
}

void make_accelerator(Reference_Data *Rd, char *xcel_FN) {
	FILE *out = fopen(xcel_FN,"wb");
	if (!out) {fprintf(stderr, "Cannot write accelerator '%s'\n",xcel_FN); exit(1);}
	char **RefSeq = Rd->RefSeq;
	uint32_t *RefIxSrt = Rd->RefIxSrt, *RefLen = Rd->RefLen, totR = Rd->totR,
		totRC = Rd->numRclumps, origTotR = Rd->origTotR;
	if (Z) fprintf(stderr,"Note: N-penalized accelerator not usable for unpenalized alignment\n");
	
	//THREADS = THREADS >> 1 ?: 1; //1; //THREADS > 4 ? 4 : THREADS;
	uint32_t *Bad = 0;
	double wtime = omp_get_wtime();
	uint32_t done = 0, badSz = 0;

	Accelerant *F = calloc(1 << (2*SCOUR_N), sizeof(*F));
	if (!F) {fputs("OOM:AccelerantF\n",stderr); exit(3);}
	const uint32_t AMLIM = Rd->skipAmbig ? 0 : 12;
	size_t fullSize = SCOUR_N > 14 ? INT32_MAX : (1 << 24); //(Rd->maxLenR*(1<<(2*AMLIM)))*16;
	//printf("Size of entity manager: %u\n",fullSize);
	uint32_t IPOW3[16] = {1,3,9,27,81,243,729,2187,6561,19683,59049,177147,
		531441,1594323,4782969,14348907};
	uint32_t IPOW4[16] = {1,4,16,61,256,1024,4096,16384,65536,262144,1048576,
		4194304,16777216,67108864,268435456,1073741824};
	uint32_t *IPOWX = Z ? IPOW3 : IPOW4;
	#pragma omp parallel num_threads(THREADS)
	{
		int tid = omp_get_thread_num(), skipAmbig = Rd->skipAmbig;
		uint8_t *WordsYN = calloc((1 << (2*SCOUR_N)) >> 3,sizeof(*WordsYN)); // 128MB
		uint32_t *WordCache = malloc(fullSize*sizeof(*WordCache)); 
		if (!WordsYN || !WordCache) {fputs("OOM:WordsYN\n",stderr); exit(3);}

		const uint32_t AMBIG = 4 + Z, RNG = SCOUR_N-1;
		#pragma omp for schedule(static,1) reduction(+:badSz)
		for (uint32_t i = 0; i < totRC; ++i) {
			if (!(done & 255)) printf("\rScanning accelerator [%u / %u]\r",done,totRC);
			uint32_t begin = i*VECSZ, end = MIN(totR, begin + VECSZ);
			uint64_t Tsum = 0;
			uint16_t doAmbig = 0;
			if (!skipAmbig) for (uint32_t z = begin; z < end; ++z) { // each seq
				uint32_t len = RefLen[RefIxSrt[z]]; 
				if (len < SCOUR_N) continue;
				char *s = RefSeq[RefIxSrt[z]]; 
				uint32_t Asum = 0;
				for (int j = 0; j < len; ++j) { // try to precompute sizes
					if (j >= RNG) {
						Tsum += IPOWX[Asum];
						if (s[j-RNG] > AMBIG) --Asum;
					}
					if (s[j] > AMBIG) ++Asum, doAmbig |= 1 << (z-begin);
					if (Tsum >= fullSize) {++badSz; goto MULTI_ACC_END;}
				}
			}
			uint32_t I = 0;
			for (uint32_t z = begin; z < end; ++z) {
				char *s = RefSeq[RefIxSrt[z]]; 
				uint32_t len = RefLen[RefIxSrt[z]];
				if (len < SCOUR_N) continue;
				if (skipAmbig) {
					for (uint32_t j = 0; j + SCOUR_N <= len; ++j) {
						uint32_t k = 0; 
						for (; k < SCOUR_N; ++k) if (s[j+k] >= 5) break;
						if (k < SCOUR_N) { j+= k; continue; }
						countAmbigScour(WordsYN, WordCache, &I, s+j, 0, 0);
					}
				}
				else if (Z) { // screen for N's first!
					for (uint32_t j = 0; j + SCOUR_N <= len; ++j) {
						uint32_t k = 0; 
						for (; k < SCOUR_N; ++k) if (s[j+k] == 5) break;
						if (k < SCOUR_N) { j+= k; continue; }
						countAmbigScour(WordsYN, WordCache, &I, s+j, 0, 0);
					}
				} else if (doAmbig << (16-(z-begin)) >> (z-begin)) 
					for (uint32_t j = 0; j + SCOUR_N <= len; ++j)
						countAmbigScour(WordsYN, WordCache, &I, s+j, 0, 0);
				else {
					uint32_t w = 0;
					for (uint32_t j = 0; j + 1 < SCOUR_N; ++j)
						w <<= 2, w |= s[j]-1;
					for (uint32_t j = SCOUR_N - 1; j < len; ++j) {
						w <<= 2, w |= s[j]-1;
						uint32_t t = (w << SCOUR_R) >> SCOUR_R;
						if (!(WordsYN[t >> 3] & (1 << (t & 7)))) 
							WordCache[I++] = t, WordsYN[t >> 3] |= (1 << (t & 7));
					}
				}
			}
			//printf("%u\t%u\n",i,I); // Clump i had I unique accelerants
			for (; I; --I) {
				#pragma omp atomic
				++F[WordCache[I-1]].cap; 
				WordsYN[WordCache[I-1] >> 3] = 0;
			}
			MULTI_ACC_END:NULL;
			#pragma omp atomic
			++done;
		}
		#pragma omp single
		{ // create data structures of the perfect size for bads and words
			printf("Finished scanning pass [%f].          \n",omp_get_wtime() - wtime);
			for (uint32_t i = 0; i < 1 << (2*SCOUR_N); ++i) if (F[i].cap) {
				F[i].Refs = malloc(F[i].cap*sizeof(*F[i].Refs));
				if (!F[i].Refs) {fputs("OOM:Accelerant.Refs\n",stderr); exit(3);}
			}
			printf("Number bad (initial) = %u\n",badSz);
			Bad = malloc(badSz*sizeof(*Bad)); 
			if (!Bad) {fputs("OOM:Bad\n",stderr); exit(3);}
			badSz = 0; // use as ix
			done = 0;
			wtime = omp_get_wtime();
		}
		#pragma omp for schedule(static,1)
		for (uint32_t i = 0; i < totRC; ++i) {
			if (!(done & 255)) printf("\rFilling accelerator [%u / %u]\r",done,totRC);
			
			uint32_t begin = i*VECSZ, end = MIN(totR, begin + VECSZ);
			uint64_t Tsum = 0;
			uint16_t doAmbig = 0;
			if (!skipAmbig) for (uint32_t z = begin; z < end; ++z) { // each seq
				uint32_t len = RefLen[RefIxSrt[z]]; 
				if (len < SCOUR_N) continue;
				char *s = RefSeq[RefIxSrt[z]]; 
				uint32_t Asum = 0;
				for (int j = 0; j < len; ++j) { // try to precompute sizes
					if (j >= RNG) {
						Tsum += IPOWX[Asum];
						if (s[j-RNG] > AMBIG) --Asum;
					}
					if (s[j] > AMBIG) ++Asum, doAmbig |= 1 << (z-begin);
					if (Tsum >= fullSize) {
						uint32_t bix;
						#pragma omp atomic capture
						bix = badSz++;
						Bad[bix] = i; 
						goto MULTI_ACC_END2;
					}
				}
			}
			uint32_t I = 0;
			for (uint32_t z = begin; z < end; ++z) { 
				char *s = RefSeq[RefIxSrt[z]]; 
				uint32_t len = RefLen[RefIxSrt[z]];
				if (len < SCOUR_N) continue;
				if (skipAmbig) {
					for (uint32_t j = 0; j + SCOUR_N <= len; ++j) {
						uint32_t k = 0; 
						for (; k < SCOUR_N; ++k) if (s[j+k] >= 5) break;
						if (k < SCOUR_N) { j+= k; continue; }
						countAmbigScour(WordsYN, WordCache, &I, s+j, 0, 0);
					}
				} else if (Z) { // screen for N's first!
					for (uint32_t j = 0; j + SCOUR_N <= len; ++j) {
						uint32_t k = 0; 
						for (; k < SCOUR_N; ++k) if (s[j+k] == 5) break;
						if (k < SCOUR_N) { j+= k; continue; }
						countAmbigScour(WordsYN, WordCache, &I, s+j, 0, 0);
					}
				} else if (doAmbig << (16-(z-begin)) >> (z-begin)) 
					for (uint32_t j = 0; j + SCOUR_N <= len; ++j)
						countAmbigScour(WordsYN, WordCache, &I, s+j, 0, 0);
				else {
					uint32_t w = 0;
					for (uint32_t j = 0; j + 1 < SCOUR_N; ++j)
						w <<= 2, w |= s[j]-1;
					for (uint32_t j = SCOUR_N - 1; j < len; ++j) {
						w <<= 2, w |= s[j]-1;
						uint32_t t = (w << SCOUR_R) >> SCOUR_R;
						if (!(WordsYN[t >> 3] & (1 << (t & 7)))) 
							WordCache[I++] = t, WordsYN[t >> 3] |= (1 << (t & 7));
					}
				}
			}
			//printf("%u\t%u\n",i,I); // Clump i had I unique accelerants
			for (; I; --I) {
				uint32_t wix;
				#pragma omp atomic capture
				wix = F[WordCache[I-1]].len++;
				F[WordCache[I-1]].Refs[wix] = i;
				WordsYN[WordCache[I-1] >> 3] = 0;
			}
			
			MULTI_ACC_END2:NULL;
			#pragma omp atomic
			++done;
		}
		free(WordsYN); free(WordCache);
	}
	printf("\nDone scanning [%f]; writing...\n",omp_get_wtime() - wtime);
	wtime = omp_get_wtime();
		
	size_t totalWords = 0;
	uint32_t num_k = 1 << (2*SCOUR_N);
	#pragma omp parallel for reduction(+:totalWords)
	for (uint32_t i = 0; i < num_k; ++i) totalWords += F[i].len;
	printf("Total accelerants stored: %llu (%u bad)\n",totalWords, badSz);

	// redux
	
	uint8_t vers = (1 << 7) | (Z << 6) | (totRC > 1048574 ? ACC_VERSION_BIG : ACC_VERSION);//vers = (1 << 7) | (Z << 6) | SCOUR_N;
	setvbuf(out,0,_IOFBF,1<<21);
	fwrite(&vers,sizeof(vers),1,out);
	fwrite(&badSz, sizeof(badSz), 1, out);
	for (uint32_t i = 0; i < num_k; ++i) 
		fwrite(&F[i].len, sizeof(F[i].len), 1, out);
	

	//for (uint32_t i = 0; i < num_k; ++i) 
	//	fwrite(F[i].Refs, sizeof(*F[i].Refs), F[i].len, out);
	//uint64_t count = 0;
	if (totRC > 16777214) {fputs("ERROR: acc error M16\n",stderr); exit(101);}
	else if (totRC > 1048574) {
		fputs(" --> [Re-Accel] Writing LARGE format acx...\n", stderr);
		//for (size_t i = 0; i < totalWords; ++i) // replace these "totalwords" things with the commented construct above
		//	fwrite(WordDump + i*sizeof(uint32_t), 3, 1, out);
		//count = totalWords;
		for (uint32_t i = 0; i < num_k; ++i)
			for (uint32_t k = 0; k < F[i].len; ++k)
				fwrite(&F[i].Refs[k],3,1,out);

	} else {
		fputs(" --> [Re-Accel] Writing SMALL format acx...\n", stderr);
		for (size_t i = 0; i < num_k; ++i) {
			for (uint32_t *P = F[i].Refs; P < F[i].Refs+F[i].len; P+=2) {
				uint64_t bay = *P;
				if (P+1 < (uint32_t*)F[i].Refs+F[i].len) {
					bay |= (uint64_t)*(P+1) << 20;
					fwrite(&bay, 1, 5, out); 
				}
				else fwrite(&bay, 1, 3, out); 
			}
		}
	}

	fwrite(Bad, sizeof(*Bad), badSz, out);
	printf("Wrote accelerator (%f).\n",omp_get_wtime()-wtime);
}

void read_accelerator(Reference_Data *Rd, char *xcel_FN) {
	FILE *in = fopen(xcel_FN,"rb");
	if (!in) {fprintf(stderr, "Cannot read accelerator '%s'\n",xcel_FN); exit(1);}
	
	uint8_t cb = fgetc(in), 
		dbVer = (uint8_t)(cb << 4) >> 4,
		didZ = (uint8_t)(cb << 1) >> 7; 
	if (didZ && !Z) {
		fprintf(stderr,"ERROR: Accelerator built without '-y'; can't use '-y'\n");
		exit(1);
	}
	int slen = 0;
	if (cb < 128 || (dbVer != ACC_VERSION && dbVer != ACC_VERSION_BIG)) {
		fprintf(stderr,"ERROR: invalid accelerator [%u:%u]\n", cb, dbVer);
		if (cb >= 128 && dbVer == SCOUR_N) 
			puts("Old accelerator version detected. Converting to .acx..."),
			slen = strlen(xcel_FN);
		else exit(1);
	}
	uint32_t szBL, hashSz = 1 << (2*SCOUR_N), *BadList = 0, 
		*Lens = malloc(hashSz*sizeof(*Lens)); 
	fread(&szBL,sizeof(szBL),1,in);
	BadList = malloc(szBL*sizeof(*BadList));
	if (!BadList || !Lens) {fputs("OOM:BadList_rd\n",stderr); exit(3);}
	printf(" --> [Accel] Accelerator found. Parsing...\n");
	fread(Lens, sizeof(*Lens), hashSz, in);
	uint32_t chunk = slen ? sizeof(uint32_t)*8 : 20;
	size_t sumBytes = 0, sumWords = 0;
	for (uint32_t i = 0; i < hashSz; ++i) sumWords += Lens[i];
	if (dbVer == ACC_VERSION) for (uint32_t i = 0; i < hashSz; ++i) 
		sumBytes += (Lens[i] / 2u) * 5 + (Lens[i] & 1) * 3;
	else if (dbVer == ACC_VERSION_BIG) sumBytes = sumWords * 3;
	else sumBytes = sumWords * 4;
	printf(" --> [Accel] Total accelerants: %llu [bytes = %llu]\n",sumWords,sumBytes);
	
	void *WordDump = malloc(sumBytes + 16);
	if (!WordDump) {fputs("OOM:WordDump_rd\n",stderr); exit(3);}
	fread(WordDump,1,sumBytes,in);
	printf(" --> [Accel] Reading %u ambiguous entries\n",szBL);
	fread(BadList,sizeof(*BadList),szBL,in);

	// pointer mapping
	void **Forest = malloc((1+hashSz)*sizeof(*Forest));
	if (!Forest) {fputs("OOM:Forest_rd\n",stderr); exit(3);}
	*Forest = WordDump;
	
	if (slen) { // convert acc to new acx format
		for (uint32_t i = 1; i <= hashSz; ++i)
			Forest[i] = Forest[i-1] + Lens[i-1] * sizeof(uint32_t);

		if (slen >= 4096) {puts("Filename error."); exit(1);}
		char newName[4096] = {0};
		strcpy(newName,xcel_FN);
		if (!strcmp(newName+slen-4,".acc")) newName[slen-1] = 'x';
		else strcpy(newName+slen,".acx");
		uint64_t max = 0;
		for (uint32_t x = 0; x < hashSz; ++x) if (Lens[x] > max) max = Lens[x];
		printf("Maximum length record in this database is %llu\n",max);
		max = 0; uint32_t *WD32 = WordDump; 
		for (uint64_t x = 0; x < sumWords; ++x) 
			if (max < WD32[x]) max = WD32[x];
		printf("Maximum pivot: %llu\n",max);
		

		FILE *out = fopen(newName,"wb");
		if (!out) {fprintf(stderr,"ERROR: Cannot write '%s'\n",newName); exit(1);}
		uint8_t vers = (1 << 7) | (didZ << 6) | (max > 1048574 ? ACC_VERSION_BIG : ACC_VERSION);
		fwrite(&vers,sizeof(vers),1,out);
		fwrite(&szBL, sizeof(szBL), 1, out);
		fwrite(Lens, sizeof(*Lens), hashSz, out);

		uint64_t count = 0;
		if (max > 16777214) {fputs("ERROR: acc error 101\n",stderr); exit(101);}
		else if (max > 1048574) {
			fputs(" --> [Re-Accel] Writing LARGE format acx...\n", stderr);
			for (size_t i = 0; i < sumWords; ++i)
				fwrite(WordDump + i*sizeof(uint32_t), 3, 1, out);
			count = sumWords;
		} else {
			fputs(" --> [Re-Accel] Writing SMALL format acx...\n", stderr);
			for (size_t i = 0; i < hashSz; ++i) {
				for (uint32_t *P = (uint32_t*)Forest[i]; P < (uint32_t*)Forest[i+1]; P+=2) {
					uint64_t bay = *P;
					bay |= (uint64_t)*(P+1) << 20;
					count += (P+1 < (uint32_t*)Forest[i+1]) ? 5 : 3;
					fwrite(&bay, 1, (P+1 < (uint32_t*)Forest[i+1]) ? 5 : 3, out); 
				}
			}
		}
		fwrite(BadList, sizeof(*BadList), szBL, out);
		printf("Re-wrote accelerator as '%s' [%llu byte table]\n",newName,count);
		puts("Re-run the program with the new ACX accelerator.");
		exit(0);
	}
	if (dbVer == ACC_VERSION) for (uint32_t i = 1; i <= hashSz; ++i)
		Forest[i] = Forest[i-1] + (Lens[i-1] / 2u) * 5 + (Lens[i-1] & 1) * 3;
	else {
		for (uint32_t i = 1; i <= hashSz; ++i)
			Forest[i] = Forest[i-1] + Lens[i-1] * 3;
		setAmbigScour = setAmbigScour24;
		setUnambigScour = setUnambigScour24;
		postScour = postScour24;
		puts(" --> [Accel] Using LARGE db mode...");
	}
	free(Lens); 
	Rd->Accelerators = Forest;
	Rd->badListSz = szBL;
	Rd->BadList = BadList;
}

static inline void do_alignments(FILE *output, Reference_Data RefDat, Query_Data QDat, int usedb) {
	#ifndef __SSSE3__
	create_sse2_profiles(&RefDat);
	#endif

	// Extract the variables
	uint32_t maxLenR = RefDat.maxLenR, totR = RefDat.totR, numRclumps = RefDat.numRclumps,
		*ClumpLen = RefDat.ClumpLen, *RefStart = RefDat.RefStart, *RefIxSrt = RefDat.RefIxSrt,
		*RefDedupIx = RefDat.RefDedupIx, *TmpRIX = RefDat.TmpRIX;
	DualCoil **ProfClump = RefDat.ProfClump, **RefClump = RefDat.RefClump;
	uint32_t totQ = QDat.totQ, maxLenQ = QDat.maxLenQ, numUniqQ = QDat.numUniqQ,
		*Offset = QDat.Offset, maxED = QDat.maxED, maxDiv = QDat.maxDiv;
	UniBin *UniBins = QDat.UniBins; 
	ShrBin *ShrBins = QDat.ShrBins;
	char **RefHead = RefDat.RefHead, **QHead = QDat.QHead, *SeqDumpRC = QDat.SeqDumpRC;
	PackaPrince Fp = QDat.FingerprintsQ, FpR = RefDat.FingerprintsR;
	Prince *Centroids = RefDat.Centroids, *RP = 0; 
	uint32_t taxa_parsed = RefDat.taxa_parsed; 
	TaxPair_t *Taxonomy = RefDat.Taxonomy;
	int taxasuppress = QDat.taxasuppress, doRC = QDat.rc, skipAmbig = QDat.skipAmbig;
	uint32_t *BadList = RefDat.BadList, szBL = RefDat.badListSz, *QBins = QDat.QBins;
	void **Forest = RefDat.Accelerators;
	int quiet = QDat.quiet;

	// Prepare dimensions, bounds
	uint32_t qdim = maxLenQ + 1, rdim = maxLenR + 1;
	++qdim; ++rdim; // add padding to right and bottom to eliminate loop logic

	uint32_t newUniqQ = numUniqQ + ((doRC && !DO_PREPASS) ? numUniqQ : 0);
	if (Xalpha && DO_FP) 
		puts("WARNING: Fingerprints are incompatible with Xalpha and will be disabled."),
		DO_FP = 0;
	if (DO_FP) {
		puts("Creating unique query fingerprints...");
		//if (!Fp.P) Fp = create_fingerprints(UniqQSeq, newUniqQ, UniqQLen, 0, 0, 1); //, NULL);
		if (!Fp.P) Fp = create_fingerprintsQ(ShrBins, UniBins, newUniqQ); 
		printf("There were %u fingerprints created for the %u unique queries.\n", Fp.nf, newUniqQ);
		uint64_t accum_fp_n = 0, accum_len = 0; // Do stats
		#pragma omp parallel for reduction(+:accum_fp_n)
		for (uint32_t i = 0; i < Fp.nf; ++i) /*printf("%u\n", Fp.N[i]),*/ accum_fp_n += Fp.N[i];
		#pragma omp parallel for reduction(+:accum_len)
		for (uint32_t i = 0; i < numUniqQ; ++i) accum_len += ShrBins[i].len;
		if (newUniqQ > numUniqQ) accum_len *= 2;
		printf("Query len: %f, FP density: %f, coverage: %f\n", (double)accum_len/newUniqQ, 
			(double)accum_fp_n/Fp.nf, ((double)accum_fp_n/Fp.nf)/((double)accum_len/newUniqQ));

		if (Z && FpR.nf) { // Ensure the individual fingerprints used are N-penalized
			#pragma omp parallel for
			for (uint32_t i = 0; i < totR; ++i) FpR.P[i] = FpR.P[FpR.Ptrs[i]];
		}
		free(FpR.Ptrs); // FPs are now ordered by the correct ambiguity
		printf("There were %u fingerprints created for the %u references.\n", FpR.nf ?: totR, totR);
		RP = FpR.P;
	} else printf("Fingerprints not enabled\n");
	
	//uint32_t *UniqQLen, *UniqQed, *UniqDiv, *NewIX, maxDiv, *RCMap, *Umap;
	//char **UniqQSeq;
	// TODO: UNCOMMENT THE IF STATEMENT!
	//if (RUNMODE == ANY || DO_PREPASS) { // preempt TmpRIX creation for inline printing
		if (!RefIxSrt && RefDedupIx) {
			RefIxSrt = malloc(totR * sizeof(*RefIxSrt));
			if (!RefIxSrt) {fputs("OOM:[DA]RefIxSrt\n",stderr); exit(3);}
			for (uint32_t i = 0; i < totR; ++i) RefIxSrt[i] = TmpRIX[RefDedupIx[i]];
		}
		else if (!RefIxSrt && !RefDedupIx) RefIxSrt = TmpRIX;
	//}
	
	Split *RefOrder = 0; // NEW - disable
	if (DO_PREPASS) {
		printf("Engaging prepass mode.\n");

		double wtime = omp_get_wtime();
		
		/*RefOrder = malloc(numRclumps*sizeof(*RefOrder));
		if (!RefOrder) {fputs("OOM:RefOrder\n",stderr); exit(3);}
		for (uint32_t i = 0; i < numRclumps; ++i) RefOrder[i] = (Split){0,i};*/

		/*Prince *p_r = FpR.P, *p_q = Fp.P;
		uint8_t *RPops = malloc(totR*sizeof(*RPops));
		uint32_t inc = reorder_only ? 16 : 1;*/
		#define ITER DO_PREPASS
		#define TOPSORTD(n,sc,ix,M,I) { \
			uint32_t i = 0, t, v, s = sc, x = ix; \
			for (; i < n; ++i) if (s > M[i]) break; \
			for (; i < n; ++i) \
				t=M[i], M[i]=s, s = t, \
				v=I[i], I[i]=x, x = v; \
		}
		// redo the main prune_ed_mat16 function!
		uint32_t doneQ = 0, quantQ = numUniqQ/100000, critQ = quantQ;
			double rUniqQ = 100.0 / numUniqQ;
		#pragma omp parallel
		{
			DualCoil *Matrices = malloc(2*rdim*sizeof(*Matrices));
			DualCoil *rclump = malloc((2+rdim)*sizeof(*rclump));
			uint16_t *Hash = calloc(numRclumps,sizeof(*Hash));
			uint32_t *Cache = calloc(numRclumps,sizeof(*Cache)),
				revlen = rdim+(15 & (16 - (rdim & 15)))+16;
			char *Rev = malloc(revlen); 
			char taxBin[4096] = {0};
			uint32_t *IXTray, *MXTray, //RefMatch[ITER << 4], RefMatchB[ITER << 4],
				FM[ITER], RM[ITER], FI[ITER], RI[ITER];
			uint8_t RefMin[ITER << 4], RefMinB[ITER << 4];
			*RM = 0;
			uint32_t attenuate = 8; // earlyterm is attenuate/ITER, so 8/16 = 1/2 for 16 iters, 1/4 for 32
			int tid = omp_get_thread_num();
			
			if (!Matrices || !rclump || !Hash || !Cache) {fputs("OOM:Prepass\n",stderr); exit(3);}
			/*
			char *BAK = "\0ACGTNKMRY S W B V H D\0"; 
			              0123456789101112131415
			char *RVC = "\0TGCANMKYR S W V B D H\0";
						  0432157698101113121514
			*/
			__m128i tr_v = _mm_setr_epi8(0,4,3,2,1,5,7,6,9,8,10,11,13,12,15,14);
			__m128i rv_v = _mm_setr_epi8(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
			
			#pragma omp for schedule(dynamic,1)
			for (uint32_t i = 0; i < numUniqQ; ++i) {
				UniBin Ub = UniBins[i]; ShrBin Sb = ShrBins[Ub.six];
				uint32_t len = Sb.len, err = Sb.ed;
				char *s = Ub.Seq; 

				uint32_t w = 0, wix = 0, n = 0; 
				//uint32_t dbgCnt = 0, letCnt = 0;
				if (postScour == postScour24) for (uint32_t j = 0; j < len; ++j) {
					if (s[j] > 4) {n = 0; continue;}
					w <<= 2, w |= s[j]-1;
					if (++n >= SCOUR_N) {
						uint32_t t = (w << SCOUR_R) >> SCOUR_R;
						for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=3) {
							uint32_t PX = *(uint32_t *)PXp & 0xFFFFFF;
							if (!Hash[PX]++) Cache[wix++] = PX; 
							//++dbgCnt;
						}
						//++letCnt;
					}
				}
				else for (uint32_t j = 0; j < len; ++j) {
					if (s[j] > 4) {n = 0; continue;}
					w <<= 2, w |= s[j]-1;
					if (++n >= SCOUR_N) {
						uint32_t t = (w << SCOUR_R) >> SCOUR_R;
						for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=5) {
							uint64_t PX = *(uint64_t *)PXp;
							uint32_t px1 = PX & 0xFFFFF, px2 = (PX >> 20) & 0xFFFFF;
							if (!Hash[px1]++) Cache[wix++] = px1; 
							if (PXp + 3 >= Forest[t+1]) break;
							if (!Hash[px2]++) Cache[wix++] = px2;
						}
					}
				}
				//printf("Query (forward dir) %u had %llu matches in %u refclumps [%llu let].\n",
				//	i, dbgCnt, wix, letCnt);
				
				for (uint32_t j = 0; j < ITER; ++j) FM[j] = 0;
				for (uint32_t j = 0; j < wix; ++j) {
					uint32_t c = Cache[j], h = Hash[c];
					TOPSORTD(ITER,h,c,FM,FI);
					Hash[c] = 0;
				}
				//printf("[%u F]: "); for (int j = 0; j < ITER; ++j) printf("%d: %u [%u], ",j,FM[j],FI[j]); printf("\n");
				
				char *rcPointer = Rev + revlen;
				if (doRC) {
					//printf("%u before: ",i); for (uint32_t j = 0; j < len; ++j) printf("%u ",s[j]); puts("");
					for (uint32_t j = 0; j < len; j+=16) {
						__m128i src = _mm_load_si128((void*)s + j);
						__m128i tr = _mm_shuffle_epi8(tr_v, src); 
						__m128i sh = _mm_shuffle_epi8(tr, rv_v); 
						_mm_store_si128((void*)(rcPointer -= 16),sh);
					}
					rcPointer += 15 & (16 - (len & 15));
					//printf("%u after:  ",i); for (uint32_t j = 0; j < len; ++j) printf("%u ",rcPointer[j]); puts("");

					w = 0, wix = 0, n = 0; 
					s = rcPointer;
					//uint32_t dbgCnt = 0, letCnt = 0;
					if (postScour == postScour24) for (uint32_t j = 0; j < len; ++j) {
						if (s[j] > 4) {n = 0; continue;}
						w <<= 2, w |= s[j]-1;
						if (++n >= SCOUR_N) {
							uint32_t t = (w << SCOUR_R) >> SCOUR_R;
							for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=3) {
								uint32_t PX = *(uint32_t *)PXp & 0xFFFFFF;
								if (!Hash[PX]++) Cache[wix++] = PX; 
								//++dbgCnt;
							}
							//++letCnt;
						}
					}
					else for (uint32_t j = 0; j < len; ++j) {
						if (s[j] > 4) {n = 0; continue;}
						w <<= 2, w |= s[j]-1;
						if (++n >= SCOUR_N) {
							uint32_t t = (w << SCOUR_R) >> SCOUR_R;
							for (void *PXp = Forest[t]; PXp < Forest[t+1]; PXp+=5) {
								uint64_t PX = *(uint64_t *)PXp;
								uint32_t px1 = PX & 0xFFFFF, px2 = (PX >> 20) & 0xFFFFF;
								if (!Hash[px1]++) Cache[wix++] = px1; 
								if (PXp + 3 >= Forest[t+1]) break;
								if (!Hash[px2]++) Cache[wix++] = px2;
							}
						}
					}
					//printf("Query (rev dir) %u had %llu matches in %u refclumps [%llu let].\n",
					//	i, dbgCnt, wix, letCnt);
					for (uint32_t j = 0; j < ITER; ++j) RM[j] = 0;
					for (uint32_t j = 0; j < wix; ++j) {
						uint32_t c = Cache[j], h = Hash[c];
						TOPSORTD(ITER,h,c,RM,RI);
						Hash[c] = 0;
					}
					//printf("[%u R] Max1 %u [%d] Max2 %u [%d]\n",i,RM[0], RI[0], RM[1], RI[1], RM[2], RI[2]);
				}
				#pragma omp atomic
				++doneQ;
				
				// Decide whether to use the RC and align accordingly
				if (!*FM && !*RM) continue;
				
				char *query = *FM >= *RM ? Ub.Seq : rcPointer;
				/* if (DO_HEUR) query = *FM >= *RM ? Ub.Seq : rcPointer;
				else if (doRC) {
					if (*FM > *RM << 1) query = Ub.Seq;
					else if (*RM > *FM << 1) query = rcPointer;
					else {
						uint32_t potF = 0, potR = 0;
						for (uint32_t j = 0; j < ITER; ++j) potF+=FM[j], potR+=RM[j];
						query = potF >= potR ? Ub.Seq : rcPointer;
					}
				} */
				
				if (query == Ub.Seq) IXTray = FI, MXTray = FM;
				else IXTray = RI, MXTray = RM;
				
				uint32_t gmin = -1, kload = err * SCOUR_N + SCOUR_N, 
					mmatch = kload < len ? len - kload : 0, 
					load = MXTray[0] * attenuate / ITER; // ITER was LATENCY
				load = MIN(MXTray[0],load);
				int p = 0; for (; p < ITER; ++p) {
					if (MXTray[p] <= mmatch || MXTray[p] < load) break;
					uint32_t ri = IXTray[p];
					DualCoil *RefSlide = RefClump[ri];
					for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
						__m128i org = _mm_stream_load_si128((void*)(RefSlide++));
						__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
						__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
						_mm_store_si128((void*)(rclump+w),ex1);
						_mm_store_si128((void*)(rclump+w+1),ex2);
					}
					uint32_t errs = len - MXTray[p] - SCOUR_N + 1;
					if (RUNMODE != FORAGE) err = MIN(gmin,err);
					errs = MIN(errs,err);
					uint8_t *MinA = RefMin + (p << 4);
					//int starred = 0;
					uint32_t min = prune_ed_mat16(rclump, query, ClumpLen[ri],
						len, rdim, Matrices, ProfClump[ri], errs, MinA);
					if (errs < err && min == -1) //printf("Retried "),
						min = prune_ed_mat16(rclump, query, ClumpLen[ri],
							len, rdim, Matrices, ProfClump[ri], err, MinA);
					gmin = gmin < min ? gmin : min;
					//printf("min = %d on clump %u\n",min,p);
					if (min == -1) _mm_store_si128((void*)MinA,_mm_set1_epi8(-1));
					else if (RUNMODE == ANY) {++p; break;}

					// TODO: both forward and reverse read scanning in all cases []
					// (normal mode becomes current -hr; non-hr becomes the above "both passes")
					// TODO: accept queries directly w/o processing, in order []
					
				}
				if (gmin == -1) {
					if (DO_HEUR || !doRC) continue;
					if (query == Ub.Seq) query = rcPointer;
						else query = Ub.Seq;
					if (query == Ub.Seq) IXTray = FI, MXTray = FM;
					else IXTray = RI, MXTray = RM;
					
					// regurgitate the upper loop here. maybe not only if gmin == -1, but in general too? 
					load = MXTray[0] * attenuate / ITER; // ITER was LATENCY
					load = MIN(MXTray[0],load);
					int p_bak = p;
					p = 0; for (; p < ITER; ++p) {
						if (MXTray[p] <= mmatch || MXTray[p] < load) break;
						uint32_t ri = IXTray[p];
						DualCoil *RefSlide = RefClump[ri];
						for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
							__m128i org = _mm_stream_load_si128((void*)(RefSlide++));
							__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
							__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
							_mm_store_si128((void*)(rclump+w),ex1);
							_mm_store_si128((void*)(rclump+w+1),ex2);
						}
						uint32_t errs = len - MXTray[p] - SCOUR_N + 1;
						if (RUNMODE != FORAGE) err = MIN(gmin,err);
						errs = MIN(errs,err);
						uint8_t *MinA = RefMin + (p << 4);
						uint32_t min = prune_ed_mat16(rclump, query, ClumpLen[ri],
							len, rdim, Matrices, ProfClump[ri], errs, MinA);
						if (errs < err && min == -1) //printf("Retried "),
							min = prune_ed_mat16(rclump, query, ClumpLen[ri],
								len, rdim, Matrices, ProfClump[ri], err, MinA);
						gmin = gmin < min ? gmin : min;
						//printf("min = %d on clump %u (kmatches %u)\n",min,p,MXTray[p]);
						if (min == -1) _mm_store_si128((void*)MinA,_mm_set1_epi8(-1));
						else if (RUNMODE == ANY) {++p; break;}
					}
					// if made for gmin more than the ==-1 case, 
					// store into the "B" container for RefMin,
					// do something with p and p_bak, and/or consolidate
					// or pad RefMin if necessary and make it double size to store these too?
					// or keep separate but do prepass to pick between them? if so, double loop for forage?
					// best: integrate into one loop. at end of loop, break if DO_HEUR otherwise swap containers
					// [swap occurs by going to upper half of containers; contiguous mins etc]
					// [print detects rc by index; higher than ITER or not?]
				}
				if (gmin == -1) continue; // still
				
				uint32_t ceil = err, k = 0;
				if (RUNMODE != FORAGE) ceil = MIN(gmin,ceil); 
				char *taxon = "";
				if (RUNMODE == CAPITALIST /* || RUNMODE == COMMUNIST */) {
					uint32_t minIX = -1, dv = 0, olen = 0;
					for (uint32_t j = 0; j < p << 4; ++j) if (RefMin[j] <= ceil) {
						uint32_t orix = (IXTray[j >> 4] << 4) + (j & 15);
						//printf("--> Matches with ref %u [%s]\n",orix,RefHead[RefIxSrt[orix]]);
						if (taxa_parsed) for (uint32_t z = RefDedupIx[orix]; z < RefDedupIx[orix+1]; ++z) {
							uint32_t rix = TmpRIX[z];
							if (minIX == -1) 
								strncpy(taxBin,taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy),4096),
								olen = strlen(taxBin);
							else {
								char *tp = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
								dv = 0;
								while (taxBin[dv] && tp[dv] && taxBin[dv]==tp[dv]) ++dv;
								taxBin[dv] = 0;
							}
						}
						if (orix < minIX) minIX = orix, k = j;
					}
					if (taxa_parsed) {
						taxon = taxBin;
						if (strlen(taxon) < olen) {
							while (dv && taxon[dv] != ';') --dv;
							taxon[dv] = 0;
						}
					}
				}
				for (; k < p << 4; ++k) if (RefMin[k] <= ceil) {
					uint32_t orix = (IXTray[k >> 4] << 4) + (k & 15);
					double fakeID = (double)((int)len-RefMin[k])/len * 100.0;
					if ((RUNMODE == FORAGE || RUNMODE == ALLPATHS) && RefDedupIx) 
						for (uint32_t z = RefDedupIx[orix]; z < RefDedupIx[orix+1]; ++z) {
							uint32_t rix = TmpRIX[z],
								stIxR = RefStart ? RefStart[rix] : 1,
								edIxR = stIxR + ClumpLen[IXTray[k >> 4]], t;
							if (taxa_parsed) taxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
							if (query != Ub.Seq) t = stIxR, stIxR = edIxR, edIxR = t;
							for (uint32_t j = Offset[Ub.six]; j < Offset[Ub.six+1]; ++j) 
								fprintf(output,"%s\t%s\t%f\t%u\t%u\t-1\t%u\t%u\t%d\t%u\t%u\t%u\t%s\n", 
									QHead[j], RefHead[rix], fakeID, 
									len + RefMin[k], RefMin[k], 1, len, stIxR, edIxR, 
									RefMin[k],j > Offset[Ub.six], taxon);
					} else {
						uint32_t rix = RefIxSrt[orix],
							stIxR = RefStart ? RefStart[rix] : 1,
							edIxR = stIxR + ClumpLen[IXTray[k >> 4]], t;
						if (taxa_parsed && taxon != taxBin) 
								taxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
						if (query != Ub.Seq) t = stIxR, stIxR = edIxR, edIxR = t;
						for (uint32_t j = Offset[Ub.six]; j < Offset[Ub.six+1]; ++j) 
							fprintf(output,"%s\t%s\t%f\t%u\t%u\t-1\t%u\t%u\t%d\t%u\t%u\t%u\t%s\n", 
								QHead[j], RefHead[rix], fakeID, 
								len + RefMin[k], RefMin[k], 1, len, stIxR, edIxR, 
								RefMin[k],j > Offset[Ub.six], taxon);
						if (RUNMODE == BEST || RUNMODE == CAPITALIST || RUNMODE == ANY) break;
					}
				}
				
				if (!tid && doneQ >= critQ) {
					critQ = doneQ + quantQ;
					if (!quiet) printf("Progress: %.3f%%\r",(double)doneQ*rUniqQ);
				}
			}
			free(Matrices);
			free(Hash);
		}
		
		printf("Time to perform prepass: %f\n",omp_get_wtime()-wtime);
		exit(101);

	}

	// Prep the main alignment data structures
	typedef struct ResultPod ResultPod;
	struct ResultPod {
		float score; // 4
		uint32_t refIx, finalPos; //8
		uint8_t numGapR, numGapQ, mismatches, rc; //4
		ResultPod *next; //8
	};
	
	ResultPod ***ThreadPods = malloc(THREADS*sizeof(*ThreadPods));
	ResultPod **BasePod = 0;

	uint32_t totDone = 0, tid = 0;
	uint64_t totSkipped = 0;

	

if (DO_ACCEL) {
	uint32_t QBUNCH = newUniqQ / (THREADS*128);
	if (QBUNCH > 16) QBUNCH = 16;
	if (!QBUNCH) QBUNCH = 1;
	printf("Setting QBUNCH to %u\n",QBUNCH);
	double totQMult = 1.0 / newUniqQ;
	
	ResultPod **Pods = calloc(newUniqQ,sizeof(*Pods)); // global container (accel only)!
	if (!Pods) {fputs("OOM:Pods\n",stderr); exit(3);}
	BasePod = Pods;
	int castUcmp(const void *a, const void *b) {
		return *(int32_t *)a - *(int32_t *)b; }
	int RefCmp(const void *a, const void *b) {
		Split *A = (Split *)a, *B = (Split *)b;
		return A->i > B->i ? -1 : B->i > A->i;
	}
	int WrdCmp2(const void *a, const void *b) {
		Split *A = (Split *)a, *B = (Split *)b;
		return A->v < B->v ? -1 : B->v < A->v ? 1 : 
			A->i < B->i ? -1 : B->i < A->i;
		// do this with if-elses and gage difference?
	}
	int WrdCmp(const void *a, const void *b) {
		Split *A = (Split *)a, *B = (Split *)b;
		if (A->v < B->v) return -1;
		else if (B->v < A->v) return 1;
		else if (A->i < B->i) return -1;
		else return B->i < A->i;
		//return A->v < B->v ? -1 : B->v < A->v ? 1 : 
		//	A->i < B->i ? -1 : B->i < A->i;
		// do this with if-elses and gage difference?
	}
	typedef struct {int32_t v, i;} SSplit;
	int WrdCmp3(const void *a, const void *b) {
		SSplit *A = (SSplit *)a, *B = (SSplit *)b;
		return A->v - B->v ?: A->i - B->i;
	}
	inline void iSort(Split *a, uint32_t n) {
		if (n > 24) qsort(a, n, sizeof(*a),RefCmp);
	    else for (uint32_t i = 1, j; i < n; ++i) {
	        Split key = a[i];
	        for (j = i; j && a[j-1].i < key.i; --j) ; //a[j] = a[j-1]; 
	        memmove(a + j + 1, a + j, sizeof(*a) * (i - j));
	        a[j] = key;
	    }
	}
	szBL = skipAmbig? 0:szBL; // experimental; skips ambiguous accelerants
	printf("Using ACCELERATOR to align %u unique queries...\n", QBins[1]);
	
	#pragma omp parallel
	{
		DualCoil *Matrices = calloc((cacheSz+2)*rdim,sizeof(*Matrices)),
				 *ScoresEX = malloc((2*rdim)*sizeof(*ScoresEX)),
				 *ShiftsEX = malloc((2*rdim)*sizeof(*ShiftsEX)),
				 *ShiftsBX = malloc((2*rdim)*sizeof(*ShiftsEX));
		uint32_t *HiBound = calloc(qdim+1,sizeof(*HiBound)),
				 *LoBound = calloc(qdim+1,sizeof(*LoBound)),
				 *StackE = malloc(qdim*sizeof(*StackE)), // num errors
				 *StackX = malloc(qdim*sizeof(*StackX)); // query ix in UniqQSeq
		if (!(Matrices && ScoresEX && ShiftsEX && ShiftsBX && HiBound && LoBound && StackE && StackX)) 
			{fputs("OOM:Aligner\n",stderr); exit(3);}
		for (int j = 0; j < (cacheSz+2); ++j) Matrices[j*rdim].v = _mm_set1_epi8(MIN(j*GAP,255));
		*LoBound = -1, LoBound[1] = 1;

		DualCoil *rclump = malloc((2+rdim)*sizeof(*rclump));

		uint16_t *Hash = calloc(numRclumps,sizeof(*Hash)), hqd = qdim;
		uint32_t *Cache = malloc(numRclumps*QBUNCH*sizeof(*Cache)); // to store all hit refs
		Split *Word_Ix = malloc((qdim+1024)*QBUNCH*SCOUR_N*sizeof(*Word_Ix)); // wordfill, sort, "dedupe", add refs in burst, delete

		//uint16_t *Srtr = calloc(numRclumps,sizeof(*Srtr));
		Split *Refs = calloc(numRclumps,sizeof(*Refs)), *RefPtr = Refs; 
		if (!rclump || !Hash || !Cache || !Refs || !Word_Ix) {fputs("OOM:Hash\n",stderr); exit(3);}
		uint32_t oldNref = 0;
		int tid = omp_get_thread_num();
		
		#pragma omp for schedule(dynamic,1)
		for (uint32_t z = 0; z < QBins[1]; z+= QBUNCH) {
			//memset(Hash,0,numRclumps*sizeof(*Hash));
			for (uint32_t j = 0; j < oldNref; ++j) Refs[j].i = 0; //Bucket[j] = 0; // Bucket memset
			uint32_t bound = MIN(z+QBUNCH,QBins[1]), nref = 0, min_mmatch = UINT32_MAX;
			uint32_t wix = 0, cix = 0, minlen = -1; // mlen NEW
			UniBins[z].div = 1; // Instead of inner loop; because it trails nothing before it
			for (uint32_t j = z; j < bound; ++j) {
				UniBin Ub = UniBins[j]; ShrBin Sb = ShrBins[Ub.six];
				uint32_t len = Sb.len, err = Sb.ed;
				minlen = len < minlen ? len : minlen; // NEW
				if (err == -1) continue;
				char *s = Ub.Seq; 
				uint32_t kload = err * SCOUR_N + SCOUR_N, 
					mmatch = kload < len ? len - kload : 0; // len - (err + 1) * SCOUR_N;
				if (DO_HEUR) DO_HEUR = (len >> 4) + 1u;
				mmatch = mmatch > DO_HEUR ? mmatch : DO_HEUR;
				if (mmatch < min_mmatch) min_mmatch = mmatch; // bank worst-case threshold
				//printf("Query %u: len = %u, err = %u, runsize = %u, mmatch = %u\n", j, len, err, (len - err)/(err+1), mmatch);
				//if (j == z) Ub.div = 1; // remove in loop below?
				if (j >= *QBins) { //storeUnambigWords(Word_Ix, s, len, &wix, j);
					uint32_t k = 0, w = s[k++] - 1;
					for (; k + 1 < SCOUR_N; ++k)
						w <<= 2, w |= s[k]-1;
					for (; k < len; ++k) {
						w <<= 2, w |= s[k]-1;
						uint32_t t = (w << SCOUR_R) >> SCOUR_R;
						Word_Ix[wix++] = (Split){t,j};
					}
				}
				else for (uint32_t k = 0; k + SCOUR_N <= len; ++k) {
					if (Z) { // N-skip case
						uint32_t p=k, e = k+SCOUR_N; 
						for (; p<e; ++p) if (s[p]==5) break;
						if (p<e) {k=p; continue;}
					}
					storeAmbigWords(Word_Ix, s + k, &wix, j, 0, 0);
				}
			}
			qsort(Word_Ix, wix, sizeof(*Word_Ix), WrdCmp3);
			postScour(Hash, Forest, Word_Ix, wix, Cache, &cix);
			for (uint32_t i = 0; i < cix; ++i) {
				uint16_t *p = Hash + Cache[i];
				//printf("Clump %u, count %u\n",Cache[i],*p);
				for (int z = 0; z < VECSZ; ++z) { 
					uint32_t refIX = Cache[i] * VECSZ + z;
					uint32_t ref = RefIxSrt[refIX];
					//printf("--->%s [%u]\n",RefHead[ref],[insert location here]);
				}
				if (*p > min_mmatch) Refs[nref++] = (Split){Cache[i],*p};
				*p = 0;
			}
			
			iSort(Refs, nref);
			
			RefPtr = Refs;
			oldNref = nref;
			
			// Run the loop 2x: 1) good R's, 2) bad R's
			for (int x = 0; x < 2; ++x) {
				for (uint32_t i = 0; i < nref; ++i) {
					uint32_t ri = RefPtr ? RefPtr[i].v : BadList[i];
					uint32_t rlen = ClumpLen[ri] + 1;
					DualCoil *pclump = ProfClump[ri];
					if (usedb) { // edx rather than raw fasta? or generate edx on the fly with fasta?
						DualCoil *RefSlide = RefClump[ri];
						for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
							__m128i org = _mm_stream_load_si128((void*)(RefSlide++));
							__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
							__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
							_mm_store_si128((void*)(rclump+w),ex1);
							_mm_store_si128((void*)(rclump+w+1),ex2);
						}
					} else rclump = RefClump[ri];
					HiBound[1] = rlen;
					//*LoBound = -1, LoBound[1] = 1; // FDR
					
					uint32_t stack_ix = 0; //thisMax = 0;
					*StackX = z, *StackE = -1; // initialize stack (worst case)
					//UniBins[z].div = 1; // force first divergence to 1 (done above)

					int8_t fp_rediv = 0;
					for (uint32_t j = z; j < bound; ++j) { 
						UniBin *Ub = UniBins + j; ShrBin *Sb = ShrBins + Ub->six;
						uint32_t thisDiv = Ub->div, Emac = Sb->ed, len = Sb->len;
						//uint32_t Emac = maxED;
						//printf("Doing q #%u [%s]: div=%u, Emac = %u\n",j,QHead[Offset[Ub->six]],thisDiv,Emac);

						// handle early k-term here
						uint32_t kload = Emac*SCOUR_N + SCOUR_N, 
							mmatch = kload < len ? len - kload : 1;
						//uint32_t mmatch = len - Emac * SCOUR_N - SCOUR_N;
						if (Emac == (uint16_t)-1 || (!x && Refs[i].i <= mmatch)) { // Emac to -1 in ANY mode to mark spent query
							fp_rediv = 1; // puts("-->Skipping (mmatch)");
							continue;
						}
						
						// handle fingerprinting here
						if (DO_FP) {
							if (Fp.N[j] - FP_intersect(Centroids[ri],Fp.P[j]) > Emac) {
								fp_rediv = 1;
								continue;
							}
							// do ALL individual FPs too.
							uint32_t trigger = 0;
							for (uint32_t k = ri*VECSZ; k < ri*VECSZ + VECSZ; ++k) {
								uint32_t t = FP_intersect(RP[k],Fp.P[j]);
								if (t > trigger) trigger = t;
							}
							if (Fp.N[j] - trigger > Emac) {fp_rediv = 1; continue;}
						}
						
						if (fp_rediv) {
							fp_rediv = 0;
							if (Emac <= StackE[stack_ix]) {
								thisDiv = 1; if (Ub->div > 1 && stack_ix) {
									uint32_t sai = StackX[stack_ix], lenD;
									lenD = MIN(ShrBins[UniBins[sai].six].len, len) - 1;
									register char *thisQ = Ub->Seq, *prevQ = UniBins[sai].Seq;
									for (uint32_t w = 0; w < lenD && thisDiv < maxDiv && 
										thisQ[w] == prevQ[w]; ++w) ++thisDiv; 
								}
							}
							//printf("--> fprediv: %u\n",thisDiv);
						}

						// if error tol exceed cur stack, increment stack, set this ptr
						if (Emac > StackE[stack_ix]) {
							while (Emac > StackE[--stack_ix]); // pop until errors <= Stack errors
							thisDiv = 1; if (Ub->div > 1 && stack_ix) {
								uint32_t sai = StackX[stack_ix], lenD;
								lenD = MIN(ShrBins[UniBins[sai].six].len, len) - 1;
								register char *thisQ = Ub->Seq, *prevQ = UniBins[sai].Seq;
								for (uint32_t w = 0; w < lenD && thisDiv < maxDiv && 
									thisQ[w] == prevQ[w]; ++w) ++thisDiv; 
							}
							//printf("--> rediv: %u\n",thisDiv);
						}
						//printf("-->stack before: ix=%u, err=%u, refix=%u ==> ",stack_ix, StackE[stack_ix],StackX[stack_ix]);
						stack_ix += Emac < StackE[stack_ix];
						StackX[stack_ix] = j;
						StackE[stack_ix] = Emac;
						//printf("ix=%u, err=%u, refix=%u\n",stack_ix, StackE[stack_ix],StackX[stack_ix]);
						
						DualCoil mins; //mins.v = _mm_set1_epi8(-1);
						uint32_t min;
						//if (Xalpha) min = aded_xalpha(rclump,UniqQSeq[j],rlen,UniqQLen[ai], rdim, 
						//	Matrices, pclump, Emac,thisDiv,LoBound,HiBound,&mins); 
						//else 
						//min = aded_mat16(rclump,Ub->Seq,rlen,len, rdim, // TODO: dip-n-dot, length restriction, 
						//	Matrices, pclump, Emac,thisDiv,LoBound,HiBound,&mins); 
						min = aded_mat16L(rclump,Ub->Seq,rlen,len, rdim, minlen,// TODO: dip-n-dot
							Matrices, pclump, Emac,thisDiv,LoBound,HiBound,&mins); 
						//printf("--> min = %u\n",min);
						if (min <= Sb->ed) { // now we get serious.
							uint32_t ai = Ub->six + (Ub->rc ? numUniqQ : 0);
							if (RUNMODE != FORAGE && RUNMODE != ANY) {
								/*if (min < Sb->ed)*/ Sb->ed = min;
								ResultPod *prv;
								if (Pods[ai] && min < Pods[ai]->mismatches) // purge
									do prv = Pods[ai]->next, free(Pods[ai]), Pods[ai] = prv; while(Pods[ai]);
							} else min = Emac; // all valid refs can be explored
							MetaPack MPK;
							reScoreM_mat16(rclump,Ub->Seq,rlen,len, rdim, ScoresEX, ShiftsEX,
								ShiftsBX, min, pclump,&MPK);
							for (int z = 0; z < VECSZ; ++z) { 
								if (mins.u8[z] > min  || ri*VECSZ + z >= totR) continue; // skip non-mins
								ResultPod *tmp = malloc(sizeof(*tmp));
								tmp->next = Pods[ai]; 
								tmp->mismatches = mins.u8[z];
								tmp->score = MPK.score[z]; //-1.f; // placeholder
								tmp->refIx = ri * VECSZ + z;
								tmp->finalPos = MPK.finalPos[z];
								tmp->numGapR = MPK.numGapR[z];
								tmp->numGapQ = MPK.numGapQ[z];
								tmp->rc = Ub->rc;
								if (RUNMODE == ANY) {
									ResultPod *rp = tmp;
									uint32_t rix = RefIxSrt[rp->refIx],
									numGap = rp->numGapR + rp->numGapQ,
									numMis = rp->mismatches - numGap,
									qlen = Sb->len,
									alLen = qlen + numGap,
									mOff = RefStart ? RefStart[rix] : 0, tmp, 
									stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
									edIxR = rp->finalPos + mOff;
									if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
									/* if (taxa_parsed) {
										char *tt = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
										if (taxasuppress) {
											uint32_t lm, s = 0; 
											strcpy(Taxon, tt);
											for (lm = 0; TAXLEVELS[lm] < rp->score; ++lm);
											if (!lm) FinalTaxon = NULLTAX;
											else for (int x = 0; Taxon[x]; ++x) {
												if (Taxon[x]==';' && ++s == lm) {
													Taxon[x] = 0; break;
												}
												FinalTaxon = Taxon;
											}
										} else FinalTaxon = tt; // i
										for (uint32_t j = Offset[Ub->six]; j < Offset[Ub->six+1]; ++j) PRINT_MATCH_TAX()
									}
									else  */
									if (Sb->ed != -1) for (uint32_t j = Offset[Ub->six]; j < Offset[Ub->six+1]; ++j) 
										fprintf(output,"%s\t%s\t%f\t%u\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\n", 
										QHead[j], RefHead[rix], rp->score * 100, 
										alLen, numMis, numGap, 1, qlen, stIxR, edIxR, 
										rp->mismatches,j > Offset[Ub->six]);
									Sb->ed = -1; 
									break;
								}
								Pods[ai] = tmp;
							}
						}
					}
					// tid = omp_get_thread_num();
					// if (!tid) printf("\rSearch Progress: [%3.2f%%]",100.0 * (double)totDone / newUniqQ);
				} // end standard
				// switch targets from good to bad
				RefPtr = 0;
				nref = szBL;
			}
			if (!quiet && !tid) printf("\rSearch Progress: [%3.2f%%]",100.0 * (double)totDone * totQMult);
			
			#pragma omp atomic
			totDone += (bound - z);
		}
		
		//ThreadPods[omp_get_thread_num()] = Pods;
		free(HiBound); free(LoBound); free(StackE); free(StackX); free(Matrices);
		free(ScoresEX); free(ShiftsEX); free(ShiftsBX); //free(Pods);
		free(Hash); 
		free(Refs);
	}
	
	// loop and consolidate thread pods in parallel
	//BasePod = *ThreadPods;
	// Fold rc into main bin
	if (newUniqQ > numUniqQ) {
		#pragma omp parallel for schedule(dynamic,1)
		for (uint32_t j = numUniqQ; j < newUniqQ; ++j) {
			if (!BasePod[j]) continue;
			// append this pod to end of first (if extant)
			if (!BasePod[j-numUniqQ]) BasePod[j-numUniqQ] = BasePod[j];
			else {
				ResultPod *t = BasePod[j-numUniqQ];
				while (t->next) t = t->next;
				t->next = BasePod[j];
			}
		}
		BasePod = realloc(BasePod, numUniqQ * sizeof(*BasePod));
	}
	
	if (!quiet) printf("\rSearch Progress: [%3.2f%%]\n",100.0 * (double)totDone * totQMult); // not necessarily 100.0%!
	free(*Forest); free(Forest); free(BadList); // NEW [mem]
} // end ACCEL

	// do remaining queries vs all refs in traditional aligner. 
	uint32_t firstQ = DO_ACCEL ? QBins[1] : 0;
	totDone = 0;
	if (firstQ == newUniqQ || (DO_ACCEL && skipAmbig)) goto EOA;
	else if (DO_ACCEL) UniBins[firstQ].div = 1; //UniqDiv[Umap[firstQ]] = 1;
	if (DO_HEUR) DO_HEUR = QDat.minLenQ / 16U + 1u; // TODO: implement in main body?
	printf("Searching best paths through %u unique queries...\n", newUniqQ-firstQ);
	#pragma omp parallel
	{
		DualCoil *Matrices = calloc((cacheSz+2)*rdim,sizeof(*Matrices)),
				 *ScoresEX = malloc((2*rdim)*sizeof(*ScoresEX)),
				 *ShiftsEX = malloc((2*rdim)*sizeof(*ShiftsEX)),
				 *ShiftsBX = malloc((2*rdim)*sizeof(*ShiftsEX));
		uint32_t *HiBound = calloc(qdim+1,sizeof(*HiBound)),
				 *LoBound = calloc(qdim+1,sizeof(*LoBound)),
				 *StackE = malloc(qdim*sizeof(*StackE)), // num errors
				 *StackX = malloc(qdim*sizeof(*StackX)); // query ix in UniqQSeq
		for (int j = 0; j < (cacheSz+2); ++j) Matrices[j*rdim].v = _mm_set1_epi8(MIN(j*GAP,255));
		*LoBound = -1, LoBound[1] = 1;
		
		ResultPod **Pods = calloc(numUniqQ,sizeof(*Pods));
		if (!Pods) {fputs("OOM:Pods\n",stderr); exit(3);}
		DualCoil *rclump = malloc((2+rdim)*sizeof(*rclump));
			
		#pragma omp for schedule(dynamic,1) 
		for (uint32_t i = 0; i < numRclumps; ++i) {
			//*LoBound = -1, LoBound[1] = 1; // FDR
			uint32_t ri = RefOrder ? RefOrder[i].i : i;
			uint32_t rlen = ClumpLen[ri] + 1;
			DualCoil *pclump = ProfClump[ri];
			if (usedb) { // edx rather than raw fasta? or generate edx on the fly with fasta?
				DualCoil *RefSlide = RefClump[ri];
				for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
					__m128i org = _mm_stream_load_si128((void*)(RefSlide++));
					__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
					__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
					_mm_store_si128((void*)(rclump+w),ex1);
					_mm_store_si128((void*)(rclump+w+1),ex2);
				}
			} else rclump = RefClump[ri];
			HiBound[1] = rlen;
			
			*StackX = /*Umap ? Umap[firstQ] : */firstQ, *StackE = -1; // initialize stack (worst case)
			uint32_t stack_ix = 0; //thisMax = 0;
			
			int8_t fp_rediv = 0;
			for (uint32_t j = firstQ; j < newUniqQ; ++j) { 
				// pre-stack variables
				UniBin *Ub = UniBins + j;
				uint32_t qi = /*Umap ? Umap[j] :*/ j, ai = Ub->six; //qi >= numUniqQ ? RCMap[qi - numUniqQ] : qi;
				ShrBin *Sb = ShrBins + ai;
				uint32_t thisDiv = Ub->div, len = Sb->len; //UniqDiv[qi]; // or convert Div to single cond stor
				
				uint32_t Emac = Sb->ed; //UniqQed[ai];
				if (Emac == (uint16_t)-1) {
					fp_rediv = 1;
					continue;
				}
				// handle fingerprinting here
				if (DO_FP) {
					if (Fp.N[qi] - FP_intersect(Centroids[ri],Fp.P[qi]) > Emac) {
						// need to trigger recalculation of divergence before next alignment.
						fp_rediv = 1;
						continue;
					}
					// do ALL individual FPs too.
					uint32_t trigger = 0;
					for (uint32_t k = ri*VECSZ; k < ri*VECSZ + VECSZ; ++k) {
						uint32_t t = FP_intersect(RP[k],Fp.P[qi]);
						if (t > trigger) trigger = t;
					}
					if (Fp.N[qi] - trigger > Emac) {fp_rediv = 1; continue;}
				}
				if (fp_rediv) {
					fp_rediv = 0;
					if (Emac <= StackE[stack_ix]) {
						thisDiv = 1; if (stack_ix) {
							uint32_t sai = StackX[stack_ix], lenD;
							lenD = MIN(len,ShrBins[UniBins[sai].six].len) - 1;
							register char *thisQ = Ub->Seq, *prevQ = UniBins[sai].Seq; 
							for (uint32_t w = 0; w < lenD && thisDiv < maxDiv && 
								thisQ[w] == prevQ[w]; ++w) ++thisDiv; 
						}
					}
				}

				// if error tol exceed cur stack, increment stack, set this ptr
				if (Emac > StackE[stack_ix]) {
					while (Emac > StackE[--stack_ix]); // pop until errors <= Stack errors
					thisDiv = 1; if (stack_ix) {
						uint32_t sai = StackX[stack_ix], lenD;
						lenD = MIN(len,ShrBins[UniBins[sai].six].len) - 1;
						register char *thisQ = Ub->Seq, *prevQ = UniBins[sai].Seq; 
						for (uint32_t w = 0; w < lenD && thisDiv < maxDiv && 
							thisQ[w] == prevQ[w]; ++w) ++thisDiv; 
					}
				}
				stack_ix += Emac < StackE[stack_ix];
				StackX[stack_ix] = qi;
				StackE[stack_ix] = Emac;
				
				//uint32_t Emac = maxED;
				DualCoil mins; //mins.v = _mm_set1_epi8(-1);
				uint32_t min;
				if (Xalpha) min = aded_xalpha(rclump,Ub->Seq,rlen,len, rdim, 
					Matrices, pclump, Emac,thisDiv,LoBound,HiBound,&mins); 
				else min = aded_mat16(rclump,Ub->Seq,rlen,len, rdim, 
					Matrices, pclump, Emac,thisDiv,LoBound,HiBound,&mins); 
				
				uint32_t fmin = min; // to allow a bulge
				if (min <= Sb->ed) { // now we get serious.
					//uint32_t ai = Ub->six;
					if (RUNMODE != FORAGE && RUNMODE != ANY) {
						//fmin -= (min? RUNMODE == COMMUNIST:0); // TODO: FINISH
						if (fmin < Sb->ed) Sb->ed = fmin;
						ResultPod *prv;
						if (Pods[ai] && fmin < Pods[ai]->mismatches) // purge
							do prv = Pods[ai]->next, free(Pods[ai]), Pods[ai] = prv; while(Pods[ai]);
					} else min = Emac; // all valid refs can be explored
					MetaPack MPK;
					if (Xalpha) reScoreM_xalpha(rclump,Ub->Seq,rlen,len, rdim, ScoresEX, ShiftsEX,
						ShiftsBX, min, pclump,&MPK);
					else reScoreM_mat16(rclump,Ub->Seq,rlen,len, rdim, ScoresEX, ShiftsEX,
						ShiftsBX, min, pclump,&MPK);
					for (int z = 0; z < VECSZ; ++z) { 
						if (mins.u8[z] > min || ri*VECSZ + z >= totR) continue; // skip non-mins
						// TODO: force quit if ri*VECSZ + z >= totR!
						ResultPod *tmp = malloc(sizeof(*tmp));
						tmp->next = Pods[ai]; 
						
						tmp->mismatches = mins.u8[z];
						tmp->score = MPK.score[z]; //-1.f; // placeholder
						tmp->refIx = ri * VECSZ + z; 

						tmp->finalPos = MPK.finalPos[z];
						tmp->numGapR = MPK.numGapR[z];
						tmp->numGapQ = MPK.numGapQ[z];
						tmp->rc = Ub->rc;
						if (RUNMODE == ANY) {
							ResultPod *rp = tmp;
							uint32_t rix = RefIxSrt[rp->refIx],
							numGap = rp->numGapR + rp->numGapQ,
							numMis = rp->mismatches - numGap,
							qlen = Sb->len,
							alLen = qlen + numGap,
							mOff = RefStart ? RefStart[rix] : 0, tmp, 
							stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
							edIxR = rp->finalPos + mOff;
							if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
							if (Sb->ed != -1) for (uint32_t j = Offset[Ub->six]; j < Offset[Ub->six+1]; ++j) 
								fprintf(output,"%s\t%s\t%f\t%u\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\n", 
								QHead[j], RefHead[rix], rp->score * 100, 
								alLen, numMis, numGap, 1, qlen, stIxR, edIxR, 
								rp->mismatches,j > Offset[Ub->six]);
							Sb->ed = -1; 
							break;
						}
						Pods[ai] = tmp;
					}
				}
			}
			#pragma omp atomic
			++totDone;
			tid = omp_get_thread_num();
			if (!quiet && !tid) printf("\rSearch Progress: [%3.2f%%]",100.0 * (double)totDone / numRclumps);
		}
		ThreadPods[omp_get_thread_num()] = Pods;
		free(HiBound); free(LoBound); free(StackE); free(StackX); free(Matrices);
		free(ScoresEX); free(ShiftsEX); free(ShiftsBX); //free(Pods);
	}
	
	// loop and consolidate thread pods in parallel
	if (!BasePod) {
		BasePod = *ThreadPods;
		#pragma omp parallel for schedule(dynamic,1)
		for (uint32_t j = 0; j < numUniqQ; ++j) {
			if (!BasePod[j]) continue;
			ResultPod *prv;
			if (BasePod[j]->mismatches > ShrBins[j].ed) 
				do prv = BasePod[j]->next, /*free(BasePod[j]),*/ BasePod[j] = prv; 
				while(BasePod[j]);
		}
	}
	
	for (uint32_t i = BasePod==*ThreadPods; i < THREADS; ++i) {
		#pragma omp parallel for schedule(dynamic,1)
		for (uint32_t j = 0; j < numUniqQ; ++j) {
			if (!ThreadPods[i][j]) continue;
			ResultPod *prv, *init;
			if (ThreadPods[i][j]->mismatches > ShrBins[j].ed) {
				do prv = ThreadPods[i][j]->next, /*free(ThreadPods[i][j]),*/ ThreadPods[i][j] = prv; 
				while(ThreadPods[i][j]);
			} else {
				init = prv = ThreadPods[i][j];
				while (prv->next) prv = prv->next;
				prv->next = BasePod[j];
				BasePod[j] = init;
			}
		}
		free(ThreadPods[i]);
	}
	if (!quiet) printf("\rSearch Progress: [100.00%%]\n");
	EOA:NULL;
	
	if (RUNMODE == ANY) return;
	
	printf("Search complete. Consolidating results...\n");
	//if (totSkipped) printf("NumSkipped = %llu (%f)\n",
	//	totSkipped,(double)totSkipped/((double)numUniqQ*numRclumps));
	if (!RefIxSrt && RefDedupIx) {
		printf("Constructing RefIxSrt from RefDedupIx...[totR %u, orig %u]\n",totR,RefDat.origTotR);
		RefIxSrt = malloc(totR * sizeof(*RefIxSrt));
		if (!RefIxSrt) {fputs("OOM:[DA]RefIxSrt\n",stderr); exit(3);}
		for (uint32_t i = 0; i < totR; ++i) RefIxSrt[i] = TmpRIX[RefDedupIx[i]];
	}
	else if (!RefIxSrt && !RefDedupIx) RefIxSrt = TmpRIX;
	
	uint32_t maxIX = 0;
	#pragma omp simd reduction(max:maxIX)
	for (uint32_t i = 0; i < totR; ++i) 
		maxIX = RefIxSrt[i] > maxIX ? RefIxSrt[i] : maxIX;
	//printf("Max IX for allocations = %u [totR = %u, origTotR = %u]\n",maxIX,totR,RefDat.origTotR);

	free(Centroids); free(FpR.initP); free(ProfClump); free(RefClump);
	free(SeqDumpRC); // free a part of the query memory
	// Remove duplicate alignments
	uint32_t *RefMap = RefDat.RefMap;
	if (!RefMap) {
		puts("Constructing RefMap...");
		RefMap = malloc((totR > RefDat.origTotR ? totR : RefDat.origTotR)*sizeof(*RefMap));
		if (!RefMap) {fputs("OOM:RefMap\n",stderr); exit(3);}
		RefDat.numRefHeads = RefDat.totR;
		for (uint64_t i = 0; i < totR; ++i) RefMap[i] = i;
	}
	
#define PRINT_MATCH() \
	fprintf(output,"%s\t%s\t%f\t%u\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\n", \
		QHead[j], RefHead[rix], rp->score * 100, \
		alLen, numMis, numGap, 1, qlen, stIxR, edIxR, \
		rp->mismatches,j > Offset[i]);
#define PRINT_MATCH_TAX() \
	fprintf(output,"%s\t%s\t%f\t%u\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\t%s\n", \
		QHead[j], RefHead[rix], rp->score * 100, \
		alLen, numMis, numGap, 1, qlen, stIxR, edIxR, \
		rp->mismatches,j > Offset[i],FinalTaxon);
#define DUPE_HUNT(ddix, R_C, S_C, rp, qlen, rix, mOff, stIxR, mapped, ql2, JMP) \
	mOff = RefStart ? RefStart[rix] : 0; \
	stIxR = rp->rc ? rp->finalPos + mOff : rp->finalPos - qlen + rp->numGapR + mOff; \
	mapped = RefMap[rix]; \
	for (uint32_t d = 0; d < ddix; ++d) \
		if (R_C[d] == mapped && S_C[d] + ql2 > stIxR \
			&& S_C[d] < stIxR + ql2) goto JMP; \
	R_C[ddix] = mapped, S_C[ddix++] = stIxR; 
#define DUPE_HUNT_S(ddix, R_C, S_C, X_C, rp, qlen, rix, mOff, stIxR, mapped, ql2, JMP) \
	mOff = RefStart ? RefStart[rix] : 0; \
	stIxR = rp->rc ? rp->finalPos + mOff : rp->finalPos - qlen + rp->numGapR + mOff; \
	mapped = RefMap[rix]; \
	for (uint32_t d = 0; d < ddix; ++d) \
		if (R_C[d] == mapped && S_C[d] + ql2 > stIxR && S_C[d] < stIxR + ql2) { \
			if (X_C[d] > rp->mismatches) X_C[d] = rp->mismatches; \
			goto JMP; \
		} \
	R_C[ddix] = mapped, X_C[ddix] = rp->mismatches, S_C[ddix++] = stIxR; 

	if (RUNMODE == ALLPATHS) {  // all on best ED path
		uint32_t *RefCache = malloc((maxIX+1)*sizeof(*RefCache)); // N L
		uint32_t *StCache = malloc((maxIX+1)*sizeof(*StCache)); // N L
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i];
			if (!rp) continue;
			uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1; // N L
			ResultPod *best = rp;
			while (rp = rp->next) if (rp->mismatches < best->mismatches) best = rp;
			rp = best;
			uint32_t bm = best->mismatches;
			if (rp->score) while (rp) {
				if (rp->mismatches == bm) { //rp->score && 
					uint32_t rix, mOff, stIxR, mapped;
					if (RefDedupIx) // handle dupes in refs
						for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
							rix = TmpRIX[k];
							DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,ALLPATHS_D_END);
							uint32_t numGap = rp->numGapR + rp->numGapQ, 
								numMis = rp->mismatches - numGap,
								alLen = qlen + numGap, 
								edIxR = rp->rc ? rp->finalPos - qlen + rp->numGapR + mOff : rp->finalPos + mOff;
							if (taxa_parsed) {
								char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
								for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
							}
							else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
							ALLPATHS_D_END:NULL;
						}
					else {
						rix = RefIxSrt[rp->refIx];
						DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,ALLPATHS_U_END);
						uint32_t numGap = rp->numGapR + rp->numGapQ, 
							numMis = rp->mismatches - numGap,
							alLen = qlen + numGap, 
							edIxR = rp->rc ? rp->finalPos - qlen + rp->numGapR + mOff : rp->finalPos + mOff;
						if (taxa_parsed) {
							char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
							for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
						}
						else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
						ALLPATHS_U_END:NULL;
					}
				}
				rp = rp->next;
			}
		}
	}
	/*else if (RUNMODE == FORAGE) { // DEBUG VERSION
		printf("Considering %u queries.\n",numUniqQ);

		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i];
			printf("Query %u:\n",i);
			while (rp) {
				// display some info about this reference: where and what it is, then dump some accelerator stats (check to see if accel is freed earlier?)
				
				//if (rp->score) {
					uint32_t rix = RefIxSrt[rp->refIx],
					numGap = rp->numGapR + rp->numGapQ,
					numMis = rp->mismatches - numGap,
					qlen = ShrBins[i].len,
					alLen = qlen + numGap,
					mOff = RefStart ? RefStart[rix] : 0, tmp, 
					stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
					edIxR = rp->finalPos + mOff;
					if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
					printf("Reference %u: \n", rix);
					if (RefDedupIx) { // handle dupes in refs
						printf("-->Range begin: %u, end: %u\n", RefDedupIx[rp->refIx], RefDedupIx[rp->refIx+1]);
						for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
							rix = TmpRIX[k];
							mOff = RefStart ? RefStart[rix] : 0;
							stIxR = rp->finalPos - qlen + rp->numGapR + mOff;
							edIxR = rp->finalPos + mOff;
							if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
							if (taxa_parsed) {
								char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
								for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
							}
							else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
							printf("  [%u] %s\n",k,RefHead[rix]);
						}
					}
					else if (taxa_parsed) {
						char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
						for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
					}
					else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
				//}
				rp = rp->next;
			}
		}
	}*/
	else if (RUNMODE == FORAGE) { // all valid alignments
		printf("Considering %u queries.\n",numUniqQ);
		uint32_t *RefCache = malloc((maxIX+1)*sizeof(*RefCache)); // N L
		uint32_t *StCache = malloc((maxIX+1)*sizeof(*StCache)); // N L
		uint8_t *ScCache = malloc((maxIX+1)*sizeof(*ScCache)); // N L
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i];
			uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1; // N L
			while (rp) {
					uint32_t rix, mOff, stIxR, mapped;
					if (RefDedupIx) // handle dupes in refs
						for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
							rix = TmpRIX[k];
							DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,FORAGE_D_END);
							// TODO: inline the dedupe with added functionality to take the _higher scoring_ of the dupes
							uint32_t numGap = rp->numGapR + rp->numGapQ, 
								numMis = rp->mismatches - numGap,
								alLen = qlen + numGap, 
								edIxR = rp->rc ? rp->finalPos - qlen + rp->numGapR + mOff : rp->finalPos + mOff;
							if (taxa_parsed) {
								char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
								for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
							}
							else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
							FORAGE_D_END:NULL; // N L
						}
					else {
						rix = RefIxSrt[rp->refIx];
						DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,FORAGE_U_END);
						uint32_t numGap = rp->numGapR + rp->numGapQ, 
							numMis = rp->mismatches - numGap,
							alLen = qlen + numGap, 
							edIxR = rp->rc ? rp->finalPos - qlen + rp->numGapR + mOff : rp->finalPos + mOff;
						if (taxa_parsed) { 
							char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
							for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
						}
						else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
						FORAGE_U_END:NULL; // N L
					}
				rp = rp->next;
			}
		}
	}
	

	else if (RUNMODE==CAPITALIST) {
		uint32_t numBins = maxIX+1;
		size_t *RefCounts = calloc(numBins,sizeof(*RefCounts)), tot = 0;
		if (!RefCounts) {fputs("OOM:RefCounts\n",stderr); exit(3);}
		uint32_t *RefCache = malloc(numBins*sizeof(*RefCache)); // N L
		uint32_t *StCache = malloc(numBins*sizeof(*StCache)); // N L
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i], *best = rp;
			if (!rp) continue;
			// Phase 1: Picking the best fund. (Set 'best' to [first] min scoring ref)
			while (rp = rp->next) if (rp->mismatches < best->mismatches) best = rp; //rp->score && 
			
			// Phase 2: Actuarial research. (Tally how much of each ref appears overall.)
			rp = best; 
			BasePod[i] = rp; // hedging derivative futures (start at the first 'best')
			
			uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1; 
			do if (rp->mismatches == best->mismatches) {
				uint32_t rix, mOff, stIxR, mapped;
				if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
					rix = TmpRIX[k];
					DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_D_END);
					++RefCounts[mapped], ++tot;
					CAPITALIST_D_END:NULL;
				}
				else {
					rix = RefIxSrt[rp->refIx];
					DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_U_END);
					++RefCounts[mapped], ++tot;
					CAPITALIST_U_END:NULL;
				}
			} while (rp = rp->next);
		}
		printf("CAPITALIST: Processed %llu investments\n", tot);
		#pragma omp parallel num_threads(MIN(THREADS,8))
		{
			int tid = omp_get_thread_num();
			uint32_t *RefCache_t, *StCache_t;
			if (!tid) RefCache_t = RefCache, StCache_t = StCache;
			else RefCache_t = malloc(numBins*sizeof(*RefCache_t)),
				StCache_t = malloc(numBins*sizeof(*StCache_t));
			if (!RefCache_t || !StCache_t) {fputs("OOM:StCache_t",stderr); exit(3);}

			char **Taxa = 0, *Taxon = 0, *FinalTaxon = 0; uint32_t *Divergence = 0; 
			if (taxa_parsed) 
				Taxa = malloc(numBins*sizeof(*Taxa)), 
				Taxon = malloc(1000000),
				Divergence = malloc(numBins*sizeof(*Divergence)), *Divergence = 0;
			
			// pass 3: deciding which investments to bank
			#pragma omp for
			for (uint32_t i = 0; i < numUniqQ; ++i) {
				if (!BasePod[i]) continue; 
				ResultPod *rp = BasePod[i];
				ResultPod *best = rp;
				uint32_t tix = 0, bestmap, bestrix;
				float best_score = -1.f;
				uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1,
					rix, mOff, stIxR, mapped;

				do {
					if (rp->mismatches < best->mismatches) continue;
					if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
						rix = TmpRIX[k];
						DUPE_HUNT(ddix,RefCache_t,StCache_t,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_D2_END);
						if (taxa_parsed) 
							Taxa[tix++] = taxa_lookup(RefHead[rix], taxa_parsed-1, Taxonomy),
							best_score = rp->score > best_score ? rp->score : best_score;
						if ( best == rp || (RefCounts[mapped] >  RefCounts[bestmap]) || 
							 (RefCounts[mapped] == RefCounts[bestmap] && mapped < bestmap) )
								best = rp, bestmap = mapped, bestrix = rix;
						CAPITALIST_D2_END:NULL;
					}
					else {
						rix = RefIxSrt[rp->refIx];
						DUPE_HUNT(ddix,RefCache_t,StCache_t,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_U2_END);
						if (taxa_parsed) 
							Taxa[tix++] = taxa_lookup(RefHead[rix], taxa_parsed-1, Taxonomy),
							best_score = rp->score > best_score ? rp->score : best_score;
						if ( best == rp || (RefCounts[mapped] >  RefCounts[bestmap]) || 
							 (RefCounts[mapped] == RefCounts[bestmap] && mapped < bestmap) )
								best = rp, bestmap = mapped, bestrix = rix;
						CAPITALIST_U2_END:NULL;
					}
				} while (rp = rp->next);

				// taxonomy interpolation pass (optional)
				if (taxa_parsed) {
					uint32_t lv = -1;
					if (tix==1) {FinalTaxon = strcpy(Taxon,*Taxa); goto END_CAP_TAX;}
					int byStr(const void *A, const void *B) {return strcmp(*(char**)A, *(char**)B);}
					qsort(Taxa, tix, sizeof(*Taxa), byStr);

					uint32_t maxDiv = 0;
					for (int z = 1, x; z < tix; ++z) {
						Divergence[z] = 0;
						for (x = 0; Taxa[z-1][x] && Taxa[z-1][x]==Taxa[z][x]; ++x) 
							Divergence[z] += Taxa[z][x] == ';';
						Divergence[z] += !Taxa[z-1][x];
						if (Divergence[z] > maxDiv) maxDiv = Divergence[z];
					}
					if (!maxDiv) {Taxon[0] = 0; FinalTaxon = Taxon; goto END_CAP_TAX;}

					// Ascend tree based on divergences
					uint32_t cutoff = tix - tix/TAXACUT; 
					//printf("    Query: %s, cutoff = %u, tix = %u\n", QHead[NewIX[Offset[i]]], cutoff, tix);
					uint32_t st = 0, ed = tix;
					for (lv = 1; lv <= maxDiv; ++lv) {
						uint32_t accum = 1;
						for (uint32_t z = st+1; z < ed; ++z) {
							if (Divergence[z] >= lv) ++accum;
							else if (accum >= cutoff) {ed = z; break;}
							else accum = 1, st = z; //, printf("reset z=%u, lv %u: %u < %u\n",z,lv,accum,cutoff);
						}
						if (accum < cutoff) break;
						cutoff = accum - accum/TAXACUT;
					}
					//for (int b = 0; b < tix; ++b) printf("%d: [%u] %s\n", b, Divergence[b], Taxa[b]); 
					
					// copy result up until lv-1 semicolon into Taxon, set FinalTaxon = Taxon;
					uint32_t s = 0;
					if (ed) --ed; --lv;
					for (st = 0; Taxa[ed][st] && (s += Taxa[ed][st] == ';') < lv; ++st) 
						Taxon[st] = Taxa[ed][st];
					Taxon[st] = 0;
					FinalTaxon = Taxon;
					//printf("--> [ed=%u, lv=%u, cut=%u] %s\n",ed, lv, cutoff, FinalTaxon);
					END_CAP_TAX:
					if (taxasuppress) {
						uint32_t lm, s = 0; 
						for (lm = 0; lm < lv && TAXLEVELS[lm] < best_score; ++lm);
						if (!lm) FinalTaxon = NULLTAX;
						else if (lm < lv) for (int x = 0; FinalTaxon[x]; ++x) 
							if (FinalTaxon[x]==';' && ++s == lm) {
								FinalTaxon[x] = 0; break;
							}
						//printf("--> Best score = %f, lim= %d/%d, newtax = %s\n",best_score,lm,lv, FinalTaxon);
					}
				}
				
				// pass 4: printing dividend reports
				rp = best; // recovering hedge fund
				rix = bestrix;
				uint32_t numGap = rp->numGapR + rp->numGapQ, 
					numMis = rp->mismatches - numGap, alLen = qlen + numGap;
				mOff = RefStart ? RefStart[rix] : 0;
				stIxR = rp->finalPos - qlen + rp->numGapR + mOff;
				uint32_t tmp, edIxR = rp->finalPos + mOff;
				if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
				
				#pragma omp critical
				if (taxa_parsed) for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
				else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
			}
		}
		//free(RefCounts), free(Taxa), free(Taxon), free(Divergence);
	}
	else if (RUNMODE==COMMUNIST) {
		typedef struct AlPod AlPod;
		struct AlPod {
			uint32_t map;
			float ID;
			uint16_t len;
			uint8_t err, gaps;
			//....
			uint32_t start, end;
			// ....
			AlPod *next;
		};
		typedef struct {
			char *Qname;
			uint32_t qlen;
			AlPod *wheel;
		} QPod;
		puts("COMMUNIST reporting mode invoked.");
		/* The idea behind this is to "smooth out" coverage of bugs hit. 
		To do this, we need to keep a tally of where queries land on their references.
		First, allocate unique best-hits. Then using these as a guide, disambiguate 
		queries based on where they would contribute most substantively to coverage.
		This objective is tricky. Want to cover barren ground, but also bias toward 
		doing so in more-covered genomes. 

		1. allocate [bytes!] arrays equal to genome length
		  1b. determine genome length...
		2. Throw in all unique best hits
		3. For the remainder, allocate according to objective above
		  3a. put this read where maximizes proportion existing coverage of genome / [min?] coverage in region
		  3b. do this in repeat until stabilized?
		4. If communist, do second pass with all available positions [redistribute reads according to same objective, include self term as if read wasn't there]
		*/

		// Determine genome sizes
		// RefMap is the ticket. For each [actual, including dupes!] shear, get its RefMap ix 
		//   and compute max given current vs [length of this shear relative to its RefStart]
		uint32_t *Osizes = calloc(RefDat.numRefHeads,sizeof(*Osizes));
		char **Oheads = malloc((size_t)RefDat.numRefHeads*sizeof(*Oheads));
		uint32_t qix = 0;
		if (!Osizes || !Oheads) {fputs("ERR:OOM Oheads\n",stderr); exit(3);}
		
		for (uint32_t i = 0; i < totR; ++i) {
			uint32_t len = ClumpLen[i/16];
			uint32_t rix, mOff, stIxR, mapped;
			if (RefDedupIx) for (uint32_t k = RefDedupIx[i]; k < RefDedupIx[i+1]; ++k) {
				rix = TmpRIX[k];
				mOff = RefStart ? RefStart[rix] : 0; 
				mapped = RefMap[rix]; 
				//printf("Mapping %u->%u->%u->%u, len=%u, mOff = %u\n",i,k,rix,mapped,len,mOff);
				Osizes[mapped] = Osizes[mapped] > mOff + len ? Osizes[mapped] : mOff + len;
				Oheads[mapped] = RefHead[rix];
			}
			else { // this branch hasn't been tested much
				rix = RefIxSrt[i];
				mOff = RefStart ? RefStart[rix] : 0; 
				mapped = RefMap[rix]; 
				Osizes[mapped] = Osizes[mapped] > mOff + len ? Osizes[mapped] : mOff + len;
				Oheads[mapped] = RefHead[rix];
			}
		}
		
		uint32_t numBins = maxIX+1;
		uint32_t *RefCache = malloc(numBins*sizeof(*RefCache)); // N L
		uint32_t *StCache = malloc(numBins*sizeof(*StCache)); // N L
		QPod **QPods = calloc(QDat.totQ,sizeof(*QPods));
		if (!RefCache || !StCache || !QPods) {fputs("ERR:OOM QPods\n",stderr); exit(3);}
		
		// now we want to run through the unique ones and increment there for coverage
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i], *best = rp;
			if (!rp) continue;
			while (rp = rp->next) if (rp->mismatches < best->mismatches) best = rp;
			rp = best; //BasePod[i]=rp; // lossy!

			uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1; 
			// enumerate number of duplicated queries 
			uint32_t copies = Offset[i+1]-Offset[i];
			for (uint32_t j = 0; j < copies; ++j) {
				QPods[qix] = malloc(sizeof(*QPods[qix]));
				*QPods[qix] = (QPod){QHead[j],qlen,0};
				++qix;
			}
			do if (rp->mismatches == best->mismatches) {
				uint32_t rix, mOff, stIxR, mapped;
				
				
					
				if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
					rix = TmpRIX[k];
					DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,SOC_D_END);
					uint32_t numGap = rp->numGapR + rp->numGapQ, 
							numMis = rp->mismatches - numGap,
							alLen = qlen + numGap, 
							edIxR = rp->rc ? rp->finalPos - qlen + rp->numGapR + mOff : rp->finalPos + mOff;
					//++RefCounts[mapped], ++tot;
					for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) {
						QPods[qix-copies--] = (QPod){};
					}
					
					SOC_D_END:NULL;
				}
				else {
					rix = RefIxSrt[rp->refIx];
					DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,SOC_U_END);
					//++RefCounts[mapped], ++tot;
					SOC_U_END:NULL;
				}
			} while (rp = rp->next);
		}
		printf("COMMUNIST: Processed %llu reallocations\n", tot);
		

		//for (uint32_t i = 0; i < RefDat.numRefHeads; ++i) {
		//	printf("%s\t%u\n",Oheads[i],Osizes[i]);
		//}
		
		

		exit(404);

		// mostly old code follows
		uint32_t numBins = maxIX+1;
		//size_t *RefCounts = calloc(numBins,sizeof(*RefCounts)), tot = 0;
		//if (!RefCounts) {fputs("OOM:RefCounts\n",stderr); exit(3);}
		uint32_t *RefCache = malloc(numBins*sizeof(*RefCache)); // N L
		uint32_t *StCache = malloc(numBins*sizeof(*StCache)); // N L
		
		/* for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i], *best = rp;
			if (!rp) continue;
			// Phase 1: Picking the best fund. (Set 'best' to [first] min scoring ref)
			while (rp = rp->next) if (rp->mismatches < best->mismatches) best = rp; //rp->score && 
			
			// Phase 2: Actuarial research. (Tally how much of each ref appears overall.)
			rp = best; 
			//BasePod[i] = rp; // hedging derivative futures (start at the first 'best')
			
			uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1; 
			do if (rp->mismatches == best->mismatches) {
				uint32_t rix, mOff, stIxR, mapped;
				if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
					rix = TmpRIX[k];
					DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_D_END);
					++RefCounts[mapped], ++tot;
					CAPITALIST_D_END:NULL;
				}
				else {
					rix = RefIxSrt[rp->refIx];
					DUPE_HUNT(ddix,RefCache,StCache,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_D_END);
					++RefCounts[mapped], ++tot;
					CAPITALIST_U_END:NULL;
				}
			} while (rp = rp->next);
		}
		printf("CAPITALIST: Processed %llu investments\n", tot);
		#pragma omp parallel num_threads(MIN(THREADS,8))
		{
			int tid = omp_get_thread_num();
			uint32_t *RefCache_t, *StCache_t;
			if (!tid) RefCache_t = RefCache, StCache_t = StCache;
			else RefCache_t = malloc(numBins*sizeof(*RefCache_t)),
				StCache_t = malloc(numBins*sizeof(*StCache_t));
			if (!RefCache_t || !StCache_t) {fputs("OOM:StCache_t",stderr); exit(3);}

			char **Taxa = 0, *Taxon = 0, *FinalTaxon = 0; uint32_t *Divergence = 0; 
			if (taxa_parsed) 
				Taxa = malloc(numBins*sizeof(*Taxa)), 
				Taxon = malloc(1000000),
				Divergence = malloc(numBins*sizeof(*Divergence)), *Divergence = 0;
			
			// pass 3: deciding which investments to bank
			#pragma omp for
			for (uint32_t i = 0; i < numUniqQ; ++i) {
				if (!BasePod[i]) continue; 
				ResultPod *rp = BasePod[i];
				ResultPod *best = rp;
				uint32_t tix = 0, bestmap, bestrix;
				float best_score = -1.f;
				uint32_t ddix = 0, qlen = ShrBins[i].len, ql2 = qlen >> 1,
					rix, mOff, stIxR, mapped;

				do {
					if (rp->mismatches < best->mismatches) continue;
					if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
						rix = TmpRIX[k];
						DUPE_HUNT(ddix,RefCache_t,StCache_t,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_D2_END);
						if (taxa_parsed) 
							Taxa[tix++] = taxa_lookup(RefHead[rix], taxa_parsed-1, Taxonomy),
							best_score = rp->score > best_score ? rp->score : best_score;
						if ( best == rp || (RefCounts[mapped] >  RefCounts[bestmap]) || 
							 (RefCounts[mapped] == RefCounts[bestmap] && mapped < bestmap) )
								best = rp, bestmap = mapped, bestrix = rix;
						CAPITALIST_D2_END:NULL;
					}
					else {
						rix = RefIxSrt[rp->refIx];
						DUPE_HUNT(ddix,RefCache_t,StCache_t,rp,qlen,rix,mOff,stIxR,mapped,ql2,CAPITALIST_U2_END);
						if (taxa_parsed) 
							Taxa[tix++] = taxa_lookup(RefHead[rix], taxa_parsed-1, Taxonomy),
							best_score = rp->score > best_score ? rp->score : best_score;
						if ( best == rp || (RefCounts[mapped] >  RefCounts[bestmap]) || 
							 (RefCounts[mapped] == RefCounts[bestmap] && mapped < bestmap) )
								best = rp, bestmap = mapped, bestrix = rix;
						CAPITALIST_U2_END:NULL;
					}
				} while (rp = rp->next);

				// taxonomy interpolation pass (optional)
				if (taxa_parsed) {
					uint32_t lv = -1;
					if (tix==1) {FinalTaxon = strcpy(Taxon,*Taxa); goto END_CAP_TAX;}
					int byStr(const void *A, const void *B) {return strcmp(*(char**)A, *(char**)B);}
					qsort(Taxa, tix, sizeof(*Taxa), byStr);

					uint32_t maxDiv = 0;
					for (int z = 1, x; z < tix; ++z) {
						Divergence[z] = 0;
						for (x = 0; Taxa[z-1][x] && Taxa[z-1][x]==Taxa[z][x]; ++x) 
							Divergence[z] += Taxa[z][x] == ';';
						Divergence[z] += !Taxa[z-1][x];
						if (Divergence[z] > maxDiv) maxDiv = Divergence[z];
					}
					if (!maxDiv) {Taxon[0] = 0; FinalTaxon = Taxon; goto END_CAP_TAX;}

					// Ascend tree based on divergences
					uint32_t cutoff = tix - tix/TAXACUT; 
					//printf("    Query: %s, cutoff = %u, tix = %u\n", QHead[NewIX[Offset[i]]], cutoff, tix);
					uint32_t st = 0, ed = tix;
					for (lv = 1; lv <= maxDiv; ++lv) {
						uint32_t accum = 1;
						for (uint32_t z = st+1; z < ed; ++z) {
							if (Divergence[z] >= lv) ++accum;
							else if (accum >= cutoff) {ed = z; break;}
							else accum = 1, st = z; //, printf("reset z=%u, lv %u: %u < %u\n",z,lv,accum,cutoff);
						}
						if (accum < cutoff) break;
						cutoff = accum - accum/TAXACUT;
					}
					//for (int b = 0; b < tix; ++b) printf("%d: [%u] %s\n", b, Divergence[b], Taxa[b]); 
					
					// copy result up until lv-1 semicolon into Taxon, set FinalTaxon = Taxon;
					uint32_t s = 0;
					if (ed) --ed; --lv;
					for (st = 0; Taxa[ed][st] && (s += Taxa[ed][st] == ';') < lv; ++st) 
						Taxon[st] = Taxa[ed][st];
					Taxon[st] = 0;
					FinalTaxon = Taxon;
					//printf("--> [ed=%u, lv=%u, cut=%u] %s\n",ed, lv, cutoff, FinalTaxon);
					END_CAP_TAX:
					if (taxasuppress) {
						uint32_t lm, s = 0; 
						for (lm = 0; lm < lv && TAXLEVELS[lm] < best_score; ++lm);
						if (!lm) FinalTaxon = NULLTAX;
						else if (lm < lv) for (int x = 0; FinalTaxon[x]; ++x) 
							if (FinalTaxon[x]==';' && ++s == lm) {
								FinalTaxon[x] = 0; break;
							}
						//printf("--> Best score = %f, lim= %d/%d, newtax = %s\n",best_score,lm,lv, FinalTaxon);
					}
				}
				
				// pass 4: printing dividend reports
				rp = best; // recovering hedge fund
				rix = bestrix;
				uint32_t numGap = rp->numGapR + rp->numGapQ, 
					numMis = rp->mismatches - numGap, alLen = qlen + numGap;
				mOff = RefStart ? RefStart[rix] : 0;
				stIxR = rp->finalPos - qlen + rp->numGapR + mOff;
				uint32_t tmp, edIxR = rp->finalPos + mOff;
				if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
				
				#pragma omp critical
				if (taxa_parsed) for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
				else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
			}
		}*/
		//free(RefCounts), free(Taxa), free(Taxon), free(Divergence);
	}
	else if (RUNMODE == BEST) {  // find first best score on min ED path
		//float e = FLT_EPSILON; int sim;
		char Taxon[UINT16_MAX] = {0}, *FinalTaxon = 0;
	 	for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i];
			if (!rp) continue;
			ResultPod *best = rp;
			while (rp = rp->next) // a bin with a dead next is also dead, see above
				if ( (rp->mismatches < best->mismatches) ||  // rp->score && (...)
					 (rp->mismatches == best->mismatches && 
					 /* !(sim=(fabsf(rp->score - best->score) <= e)) && */ rp->score > best->score) ||
					 (rp->mismatches == best->mismatches && rp->score == best->score && /* sim &&  */
					  RefIxSrt[rp->refIx] < RefIxSrt[best->refIx]) )
						best = rp;
			rp = best;
			//if (rp->score) {
				uint32_t rix = RefIxSrt[rp->refIx],
				numGap = rp->numGapR + rp->numGapQ,
				numMis = rp->mismatches - numGap,
				qlen = ShrBins[i].len,
				alLen = qlen + numGap,
				mOff = RefStart ? RefStart[rix] : 0, tmp, 
				stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
				edIxR = rp->finalPos + mOff;
				if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
				if (taxa_parsed) {
					char *tt = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
					if (taxasuppress) {
						uint32_t lm, s = 0; 
						strcpy(Taxon, tt);
						for (lm = 0; TAXLEVELS[lm] < rp->score; ++lm);
						if (!lm) FinalTaxon = NULLTAX;
						else for (int x = 0; Taxon[x]; ++x) {
							if (Taxon[x]==';' && ++s == lm) {
								Taxon[x] = 0; break;
							}
							FinalTaxon = Taxon;
						}
					} else FinalTaxon = tt;
					for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
				}
				else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
			//}
		}
	}
}
int isRefEDB(char *ref_FN) {
	FILE *temp_handle = fopen(ref_FN, "rb");
	if (!temp_handle) {fputs("ERROR: invalid input file.\n",stderr); exit(1);}
	int control = fgetc(temp_handle);
	fclose(temp_handle);
	if (control == EOF) {fputs("ERROR: invalid input file.\n",stderr); exit(1);}
	return (uint8_t)control >> 7;
}
int main( int argc, char *argv[] ) {
	THREADS = getNumberOfCores();
	Reference_Data RefDat = (Reference_Data){0};
	Query_Data QDat = (Query_Data){0};
	int makedb = 0, doDedupe = 0; // clustradius = 0, incl_whitespace = 0; 
	char *ref_FN = 0, *query_FN = 0, *output_FN = 0, *xcel_FN = 0, *tax_FN = 0;
	DBType dbType = QUICK;
	RefDat.dbType = dbType;
	printf("This is BURST ["VER"]\n");
	if (argc < 2) PRINT_USAGE()
	for (int i = 1; i < argc; ++i) {
		if (!strcmp(argv[i],"--references") || !strcmp(argv[i],"-r")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --references requires filename argument"); exit(1);}
			ref_FN = argv[i];
		}
		else if (!strcmp(argv[i],"--queries") || !strcmp(argv[i],"-q")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --queries requires filename argument"); exit(1);}
			query_FN = argv[i];
		}
		else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --output requires filename argument"); exit(1);}
			output_FN = argv[i];
		}
		//puts("--forwardreverse (-fr): also search the reverse complement of queries"); 
		else if (!strcmp(argv[i],"--forwardreverse") || !strcmp(argv[i],"-fr")) {
			QDat.rc = 1; 
			printf(" --> Also considering the reverse complement of reads\n");
		}
		else if (!strcmp(argv[i],"--whitespace") || !strcmp(argv[i],"-w")) {
			QDat.incl_whitespace = 1; 
			printf(" --> Allowing whitespace in query name output\n");
		}
		else if (!strcmp(argv[i],"--npenalize") || !strcmp(argv[i],"-n")) {
			Z = 1; // global
			printf(" --> Setting N penalty (ref N vs query A/C/G/T)\n");
		}
		else if (!strcmp(argv[i],"--nwildcard") || !strcmp(argv[i],"-y")) {
			Z = 0; // global
			printf(" --> Setting N's and X's to wildcards (match anything)\n");
		}
		else if (!strcmp(argv[i],"--xalphabet") || !strcmp(argv[i],"-x")) {
			Xalpha = 1; // global
			printf(" --> Allowing any alphabet (unambiguous ID matching)\n");
		}
		else if (!strcmp(argv[i],"--taxonomy") || !strcmp(argv[i],"-b")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --taxonomy requires filename argument"); exit(1);}
			tax_FN = argv[i];
			printf(" --> Assigning taxonomy based on mapping file: %s\n",tax_FN);
		}
		else if (!strcmp(argv[i],"--mode") || !strcmp(argv[i],"-m")) {
			if (++i == argc || argv[i][0] == '-')  
				{ puts("ERROR: --mode requires an argument (see -h)"); exit(1); }
			if (!strcmp(argv[i],"BEST")) RUNMODE = BEST;
			else if (!strcmp(argv[i],"ALLPATHS")) RUNMODE = ALLPATHS;
			else if (!strcmp(argv[i],"CAPITALIST")) RUNMODE = CAPITALIST;
			else if (!strcmp(argv[i],"COMMUNIST")) RUNMODE = COMMUNIST;
			else if (!strcmp(argv[i],"ANY")) RUNMODE = ANY;
			else if (!strcmp(argv[i],"MATRIX")) 
				{fputs("ERROR: Matrix mode is no longer supported\n",stderr); exit(1);}
			else if (!strcmp(argv[i],"FORAGE")) RUNMODE = FORAGE;
			else {printf("Unsupported run mode '%s'\n",argv[i]); exit(1);}
			printf(" --> Setting run mode to %s\n",argv[i]);
		}
		else if (!strcmp(argv[i],"--makedb") || !strcmp(argv[i],"-d")) {
			makedb = 1; char *dbsel = "QUICK"; // make this read the default
			if (i + 1 != argc && argv[i+1][0] != '-' && !atol(argv[i+1])) { // non-numeric arg
				if (!strcmp(argv[++i],"DNA")) dbType = DNA_16;
				else if (!strcmp(argv[i],"RNA")) dbType = DNA_16;
				/*else if (!strcmp(argv[i],"PROTEIN")) dbType = PROT;*/
				else if (!strcmp(argv[i],"QUICK")) dbType = QUICK;
				else {printf("Unsupported makedb mode '%s'\n",argv[i]); exit(1);};
				dbsel = argv[i];
			}
			if (i + 1 != argc && argv[i+1][0] != '-') { // numeric arg provided
				DB_QLEN = atol(argv[++i]);
				if (DB_QLEN <= 0) {fprintf(stderr,"ERROR: bad max query length '%s'\n",argv[i]); exit(1);}
				if (DB_QLEN < 8) printf("WARNING: query length very short (%ld)\n",DB_QLEN);
			}
			printf(" --> Creating %s database (assuming max query length %ld)\n",dbsel, DB_QLEN);
		}
		else if (!strcmp(argv[i],"--dbpartition") || !strcmp(argv[i],"-dp")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --dbpartition requires integer argument"); exit(1);}
			int temp = atoi(argv[i]);
			if (temp < 0) {fputs("ERROR: numb partitions must be >= 1\n",stderr); exit(1);}
			RefDat.cparts = temp;
			printf(" --> Partitioning database into %d slices\n",temp);
		}
		else if (!strcmp(argv[i],"--accelerator") || !strcmp(argv[i],"-a")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --accelerator requires filename argument"); exit(1);}
			xcel_FN = argv[i];
			DO_ACCEL = 1;
			printf(" --> Using accelerator file %s\n",xcel_FN);
		}
		else if (!strcmp(argv[i],"--taxacut") || !strcmp(argv[i],"-bc")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --taxacut requires numeric argument"); exit(1);}
			int temp = atoi(argv[i]);
			if (temp < 2) { // interpret as float
				double fl = 1.0 / (1.0 - atof(argv[i]));
				temp = (int)(fl+0.5);
				printf(" --> Taxacut: converting %s to %d...\n",argv[i],temp);
			}
			if (temp < 2) {fputs("ERROR: taxacut must be >= 2\n",stderr); exit(1);}
				
			TAXACUT = temp;
			printf(" --> Ignoring 1/%d disagreeing taxonomy calls\n",TAXACUT);
		}
		else if (!strcmp(argv[i],"--taxa_ncbi") || !strcmp(argv[i],"-bn")) {
			taxa_lookup = taxa_lookup_ncbi;
			printf(" --> Using NCBI header formatting for taxonomy lookups\n");
		}
		else if (!strcmp(argv[i],"--skipambig") || !strcmp(argv[i],"-sa")) {
			QDat.skipAmbig = RefDat.skipAmbig = 1;
			printf(" --> Skipping highly ambiguous sequences\n");
		}
		else if (!strcmp(argv[i],"--taxasuppress") || !strcmp(argv[i],"-bs")) {
			QDat.taxasuppress = 1;
			if (i + 1 != argc && argv[i+1][0] != '-') { // arg provided
				if (!strcmp(argv[++i],"STRICT")) TAXLEVELS = TAXLEVELS_STRICT;
				else {fprintf(stderr,"ERROR: Unrecognized taxasuppress '%s'\n",argv[i]); exit(1);}
			}
			printf(" --> Surpressing taxonomic specificity by alignment identity%s\n",
				TAXLEVELS == TAXLEVELS_STRICT ? " [STRICT]" : "");
		}
		else if (!strcmp(argv[i],"--id") || !strcmp(argv[i],"-i")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --id requires decimal argument"); exit(1);}
			THRES = atof(argv[i]);
			if (THRES > 1.f || THRES < 0.f) {puts("Invalid id range [0-1]"); exit(1);}
			if (THRES < 0.01f) THRES = 0.01f;
			printf(" --> Setting identity threshold to %f\n",THRES);
		}
		else if (!strcmp(argv[i],"--threads") || !strcmp(argv[i],"-t")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --threads requires integer argument"); exit(1);}
			THREADS = HAS_OMP ? atoi(argv[i]) : 1;
			printf(" --> Setting threads to %d%s\n",THREADS, HAS_OMP ? "" : " (no MT!)");
		}
		else if (!strcmp(argv[i],"--shear") || !strcmp(argv[i],"-s")) {
			REBASE = 1; // global
			if (i + 1 != argc && argv[i+1][0] != '-') { // arg provided
				REBASE_AMT = atol(argv[++i]);
				if (REBASE_AMT < 0) {printf("ERROR: bad shear length '%s'\n",argv[i]); exit(1);}
			}
			if (!REBASE_AMT) REBASE = 0, puts("DISABLING shearing");
			else printf(" --> Shearing references longer than %ld\n",REBASE_AMT);
		}
		else if (!strcmp(argv[i],"--unique") || !strcmp(argv[i],"-u")) {
			doDedupe = 1; // global
			printf(" --> Preprocessing references (dereplicating)\n");
		}
		else if (!strcmp(argv[i],"--fingerprint") || !strcmp(argv[i],"-f")) {
			DO_FP = 1; // global
			printf(" --> Using fingerprint profiling\n");
		}
		else if (!strcmp(argv[i],"--prepass") || !strcmp(argv[i],"-p")) {
			DO_PREPASS = 16;
			if (i + 1 != argc && argv[i+1][0] != '-')
				DO_PREPASS = atoi(argv[++i]);
			printf(" --> Using ultra-heuristic prepass with effort %u%s\n",
				DO_PREPASS,DO_PREPASS?"":" [DISABLED]");
		}
		else if (!strcmp(argv[i],"--heuristic") || !strcmp(argv[i],"-hr")) {
			DO_HEUR = 1; // global
			printf(" --> WARNING: Heuristic mode set; optimality not guaranteed at low ids\n");
		}
		//
		else if (!strcmp(argv[i],"--noprogress")) {
			QDat.quiet = 1;
			printf(" --> Surpressing progress indicator\n");
		}
		else if (!strcmp(argv[i],"--cache") || !strcmp(argv[i],"-c")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --cache requires integer argument"); exit(1);}
			cacheSz = atoi(argv[i]); // global
			printf(" --> Setting number of cached lines in matrix to %d\n",cacheSz);
		}
		else if (!strcmp(argv[i],"--latency") || !strcmp(argv[i],"-l")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --latency requires integer argument"); exit(1);}
			LATENCY = atoi(argv[i]); // global
			printf(" --> Setting clump formation latency to %d bases\n",LATENCY);
		}
		else if (!strcmp(argv[i],"--clustradius") || !strcmp(argv[i],"-cr")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --clustradius requires integer argument"); exit(1);}
			RefDat.clustradius = atoi(argv[i]); // Reference_data member
			printf(" --> Setting FP cluster search radius to %d members\n",RefDat.clustradius);
		}
		else if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) PRINT_USAGE()
		else {
			printf("ERROR: Unrecognized command-line option: %s\n",argv[i]);
			puts("See help by running with just '-h'");
			exit(1);
		}
	}
	FILE *output = fopen(output_FN,"wb");
	if (!output) {fprintf(stderr,"ERROR: Cannot open output: %s\n",output_FN); exit(2);}

#ifdef _OPENMP
	omp_set_num_threads(THREADS);
#else
	THREADS = 1;
	puts("WARNING: Multi-threading not enabled in this build.");
#endif
	printf("Using up to "ASMVER" with %d threads.\n",THREADS);
	double start = omp_get_wtime();
	setScore(); // Prepare the scoring tables

	// Determine what to do
	if (makedb) {
		RefDat.dbType = dbType;
		puts("");
		if (!REBASE) DB_QLEN = 0;
		if (isRefEDB(ref_FN)) {fputs("ERROR: DBs can't make DBs.\n",stderr); exit(1);}
		process_references(ref_FN, &RefDat, DB_QLEN, 2); // DO_FP is integrated
		puts("Writing database...");
		dump_edb(output, &RefDat);
		puts("Database written.");
		if (DO_ACCEL) { // NEW [mem]
			printf("Generating accelerator '%s'\n",xcel_FN);
			free(RefDat.FingerprintsR.initP); free(RefDat.FingerprintsR.N);
			free(RefDat.Centroids); free(RefDat.RefClump); free(RefDat.ProfClump);
			make_accelerator(&RefDat, xcel_FN);
		}
		exit(0);
	}
	else { // alignment
		// Determine whether we're using a db or not and proceed to alignment	
		uint32_t dShear = 0, usedb;
		if (DO_ACCEL) read_accelerator(&RefDat, xcel_FN);
		if (usedb = isRefEDB(ref_FN)) 
			puts("\nEDB database provided. Parsing..."),
			dShear = read_edb(ref_FN, &RefDat);
		if (tax_FN) {
			size_t taxa_parsed = parse_taxonomy(tax_FN, &RefDat.Taxonomy);
			if (!taxa_parsed) {fputs("ERROR: invalid taxonomy\n",stderr); exit(1);}
			RefDat.taxa_parsed = taxa_parsed;
			int srtFuncStr(const void *a, const void *b) {
				return strcmp(((TaxPair_t *)a)->Head, ((TaxPair_t *)b)->Head); }
			qsort(RefDat.Taxonomy, taxa_parsed, sizeof(*RefDat.Taxonomy), srtFuncStr);
		}
		process_queries(query_FN, &QDat);
		if (!usedb) process_references(ref_FN, &RefDat, QDat.maxLenQ, doDedupe);
		else if (dShear && (uint32_t)(QDat.maxLenQ / THRES) > dShear) {
			fputs("ERROR: DB incompatible with selected queries/identity.\n",stderr);
			if (!DO_PREPASS && !DO_HEUR) exit(1);
			fputs("!!! WARNING: Error overridden by use of heuristic mode!\n",stderr);
		}
		REBASE_AMT = dShear;
		//alignNVU(RefDat.RefSeq[0], QDat.QSeq[0], RefDat.RefLen[0], QDat.QLen[0]);
		do_alignments(output, RefDat, QDat, usedb);
	}
	
	printf("\nAlignment time: %f seconds\n", omp_get_wtime() - start);
	return 0;
}
