/* EMBALMER aligner -- fast optimal aligner by Gabe. 
Copyright (C) 2015-2017 Knights Lab, Regents of the University of Minnesota.
This software is released under the GNU Affero General Public License (AGPL) v3.0.
*/
#define _LARGEFILE_SOURCE_
#define FILE_OFFSET_BITS 64
#include <stdio.h>
#include <inttypes.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#ifdef _WIN32
	#include <windows.h>
	#define fseek _fseeki64
	#define ftell _ftelli64
#endif
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
#include <immintrin.h>
#ifdef __AVX__
	#define ASMVER "AVX-128"
#elif __SSE4_1__
	#define ASMVER "SSE4.1"
#else
	#define _mm_popcnt_u64 __builtin_popcountll
	#define _mm_blendv_epi8(f,t,m) \
		_mm_xor_si128(_mm_and_si128(_mm_xor_si128(f,t),m),f) 
	#define _mm_blendv_ps(f,t,m) \
		_mm_xor_ps(_mm_and_ps(_mm_xor_ps(f,t),m),f)
	#define _mm_cvtepu8_epi32(x) \
		_mm_unpacklo_epi16(_mm_unpacklo_epi8(x,_mm_setzero_si128()),_mm_setzero_si128())
	#define _mm_extract_epi8(a,x) \
		(0xff & (_mm_extract_epi16(a,((x)>>1))>>(8*((x) & 0x01))))
	#ifdef __SSSE3__
		#define ASMVER "SSSE3"
	#else
		#define ASMVER "SSE2"
		#define _mm_cmple_epu16(x,y) \
			_mm_cmpeq_epi16(_mm_subs_epu16(x, y), _mm_setzero_si128())
		#define _mm_blendv_si128 (x, y, mask) \
			_mm_or_si128(_mm_andnot_si128(mask, x), _mm_and_si128(mask, y))
		#define _mm_min_epu16(x,y) \
			_mm_blendv_si128(y, x, _mm_cmple_epu16(x, y))
	#endif
#endif
#define PADDING 0

#include "math.h"
#define LORAM 1

typedef enum {
	FORAGE, BEST, ALLPATHS, CAPITALIST, 
	MATRIX, PAIRED, MAPPED, INLINE
} Mode;
typedef enum {DNA_16, DNA_8, DNA_4, PROT, DATA, QUICK} DBType;
typedef enum {NONE = 0, REORDER = 1, FAST = 2, FULL = 3, AUTO = 4} Prepass_t;
Mode RUNMODE = BEST;
Prepass_t DO_PREPASS = NONE; // don't turn it on by default
#define PP_DEF_ST "auto"
Prepass_t PP_DEF_ON = AUTO; // but default to auto if turned on
int LATENCY = 16;
int THREADS = 1;
int cacheSz = 150;
int Xalpha = 0;
int REBASE = 0;
int DO_FP = 0;
int DO_ACCEL = 0;
uint32_t TAXACUT = 10;
float THRES = 0.97f;
long REBASE_AMT = 500, DB_QLEN = 500;
#define VECSZ 16
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define PRINT_USAGE() { \
	puts("\nEMBALMER aligner (IMB_All-mer v0.99.3b) by Gabe."); \
	printf("Compiled with " ASMVER " and%s multithreading\n", !HAS_OMP ? " WITHOUT" : "");\
	puts("\nRequired parameters:");\
	puts("--references (-r) <name>: FASTA file of reference sequences to search");\
	puts("--queries (-q) <name>: FASTA file of queries to match against references");\
	puts("--output (-o) <name>: Text file for writing results (BLAST-like tdf)");\
	puts("\nBehavior parameters:"); \
	puts("--forwardreverse (-fr): also search the reverse complement of queries"); \
	puts("--whitespace (-w): write full query names in output (incl. whitespace)"); \
	puts("--npenalize (-n): Make A,C,G,T,U in query not match X,N in reference"); \
	puts("--xalphabet (-x): Allow any alphabet and disable ambiguity matching");\
	puts("--taxonomy (-b) <name>: assign taxonomy w/file. Interpolate with -m CAPITALIST");\
	puts("--mode (-m) <name>: Operating mode. [BEST] Pick one of these:");\
	puts("  BEST (reports first best match by hybrid BLAST id)"); \
	puts("  ALLPATHS (reports all ties with same error profile)"); \
	puts("  CAPITALIST (reports minimal set of references among ties)"); \
	puts("  FORAGE (reports all matches above specified threshold"); \
	puts("  [not enabled]: PAIRED, MAPPED, INLINE"); \
	puts("--makedb (-d) [name qLen]: Create a database from input references");\
	printf("  [name]: Optional. Can be PROT, DNA, RNA, or QUICK [QUICK]\n");\
	printf("     NOTE: Only QUICK supported in this release!\n"); \
	printf("  [qLen]: Optional, reqs [name]. Max query length to search in DB [%d]\n",DB_QLEN); \
	printf("--accelerator (-a) <name>: Create or use a helper DB\n"); \
	printf("\nPerformance parameters:\n"); \
	printf("--taxacut (-bc) <int>: ignore 1/[%u] disagreeing taxonomy calls\n",TAXACUT); \
	printf("--taxa_ncbi (-bn): Assume NCBI header '>xxx|accsn...' for taxonomy\n"); \
	printf("--skipambig (-sa): Do not consider highly ambiguous queries\n"); \
	printf("--taxasuppress (-bs) [STRICT]: Adjust taxonomy calls by %%ID\n"); \
	printf("--id (-i) <decimal>: similarity (range 0-1) needed to report match [%.2f]\n",THRES);\
	printf("--threads (-t) <int>: How many processors to use [%u]\n",THREADS);\
	printf("--shear (-s) [len]: Shear references longer than len bases [%ld]\n",REBASE_AMT);\
	printf("--unique (-u): Dereplicate references (lossless preprocessing)\n"); \
	printf("--fingerprint (-f): Use sketch fingerprinting to precheck matches\n"); \
	printf("--prepass (-p) [speed]: use fingerprints to pre-filter matches ["PP_DEF_ST"]\n"); \
	printf("  [speed]: Optional. Can be 'auto', full', 'fast', 'reorder', 'none'\n"); \
	printf("--cache (-c) <int>: Performance tweaking parameter [%d]\n",cacheSz); \
	printf("--latency (-l) <int>: Performance tweaking parameter [%d]\n",LATENCY); \
	printf("--clustradius (-cr) <int>: Performance tweaking parameter [auto]\n"); \
	puts("\n--help (-h): Shows this help screen with version info"); \
	puts("Example: embalm -r myRefs.fasta -q myQs.fasta -o outputs.txt -i 0.98");\
	exit(1);\
}

#define GAP 1
#define SCD 16
#define EDB_VERSION 2
#define EDX_VERSION 3
#define ACC_VERSION 0
#define ACC_VERSION_BIG 1
char Z = 0; // N mismatch penalty
char *BAK = "\0ACGTNKMRYSWBVHD\0"; 
char *RVC = "\0TGCANMKYRSWVBDH\0";
char RVT[] = {0,4,3,2,1,5,7,6,9,8,10,11,13,12,15,14};
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
	// #ifdef __AVX2__
	// __m256i v;
	// #elif __AVX__
	// __m256d v;
	// #endif
} Prince;

typedef struct {
	Prince *P;   // Prince
	uint8_t *N;  // Pops
	void *initP; 
	uint32_t nf, *Ptrs;
} PackaPrince;

typedef struct {
	uint32_t len, cap;
	uint32_t *Refs;
} Accelerant;
typedef struct AccelNode AccelNode;
struct AccelNode {
	AccelNode *next;
	uint32_t ref;
};

typedef struct {char *Head, *Tax;} TaxPair_t;
                           //K     P     C     O     F     G     S     SS+
float TAXLEVELS_STRICT[] = {.65f, .75f, .78f, .82f, .86f, .94f, .98f, .995f},
	 TAXLEVELS_LENIENT[] = {.55f, .70f, .75f, .80f, .84f, .93f, .97f, .985f},
	 *TAXLEVELS = TAXLEVELS_LENIENT;

typedef struct {
	char *Seq;    //8
	uint16_t div; //2
	uint32_t six; //4
	uint8_t rc;   //1
} UniBin; //15

typedef struct {
	uint32_t len; //4
	uint16_t ed;  //2
} ShrBin; //6

typedef struct {
	char **RefHead, **RefSeq;
	uint32_t *RefLen, *ClumpLen, *RefStart, *RefIxSrt, *TmpRIX, *RefDedupIx;
	uint32_t totR, origTotR, numRclumps, maxLenR;
	DualCoil **RefClump, **ProfClump;
	
	Prince *Centroids;
	PackaPrince FingerprintsR;
	TaxPair_t *Taxonomy;
	uint32_t taxa_parsed;
	uint32_t clustradius;
	uint32_t numRefHeads, *RefMap;

	uint32_t *BadList, badListSz;
	void **Accelerators;
} Reference_Data;

typedef struct {
	char **QHead, *SeqDumpRC; 
	uint32_t totQ, numUniqQ, maxLenQ, *Offset, *QBins, maxED, maxDiv;
	int incl_whitespace, taxasuppress, rc, skipAmbig;

	PackaPrince FingerprintsQ;
	UniBin *UniBins;
	ShrBin *ShrBins;
} Query_Data;

void * malloc_a(size_t algn, size_t size, void **oldPtr) {
    uintptr_t mask = ~(uintptr_t)(algn - 1);
	*oldPtr = malloc(size+algn-1);
    return (void *)(((uintptr_t)*oldPtr+algn-1) & mask);
}
void * calloc_a(size_t algn, size_t size, void **oldPtr) {
    uintptr_t mask = ~(uintptr_t)(algn - 1);
	*oldPtr = calloc(size+algn-1,1);
    return (void *)(((uintptr_t)*oldPtr+algn-1) & mask);
}
int cmpStrIx(const void *a, const void *b) {
	return strcmp(((StrIxPair*)a)->s,((StrIxPair*)b)->s); }

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
#define NIB4 s[0] << 12 | s[1] << 8 | s[2] << 4 | s[3]
#define NIB5 s[0] << 16 | s[1] << 12 | s[2] << 8 | s[3] << 4 | s[4]
static inline void parallel_sort_strpack(StrIxPair *Pack, uint32_t len) 
	PARALLEL_SORT_PROTOTYPE(StrIxPair, s, 5, NIB5)
static inline void parallel_sort_unibin(UniBin *Pack, uint32_t len) 
	PARALLEL_SORT_PROTOTYPE(UniBin, Seq, 5, NIB5)

// Taxonomy handling
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

char * (*taxa_lookup)(char *, uint32_t, TaxPair_t *) = taxa_lookup_generic;

size_t parse_taxonomy(char *filename, TaxPair_t **Obj) {
	static const size_t linelen = 10000000; //1k entries, 10 meg line lines
	size_t cur_sz = 1000, ns = -1;
	FILE *file = fopen(filename,"rb");
	if (!file) 
		{fprintf(stderr,"Cannot open TAXONOMY file: %s.\n",filename); exit(2);}
	TaxPair_t *T = malloc(cur_sz * sizeof(*T));
	char *line = calloc(linelen,1), *lineO = line; 
	while (line = fgets(line, linelen, file)) {
		if (++ns == cur_sz - 1) { // double all data structures
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
	}
	free(lineO);
	*Obj = realloc(T, ++ns*sizeof(*T)); // shrink T
	fclose(file);
	return ns;
}

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
				else Seqs[ns] = memcpy(malloc(len16+1),Seqs[ns],Lengths[ns]); 
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
// hasty, less robust parser; may be useful for parsing queries
size_t parse_tlY_fasta(char * filename, char ***HeadersP, char ***SeqsP, uint32_t **LengthsP) {
	static const size_t linelen = INT32_MAX; //1k entries, 2-gig lines
	FILE *file = fopen(filename,"rb");
	if (!file) { 
		fprintf(stderr,"Cannot open FASTA file: %s.\n",filename); exit(2); }
	size_t cur_sz = 10000, ns = 0;
	char **Headers = malloc(cur_sz * sizeof(*Headers)),
		 **Seqs = malloc(cur_sz * sizeof(*Seqs));
	uint32_t *Lengths = malloc(cur_sz * sizeof(*Lengths));
	char *line = malloc(linelen), *lineO = line; 
	if (!Headers || !Seqs || !Lengths || !line) {
		fputs("OOM:parse_tl_fasta\n",stderr); exit(3); }
	while (line = fgets(line, linelen, file)) {
		char *end = strchr(line,'\n');
		if (!end) {puts("End of FASTA reached"); break;}
		*end = 0;
		if (ns == cur_sz) {
			cur_sz += cur_sz;
			Headers = realloc(Headers, cur_sz*sizeof(*Headers));
			Seqs = realloc(Seqs, cur_sz*sizeof(*Seqs));
			Lengths = realloc(Lengths, cur_sz*sizeof(*Lengths));
			if (!Headers || !Seqs || !Lengths) {
				fputs("OOM:parse_tl_fasta\n",stderr); exit(3); }
		}
		Headers[ns] = memcpy(malloc(end-line),line+1,end-line);
		line = fgets(line, linelen, file);
		if (!line || !(end = strchr(line,'\n'))) {
			puts("WARNING: Premature end of FASTA reached"); break;}
		*end = 0;
		Seqs[ns] = memcpy(calloc(end-line+15,1),line,end-line+1);
		Lengths[ns] = end-line;
		++ns;
	}
	free(lineO);
	Headers = realloc(Headers,ns*sizeof(*Headers)), Seqs = realloc(Seqs,ns*sizeof(*Seqs));
	Lengths = realloc(Lengths,ns*sizeof(*Lengths));
	*HeadersP = Headers; *SeqsP = Seqs; *LengthsP = Lengths;
	if (ns >= UINT32_MAX) puts("WARNING: >4 billion sequences processed.");
	return ns;
}

static inline char * findNL(char *a) {
	__m128i nil = _mm_set1_epi8(0), nl = _mm_set1_epi8('\n');
	while (1) {
		__m128i v = _mm_lddqu_si128((void*)a);
		__m128i n = _mm_cmpeq_epi8(v,nil);
		__m128i e = _mm_cmpeq_epi8(v,nl);
		if (_mm_movemask_epi8(n)) return 0;
		uint16_t x = _mm_movemask_epi8(e);
		if (x) return a + __builtin_ctz(x);
		a += 16;
	}
}
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

// hastier, even less robust parser; may still be useful for parsing queries
size_t parse_tlX_fasta(char * filename, char ***HeadersP, char ***SeqsP, uint32_t **LengthsP) {
	static const size_t linelen = 1<<24; 
	FILE *file = fopen(filename,"rb");
	if (!file) { 
		fprintf(stderr,"Cannot open FASTA file: %s.\n",filename); exit(2); }
	

	size_t cur_sz = 10000, ns = 0;
	char **Headers = malloc(cur_sz * sizeof(*Headers)),
		 **Seqs = malloc(cur_sz * sizeof(*Seqs));
	uint32_t *Lengths = malloc(cur_sz * sizeof(*Lengths));
	char *buf1 = calloc(linelen+16,1); 
	if (!Headers || !Seqs || !Lengths || !buf1) {
		fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
	int last = 0; // denotes this is the last buffer in the file
	size_t numRead = fread(buf1,1,linelen,file);
	if (numRead < linelen) last = 1;
	char *begin, *end = buf1 - 1;
	while (1) {
		begin = end + 1, end = findNL(begin); //strchr(begin,'\n');
		if (!end) {
			if (!last) {
				uint64_t rem = linelen + buf1 - begin;
				buf1 = memmove(buf1, begin, rem);
				numRead = fread(buf1+rem, 1, linelen - rem, file);
				if (numRead < linelen - rem) last = 1;
				begin = buf1, end = findNL(begin); //strchr(begin, '\n');
			}
			if (!end) {printf("End of FASTA reached [H %u]\n",ns); break;} 
		} 
		*end = 0;
		if (ns == cur_sz) {
			cur_sz += cur_sz;
			Headers = realloc(Headers, cur_sz*sizeof(*Headers));
			Seqs = realloc(Seqs, cur_sz*sizeof(*Seqs));
			Lengths = realloc(Lengths, cur_sz*sizeof(*Lengths));
			if (!Headers || !Seqs || !Lengths) {
				fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
		}
		Headers[ns] = memcpy(malloc(end-begin),begin+1,end-begin);
		//printf(">%s\n",Headers[ns]);

		begin = end + 1, end = findNL(begin); //strchr(begin,'\n');
		if (!end) {
			if (!last) {
				uint64_t rem = linelen + buf1 - begin;
				buf1 = memmove(buf1, begin, rem);
				numRead = fread(buf1+rem, 1, linelen - rem, file);
				if (numRead < linelen - rem) last = 1;
				begin = buf1, end = findNL(begin); //strchr(begin, '\n');
			}
			if (!end) {
				printf("End of FASTA reached [SE %u]\n",ns); 
				end = strchr(begin,0);
				Seqs[ns] = memcpy(calloc(end-begin+15,1),begin,end-begin+1);
				//puts(Seqs[ns]);
				Lengths[ns++] = end-begin;
				break;
			}
		}
		*end = 0;
		Seqs[ns] = memcpy(calloc(end-begin+15,1),begin,end-begin+1);
		//printf("%s\n",Seqs[ns]);
		Lengths[ns] = end-begin;
		++ns;
	}
	free(buf1);
	Headers = realloc(Headers,ns*sizeof(*Headers)), Seqs = realloc(Seqs,ns*sizeof(*Seqs));
	Lengths = realloc(Lengths,ns*sizeof(*Lengths));
	*HeadersP = Headers; *SeqsP = Seqs; *LengthsP = Lengths;
	if (ns >= UINT32_MAX) puts("WARNING: >4 billion sequences processed.");
	return ns;
}

// fastest cached fasta file parser in the upper midwest
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
	char *nl = findNL_abs(Headers[0]); //,'\n'); 
	*nl = 0;
	Seqs[0] = nl + 1;
	nl = findNL_or_eol(Seqs[0]); //,'\n');
	//if (!nl) nl = strchr(Seqs[0],0);
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
		char *nl = findNL_abs(s); //, '\n');
		*nl = 0;
		s = nl + 1;
		Seqs[tix] = s;
		nl = findNL_or_eol(s); //,'\n');
		//if (!nl) nl = strchr(s,0);
		*nl = 0;
		Lengths[tix] = nl - s;
	}
	if (ix > UINT32_MAX) puts("WARNING: greater than 4 billion queries");
	*HeadersP = Headers; *SeqsP = Seqs; *LengthsP = Lengths;
	return ix;
}

// MT version of the hasty original
size_t parse_tl_fastaMT(char * filename, char ***HeadersP, char ***SeqsP, uint32_t **LengthsP) {
	static const size_t linelen = INT32_MAX; 
	FILE *file = fopen(filename,"rb");
	if (!file) { 
		fprintf(stderr,"Cannot open FASTA file: %s.\n",filename); exit(2); }
	size_t cur_sz = 10000, ns = 0;
	char **Headers = malloc(cur_sz * sizeof(*Headers)),
		 **Seqs = malloc(cur_sz * sizeof(*Seqs));
	uint32_t *Lengths = malloc(cur_sz * sizeof(*Lengths));
	if (!Headers || !Seqs || !Lengths) {
		fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
	int num_th = MIN(THREADS,2);
	char *buffer = malloc(linelen+16), *lastStop = 0; 
	memset(buffer+linelen-1,0,16);
	uint32_t leftover = 0, moreLeft = 1, toGetThru, tTot = 0;
	#pragma omp parallel num_threads(num_th)
	{
		uint32_t tsz = 10000, tix = 0, tid = omp_get_thread_num();
		uint32_t *LenT = malloc(tsz*sizeof(*LenT));
		char **HeadT = malloc(tsz*sizeof(*HeadT)), **SeqT = malloc(tsz*sizeof(*SeqT));
		
		while (moreLeft) {
			tix = 0; // reset thread-local compartment index; keep old size intact
			//puts("Beginning moreLeft loop");
			#pragma omp barrier
			#pragma omp master
			{
				puts("Getting a new chunk.");
				if (leftover) memmove(buffer,lastStop,leftover);
				size_t numRead = fread(buffer+leftover,1,linelen-leftover,file);
				toGetThru = numRead + leftover;
				if (numRead < linelen-leftover) moreLeft = 0, leftover = 0;
				else {
					lastStop = buffer + linelen;
					while (lastStop > buffer) if (*--lastStop=='>') break;  
					leftover = buffer + linelen - lastStop;
				}
				printf("Numread: %u\n",numRead);
			}
		#pragma omp barrier
		//#pragma omp parallel num_threads(num_th)
		//{
			//uint32_t tid = omp_get_thread_num();
			uint32_t radius = toGetThru / num_th; //, ix;
			if (radius < (1 << 24)) radius = MIN((1 << 24), toGetThru);

			size_t stIx = tid * radius;
			if (stIx <= toGetThru) {
				if (stIx + radius > toGetThru) radius = MIN(toGetThru, 2*radius);

				//stIx += ns; radius += stIx; // so it becomes start to end
				//radius += stIx; // start to end of buffer
				char *st = buffer + stIx, *bounds = MIN(lastStop,st + radius);
				while (st < bounds) {
					st = strchr(st,'>');
					if (!st) {printf("Bollocks. %u\n",ns); exit(11);}
					//uint32_t hasNL = 1;
					char *ed = strchr(st,'\n');
					if (!ed) ed = strchr(st,0); //, hasNL = 0;
					if (!ed) {printf("OAW ballocks.\n"); continue;}
					
					
					
					char *head = memcpy(malloc(ed-st),st+1,ed-st-1);
					head[ed-st-1] = 0;
					//Headers[ix][ed-st-1] = 0;
					st = ed + 1;
					ed = strchr(st,'\n'); //findNL_al(line2);
					if (!ed) {
						//hasNL = 0;
						puts("End of FASTA reached (sans newline)");
						ed = strchr(st,0);
						if (ed >= bounds) {printf("I go HARD.\n"); continue;}
						printf("Bounds, ed now: %llu,%llu\n",bounds,ed);
						if (!ed) {puts("ERROR1"); exit(10);}
					} //else hasNL = 1;
					//memset(end,0,16);
					//printf("%s\n\n",line2);
					uint32_t sa = ed - st + 1; // + hasNL;
					sa += 15 & (16 - (sa & 15));
					// char *dest = malloc(sa);
					// for (uint32_t i = 0; i < sa; i+=16) 
					// 	_mm_store_si128((void*)dest+i,_mm_load_si128((void*)line2+i));
					// Seqs[ix] = dest;
					if (tix >= tsz) {
						//printf("Tid %u Current size: %u, ns: %u\n",tid, tsz, tix);
						tsz += tsz;
						HeadT = realloc(HeadT, tsz*sizeof(*HeadT));
						SeqT = realloc(SeqT, tsz*sizeof(*SeqT));
						LenT = realloc(LenT, tsz*sizeof(*LenT));
						if (!HeadT || !SeqT || !LenT) {
							fputs("OOM:parse_MT_fasta\n",stderr); exit(3); }
					}
					HeadT[tix] = head;
					SeqT[tix] = memcpy(malloc(sa),st,ed - st); //realloc(line2,end-line2);
					memset(SeqT[tix]+(ed-st),0,sa-(ed-st));
					LenT[tix] = ed-st;
					/* 
					#pragma omp critical
					if (ns + num_th + num_th >= cur_sz) {
						printf("Current size: %u, ns: %u\n", cur_sz, ns);
						cur_sz += cur_sz;
						Headers = realloc(Headers, cur_sz*sizeof(*Headers));
						Seqs = realloc(Seqs, cur_sz*sizeof(*Seqs));
						Lengths = realloc(Lengths, cur_sz*sizeof(*Lengths));
						if (!Headers || !Seqs || !Lengths) {
							fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
					}
					if (ed-st) {  
						#pragma omp atomic capture
						ix = ns++;
						Headers[ix] = head;
						Seqs[ix] = memcpy(malloc(sa),st,ed - st); //realloc(line2,end-line2);
						memset(Seqs[ix]+(ed-st),0,sa-(ed-st));
						Lengths[ix] = ed-st;
					} else printf("If I don't go all the way, it'll be a complete waste of time\nAnd EnerGAAAAY\nand Money\nand MONEY\nmoney\nmoney\nmoney\n");
					 */
					st = ed + 1;
					++tix;
				}
				printf("Thread %u committed %u entries.\n",tid,tix);
				
				// commit to global structures
				//#pragma omp barrier
				#pragma omp atomic
				tTot += tix;
				#pragma omp barrier
				#pragma omp master
				{
					printf("GLOBAL new size: %u\n", tTot);
					//cur_sz += cur_sz;
					Headers = realloc(Headers, tTot*sizeof(*Headers));
					Seqs = realloc(Seqs, tTot*sizeof(*Seqs));
					Lengths = realloc(Lengths, tTot*sizeof(*Lengths));
					if (!Headers || !Seqs || !Lengths) {
						fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
				}
				/* if (tTot >= cur_sz) {
					printf("GLOBAL resizing from %u [current combined num %u]\n",cur_sz,tTot);
					cur_sz = tTot + tTot;
					Headers = realloc(Headers, cur_sz*sizeof(*Headers));
					Seqs = realloc(Seqs, cur_sz*sizeof(*Seqs));
					Lengths = realloc(Lengths, cur_sz*sizeof(*Lengths));
					if (!Headers || !Seqs || !Lengths) {
						fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
				} */
				#pragma omp barrier
				//#pragma omp critical
				for (uint32_t i = 0; i < tix; ++i) {
					uint32_t ix;
					#pragma omp atomic capture
					ix = ns++;
					//printf("[%u] Writing entry %u of %u into %u...\n",tid,i,tix,ix);
					/* if (ns == cur_sz) {
						printf("GLOBAL Current size: %u, ns: %u\n", cur_sz, ns);
						cur_sz += cur_sz;
						Headers = realloc(Headers, cur_sz*sizeof(*Headers));
						Seqs = realloc(Seqs, cur_sz*sizeof(*Seqs));
						Lengths = realloc(Lengths, cur_sz*sizeof(*Lengths));
						if (!Headers || !Seqs || !Lengths) {
							fputs("OOM:parse_tlX_fasta\n",stderr); exit(3); }
					} */
					Headers[ix] = HeadT[i];
					Seqs[ix] = SeqT[i];
					Lengths[ix] = LenT[i];
				}
			}
			printf("I'm thread %u. More to go? %u\n",tid,moreLeft);
		}
		
		free(LenT), free(HeadT), free(SeqT);
		puts("Freed variables; ending parallel section");
	}
	
	Headers = realloc(Headers,ns*sizeof(*Headers)), Seqs = realloc(Seqs,ns*sizeof(*Seqs));
	Lengths = realloc(Lengths,ns*sizeof(*Lengths));
	*HeadersP = Headers; *SeqsP = Seqs; *LengthsP = Lengths;
	if (ns >= UINT32_MAX) puts("WARNING: >4 billion sequences processed.");
	return ns;
}



// universal non-vectorized aligner; needs nothing but ref, q, and lengths
inline float alignNVU(char *ref, char *query, unsigned rlen, unsigned qlen) {
	static const int DIAG = 0, UP = 1, LEFT = 2;
	++rlen, ++qlen;
	unsigned *scores = calloc(rlen*qlen,sizeof(*scores)); //malloc
	unsigned *shifts = calloc(rlen*qlen,sizeof(*shifts));
	unsigned *traces = calloc(rlen*qlen,sizeof(*traces));
	for (int i = 0; i < qlen; ++i) scores[i*rlen] = i*GAP; // for edit distance
	//for (int j = 0; j < rlen; ++j) matrix[j] = 0; //-j; // fill left to right
	for (int y = 1; y < qlen; ++y) {
		for (int x = 1; x < rlen; ++x) {
			int scoreD = scores[(y-1)*rlen + x-1] + SCORENVedN[ref[x-1]*SCD + query[y-1]];
			int scoreU = scores[(y-1)*rlen + x] + GAP;
			int scoreL = scores[y*rlen + x-1] + GAP;
			int shiftD = shifts[(y-1)*rlen + x-1];
			int shiftU = shifts[(y-1)*rlen + x];
			int shiftL = shifts[y*rlen + x-1] + 1;
			int score, shift, trace;
			
			#define GOLEFT {score = scoreL, shift = shiftL, trace = LEFT;}
			#define GOUP {score = scoreU, shift = shiftU, trace = UP;}
			#define GODIAG {score = scoreD, shift = shiftD, trace = DIAG;}
			GODIAG
			
			// Edit-distance only
			if (scoreU < score) GOUP
			if (scoreL < score) GOLEFT
			
			// Hybrid distance
			//if ((scoreU < score) || (scoreU==score && shiftU > shift)) GOUP
			//if ((scoreL < score) || (scoreL==score && shiftL > shift)) GOLEFT
			
			// ID only
			//if ((float)scoreU/(y+shiftU) < (float)score/(y+shift)) GOUP
			//if ((float)scoreL/(y+shiftL) < (float)score/(y+shift)) GOLEFT
			
			scores[y*rlen + x] = score;
			shifts[y*rlen + x] = shift;
			traces[y*rlen + x] = trace;
		}
	}
	int score = INT_MAX, lasty = (qlen-1)*rlen, l=0;
	//for (int i = 1; i < rlen; ++i) if (scores[lasty + i] < score) l=i, score = scores[lasty + i];
	//printf("\nMin score found at %d = %d [%d]\n", l, (int)score, (int)shifts[lasty + l]);
	
	float c = 0, t; int ti, last=INT_MAX, lv=0;
	for (int i = 1; i < rlen; ++i) 
		//if ((t=(scores[lasty+i]+(qlen-1))/2+(qlen-1+shifts[lasty+i])) > c) lv=i, c = t;
		if ((t=1.f-(float)scores[lasty+i]/(qlen-1+shifts[lasty+i])) > c) lv=i, c = t;
		//{ti=scores[lasty+i]; if (ti < last || (ti==last && shifts[lasty+i] > shifts[lasty+lv])) lv=i, last=ti;}
	//c = 1.f - (float)last/(qlen-1+shifts[lasty+lv]);
	
	///////////////// To skip output //////////////////
			//free(scores); free(shifts); free(traces);
			//return c; 
	///////////////////////////////////////////////////
	
	printf("composite at %d = {%f} %d [%d]\n", lv, c, (int)scores[lasty+lv], 
		(int)shifts[lasty+lv]); //, (int)gaps[lasty + lv]);
	
	// print output
	int tb = qlen-1+shifts[lasty + lv], curX = lv, curY = qlen-1; 
	char refString[tb+1], qString[tb+1], gString[tb+1], dString[tb+1];
	for (int i=0; i<tb;++i) refString[i]='P', qString[i]='P',gString[i]='P',dString[i]='P';
	refString[tb] = 0; qString[tb] = 0; gString[tb] = 0, dString[tb] = 0;
	printf("tb = %d, curX = %d, curY = %d\n",tb,curX,curY);
	
	printf("Length of alignment = %d\n",tb);
	--tb;
	do {
		if (traces[(curY)*rlen + curX] == DIAG) {
			refString[tb] = curX ? BAK[ref[curX-1]] : '-';
			qString[tb] = curY ? BAK[query[curY-1]] : '-';
			gString[tb] = SCORENVedN[ref[curX-1]*SCD + query[curY-1]] ? ' ' : '|';
			--curX; --curY;
			dString[tb] = '='; //printf("=");
		}
		else if (traces[(curY)*rlen + curX] == UP) { // U 
			refString[tb] = '-';
			qString[tb] = curY ? BAK[query[curY-1]] : '-';
			gString[tb] = ' ';
			--curY;
			dString[tb] = '<'; //printf("<");
		}
		else { // LEFT
			refString[tb] = curX ? BAK[ref[curX-1]] : '-';
			qString[tb] = '-'; //BAK[query[curY-1]];
			gString[tb] = ' ';
			--curX;
			dString[tb] = '>'; //printf(">");
		}
	} while (tb--);
	puts("");
	printf("%s\n%s\n%s\n%s\n",dString, refString, gString, qString);
	free(scores); free(shifts); free(traces);
	return c; // or return score;
}

// Must set HiBound[1]=1 at init and LoBound[1]=rwidth in loop
			/* if (Xalpha) score = _mm_add_epi8(_mm_cmpeq_epi8(_mm_set1_epi8(qLet), 
				rChunk),_mm_set1_epi8(1));
			else
			#ifdef __SSSE3__
				score = _mm_shuffle_epi8(SCOREFAST[qLet], rChunk); 
			#else
				score = _mm_load_si128((void*)(profile + (x-1)*SCD + qLet));
			#endif 
			if (Xalpha) diagSc = _mm_add_epi8(_mm_cmpeq_epi8(_mm_set1_epi8(qLet), 
				rChunk),_mm_set1_epi8(1));
			else
			#ifdef __SSSE3__
				diagSc = _mm_shuffle_epi8(SCOREFAST[qLet], rChunk); 
			#else
				diagSc = _mm_load_si128((void*)(profile + (x-1)*SCD + qLet));
			#endif*/

#define DIAGSC_XALPHA _mm_add_epi8(_mm_cmpeq_epi8(_mm_set1_epi8(qLet), \
	rChunk),_mm_set1_epi8(1))
#ifdef __SSSE3__
	#define DIAGSC_MAT16 _mm_shuffle_epi8(SCOREFAST[qLet], rChunk)
#else
	#define DIAGSC_MAT16 _mm_load_si128((void*)(profile + (x-1)*SCD + qLet))
#endif

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

inline void reScoreM_mat16(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, 
 DualCoil *Matrix, DualCoil *Shifts, DualCoil *ShiftR, uint32_t maxED, DualCoil *profile, MetaPack *M16) 
	RESCOREM_PROTYPE(DIAGSC_MAT16)

inline void reScoreM_xalpha(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, 
 DualCoil *Matrix, DualCoil *Shifts, DualCoil *ShiftR, uint32_t maxED, DualCoil *profile, MetaPack *M16) 
	RESCOREM_PROTYPE(DIAGSC_XALPHA)

// No LoBound, StartQ, HiBound, MinA
#define PRUNE_ED_PROTOTYPE(DIAG_FUNC) {\
	uint32_t y, x; \
	__m128i maxEDv = _mm_set1_epi8(maxED+1); /* < 255 ? maxED+1 : 255); /* BAD if score >= maxED+1 */ \
	--query, --ref, ++qlen, ++rwidth; \
	maxED = maxED < qlen ? maxED : qlen; \
	DualCoil *restrict prev = Matrix + width, *restrict cur = Matrix; \
	_mm_store_si128((void*)(cur),_mm_set1_epi8(1)); \
	char qLet = query[1]; \
	for (x = 1; x < rwidth; ++x) { \
		__m128i rChunk = _mm_load_si128((void*)(ref+x)); \
		__m128i score = DIAG_FUNC; \
		_mm_store_si128((void*)(cur+x),score); \
	} \
	for (y = 2; y <= maxED; ++y) { \
		char qLet = query[y]; \
		DualCoil *temp = cur; cur = prev, prev = temp; \
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); /* column 0 */ \
		for (x = 1; x < rwidth; ++x) { \
			__m128i rChunk = _mm_load_si128((void*)(ref+x)); \
			__m128i prevRow_x = _mm_load_si128((void*)(prev+x)); \
			__m128i prevRow_x_1 = _mm_load_si128((void*)(prev+(x-1))); \
			__m128i curRow_x_1 = _mm_load_si128((void*)(cur+(x-1))); \
			__m128i diagSc = DIAG_FUNC; \
			__m128i score = _mm_adds_epu8(prevRow_x_1, diagSc); \
			__m128i scoreU = _mm_adds_epu8(prevRow_x, _mm_set1_epi8(GAP)); \
			score = _mm_min_epu8(scoreU,score); \
			__m128i scoreL = _mm_adds_epu8(curRow_x_1,_mm_set1_epi8(GAP)); \
			score = _mm_min_epu8(scoreL,score); \
			_mm_store_si128((void*)(cur+x),score); \
		} \
	} \
	uint32_t LB = 1, HB = rwidth, LBN = 1, HBN = rwidth; \
	for (; y < qlen; ++y) { \
		LB = LBN, HB = HBN; \
		LBN = 0; \
		char qLet = query[y]; \
		DualCoil *temp = cur; cur = prev, prev = temp; \
		_mm_store_si128((void*)(cur),_mm_set1_epi8(MIN(y,255))); \
		for (x = LB; x < HB; ++x) { \
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
				if (!LBN) LBN = x; \
				HBN = x; \
			} \
			_mm_store_si128((void*)(cur+x),score); \
		} \
		\
		if (!LBN) return -1; \
		++LBN; /* we guarantee the bounds will close in */ \
		++HBN; /* since it'll bound the 'for (...x < rightBound; ++x)' */ \
		_mm_store_si128((void*)(cur+HBN),_mm_set1_epi8(-1)); /* kill the right block */ \
		/*if (LBN > 1)*/ \
		_mm_store_si128((void*)(prev+LBN - 1), _mm_set1_epi8(-1)); /* kill left of new */ \
		HBN += HBN < rwidth; \
	} \
	/* Vectorized min reduction */ \
	__m128i minIntV = _mm_set1_epi8(-1); \
	for (uint32_t i = LB; i < HB; ++i) \
		minIntV = _mm_min_epu8(minIntV,_mm_load_si128((void*)(cur+i))); \
	__m128i c = _mm_srli_si128(minIntV,8); \
	minIntV = _mm_min_epu8(minIntV,c); \
	c = _mm_srli_si128(minIntV,4); \
	minIntV = _mm_min_epu8(minIntV,c); \
	c = _mm_srli_si128(minIntV,2); \
	minIntV = _mm_min_epu8(minIntV,c); \
	c = _mm_srli_si128(minIntV,1); \
	minIntV = _mm_min_epu8(minIntV,c); \
	return _mm_extract_epi8(minIntV,0); /* save min of mins as new bound */ \
}

inline uint32_t prune_ed_mat16(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, 
	uint32_t width, DualCoil *Matrix, DualCoil *profile, uint32_t maxED) PRUNE_ED_PROTOTYPE(DIAGSC_MAT16)


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
inline uint32_t aded_mat16(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, DualCoil *Matrix,
 DualCoil *profile, uint32_t maxED, uint32_t startQ, uint32_t *LoBound, uint32_t *HiBound, DualCoil *MinA) 
 ADED_PROTOTYPE(DIAGSC_MAT16)
inline uint32_t aded_xalpha(DualCoil *ref, char *query, uint32_t rwidth, uint32_t qlen, uint32_t width, DualCoil *Matrix,
 DualCoil *profile, uint32_t maxED, uint32_t startQ, uint32_t *LoBound, uint32_t *HiBound, DualCoil *MinA) 
 ADED_PROTOTYPE(DIAGSC_XALPHA)


inline void translateNV(char* string, size_t len) { // makes string into nums
	for (size_t i = 0; i < len; ++i) string[i] = CHAR2NUM[string[i]]; }

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

void setScore() {
/*{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //. [0]
	//. A C G T N K M R Y S W B V H D
	 -1,0,1,1,1,Z,1,0,0,1,1,0,1,0,0,0, //A [1]
	 -1,1,0,1,1,Z,1,0,1,0,0,1,0,0,0,1, //C [2]
	 -1,1,1,0,1,Z,0,1,0,1,0,1,0,0,1,0, //G [3]
	 -1,1,1,1,0,Z,0,1,1,0,1,0,0,1,0,0, //T/U [4]
	 -1,0,0,0,0,Z,0,0,0,0,0,0,0,0,0,0, //N/X [5]
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
	SCOREFAST[5]  = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,Z,0,0,0,0,-1); //N or X
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

typedef struct {
	char *seq;
	uint32_t len, ix;
} Tuxedo;
static int cpcmp(const void *a, const void *b) 
	{ return strcmp(**(char ***)a, **(char ***)b); }
static int cmpPackLen(const void *first, const void *second) {
	Tuxedo *a = (Tuxedo *)first, *b = (Tuxedo *)second;
	return a->len < b->len ? -1 : a->len > b->len;
}
static int cmpPackSeq(const void *first, const void *second) {
	Tuxedo *a = (Tuxedo *)first, *b = (Tuxedo *)second;
	return strcmp(a->seq, b->seq);
}

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

// make_db_fasta() GOES HERE

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

// curate is 0 (no dedupe), 1 (dedupe), 2 (make db so skip sse2 profiles)
static inline void process_references(char *ref_FN, Reference_Data *Rd, uint32_t maxLenQ, int curate) {
	Rd->totR = parse_tl_fasta(ref_FN, &Rd->RefHead, &Rd->RefSeq, &Rd->RefLen);
	printf("Parsed %u references.\n",Rd->totR);

	// Translate nucleotides into parallel register lookups
	if (!Xalpha) for (uint32_t i = 0; i < Rd->totR; ++i) 
		translate16aln(Rd->RefSeq[i],Rd->RefLen[i]);

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

		//Now sort refs lexicographically within tolerance range. [TODO: MT]
		uint32_t curTol = SeqLenPair[0].len, prev_ix = 0, lat = DO_FP ? 0 : LATENCY;
		for (uint32_t i = 1; i < Rd->totR; ++i) {
			if (SeqLenPair[i].len > curTol + lat) {
				curTol = SeqLenPair[i].len;
				qsort(SeqLenPair + prev_ix, i - prev_ix, sizeof(*SeqLenPair),cmpPackSeq);
				prev_ix = i;
			}
		}
		if (prev_ix < Rd->totR-1)
			qsort(SeqLenPair + prev_ix, Rd->totR - prev_ix, sizeof(*SeqLenPair),cmpPackSeq);
		
		for (uint32_t i = 0; i < Rd->totR; ++i) Rd->RefIxSrt[i] = SeqLenPair[i].ix;
		free(SeqLenPair);
		
	} else for (uint32_t i = 0; i < Rd->totR; ++i) 
		Rd->maxLenR = Rd->RefLen[i] > Rd->maxLenR ? Rd->RefLen[i] : Rd->maxLenR,
		Rd->RefIxSrt[i] = i;
	Rd->origTotR = Rd->totR, Rd->TmpRIX = Rd->RefIxSrt;

	if (curate) {   // dedupe references (redefine RefIxSrt using RefDedupIx map)
		uint64_t tot_divR = 0, rdupes = 0;
		uint32_t uix = 0; // dupe-map, map[uniq_ix] -> orig_ix
		Rd->RefDedupIx = calloc(Rd->totR+1,sizeof(*Rd->RefDedupIx));
		for (uint32_t i = 1, j; i < Rd->totR; ++i) {
			for (j = 0; j < Rd->RefLen[Rd->RefIxSrt[i]]; ++j) 
				if (Rd->RefSeq[Rd->RefIxSrt[i]][j] != Rd->RefSeq[Rd->RefIxSrt[i-1]][j]) break;
			tot_divR += j;
			if (j == Rd->RefLen[Rd->RefIxSrt[i]] && j == Rd->RefLen[Rd->RefIxSrt[i-1]]) ++rdupes; 
			else Rd->RefDedupIx[++uix] = i;
		}
		Rd->RefDedupIx[++uix] = Rd->totR; // add an end-cap for closure
		printf("Avg. R Divergence: %f (%u dupes, %u uniq)\n",(double)tot_divR/Rd->totR,rdupes,uix);

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
			puts("WARNING: penalizing N's in fingerprints may slow unpenalized alignments.");
			Prince t; // New: swap to error-free fingerprints anyway (higher quality clusters for N-case and without)
			for (uint32_t i = 0; i < Rd->totR; ++i) t = p[i], p[i] = p[Fp.Ptrs[i]], p[Fp.Ptrs[i]] = t; // swap FPs
		}
		
		// prelim stats: individual coverage.
		uint32_t *RefIxSrt = Rd->RefIxSrt, *RefLen = Rd->RefLen, totR = Rd->totR;
		
		uint64_t tot_cl = 0;
		double similarity = 0.06;
		#pragma omp parallel for reduction(+:tot_cl)
		for (uint32_t i = 0; i < totR; ++i) tot_cl += FP_pop(p+i);
		uint32_t sig_thres = Rd->clustradius ?: 1 + similarity * (double)tot_cl/totR;
		if (!Rd->clustradius && sig_thres < 5) sig_thres = 5;
		else if (sig_thres > 999) similarity = (double)(sig_thres-1000)/1000, sig_thres = 0;
		printf("Average coverage (atomic) = %f, cluster radius: %u\n",(double)tot_cl/totR, sig_thres);
		#define PTHRES 240
		#define PTHRES2 240
		uint32_t IXBIN[THREADS];
		
		
		// The direct fingerprint way
		double wtime = omp_get_wtime();
		#define BANDED_FP 1000000
		// Heuristic: sort FPs by first 16 bits (Prince.h), then iterate over bounded range
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
		}


		Split **CC = calloc((totR-0),sizeof(*CC));
		if (!CC) {fputs("OOM:FP:CC\n",stderr); exit(3);}
		uint64_t numLinks = 0;
		#pragma omp parallel
		{
			Split *CCL = malloc(totR*sizeof(*CCL));
			if (!CCL) {fputs("OOM:CCL\n",stderr); exit(3);}
			#pragma omp for schedule(dynamic,1) reduction(+:numLinks)
			for (uint32_t i = 0; i < totR; ++i) {
				// Get all the acceptable unions for this reference
				uint32_t cpop = FP_pop(p+i), x = 0,
					st = cpop + (sig_thres ?: similarity * cpop);
				if (cpop > PTHRES) continue; // otherwise they'll match mostly anything
				st = MIN(st,PTHRES2); // FPs aren't useful past this range
				uint32_t bounds = MIN(totR, (i+1)+BANDED_FP);
				for (uint32_t j = i + 1, t; j < bounds; ++j) 
					if ((t=FP_union(p[i],p[j])) < st) CCL[x++] = (Split){t,j};
				// or do something like MIN16 instead of arbitrary threshold?
				
				// Store results. First bin is {length,max_ix}
				if (x) { // make a new bin, copy these in. 
					CC[i] = malloc((x+1)*sizeof(*CC[i]));
					for (uint32_t j = 0; j < x; ++j) {
						CC[i][j+1] = CCL[j];
					}
					CC[i][0] = (Split){x,0}; //mix}; // store num bins for easy lookup
					numLinks += x;
				}
			}
			free(CCL);
		}
		// To cluster: sort all pairs within threshold cost of fusion
		// Pick off the lowest available, use as centroid if val not -1.
		// Search through everything in 32u_1 and 32u_2's bins for best match, combine, finalize cluster
		// Set each used centroid and all matches to {-1,-1} to indicate it's spent

		uint64_t numEmpty = 0;
		for (uint32_t i = 0; i < totR; ++i) numEmpty+=!CC[i];
		printf("There are %llu links (%f) and %llu empties (%f)\n",
			numLinks,(double)numLinks/totR,numEmpty,(double)numEmpty/totR);
		if (!numLinks) {fputs("ERROR: no qualifying clusters.\n",stderr); exit(4);}

		Split *SrtIx_cent = malloc(numLinks*sizeof(*SrtIx_cent));
		uint16_t *Connectivity = malloc(totR*sizeof(*Connectivity));
		for (uint32_t i = 0; i < totR; ++i) Connectivity[i] = 1;
		
		if (!SrtIx_cent) // || !SrtIx_y || !Y_bins) 
			{fputs("OOM:SrtIx_cent\n",stderr); exit(3);}
		
		printf("Time to make clusters: %f\n",omp_get_wtime()-wtime);
		wtime = omp_get_wtime();
		#pragma omp parallel sections
		{
			#pragma omp section
			{	// Create population queue
				uint64_t Popper[258] = {0}, *P_bins = Popper + 1;
				for (uint32_t i = 0; i < totR; ++i) if (CC[i]) {
					for (uint32_t j = 1; j <= CC[i][0].v; ++j)
						++P_bins[CC[i][j].v];
				}
				--P_bins;
				for (uint32_t i = 1; i <= 256; ++i) 
					P_bins[i] += P_bins[i-1];
				for (uint32_t i = 0; i < totR; ++i) if (CC[i])
					for (uint32_t j = 1; j <= CC[i][0].v; ++j)
						SrtIx_cent[P_bins[CC[i][j].v]++] = (Split){i,j};
			}
			#pragma omp section
			{
				for (uint32_t i = 0; i < totR; ++i) if (CC[i]) {
					uint64_t pop_tot = 0;
					for (uint32_t j = 1; j <= CC[i][0].v; ++j) {
						pop_tot += CC[i][j].v;
						Connectivity[CC[i][j].i] += Connectivity[CC[i][j].i] < UINT16_MAX;
					}
					Connectivity[i] = MIN(UINT16_MAX, (uint64_t)Connectivity[i] + CC[i][0].v);
				}
			}
		}
		printf("Time to sort: %f\n",omp_get_wtime()-wtime);
		// do the clustering
		wtime = omp_get_wtime();
		uint64_t c_ix = 0;
		uint32_t fp_ix = 0, *ClusIX = malloc(totR*sizeof(*ClusIX)); // to store the re-ordered references
		Prince *UnionPrince = malloc(totRC*sizeof(*UnionPrince));
		uint64_t totalPop = 0;
		for (uint64_t i = 0; i < numLinks; ++i) {
			uint32_t c1 = SrtIx_cent[i].v, y = SrtIx_cent[i].i;
			if (!Connectivity[c1]) continue; // first member is dead
			uint32_t c2 = CC[c1][y].i;
			if (!Connectivity[c2]) continue;
			Prince centroid = FP_combine(p[c1],p[c2]); // and now its pop is in "cost"
			ClusIX[c_ix++] = c1;
			ClusIX[c_ix++] = c2;
			Connectivity[c1] = Connectivity[c2] = 0;
			
			uint32_t min, mix;
			uint32_t start, end;
			if (c2 > BANDED_FP) start = c2 - BANDED_FP,	end = c2;
			else start = 0,	end = MIN(totR, BANDED_FP);
			for (int z = 0; z < (VECSZ-2); ++z) {
				min = 257; //, mix = -1
				#pragma omp parallel
				{
					uint32_t mp = -1, mi, bp = -1;
					#pragma omp for reduction(min:min)
					for (uint32_t j = start; j < end; ++j) if (Connectivity[j]) { 
						uint32_t t = FP_union(centroid, p[j]), tp;
						if (t < min) min = mp = t, mi = j, bp = FP_dist(p[j],centroid);
						else if (t == min && Connectivity[j] < Connectivity[mi]) 
							mi = j, bp = FP_dist(p[j],centroid);
						else if (t == min && Connectivity[j] == Connectivity[mi] && (tp=FP_dist(p[j],centroid)) < bp)
							mi = j, bp = tp;
					}
					IXBIN[omp_get_thread_num()] = mp == min ? mi : -1;
				}
				uint32_t bestC = -1, bestP = -1, t;
				for (int j = 0; j < THREADS; ++j) {
					if (IXBIN[j] != -1) {
						if (Connectivity[IXBIN[j]] < bestC) //(t=FP_dist(p[IXBIN[j]],centroid)) < bestP)  
							mix = IXBIN[j], bestC = Connectivity[mix], bestP = FP_dist(p[IXBIN[j]],centroid); 
						else if (Connectivity[IXBIN[j]] == bestC && (t=FP_dist(p[IXBIN[j]],centroid)) < bestP) 
							mix = IXBIN[j], bestP = t; 
						else if (Connectivity[IXBIN[j]] == bestC && FP_dist(p[IXBIN[j]],centroid) == bestP && IXBIN[j] < mix)
							mix = IXBIN[j]; 
					}
				}
				if (min > 256) { // Expand the reach -- backwards first	
					uint32_t j = start; 
					if (j--) do if (Connectivity[j]) {
						centroid = FP_combine(centroid, p[j]);
						ClusIX[c_ix++] = j;
						Connectivity[j] = 0;
						if (++z == VECSZ-2) break;
					} while (j--);
					if (j==-1) for (j = end+1; j < totR; ++j) if (Connectivity[j]) {
						centroid = FP_combine(centroid, p[j]);
						ClusIX[c_ix++] = j;
						Connectivity[j] = 0;
						if (++z == VECSZ-2) break;
					}
					if (z != VECSZ-2) {
						printf("\nCluster pool depleted"); 
						goto FP_ENDGAME; 
					}
					min = FP_pop(&centroid); 
					continue;
				}
				centroid = FP_combine(centroid, p[mix]);
				ClusIX[c_ix++] = mix;
				Connectivity[mix] = 0;
			}
			UnionPrince[fp_ix++] = centroid;
			totalPop += min;
			printf("\r[%f]: Finalizing cluster %u", (double)totalPop/fp_ix, fp_ix);
		}
		FP_ENDGAME:NULL;
		printf("\nTime to cluster: %f\n",omp_get_wtime()-wtime);
		printf("\nDone. Reads safely tucked into clusters: %u\n",c_ix);
		printf("Considering all clusters: %f\n",(double)(totalPop+256*(totRC-fp_ix))/totRC);

		if (c_ix < totR) for (uint32_t i = 0; i < totR; ++i) if (Connectivity[i]) ClusIX[c_ix++] = i;
		Prince badP; badP.w[0] = -1, badP.w[1] = -1, badP.w[2] = -1, badP.w[3] = -1;
		for (uint32_t i = fp_ix; i < totRC; ++i) UnionPrince[i] = badP;
		totalPop = 0;
		for (uint32_t i = 0; i < totRC; ++i) totalPop += FP_pop(UnionPrince + i);

		printf("Added remainder to the bin. Total avg: %f\n",(double)totalPop/totRC);
		for (uint32_t i = 0; i < totR; ++i) free(CC[i]);
		free(CC);
		free(Connectivity); // reuse Connectivity for new RefIxSrt, free old RefIxSrt
		free(SrtIx_cent);
		//#pragma omp sections
		{
			//#pragma omp section
			if (curate) { // generate new RefIxSrt, RefDedupIx, Original (TmpRIX)
				uint32_t *TmpRIX = Rd->TmpRIX, *NewOrig = malloc(Rd->origTotR*sizeof(*NewOrig)),
					*NewDedup = Rd->RefIxSrt, *RefDedupIx = Rd->RefDedupIx;
				if (!NewOrig) {fputs("OOM:NewOrig\n",stderr); exit(3);}
				uint32_t j = 0; for (uint32_t i = 0; i < totR; ++i) {
					NewDedup[i] = j;
					//printf("clusIX[%u] = %u\n",i, ClusIX[i]);
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
	if (!DO_FP) for (uint32_t i = 0; i < origR; ++i) free(origRefSeq[i]);
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
	
	for (uint32_t i = 0; i < Rd->origTotR; ++i) 
		totRefHeadLen += fprintf(output,"%s",Rd->RefHead[i]) + 1,
		fputc(0, output); 
	if (REBASE) 
		fwrite(Rd->RefStart, sizeof(*Rd->RefStart), Rd->origTotR, output);
	if (Rd->totR != Rd->origTotR) //puts("Dupewriting"), 
		fwrite(Rd->RefDedupIx, sizeof(*Rd->RefDedupIx), Rd->totR+1, output);

	// Write (sheared/deduped) clumps and mask into originals
	fwrite(Rd->TmpRIX, sizeof(*Rd->TmpRIX), Rd->origTotR, output); // can infer RefIxSrt
	fwrite(Rd->ClumpLen, sizeof(*Rd->ClumpLen), Rd->numRclumps, output);
	for (uint32_t i = 0; i < Rd->numRclumps; ++i)
		fwrite(Rd->RefClump[i], sizeof(*Rd->RefClump[i]), Rd->ClumpLen[i], output);
	
	if (DO_FP) { // write centroid and nf prints
		uint32_t nf = Rd->FingerprintsR.nf;
		fwrite(Rd->Centroids, sizeof(*Rd->Centroids), Rd->numRclumps, output);
		fwrite(&nf, sizeof(nf), 1, output);
		if (nf) fwrite(Rd->FingerprintsR.Ptrs, sizeof(*Rd->FingerprintsR.Ptrs), Rd->totR, output);
		else nf = Rd->totR; // if nf == 0, it means the FPs are unambiguous (N-penalized)
		fwrite(Rd->FingerprintsR.P, sizeof(*Rd->FingerprintsR.P), nf, output);
	}
	rewind(output);
	fputc(1 << 7 | REBASE << 6 | DO_FP << 5 | Xalpha << 4 | 
		(EDB_VERSION - Xalpha), output); // 5: xalpha, 6: ?
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
		fputs(" --> EDB: Parsing compressed headers\n",stderr);
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
#ifndef SCOUR_N
#define SCOUR_N 15
#endif
#define SCOUR_R (32 - 2*SCOUR_N)
#define SCOUR_L 4

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
	numUniqQ_rc = numUniqQ + (Qd->rc ? numUniqQ : 0);
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
	if (Qd->rc) {
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

	uint32_t *QBins = calloc(4,sizeof(*QBins));
	uint32_t  nq_clear, nq_ambig, nq_bad;
	if (DO_ACCEL) {
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
			if (len < SCOUR_N || Sb.ed >= len/SCOUR_N) QStat[i] = 2;
		
			else for (uint32_t j = 0; j < len; ++j) if (s[j] > 4) {
				QStat[i] = 0; // set as ambig-containing
				if (++totN > 5 || ++rowN > 3) {QStat[i] = 2; break;}
			} else rowN = 0;
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
	if (Qd->rc || (DO_ACCEL && nq_clear != numUniqQ_rc)) {
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
	Qd->totQ = totQ, Qd->numUniqQ = numUniqQ, Qd->maxLenQ = maxLenQ;
	Qd->UniBins = UniBins; Qd->ShrBins = ShrBins; Qd->maxED = maxED;
	Qd->QBins = QBins; Qd->maxDiv = maxDiv;
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
	//printf("Call: r = %u, w = %u, ix = %u\n", r, w, ix);
	if (ix == SCOUR_N) {
		if (F[w].len == F[w].cap) // check on success of realloc here?
			F[w].Refs = realloc(F[w].Refs,(F[w].cap+=F[w].cap+1)*sizeof(*F[w].Refs));
		if (!F[w].len || F[w].Refs[F[w].len-1] != r) F[w].Refs[F[w].len++] = r;
	}
	else for (int8_t i = 0; AMBIGS[S[ix]][i] != UINT8_MAX; ++i) 
		addAmbigScour(F, r, S, w << 2 | AMBIGS[S[ix]][i], ix + 1);
}
void make_accelerator(Reference_Data *Rd, char *xcel_FN) {
	FILE *out = fopen(xcel_FN,"wb");
	if (!out) {fprintf(stderr, "Cannot write accelerator '%s'\n",xcel_FN); exit(1);}
	char **RefSeq = Rd->RefSeq;
	uint32_t *RefIxSrt = Rd->RefIxSrt, *RefLen = Rd->RefLen, totR = Rd->totR,
		totRC = Rd->numRclumps;

	//AccelNode **Forest = calloc(1 << (2*SCOUR_N), sizeof(*Forest));
	Accelerant *F = calloc(1 << (2*SCOUR_N), sizeof(*F));
	uint32_t szBL = 1000, badIX = 0, *BadList = malloc(szBL*sizeof(*BadList));
	if (!F || !BadList) {fputs("OOM:BadList\n",stderr); exit(3);}
	double wtime = omp_get_wtime();
	if (Z) fprintf(stderr,"WARNING: N-penalized accelerator not usable for unpenalized alignment\n");
	
	for (uint32_t i = 0; i < totRC; ++i) {
		printf("Generating accelerator [%u / %u]\r",i,totRC);
		
		uint32_t begin = i*VECSZ, end = MIN(totR,begin + VECSZ); 
		uint16_t doAmbig = 0;

		for (uint32_t z = begin; z < end; ++z) { // specific sequence
			char *s = RefSeq[RefIxSrt[z]]; 
			uint32_t len = RefLen[RefIxSrt[z]], rowN = 0;
			if (len < SCOUR_N) continue;
			
			// Perfect candidate for vectorization?
			for (uint32_t j = 0; j < len; ++j) if (s[j] > 4 + Z) {
				if (s[j]==5) {
					if (++rowN > SCOUR_L) {
						if (badIX == szBL) {
							BadList = realloc(BadList, (szBL*=2)*sizeof(*BadList));
							if (!BadList) {fputs("OOM:BadList2\n",stderr); exit(3);}
						}
						//if (!badIX || BadList[badIX-1] != i) 
						BadList[badIX++] = i;
						goto END_ACCEL_LP;
					}
				}
				else rowN = 0;
				doAmbig |= 1 << (z-begin); // per-sequence control of ambigging
			} else rowN = 0;
		}
		
		for (uint32_t z = begin; z < end; ++z) { 
			char *s = RefSeq[RefIxSrt[z]]; 
			uint32_t len = RefLen[RefIxSrt[z]];
			if (len < SCOUR_N) continue;
			if (Z) { // screen for N's first!
				for (uint32_t j = 0; j + SCOUR_N <= len; ++j) {
					uint32_t k = 0; 
					for (; k < SCOUR_N; ++k) if (s[j+k] == 5) break;
					if (k < SCOUR_N) { j+= k; continue; }
					addAmbigScour(F, i, s + j, 0, 0);
				}
			}
			else if (doAmbig << (16-(z-begin)) >> (z-begin)) //addAmbigScourLong (F, i, s, 0, 0, len);
				for (uint32_t j = 0; j + SCOUR_N <= len; ++j)
					addAmbigScour(F, i, s + j, 0, 0);
			else {
				uint32_t w = 0;
				for (uint32_t j = 0; j + 1 < SCOUR_N; ++j)
					w <<= 2, w |= s[j]-1;
				for (uint32_t j = SCOUR_N - 1; j < len; ++j) {
					w <<= 2, w |= s[j]-1;
					uint32_t t = (w << SCOUR_R) >> SCOUR_R;
					//printf("Word = %u [in %u]\n",t,i);
					if (F[t].len == F[t].cap) // check on success of realloc here?
						F[t].Refs = realloc(F[t].Refs,(F[t].cap+=F[t].cap+1)*sizeof(*F[t].Refs));
					if (!F[t].len || F[t].Refs[F[t].len-1] != i) F[t].Refs[F[t].len++] = i;
				}
			}
		}
		END_ACCEL_LP:NULL;
	}
	printf("Done calculating accelerator (%f s), %u bad\n",
		omp_get_wtime()-wtime, badIX);
	wtime = omp_get_wtime();
	size_t totalWords = 0;
	for (uint32_t i = 0; i < 1 << (2*SCOUR_N); ++i) 
		totalWords += F[i].len;
	printf("Total accelerants stored: %llu\n",totalWords);
	//BadList = realloc(BadList,badIX*sizeof(*BadList));
	uint8_t vers = (1 << 7) | (Z << 6) | SCOUR_N;
	fwrite(&vers,sizeof(vers),1,out);
	fwrite(&badIX, sizeof(badIX), 1, out);
	for (uint32_t i = 0; i < 1 << (2*SCOUR_N); ++i) 
		fwrite(&F[i].len, sizeof(F[i].len), 1, out);
	for (uint32_t i = 0; i < 1 << (2*SCOUR_N); ++i) 
		fwrite(F[i].Refs, sizeof(*F[i].Refs), F[i].len, out);
	fwrite(BadList, sizeof(*BadList), badIX, out);
	printf("Wrote accelerator (%f).\n",omp_get_wtime()-wtime);
	
	//for (uint32_t i = 0; i < 1 << (2*SCOUR_N); ++i) free(F[i].Refs);
	free(F), free(BadList);
}

void read_accelerator(Reference_Data *Rd, char *xcel_FN) {
	FILE *in = fopen(xcel_FN,"rb");
	if (!in) {fprintf(stderr, "Cannot read accelerator '%s'\n",xcel_FN); exit(1);}
	
	uint8_t cb = fgetc(in), 
		dbVer = (uint8_t)(cb << 4) >> 4,
		didZ = (uint8_t)(cb << 1) >> 7; 
	if (didZ && !Z) {
		fprintf(stderr,"ERROR: Accelerator built with '-n'; must use '-n'\n");
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

static inline void do_alignments(FILE *output, Reference_Data RefDat, Query_Data QDat) {
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

	// Prepare dimensions, bounds
	uint32_t qdim = maxLenQ + 1, rdim = maxLenR + 1;
	++qdim; ++rdim; // add padding to right and bottom to eliminate loop logic

	uint32_t newUniqQ = numUniqQ + (doRC ? numUniqQ : 0);
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
	
	Split *RefOrder = 0;
	if (DO_PREPASS == AUTO && THRES >= 0.99f) DO_PREPASS = NONE, 
		printf("Auto-prepass: DISABLING prepass.\n");
	if (DO_FP && DO_PREPASS && maxED) {
		int reorder_only = DO_PREPASS == REORDER, fast_prepass = DO_PREPASS == FAST;
		if (DO_PREPASS == AUTO) fast_prepass = THRES >= 0.93f, reorder_only = THRES >= 0.95f,
			printf("Auto-prepass: enabling '%s' prepass.\n", 
				reorder_only ? "reorder-only" : fast_prepass ? "fast" : "full");
		uint32_t curtail = (maxED*1/(1+fast_prepass)) < 32 ? (maxED*1/(1+fast_prepass)) : 32,
			curtail2 = fast_prepass ? MIN(2,curtail/2) : curtail/2, curtail3 = curtail2/2;
		uint32_t num_ties_prepass = fast_prepass ? curtail : 128;

		double wtime = omp_get_wtime();
		
		RefOrder = malloc(numRclumps*sizeof(*RefOrder));
		if (!RefOrder) {fputs("OOM:RefOrder\n",stderr); exit(3);}
		for (uint32_t i = 0; i < numRclumps; ++i) RefOrder[i] = (Split){0,i};

		Prince *p_r = FpR.P, *p_q = Fp.P;
		uint8_t *RPops = malloc(totR*sizeof(*RPops));
		uint32_t inc = reorder_only ? 16 : 1;
		#pragma omp parallel
		{
			DualCoil *Matrices = malloc(2*rdim*sizeof(*Matrices));
			uint32_t *RealQBest = malloc(numRclumps*sizeof(*RealQBest));
			if (!Matrices || !RealQBest) {fputs("OOM:Prepass\n",stderr); exit(3);}
			DualCoil *rclump = malloc((2+rdim)*sizeof(*rclump));
			uint32_t MinBin[num_ties_prepass];
			#pragma omp for 
			for (uint32_t i = 0; i < totR; ++i)
				RPops[i] = MIN(255,FP_pop(p_r + i));
			
			#pragma omp for schedule(dynamic,1)
			for (uint32_t i = 0; i < newUniqQ; i+=inc) {
				uint32_t qmax = 0; //, pmin = 257, ri = -1, rp;
				uint32_t ties = 0;
				for (uint32_t j = 0; j < totR; ++j) if (RPops[j] < 250) {
					uint32_t t = FP_intersect(p_r[j],p_q[i]);
					if (t > qmax) qmax = t, ties = 1, RealQBest[0] = j/16;
					else if (t == qmax && qmax && RealQBest[ties-1] != j/16) RealQBest[ties++] = j/16;
				}
				// Sort first by XOR (or pop?)
				if (reorder_only) for (uint32_t j = 0; j < ties; ++j)
					#pragma omp atomic
					++RefOrder[RealQBest[j]].v;
				else {
					uint32_t Pops[258] = {0}, *P_bins = Pops + 1; 
					for (uint32_t j = 0; j < ties; ++j) {
						#pragma omp atomic
						++RefOrder[RealQBest[j]].v;
						++P_bins[FP_dist(p_q[i],p_r[RealQBest[j]])];
					}
					--P_bins;
					for (uint32_t j = 1; j < 258; ++j) P_bins[j] += P_bins[j-1];
					for (uint32_t j = 0; j < ties; ++j) {
						uint32_t bin = FP_dist(p_q[i],p_r[RealQBest[j]]); 
						if (P_bins[bin] < num_ties_prepass) MinBin[P_bins[bin]++] = RealQBest[j];
					}
					ties = MIN(ties, num_ties_prepass); // limit the alignments per query to this number
					uint32_t fp_err = Fp.N[i] - qmax;
					uint32_t lastChange = 0, stop = curtail; //, stop2 = curtail2*2; ///2;
					UniBin Ub = UniBins[i]; ShrBin Sb = ShrBins[Ub.six];
					if (Sb.ed > fp_err) for (uint32_t j = 0; j < ties; ++j) {
						uint32_t ri = MinBin[j]; // or RealQBest[j] if no min-binning
						DualCoil *RefSlide = RefClump[ri];
						for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
							__m128i org = _mm_stream_load_si128((void*)(RefSlide++));
							__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
							__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
							_mm_store_si128((void*)(rclump+w),ex1);
							_mm_store_si128((void*)(rclump+w+1),ex2);
						}
						uint32_t min = prune_ed_mat16(rclump, Ub.Seq, ClumpLen[ri],
							Sb.len, rdim, Matrices, ProfClump[ri], Sb.ed);
						if (min < Sb.ed) {
							Sb.ed = min;
							lastChange = j;
							stop = stop==curtail2 ? curtail3 : curtail2;
						}
						if (j - lastChange >= stop) break; 
					}
				}
			}
			free(RealQBest);
			free(Matrices);
		}
		free(RPops);

		void pivSortPopD(Split *A, uint64_t len, uint64_t min, uint64_t max, uint32_t depth) { // WARNING: len MUST BE >= 1
			uint64_t lp = 0, hp = len-1, pivot = (max - min)/2 + min;
			Split temp;
			while (lp < hp) if (A[lp].v > pivot) ++lp;
				else temp = A[lp], A[lp] = A[hp], A[hp--] = temp;
			lp += A[lp].v > pivot;
			if (lp) 
				#pragma omp task final(depth > 24) mergeable
				pivSortPopD(A, lp,  pivot+1, max, depth+1);
			if (lp < len-1 && pivot > min)
				pivSortPopD(A+lp, len-lp, min, pivot, depth+1);
		}
		if (DO_ACCEL) free(RefOrder), RefOrder = 0; 
		else pivSortPopD(RefOrder, numRclumps, 0, newUniqQ, 0);

		printf("Time to perform prepass: %f\n",omp_get_wtime()-wtime);
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

	//uint32_t *UniqQLen, *UniqQed, *UniqDiv, *NewIX, maxDiv, *RCMap, *Umap;
	//char **UniqQSeq;

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
	
	printf("Using ACCELERATOR to align %u unique queries...\n", QBins[1]);
	#pragma omp parallel
	{
		DualCoil *Matrices = calloc(PADDING+(cacheSz+2)*rdim+PADDING,sizeof(*Matrices)) + PADDING,
				 *ScoresEX = malloc((PADDING+2*rdim+PADDING)*sizeof(*ScoresEX)) + PADDING,
				 *ShiftsEX = malloc((PADDING+2*rdim+PADDING)*sizeof(*ShiftsEX)) + PADDING,
				 *ShiftsBX = malloc((PADDING+2*rdim+PADDING)*sizeof(*ShiftsEX)) + PADDING;
		uint32_t *HiBound = calloc(PADDING+qdim+PADDING,sizeof(*HiBound)) + PADDING,
				 *LoBound = calloc(PADDING+qdim+PADDING,sizeof(*LoBound)) + PADDING,
				 *StackE = malloc(qdim*sizeof(*StackE)), // num errors
				 *StackX = malloc(qdim*sizeof(*StackX)); // query ix in UniqQSeq
		if (!(Matrices && ScoresEX && ShiftsEX && ShiftsBX && HiBound && LoBound && StackE && StackX)) 
			{fputs("OOM:Aligner\n",stderr); exit(3);}
		for (int j = 0; j < (cacheSz+2); ++j) Matrices[j*rdim].v = _mm_set1_epi8(MIN(j*GAP,255));
		*LoBound = -1, LoBound[1] = 1;
		
		//ResultPod **Pods = calloc(numUniqQ,sizeof(*Pods));
		//if (!Pods) {fputs("OOM:Pods\n",stderr); exit(3);}
		
		//for (size_t i = 0; i < numUniqQ; ++i) 
		//	Pods[i] = &PodBase[i],
		//	PodBase[i].mismatches = -1;

		DualCoil *rclump = malloc((2+rdim)*sizeof(*rclump));

		uint16_t *Hash = calloc(numRclumps,sizeof(*Hash)), hqd = qdim;
		uint32_t *Cache = malloc(numRclumps*QBUNCH*sizeof(*Cache)); // to store all hit refs
		Split *Word_Ix = malloc((qdim+1024)*QBUNCH*sizeof(*Word_Ix)); // wordfill, sort, "dedupe", add refs in burst, delete

		//uint16_t *Srtr = calloc(numRclumps,sizeof(*Srtr));
		Split *Refs = calloc(numRclumps,sizeof(*Refs)), *RefPtr = Refs; 
		if (!Hash || !Refs) {fputs("OOM:Hash\n",stderr); exit(3);}
		uint32_t oldNref = 0;
		int tid = omp_get_thread_num();
		
		#pragma omp for schedule(dynamic,1)
		for (uint32_t z = 0; z < QBins[1]; z+= QBUNCH) {
			//memset(Hash,0,numRclumps*sizeof(*Hash));
			for (uint32_t j = 0; j < oldNref; ++j) Refs[j].i = 0; //Bucket[j] = 0; // Bucket memset
			uint32_t bound = MIN(z+QBUNCH,QBins[1]), nref = 0, min_mmatch = UINT32_MAX;
			uint32_t wix = 0, cix = 0;
			for (uint32_t j = z; j < bound; ++j) {
				UniBin Ub = UniBins[j]; ShrBin Sb = ShrBins[Ub.six];
				uint32_t len = Sb.len, err = Sb.ed;
				char *s = Ub.Seq; 
				uint32_t mmatch = len - SCOUR_N * err - SCOUR_N; // len - (err + 1) * SCOUR_N;
				if (mmatch < min_mmatch) min_mmatch = mmatch; // bank worst-case threshold
				//printf("Query %u: len = %u, err = %u, runsize = %u, mmatch = %u\n", j, len, err, (len - err)/(err+1), mmatch);
				if (j >= *QBins) {//storeUnambigWords(Word_Ix, s, len, &wix, j);
					uint32_t w = *s - 1;
					for (uint32_t k = 1; k + 1 < SCOUR_N; ++k)
						w <<= 2, w |= s[k]-1;
					for (uint32_t k = SCOUR_N-1; k < len; ++k) {
						w <<= 2, w |= s[k]-1;
						uint32_t t = (w << SCOUR_R) >> SCOUR_R;
						Word_Ix[wix++] = (Split){t,j};
					}
				}
				else for (uint32_t k = 0; k + SCOUR_N <= len; ++k) 
					storeAmbigWords(Word_Ix, s + k, &wix, j, 0, 0);
				/* if (j >= *QBins) 
					setUnambigScour(Hash, Refs, &nref, mmatch, Forest, s, qdim, len, Cache, &cix);
				else for (uint32_t k = 0; k + SCOUR_N <= len; ++k)
					setAmbigScour(Hash, Refs, &nref, mmatch, Forest, s + k, qdim, 0, 0, Cache, &cix);
				for (uint32_t i = 0; i < nref; ++i) // bank worst-case limits
					Refs[i].i = Refs[i].i >= Hash[Refs[i].v] ? Refs[i].i : Hash[Refs[i].v];
				//memset(Hash,0,numRclumps*sizeof(*Hash));
				//if (j+1 < bound) for (uint32_t i = 0; i < nref; ++i) Hash[Refs[i].v] = hqd;
				//if (j >= *QBins) 
					for (uint32_t i = cix; i; --i) Hash[Cache[i-1]] = 0;
				//else memset(Hash,0,numRclumps*sizeof(*Hash));
				if (j+1 < bound) for (uint32_t i = 0; i < nref; ++i) Hash[Refs[i].v] = hqd; */
			}
			//if (wix >= (qdim+1024)*QBUNCH) {puts("\n\nERROR: WIX"); exit(4);}
			//sort the Word_Ix pairs, then feed to big reference hash table
			qsort(Word_Ix, wix, sizeof(*Word_Ix), WrdCmp3);
			/* #pragma omp critical
			{
				printf("Words and which q:\n");
				for (uint32_t i = 0; i < wix; ++i) {
					printf(" [%u] %u (%u)\n",i,Word_Ix[i].v, Word_Ix[i].i); 
				}
			} */
			postScour(Hash, Forest, Word_Ix, wix, Cache, &cix);
			for (uint32_t i = 0; i < cix; ++i) {
				uint16_t *p = Hash + Cache[i];
				if (*p > min_mmatch) Refs[nref++] = (Split){Cache[i],*p};
				*p = 0;
			}
			//for (uint32_t i = 0; i < cix; ++i) Hash[Cache[i]] = 0;
			//memset(Hash,0,numRclumps*sizeof(*Hash));
			
			
			//for (uint32_t i = 0; i < cix; ++i) Hash[Cache[i]] = 0; // New memset replacement
			//memset(Hash,0,numRclumps*sizeof(*Hash));
			         //uint16_t maxK = 0; uint32_t maxI=0, tmp, minI=-1;
			//for (uint32_t i = 0; i < nref; ++i) Refs[i].i -= qdim; 
			
			
			iSort(Refs, nref);
			/* #pragma omp critical
			{
				printf("Words and which q:\n");
				for (uint32_t i = 0; i < nref; ++i) {
					printf(" [%u] %u (%u)\n",i,Refs[i].v, Refs[i].i); 
				}
			} */
			
			RefPtr = Refs;
			oldNref = nref;
			
			// Run the loop 2x: 1) good R's, 2) bad R's
			for (int x = 0; x < 2; ++x) {
				for (uint32_t i = 0; i < nref; ++i) {
					uint32_t ri = RefPtr ? RefPtr[i].v : BadList[i];
					uint32_t rlen = ClumpLen[ri] + 1;
					DualCoil *pclump = ProfClump[ri];
					DualCoil *RefSlide = RefClump[ri];
					for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
						__m128i org = _mm_stream_load_si128((void*)(RefSlide++));//_mm_load_si128((void*)(RefSlide++));
						__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
						__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
						_mm_store_si128((void*)(rclump+w),ex1);
						_mm_store_si128((void*)(rclump+w+1),ex2);
					}
					HiBound[1] = rlen;
					//*LoBound = -1, LoBound[1] = 1; // FDR
					
					uint32_t stack_ix = 0; //thisMax = 0;
					*StackX = z, *StackE = -1; // initialize stack (worst case)
					UniBins[z].div = 1; // force first divergence to 1

					int8_t fp_rediv = 0;
					for (uint32_t j = z; j < bound; ++j) { 
						UniBin *Ub = UniBins + j; ShrBin *Sb = ShrBins + Ub->six;
						uint32_t thisDiv = Ub->div, Emac = Sb->ed, len = Sb->len;
						//uint32_t Emac = maxED;
						//printf("Doing q #%u [%s]: div=%u, Emac = %u\n",j,QHead[Offset[Ub->six]],thisDiv,Emac);

						// handle early k-term here
						uint32_t mmatch = len - Emac * SCOUR_N - SCOUR_N;
						if (!x && Refs[i].i <= mmatch) {
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
						min = aded_mat16(rclump,Ub->Seq,rlen,len, rdim, 
							Matrices, pclump, Emac,thisDiv,LoBound,HiBound,&mins); 
						//printf("--> min = %u\n",min);
						if (min <= Sb->ed) { // now we get serious.
							uint32_t ai = Ub->six + (Ub->rc ? numUniqQ : 0);
							if (RUNMODE != FORAGE) {
								/*if (min < Sb->ed)*/ Sb->ed = min;
								ResultPod *prv;
								if (Pods[ai] && min < Pods[ai]->mismatches) // purge
									do prv = Pods[ai]->next, free(Pods[ai]), Pods[ai] = prv; while(Pods[ai]);
							} else min = Emac; // all valid refs can be explored
							MetaPack MPK;
							reScoreM_mat16(rclump,Ub->Seq,rlen,len, rdim, ScoresEX, ShiftsEX,
								ShiftsBX, min, pclump,&MPK);
							for (int z = 0; z < VECSZ; ++z) { 
								if (mins.u8[z] > min) continue; // skip non-mins
								ResultPod *tmp = malloc(sizeof(*tmp));
								tmp->next = Pods[ai]; Pods[ai] = tmp;
								
								Pods[ai]->mismatches = mins.u8[z];
								Pods[ai]->score = MPK.score[z]; //-1.f; // placeholder
								Pods[ai]->refIx = ri * VECSZ + z;
								Pods[ai]->finalPos = MPK.finalPos[z];
								Pods[ai]->numGapR = MPK.numGapR[z];
								Pods[ai]->numGapQ = MPK.numGapQ[z];
								Pods[ai]->rc = Ub->rc;
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
			if (!tid) printf("\rSearch Progress: [%3.2f%%]",100.0 * (double)totDone * totQMult);
			
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
	
	// this needs a re-write. I need to go through the whole list and prune pods with more errors than the thres.
	
	/* #pragma omp parallel for schedule(dynamic,1)
	for (uint32_t j = 0; j < numUniqQ; ++j) {
		//if (!BasePod[j]) continue;
		ResultPod *prv;
		while (BasePod[j] && BasePod[j]->mismatches > ShrBins[j].ed) 
			BasePod[j] = BasePod[j]->next; // freeing version below \/
			//prv = BasePod[j]->next, free(BasePod[j]), BasePod[j] = prv; 
			//while(BasePod[j]);
	} */
	/*for (uint32_t i = 1; i < THREADS; ++i) {
		#pragma omp parallel for schedule(dynamic,1)
		for (uint32_t j = 0; j < numUniqQ; ++j) {
			if (!ThreadPods[i][j]) continue;
			ResultPod *prv, *init;
			if (ThreadPods[i][j]->mismatches > ShrBins[j].ed) {
				do prv = ThreadPods[i][j]->next, free(ThreadPods[i][j]), ThreadPods[i][j] = prv; 
				while(ThreadPods[i][j]);
			} else {
				init = prv = ThreadPods[i][j];
				while (prv->next) prv = prv->next;
				prv->next = BasePod[j];
				BasePod[j] = init;
			}
		}
		free(ThreadPods[i]);
	}*/
	//printf("\rSearch Progress: [100.00%%]\n");
	printf("\rSearch Progress: [%3.2f%%]\n",100.0 * (double)totDone * totQMult); // not necessarily 100.0%!
	free(*Forest); free(Forest); free(BadList); // NEW [mem]
} // end ACCEL

	// do remaining queries vs all refs in traditional aligner. 
	uint32_t firstQ = DO_ACCEL ? QBins[1] : 0;
	totDone = 0;
	if (firstQ == newUniqQ || (DO_ACCEL && skipAmbig)) goto EOA;
	else if (DO_ACCEL) UniBins[firstQ].div = 1; //UniqDiv[Umap[firstQ]] = 1;
	printf("Searching best paths through %u unique queries...\n", newUniqQ-firstQ);
	#pragma omp parallel
	{
		DualCoil *Matrices = calloc(PADDING+(cacheSz+2)*rdim+PADDING,sizeof(*Matrices)) + PADDING,
				 *ScoresEX = malloc((PADDING+2*rdim+PADDING)*sizeof(*ScoresEX)) + PADDING,
				 *ShiftsEX = malloc((PADDING+2*rdim+PADDING)*sizeof(*ShiftsEX)) + PADDING,
				 *ShiftsBX = malloc((PADDING+2*rdim+PADDING)*sizeof(*ShiftsEX)) + PADDING;
		uint32_t *HiBound = calloc(PADDING+qdim+PADDING,sizeof(*HiBound)) + PADDING,
				 *LoBound = calloc(PADDING+qdim+PADDING,sizeof(*LoBound)) + PADDING,
				 *StackE = malloc(qdim*sizeof(*StackE)), // num errors
				 *StackX = malloc(qdim*sizeof(*StackX)); // query ix in UniqQSeq
		for (int j = 0; j < (cacheSz+2); ++j) Matrices[j*rdim].v = _mm_set1_epi8(MIN(j*GAP,255));
		*LoBound = -1, LoBound[1] = 1;
		
		ResultPod **Pods = calloc(numUniqQ,sizeof(*Pods));
		if (!Pods) {fputs("OOM:Pods\n",stderr); exit(3);}
		/* ResultPod *PodBase = calloc(numUniqQ,sizeof(*PodBase));
		ResultPod **Pods = malloc(numUniqQ*sizeof(*Pods));
		for (size_t i = 0; i < numUniqQ; ++i) 
			Pods[i] = &PodBase[i],
			PodBase[i].mismatches = -1; */
		DualCoil *rclump = malloc((2+rdim)*sizeof(*rclump));
			
		#pragma omp for schedule(dynamic,1) 
		for (uint32_t i = 0; i < numRclumps; ++i) {
			//*LoBound = -1, LoBound[1] = 1; // FDR
			uint32_t ri = RefOrder ? RefOrder[i].i : i;
			uint32_t rlen = ClumpLen[ri] + 1;
			DualCoil *pclump = ProfClump[ri];
			DualCoil *RefSlide = RefClump[ri];
			for (uint32_t w = 0; w < ClumpLen[ri]; w += 2) {
				__m128i org = _mm_stream_load_si128((void*)(RefSlide++));
				__m128i ex1 = _mm_and_si128(org,_mm_set1_epi8(0xF));
				__m128i ex2 = _mm_and_si128(_mm_srli_epi16(org,4),_mm_set1_epi8(0xF));
				_mm_store_si128((void*)(rclump+w),ex1);
				_mm_store_si128((void*)(rclump+w+1),ex2);
			}
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
				//printf("[%u,qi=%u] UniqDiv = %u, err = %u\n",j,qi,UniqDiv[qi], UniqQed[qi]);
				
				uint32_t Emac = Sb->ed; //UniqQed[ai];

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

					if (fp_rediv) {
						fp_rediv = 0;
						if (Emac <= StackE[stack_ix]) {
							thisDiv = 1; if (stack_ix) {
								uint32_t sai = StackX[stack_ix], lenD;
								//sai = sai >= numUniqQ ? RCMap[sai-numUniqQ] : sai;
								lenD = MIN(len,ShrBins[UniBins[sai].six].len) - 1;
								register char *thisQ = Ub->Seq, *prevQ = UniBins[sai].Seq; //UniqQSeq[StackX[stack_ix]];
								for (uint32_t w = 0; w < lenD && thisDiv < maxDiv && 
									thisQ[w] == prevQ[w]; ++w) ++thisDiv; 
							}
						}
						//Emac = StackE[stack_ix ? stack_ix - 1 : stack_ix];
					}
				}

				// if error tol exceed cur stack, increment stack, set this ptr
				if (Emac > StackE[stack_ix]) {
					while (Emac > StackE[--stack_ix]); // pop until errors <= Stack errors
					thisDiv = 1; if (stack_ix) {
						uint32_t sai = StackX[stack_ix], lenD;
						//sai = sai >= numUniqQ ? RCMap[sai-numUniqQ] : sai;
						lenD = MIN(len,ShrBins[UniBins[sai].six].len) - 1;
						register char *thisQ = Ub->Seq, *prevQ = UniBins[sai].Seq; //UniqQSeq[StackX[stack_ix]];
						for (uint32_t w = 0; w < lenD && thisDiv < maxDiv && 
							thisQ[w] == prevQ[w]; ++w) ++thisDiv; 
					}
					/* register char *thisQ = UniqQSeq[qi], *prevQ = UniqQSeq[StackX[stack_ix]];
					//thisDiv = 1; if (stack_ix) while (*thisQ++ == *prevQ++ && ++thisDiv < maxDiv);
					thisDiv = 1; if (stack_ix) 
						for (uint32_t w = 0; w < UniqQLen[qi] && w < UniqQLen[StackX[stack_ix]] && 
							thisDiv < maxDiv && thisQ[w] == prevQ[w]; ++w) ++thisDiv;  */
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
				
				if (min <= Sb->ed) { // now we get serious.
					//uint32_t ai = Ub->six;
					if (RUNMODE != FORAGE) {
						if (min < Sb->ed) Sb->ed = min;
						ResultPod *prv;
						if (Pods[ai] && min < Pods[ai]->mismatches) // purge
							do prv = Pods[ai]->next, free(Pods[ai]), Pods[ai] = prv; while(Pods[ai]);
					} else min = Emac; // all valid refs can be explored
					MetaPack MPK;
					if (Xalpha) reScoreM_xalpha(rclump,Ub->Seq,rlen,len, rdim, ScoresEX, ShiftsEX,
						ShiftsBX, min, pclump,&MPK);
					else reScoreM_mat16(rclump,Ub->Seq,rlen,len, rdim, ScoresEX, ShiftsEX,
						ShiftsBX, min, pclump,&MPK);
					for (int z = 0; z < VECSZ; ++z) { 
						if (mins.u8[z] > min) continue; // skip non-mins
						ResultPod *tmp = malloc(sizeof(*tmp));
						tmp->next = Pods[ai]; Pods[ai] = tmp;
						
						Pods[ai]->mismatches = mins.u8[z];
						Pods[ai]->score = MPK.score[z]; //-1.f; // placeholder
						Pods[ai]->refIx = ri * VECSZ + z;
						Pods[ai]->finalPos = MPK.finalPos[z];
						Pods[ai]->numGapR = MPK.numGapR[z];
						Pods[ai]->numGapQ = MPK.numGapQ[z];
						Pods[ai]->rc = Ub->rc;
					}
				}
				/*if (min <= UniqQed[ai]) { // now we get serious.
					if (RUNMODE != FORAGE) {
						if (min < UniqQed[ai]) UniqQed[ai] = min;
						ResultPod *prv;
						if (Pods[ai] && min < Pods[ai]->mismatches) // purge
							do prv = Pods[ai]->next, free(Pods[ai]), Pods[ai] = prv; while(Pods[ai]);
					} else min = Emac; // all valid refs can be explored
					MetaPack MPK;
					reScoreM(rclump,UniqQSeq[qi],rlen,UniqQLen[ai], rdim, ScoresEX, ShiftsEX,
						ShiftsBX, min, pclump,&MPK);
					for (int z = 0; z < VECSZ; ++z) { 
						if (mins.u8[z] > min) continue; 
						ResultPod *tmp = malloc(sizeof(*tmp));
						tmp->next = Pods[ai]; Pods[ai] = tmp;
						
						Pods[ai]->mismatches = mins.u8[z];
						Pods[ai]->score = MPK.score[z];
						Pods[ai]->refIx = ri * VECSZ + z;
						Pods[ai]->finalPos = MPK.finalPos[z];
						Pods[ai]->numGapR = MPK.numGapR[z];
						Pods[ai]->numGapQ = MPK.numGapQ[z];
						Pods[ai]->rc = ai != qi;
					}
				}*/
			}
			#pragma omp atomic
			++totDone;
			tid = omp_get_thread_num();
			if (!tid) printf("\rSearch Progress: [%3.2f%%]",100.0 * (double)totDone / numRclumps);
		}
		/* #pragma omp master
		printf("\rSearch Progress: [100.00%%]\n");
		#pragma omp critical
		{
			//#pragma omp parallel for
			for (uint32_t i = 0; i < numUniqQ; ++i) {
				Pods[i]->next = FinalPod[i];
				FinalPod[i] = &PodBase[i];
			}
		} */
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
	printf("\rSearch Progress: [100.00%%]\n");
	EOA:NULL;
	printf("Search complete. Consolidating results...\n");

	//if (totSkipped) printf("NumSkipped = %llu (%f)\n",
	//	totSkipped,(double)totSkipped/((double)numUniqQ*numRclumps));
	if (!RefIxSrt && RefDedupIx) {
		RefIxSrt = malloc(totR * sizeof(*RefIxSrt));
		if (!RefIxSrt) {fputs("OOM:[DA]RefIxSrt\n",stderr); exit(3);}
		for (uint32_t i = 0; i < totR; ++i) RefIxSrt[i] = TmpRIX[RefDedupIx[i]];
	}
	else if (!RefIxSrt && !RefDedupIx) RefIxSrt = TmpRIX;

	free(Centroids); free(FpR.initP); free(ProfClump); free(RefClump);
	free(SeqDumpRC); // free a part of the query memory

	// Remove duplicate alignments
	uint32_t *RefMap = RefDat.RefMap;
	if (RUNMODE != BEST) {
		typedef struct PosP PosP;
		struct PosP {uint32_t i; PosP *p;};
		#pragma omp parallel
		{
			PosP **Grid = calloc(RefDat.numRefHeads,sizeof(*Grid));
			#pragma omp for schedule(dynamic,1)
			for (uint32_t i = 0; i < numUniqQ; ++i) {
				ResultPod *rp = BasePod[i], *prev = 0, *next;
				while (rp) {
					uint32_t rix = RefIxSrt[rp->refIx],
					map = RefMap[rix],
					mOff = RefStart ? RefStart[rix] : 0,
					stIxR = rp->finalPos - ShrBins[i].len + rp->numGapR + mOff;
					ResultPod *next = rp->next;

					// if (map >= RefDat.numRefHeads) {
					// 	printf("\nERROR: map is out of range: query %u; map = %u; rng = %u\n",i,map,RefDat.numRefHeads);
					// 	printf("-->rix = %u, Head = %s\n",rix,RefHead[rix]);
					// 	rp = next; continue;
					// 	//exit(4);
					// }
					PosP *pt = Grid[map], *pto = 0;
					while (pt) {
						if (stIxR + (DB_QLEN>>1) > pt->i && stIxR < pt->i + (DB_QLEN>>1)) {
							prev->next = next;
							free(rp); rp = 0;
							break;
						}
						pto = pt, pt = pt->p;
					} 
					if (!rp) { rp = next; continue; }
					PosP *newPos = malloc(sizeof(*newPos));
					*newPos = (PosP){stIxR,0};
					if (pto) pto->p = newPos;
					else Grid[map] = newPos;
					prev = rp; rp = next;
				}
				rp = BasePod[i];
				while (rp) {
					uint32_t map = RefMap[RefIxSrt[rp->refIx]];
					PosP *tmp;
					while (Grid[map]) 
						tmp = Grid[map]->p,
						free(Grid[map]),
						Grid[map] = tmp;
					rp = rp->next;
				}
				//memset(Grid,0,RefDat.numRefHeads*sizeof(*Grid));
			}
			free(Grid);
		}
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

#define SPRINT_MATCH() \
	snprintf(LineDmp+sptr,LNL,"%s\t%s\t%f\t%u\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\n", \
		QHead[j], RefHead[rix], rp->score * 100, \
		alLen, numMis, numGap, 1, qlen, stIxR, edIxR, \
		rp->mismatches,j > Offset[i]);
#define SPRINT_MATCH_TAX() \
	snprintf(LineDmp+sptr,LNL,"%s\t%s\t%f\t%u\t%u\t%u\t%u\t%u\t%d\t%u\t%u\t%u\t%s\n", \
		QHead[j], RefHead[rix], rp->score * 100, \
		alLen, numMis, numGap, 1, qlen, stIxR, edIxR, \
		rp->mismatches,j > Offset[i],FinalTaxon);

	if (RUNMODE == ALLPATHS) {  // all on best ED path
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i];
			if (!rp) continue;
			ResultPod *best = rp;
			while (rp = rp->next) 
				if (rp->mismatches < best->mismatches) best = rp;
			rp = best;
			uint32_t bm = best->mismatches;
			if (rp->score) while (rp) {
				if (rp->mismatches == bm) { //rp->score && 
					uint32_t rix = RefIxSrt[rp->refIx],
					numGap = rp->numGapR + rp->numGapQ,
					numMis = rp->mismatches - numGap,
					qlen = ShrBins[i].len,
					alLen = qlen + numGap,
					mOff = RefStart ? RefStart[rix] : 0, tmp, 
					stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
					edIxR = rp->finalPos + mOff;
					if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
					if (RefDedupIx) // handle dupes in refs
						for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) {
							rix = TmpRIX[k];
							mOff = RefStart ? RefStart[rix] : 0;
							stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
							edIxR = rp->finalPos + mOff, tmp;
							if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
							if (taxa_parsed) {
								char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
								for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
							}
							else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
						}
					else if (taxa_parsed) {
						char *FinalTaxon = taxa_lookup(RefHead[rix],taxa_parsed-1,Taxonomy);
						for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
					}
					else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
				}
				rp = rp->next;
			}
		}
	}
	else if (RUNMODE == FORAGE) { // all valid alignments
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i];
			while (rp) {
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
					if (RefDedupIx) // handle dupes in refs
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
	}
	else if (RUNMODE==CAPITALIST) {
		uint32_t numBins = RefDat.numRefHeads; // only unique headers get separate entries
		size_t *RefCounts = calloc(numBins,sizeof(*RefCounts)), tot = 0;
		for (uint32_t i = 0; i < numUniqQ; ++i) {
			ResultPod *rp = BasePod[i], *best = rp;
			if (!rp) continue;
			while (rp = rp->next) if (rp->mismatches < best->mismatches) best = rp; //rp->score && 
			// pass 2: actuarial research
			rp = best; 
			BasePod[i] = rp; // hedging derivative futures

			while (rp) {
				if (rp->mismatches == best->mismatches) {
					//if (RefDedupIx) 
					//	for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k)
					//		++RefCounts[RefMap[TmpRIX[k]]], ++tot;
					//else 
					++RefCounts[RefMap[RefIxSrt[rp->refIx]]], ++tot;
				}
				rp = rp->next;
			}
		}
		printf("CAPITALIST: Processed %llu investments\n", tot);
		//size_t LNL = 2000000; //, LNL_T = LNL/2; // arbitrary
		#pragma omp parallel
		{
			char **Taxa = 0, *Taxon = 0, *FinalTaxon = 0; uint32_t *Divergence = 0; //, *Dp = 0; 
			//char *LineDmp = malloc(LNL); 
			if (taxa_parsed) 
				Taxa = malloc(numBins*sizeof(*Taxa)), 
				Taxon = malloc(1000000),
				/*Dp = malloc(numBins*sizeof(*Dp)),*/
				Divergence = malloc(numBins*sizeof(*Divergence)), *Divergence = 0;
			
			// pass 3: deciding which investments to bank
			#pragma omp for
			for (uint32_t i = 0; i < numUniqQ; ++i) {
				if (!BasePod[i]) continue; 
				ResultPod *rp = BasePod[i];
				ResultPod *best = rp;
				uint32_t tix = 0;
				float best_score;

				if (taxa_parsed) {
					if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) 
						Taxa[tix++] = taxa_lookup(RefHead[TmpRIX[k]], taxa_parsed-1, Taxonomy);
					else Taxa[tix++] = taxa_lookup(RefHead[RefIxSrt[rp->refIx]], taxa_parsed-1, Taxonomy);
					best_score = rp->score;
				}
					
				while (rp = rp->next) if (rp->mismatches == best->mismatches) {
					if (taxa_parsed) {
						if (RefDedupIx) for (uint32_t k = RefDedupIx[rp->refIx]; k < RefDedupIx[rp->refIx+1]; ++k) 
							Taxa[tix++] = taxa_lookup(RefHead[TmpRIX[k]], taxa_parsed-1, Taxonomy);
						else Taxa[tix++] = taxa_lookup(RefHead[RefIxSrt[rp->refIx]], taxa_parsed-1, Taxonomy);
						if (rp->score > best_score) best_score = rp->score;
					}
					if ( (RefCounts[RefMap[RefIxSrt[rp->refIx]]] > RefCounts[RefMap[RefIxSrt[best->refIx]]]) || 
						(RefCounts[RefMap[RefIxSrt[rp->refIx]]] == RefCounts[RefMap[RefIxSrt[best->refIx]]] && 
							RefIxSrt[rp->refIx] < RefIxSrt[best->refIx]) )
								best = rp;
				} 
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
					uint32_t cutoff = tix - tix/TAXACUT; // need 90% support?
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
					--ed, --lv;
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
				uint32_t rix = RefIxSrt[rp->refIx], 
				numGap = rp->numGapR + rp->numGapQ,
				numMis = rp->mismatches - numGap,
				qlen = ShrBins[i].len,
				alLen = qlen + numGap,
				mOff = RefStart ? RefStart[rix] : 0, tmp, 
				stIxR = rp->finalPos - qlen + rp->numGapR + mOff,
				edIxR = rp->finalPos + mOff;
				if (rp->rc) tmp = stIxR, stIxR = edIxR, edIxR = tmp;
				
				#pragma omp critical
				if (taxa_parsed) for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH_TAX()
				else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) PRINT_MATCH()
				/*size_t sptr = 0;
				if (taxa_parsed) for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) {
					sptr+=SPRINT_MATCH_TAX()
					if (sptr >= LNL >> 1) {
						LineDmp = realloc(LineDmp,2*sptr*sizeof(*LineDmp));
						if (!LineDmp) {fputs("OOM:LineDmp\n",stderr); exit(3);}
					}
				}
				else for (uint32_t j = Offset[i]; j < Offset[i+1]; ++j) {
					sptr+=SPRINT_MATCH()
					if (sptr >= LNL >> 1) {
						LineDmp = realloc(LineDmp,2*sptr*sizeof(*LineDmp));
						if (!LineDmp) {fputs("OOM:LineDmp\n",stderr); exit(3);}
					}
				}
				#pragma omp critical
				fwrite(LineDmp,1,sptr,output);*/
			}
		}
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
	if (argc < 2) PRINT_USAGE()
	for (int i = 1; i < argc; ++i) {
		if (!strcmp(argv[i],"--references") || !strcmp(argv[i],"-r")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --references requires filename argument."); exit(1);}
			ref_FN = argv[i];
		}
		else if (!strcmp(argv[i],"--queries") || !strcmp(argv[i],"-q")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --queries requires filename argument."); exit(1);}
			query_FN = argv[i];
		}
		else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --output requires filename argument."); exit(1);}
			output_FN = argv[i];
		}
		//puts("--forwardreverse (-fr): also search the reverse complement of queries"); 
		else if (!strcmp(argv[i],"--forwardreverse") || !strcmp(argv[i],"-fr")) {
			QDat.rc = 1; 
			printf(" --> Also considering the reverse complement of reads.\n");
		}
		else if (!strcmp(argv[i],"--whitespace") || !strcmp(argv[i],"-w")) {
			QDat.incl_whitespace = 1; 
			printf(" --> Allowing whitespace in query name output.\n");
		}
		else if (!strcmp(argv[i],"--npenalize") || !strcmp(argv[i],"-n")) {
			Z = 1; // global
			printf(" --> Setting N penalty (ref N vs query A/C/G/T).\n");
		}
		else if (!strcmp(argv[i],"--xalphabet") || !strcmp(argv[i],"-x")) {
			Xalpha = 1; // global
			printf(" --> Allowing any alphabet (unambiguous ID matching).\n");
		}
		else if (!strcmp(argv[i],"--taxonomy") || !strcmp(argv[i],"-b")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --taxonomy requires filename argument."); exit(1);}
			tax_FN = argv[i];
			printf(" --> Assigning taxonomy based on input file %s.\n",tax_FN);
		}
		else if (!strcmp(argv[i],"--mode") || !strcmp(argv[i],"-m")) {
			if (++i == argc || argv[i][0] == '-')  
				{ puts("ERROR: --mode requires an argument (see -h)"); exit(1); }
			if (!strcmp(argv[i],"BEST")) RUNMODE = BEST;
			else if (!strcmp(argv[i],"ALLPATHS")) RUNMODE = ALLPATHS;
			else if (!strcmp(argv[i],"CAPITALIST")) RUNMODE = CAPITALIST;
			else if (!strcmp(argv[i],"MATRIX")) 
				{fputs("ERROR: Matrix mode is no longer supported\n",stderr); exit(1);}
			else if (!strcmp(argv[i],"FORAGE")) RUNMODE = FORAGE;
			else {printf("Unsupported run mode '%s'\n",argv[i]); exit(1);}
			printf(" --> Setting run mode to %s.\n",argv[i]);
		}
		else if (!strcmp(argv[i],"--makedb") || !strcmp(argv[i],"-d")) {
			makedb = 1; char *dbsel = "QUICK"; // make this read the default
			if (i + 1 != argc && argv[i+1][0] != '-') { // arg provided
				if (!strcmp(argv[++i],"DNA")) dbType = DNA_16;
				else if (!strcmp(argv[i],"RNA")) dbType = DNA_16;
				else if (!strcmp(argv[i],"PROTEIN")) dbType = PROT;
				else if (!strcmp(argv[i],"QUICK")) dbType = QUICK;
				else {printf("Unsupported makedb mode '%s'\n",argv[i]); exit(1);};
				dbsel = argv[i];
			}
			if (i + 1 != argc && argv[i+1][0] != '-') { // another arg provided
				DB_QLEN = atol(argv[++i]);
				if (DB_QLEN < 0) {fprintf(stderr,"ERROR: bad max query length '%s'\n",argv[i]); exit(1);}
				if (DB_QLEN < 5) printf("WARNING: query length very short (%ld)\n",DB_QLEN);
			}
			printf(" --> Creating %s database (assuming max query length %ld)\n",dbsel, DB_QLEN);
		}
		else if (!strcmp(argv[i],"--accelerator") || !strcmp(argv[i],"-a")) {
			//puts("ERROR: Accelerator not currently implemented [target 0.99.2]"); exit(1);
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --accelerator requires filename argument."); exit(1);}
			xcel_FN = argv[i];
			DO_ACCEL = 1;
			printf(" --> Using accelerator file %s\n",xcel_FN);
		}
		else if (!strcmp(argv[i],"--taxacut") || !strcmp(argv[i],"-bc")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --taxacut requires integer argument."); exit(1);}
			int temp = atoi(argv[i]);
			if (temp < 2) {fputs("ERROR: taxacut must be >= 2\n",stderr); exit(1);}
			TAXACUT = temp;
			printf(" --> Ignoring 1/%d disagreeing taxonomy calls\n",TAXACUT);
		}
		else if (!strcmp(argv[i],"--taxa_ncbi") || !strcmp(argv[i],"-bn")) {
			taxa_lookup = taxa_lookup_ncbi;
			printf(" --> Using NCBI header formatting for taxonomy lookups.\n");
		}
		else if (!strcmp(argv[i],"--skipambig") || !strcmp(argv[i],"-sa")) {
			QDat.skipAmbig = 1;
			printf(" --> Skipping highly ambiguous queries (4+ consecutive ambigs).\n");
		}
		else if (!strcmp(argv[i],"--taxasuppress") || !strcmp(argv[i],"-bs")) {
			QDat.taxasuppress = 1;
			if (i + 1 != argc && argv[i+1][0] != '-') { // arg provided
				if (!strcmp(argv[++i],"STRICT")) TAXLEVELS = TAXLEVELS_STRICT;
				else {fprintf(stderr,"ERROR: Unrecognized taxasuppress '%s'\n",argv[i]); exit(1);}
			}
			printf(" --> Surpressing taxonomic specificity by alignment identity%s.\n",
				TAXLEVELS == TAXLEVELS_STRICT ? " [STRICT]" : "");
		}
		else if (!strcmp(argv[i],"--id") || !strcmp(argv[i],"-i")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --id requires decimal argument."); exit(1);}
			THRES = atof(argv[i]);
			if (THRES > 1.f || THRES < 0.f) {puts("Invalid id range [0-1]"); exit(1);}
			if (THRES < 0.01f) THRES = 0.01f;
			printf(" --> Setting identity threshold to %f.\n",THRES);
		}
		else if (!strcmp(argv[i],"--threads") || !strcmp(argv[i],"-t")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --threads requires integer argument."); exit(1);}
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
			printf(" --> Preprocessing references (dereplicating).\n");
		}
		else if (!strcmp(argv[i],"--fingerprint") || !strcmp(argv[i],"-f")) {
			DO_FP = 1; // global
			printf(" --> Using fingerprint profiling.\n");
		}
		else if (!strcmp(argv[i],"--prepass") || !strcmp(argv[i],"-p")) {
			printf(" --> Using fingerprint pre-checks ");
			if (i + 1 != argc && argv[i+1][0] != '-') { // arg provided
				if (!strcmp(argv[++i],"reorder")) DO_PREPASS = REORDER, printf("[reorder only]\n");
				else if (!strcmp(argv[i],"fast")) DO_PREPASS = FAST, printf("[fast]\n");
				else if (!strcmp(argv[i],"full")) DO_PREPASS = FULL, printf("[full]\n");
				else if (!strcmp(argv[i],"auto")) DO_PREPASS = AUTO, printf("[auto]\n");
				else if (!strcmp(argv[i],"off")) DO_PREPASS = NONE, printf("is DISABLED\n");
				else {printf("\nERROR: unsupported prepass mode '%s'\n",argv[i]); exit(1);}
			}
			else DO_PREPASS = PP_DEF_ON, printf("["PP_DEF_ST"]\n");
		}
		else if (!strcmp(argv[i],"--cache") || !strcmp(argv[i],"-c")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --cache requires integer argument."); exit(1);}
			cacheSz = atoi(argv[i]); // global
			printf(" --> Setting number of cached lines in matrix to %d.\n",cacheSz);
		}
		else if (!strcmp(argv[i],"--latency") || !strcmp(argv[i],"-l")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --latency requires integer argument."); exit(1);}
			LATENCY = atoi(argv[i]); // global
			printf(" --> Setting clump formation latency to %d bases.\n",LATENCY);
		}
		else if (!strcmp(argv[i],"--clustradius") || !strcmp(argv[i],"-cr")) {
			if (++i == argc || argv[i][0] == '-') 
				{ puts("ERROR: --clustradius requires integer argument."); exit(1);}
			RefDat.clustradius = atoi(argv[i]); // Reference_data member
			printf(" --> Setting FP cluster search radius to %d members.\n",RefDat.clustradius);
		}
		else if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) PRINT_USAGE()
		else {
			printf("ERROR: Unrecognized command-line option: %s.\n",argv[i]);
			puts("See help by running with just '-h'");
			exit(1);
		}
	}
	FILE *output = fopen(output_FN,"wb");
	if (!output) {fprintf(stderr,"ERROR: Cannot open output: %s",output_FN); exit(2);}

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
	if (makedb && dbType != QUICK) { 
		puts("DB mode not implemented in this release. Use QUICK instead.");
		exit(0);
	} 
	else if (makedb) {
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
			exit(1);
		}
		REBASE_AMT = dShear;
		//alignNVU(RefDat.RefSeq[0], QDat.QSeq[0], RefDat.RefLen[0], QDat.QLen[0]);
		do_alignments(output, RefDat, QDat);
	}
	
	printf("\nAlignment time: %f seconds\n", omp_get_wtime() - start);
	return 0;
}
