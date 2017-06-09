// assumptions: <256 levels, >25M tids, <64k linlen, rs82 dmp format, virus tid = 10239
#define HOWTO "Usage: a2gg_make in.tid2gg in.acc2tid out.acc2gg [threads]"
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define DBVER -1

char NULLTAX[1] = {0}, BAY[1<<16] = {0};
typedef struct {char *key; uint32_t val;} Tuple;

static int tupSrt(const void *a, const void *b) {
	return strcmp(((Tuple*)a)->key,((Tuple*)b)->key); }
#define SYM_BYTES 7
// To use, pass in Pack, len, and define structure type, num letters per block, and nib function
#define PARALLEL_SORT_PROTOTYPE(STRUCT_TYPE, STR_MEMBER, NUMLET, NIBFUNC) { \
	static const uint32_t NLB = 1 << (NUMLET*SYM_BYTES); \
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
#define LET4 s[0] << 21 | s[1] << 14 | s[2] << 7 | s[3]
static inline void parallel_sort_tuple(Tuple *Pack, uint32_t len) 
	PARALLEL_SORT_PROTOTYPE(Tuple, key, 4, LET4)

static uint32_t crBST(char *key, uint32_t sz, Tuple *Dict) {
	Tuple *p = Dict;
	while (sz) {
		uint32_t w = sz >> 1; 
		char *ref_s = p[w+1].key, *key_s = key;
		while (*ref_s == *key_s++) if (!*ref_s++) return p[w+1].val; 
		if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	char *ref_s = p->key, *key_s = key;
	while (*ref_s == *key_s++) if (!*ref_s++) return p->val;
	return 0;
}

void main(int argc, char *argv[]) {
	if (argc < 4) {puts(HOWTO); exit(1);}
	FILE *gg = fopen(argv[1],"rb"), *map = fopen(argv[2],"rb");
	FILE *out = fopen(argv[3],"wb");
	if (!gg || !map || !out) {puts("Can't open file(s)"); exit(1);}
	int threads = omp_get_max_threads();
	if (argc > 4) threads = atoi(argv[4]);
	omp_set_num_threads(threads);
	
	uint64_t nm = 0, nl = 0, sz = 1000, lnSz = ((uint64_t)1 << 31)-1; // 30: giglines
	char *line = malloc(lnSz+1), *lineO = line; line[lnSz] = 0;
	Tuple *Entries = malloc(sz * sizeof(*Entries));
	char *begin, *end;
	line = fgets(line, lnSz, map); // discard the header
	uint64_t biggest = 0, whence = 0, highest = 0;
	while (line = fgets(line, lnSz, map)) {
		end = strchr(begin = line,'\t'); 
		if (!begin || !end) {printf("Error on map line %llu\n",nm+1); break;}
		*end = 0; 
		Entries[nm].key = strcpy(malloc(end-begin+1),begin);
		if (end-begin > biggest) biggest = end-begin, whence = nm;
		Entries[nm].val = atol(strchr(end + 1,'\t')); 
		//printf("[%llu] %s -> %u\n", nm, Entries[nm].key, Entries[nm].val);
		if (++nm == sz) Entries = realloc(Entries, (sz*=2)*sizeof(*Entries));
	}
	free(lineO);
	printf("Parsed %llu accessions. Largest: %llu at ln %llu [%s]\n",nm,biggest,whence,Entries[whence].key);
	Entries = realloc(Entries, nm*sizeof(*Entries));
	parallel_sort_tuple(Entries,nm);
	
	// read the gg index and make an index of tid -> line number
	sz = 1000, biggest = 0; Tuple *GG = malloc(sz * sizeof(*GG));
	line = malloc(lnSz+1), lineO = line; line[lnSz] = 0;
	while (line = fgets(line, lnSz, gg)) {
		GG[nl].val = atoi(line);
		if (GG[nl].val > highest) highest = GG[nl].val;
		begin = strchr(line,'\t') + 1;
		end = strchr(begin,'\n'); 
		if (begin==(void*)1 || !end) {printf("Error on gg line %llu\n",nl+1); break;}
		*end = 0; 
		GG[nl].key = strcpy(malloc(end-begin+1),begin);
		if (end-begin > biggest) biggest = end-begin, whence = nl;
		//printf("[%llu] %s -> %u\n", nl, GG[nl].key, GG[nl].val);
		if (++nl == sz) GG = realloc(GG, (sz*=2)*sizeof(*GG));
	}
	free(lineO);
	printf("Parsed %llu tids. Largest: %llu at ln %llu [%s] / %llu\n",nl,biggest,whence,GG[whence].key, highest);
	GG = realloc(GG, nl*sizeof(*GG));
	uint32_t *RevMap = calloc(highest,sizeof(*RevMap));
	for (uint32_t i = 0; i < nl; ++i) RevMap[GG[i].val] = i+1;
	++nl;
	
	// Write database
	printf("\nCreating output files...\n");
	uint64_t charsInAcc = 0, charsInGG = 0;
	char control = DBVER; // also acts as version
	fputc(control, out);
	fwrite(&charsInAcc, 8, 1, out);
	fwrite(&charsInGG, 8, 1, out);
	fwrite(&nm, 8, 1, out);
	fwrite(&nl, 8, 1, out);
	for (uint64_t i = 0; i < nm; ++i) 
		charsInAcc += fprintf(out, "%s", Entries[i].key) + 1, 
		fputc(0,out);
	for (uint64_t i = 0; i < nm; ++i) {
		uint32_t lnIx = RevMap[Entries[i].val];
		fwrite(&lnIx, 4, 1, out);
	}
	charsInGG += fprintf(out, "UNKNOWN") + 1, fputc(0,out);
	for (uint64_t i = 0; i < nl; ++i)
		charsInGG += fprintf(out, "%s", GG[i].key) + 1,
		fputc(0,out);
	//for (uint64_t i = 0; i < nl; ++i) fwrite(&GG[i].val, 4, 1, out);
	rewind(out);
	fputc(control, out);
	fwrite(&charsInAcc, 8, 1, out);
	fwrite(&charsInGG, 8, 1, out);
	printf("Database written! [nm %llu, nl %llu, ACC %llu, GG %llu]\n",
		nm, nl, charsInAcc, charsInGG);
}