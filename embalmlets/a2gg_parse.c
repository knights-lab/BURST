#define HOWTO "Usage: a2gg_parse in.fasta in.mapDB outPrefix [d] [e] [FULL] [GUESS]"
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#define DBVER -1

char BAY[1<<16] = {0};
typedef struct {char *key, *val;} Tuple;
static inline uint64_t crBST(char *key, uint64_t sz, char **Dict) {
	char **p = Dict;
	while (sz) {
		uint64_t w = sz >> 1; 
		char *ref_s = p[w+1], *key_s = key;
		while (*ref_s == *key_s++) if (!*ref_s++) return p+(w+1)-Dict; 
		if (*ref_s < *(key_s-1)) p+=w+1, sz-=w+1; 
		else sz = w;
	}
	char *ref_s = *p, *key_s = key;
	while (*ref_s == *key_s++) if (!*ref_s++) return p-Dict;
	return 0;
}

void main(int argc, char *argv[]) {
	if (argc < 4) {puts(HOWTO); exit(1);}
	FILE *in = fopen(argv[1],"rb"), *map = fopen(argv[2],"rb");
	sprintf(BAY,"%s.fna",argv[3]); FILE *outF = fopen(BAY,"wb");
	sprintf(BAY,"%s.tax",argv[3]); FILE *outT = fopen(BAY,"wb");
	if (!in || !map || !outT || !outF) {puts("Can't open file(s)"); exit(1);}
	int guess = !strcmp(argv[argc-1],"GUESS"); argc -= guess != 0;
	int full = !strcmp(argv[argc-1],"FULL"); argc -= full != 0;
	char begD = argc >=5 ? *argv[4] : '>', endD = argc >= 6 ? *argv[5] : '\n';
	printf("Starting after char: '%c', ending at: '%c'\n",begD,endD);
	
	// Database import
	char control = fgetc(map);
	if (control != DBVER) {fputs("ERROR: Incompatible DB. Use a2gg_make\n",stderr); exit(2);}
	uint64_t charsInAcc, charsInGG, nm, nl, token;
	fread(&charsInAcc, 8, 1, map);
	fread(&charsInGG, 8, 1, map);
	fread(&nm, 8, 1, map);
	fread(&nl, 8, 1, map);
	printf("--> MapDB: nm %llu, nl %llu, ACC %llu, GG %llu\n",
		nm, nl, charsInAcc, charsInGG);
	char *AccDmp = malloc(charsInAcc), *GGDmp = malloc(charsInGG),
		**AccDict = malloc(nm*sizeof(*AccDict)), **GGDict = malloc(nl*sizeof(*GGDict));
	uint32_t *AccLn = malloc(nm*sizeof(*AccLn)), *GGLn = malloc(nl*sizeof(*GGLn));
	if (!(AccDmp && GGDmp && AccDict && GGDict && AccLn && GGLn)) {
		fputs("ERROR: Out of memory\n",stderr); exit(3); }
	token = fread(AccDmp,1,charsInAcc,map);
	if (token != charsInAcc) {fputs("ERR:AccDmp\n",stderr); exit(3);}
	token = fread(AccLn,4,nm,map);
	if (token != nm) {fputs("ERR:AccLn\n",stderr); exit(3);}
	token = fread(GGDmp,1,charsInGG,map);
	if (token != charsInGG) {fputs("ERR:GGDmp\n",stderr); exit(3);}
	//token = fread(GGLn,4,nl,map);
	//if (token != nl) {fputs("ERR:GGLn\n",stderr); exit(3);}
	
	// Now enumerate the pointers
	char *DmpP = AccDmp;
	*AccDict = AccDmp;
	for (uint64_t i = 1; i < nm; ++i) {
		while (*DmpP++);
		AccDict[i] = DmpP;
	}
	DmpP = GGDmp;
	*GGDict = GGDmp;
	for (uint64_t i = 1; i < nl; ++i) {
		while (*DmpP++);
		GGDict[i] = DmpP;
	}
	puts("--> MapDB: All data read successfully.");
	
	// Do the sequence mapping
	printf("\nCreating output files...\n");
	uint32_t lnSz = INT32_MAX;
	char *line = malloc(lnSz), *lineO = line; line[lnSz] = 0;
	if (!line) {fputs("OOM:fasta\n",stderr); exit(3);}
	--nm;  // for mapper
	char *taxon = "UNKNOWN", *begin, *end;
	while (++nl, line = fgets(line, lnSz, in)) {
		begin = strchr(line, begD);
		if (!begin) {
			if (*line == '>') {
				printf("Ln %llu: '%s' *DELIM* not found\n",nl,begin); 
				line = fgets(line, lnSz, in);
				continue;
			}
			else {printf("End of fasta reached (ln %llu)\n",nl); break;}
		}
		end = strchr(++begin, endD);
		if (!end) {
			printf("Ln %llu: '%s' *END* not found\n",nl,begin); 
			line = fgets(line, lnSz, in);
			continue;
		}
		*end = 0;
		uint64_t ix = crBST(begin,nm,AccDict); 
		//printf("Found %s at %s [ix %llu, maps to line %llu in tax]\n",begin,AccDict[ix],ix,AccLn[ix]);
		char *taxonTest = ix ? GGDict[AccLn[ix]] : 0;
		//char *taxonTest = crBST(begin,nm,Entries);
		if (!taxonTest && !guess) {
			printf("Ln %llu: '%s' *TAXON* not found [%llu]\n",nl,begin,ix); 
			line=fgets(line, lnSz, in);
		} else if (full) {
			taxon = taxonTest ?: taxon;
			if (!taxonTest) printf("Ln %llu: '%s' *TAXON* interpolation: %s\n",nl,begin,taxon); 
			if (endD != '\n') {
				char *p = strchr(end+1,'\n');
				if (p) *p = 0;
				fprintf(outT,"%s%c%s\t%s\n",line+1,endD,end+1,taxon);
				fprintf(outF,"%s%c%s\n",line,endD,end+1);
			} else fprintf(outT,"%s\t%s\n",line,taxon),
				fprintf(outF,"%s\n",line);
			fputs(line=fgets(line, lnSz, in),outF);
		} else {
			taxon = taxonTest ?: taxon;
			if (!taxonTest) printf("Ln %llu: '%s' *TAXON* interpolation: %s\n",nl,begin,taxon); 
			if (endD != '\n') 
				fprintf(outT,"%s\t%s\n",begin,taxon),
				fprintf(outF,"%s%c%s",line,endD,end+1);
			else fprintf(outT,"%s\t%s\n",begin,taxon),
				fprintf(outF,"%s\n",line);
			fputs(line=fgets(line, lnSz, in),outF);
		}
	}
}