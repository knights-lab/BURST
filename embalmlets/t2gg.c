// assumptions: <256 levels, >25M tids, <64k linlen, rs82 dmp format, virus tid = 10239
#define HOWTO "Usage: t2gg nodes.dmp names.dmp tid2gg.txt"
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
char *D[] = {
"species\0\0\0\0\0\0\0\0\0\0",
"genus\0\0\0\0\0\0\0\0\0\0\0",
"family\0\0\0\0\0\0\0\0\0\0",
"order\0\0\0\0\0\0\0\0\0\0\0",
"class\0\0\0\0\0\0\0\0\0\0\0",
"phylum\0\0\0\0\0\0\0\0\0\0",
"kingdom\0\0\0\0\0\0\0\0\0\0",
"superkingdom\0\0\0\0"}; 
char L[] = "sgfocpkk"; //"kkpcofgs";
char O[] = "kpcofgst";
typedef struct { uint32_t p; char *n; char r; } Branch;

void main(int argc, char *argv[]) {
	if (argc != 4) {puts(HOWTO); exit(1);}
	FILE *nd = fopen(argv[1],"rb"), *nm = fopen(argv[2],"rb"), 
		*out = fopen(argv[3],"wb");
	if (!nd || !nm || !out) {fputs("I/O error\n",stderr); exit(2);}
	char *line = malloc(UINT16_MAX), *lineO = line;
	Branch *Tree = calloc(25000000,sizeof(*Tree));
	if (!line || !Tree) {fputs("ERROR:OOM:TREE\n",stderr); exit(3); }
	
	// make the ncbi taxonomy tree
	uint32_t nl = 0, lastIx = 0;
	while (line = fgets(line,UINT16_MAX,nd)) {
		uint32_t ix = atol(line);
		char *atParent = strchr(line,'|');
		if (!atParent) {fputs("\nBad tree\n",stderr); exit(2);}
		Tree[ix].p = atol(atParent+1);
		uint32_t len = strlen(line);
		char *atRank = strchr(atParent + 1, '|');
		if (!atRank) {fputs("\nBad tree\n",stderr); exit(2);}
		atRank += 2; 
		char *end = strchr(atRank,'\t');
		if (!end) {fputs("\nBad tree\n",stderr); exit(2);}
		*end = 0;
		Tree[ix].r = 'x';
		for (int i = 0; i < 8; ++i) if (!strcmp(atRank,D[i])) {
			Tree[ix].r = L[i]; 
			break;
		}
		if (ix > lastIx) lastIx = ix;
		++nl;
	}
	line = lineO;
	Tree = realloc(Tree, ++lastIx*sizeof(*Tree));
	fputs("Done with node parse\n",stderr);
	
	// massage the taxonomy labels
	for (uint32_t i = 0; i < lastIx; ++i) {
		if (Tree[i].p == 10239) Tree[i].r = 'p';
		else if (Tree[Tree[i].p].p == 10239) Tree[i].r = 'c';
		else if (Tree[Tree[i].p].r == 's') Tree[i].r = 't';
		else if (Tree[i].r == 'k') Tree[i].p = 1; 
	}
	fputs("Done with taxonomizing\n",stderr);
	
	// grab the names
	uint32_t lix = 0;
	while (line = fgets(line,UINT16_MAX,nm)) {
		uint32_t ix = atol(line);
		if (ix == lix || Tree[ix].r == 'x') {lix = ix; continue;}
		char *begin = strchr(line,'|');
		if (!begin) {fputs("\nBad names\n",stderr); exit(2);}
		begin += 2;
		char *end = strchr(begin,'\t');
		if (!end) {fputs("\nBad names\n",stderr); exit(2);}
		*end = 0;
		char *nameB = strchr(end+2,'|');
		if (!nameB) {fputs("\nBad names\n",stderr); exit(2);}
		nameB += 2;
		if (*nameB == 's' && *(nameB+1) == 'c') {
			lix = ix; 
			char *name = strcpy(malloc(end - begin + 1),begin);
			name[end-begin] = 0;
			Tree[ix].n = name;
		}
	}
	free(lineO);
	fputs("Done with name assignment\n",stderr);
	
	// assign each tid a GG taxonomy string
	uint32_t Composer[256] = {0}; 
	for (uint32_t i = 2; i < lastIx; ++i) {
		if (!Tree[i].p) continue;
		uint32_t node = i, lv = 0;
		while (node > 1) 
			Composer[++lv] = node,
			node = Tree[node].p;
		char cur = 0;
		fprintf(out,"%u\t",i);
		for (uint32_t j = lv; j; --j) {
			Branch T = Tree[Composer[j]];
			if (T.r == 'x') continue;
			for (; O[cur] != T.r; ++cur) // fill gap
				fprintf(out,"%c__;",O[cur]);
			fprintf(out,"%c__%s%s",O[cur],T.n,cur < 7 ? ";" : "");
			++cur;
		}
		for (; cur < 8; ++cur) fprintf(out,"%c__%s", O[cur], cur < 7 ? ";" : "");
		fprintf(out,"\n");
	}
	fputs("Done with file writing!\n",stderr);
}