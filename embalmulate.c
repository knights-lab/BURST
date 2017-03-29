#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
// Version 0.2
typedef struct T T;
struct T { char *s; uint64_t i; T *left, *right; };
static inline T* T_add(T *t, char *s, uint64_t l, uint64_t *i) {
	if (!t->s) {
		*t = (T){strcpy(malloc(l+1),s),(*i)++,0,0};
		return t;
	}
	int cmp = strcmp(s,t->s);
	do {
		if (cmp > 0) {
			if (!t->right) {
				t->right = malloc(sizeof(T));
				*t->right = (T){strcpy(malloc(l+1),s),(*i)++,0,0};
				return t->right;
			} t = t->right;
		} else if (cmp < 0) {
			if (!t->left) {
				t->left = malloc(sizeof(T)); 
				*t->left = (T){strcpy(malloc(l+1),s),(*i)++,0,0};
				return t->left;
			} t = t->left;
		}
	} while (cmp = strcmp(s,t->s));
	return t;
}
typedef struct L L;
struct L { uint64_t ix, cnt; L *left, *right; };
static inline void L_add(L *t, uint64_t ix) {
	do {
		if (ix > t->ix) {
			if (!t->right) {
				t->right = malloc(sizeof(L));
				*t->right = (L){ix,1,0,0};
				return;
			} t = t->right;
		} else if (ix < t->ix) {
			if (!t->left) {
				t->left = malloc(sizeof(L));
				*t->left = (L){ix,1,0,0};
				return;
			} t = t->left;
		}
	} while (t->ix != ix);
	++t->cnt;
}
typedef struct XT XT;
struct XT { char *s; L* ls; XT *left, *right; };
static inline uint64_t XT_add(XT *t, char *s, uint64_t l, uint64_t ix) {
	if (!t->s) {
		L *nl = malloc(sizeof(L));
		*nl = (L){ix, 1, 0, 0};
		*t = (XT){strcpy(malloc(l+1),s),nl,0,0};
		return 1;
	}
	int cmp = strcmp(s,t->s);
	do {
		if (cmp > 0) {
			if (!t->right) {
				t->right = malloc(sizeof(XT));
				L *nl = malloc(sizeof(L));
				*nl = (L){ix, 1, 0, 0};
				*t->right = (XT){strcpy(malloc(l+1),s),nl,0,0};
				return 1;
			} t = t->right;
		} else if (cmp < 0) {
			if (!t->left) {
				t->left = malloc(sizeof(XT));
				L *nl = malloc(sizeof(L));
				*nl = (L){ix, 1, 0, 0};
				*t->left = (XT){strcpy(malloc(l+1),s),nl,0,0};
				return 1;
			} t = t->left;
		}
	} while (cmp = strcmp(s,t->s));
	L_add(t->ls,ix); 
	return 0;
}
static inline void L_write(L *t, uint32_t *row) {
	row[t->ix] = t->cnt;
	if (t->left) L_write(t->left, row);
	if (t->right) L_write(t->right, row);
}
static inline void XT_dump(XT *t, uint32_t *row, uint32_t l, FILE *out) {
	memset(row,'\0',l*sizeof(*row));
	L_write(t->ls, row);
	fprintf(out,"\n%s",t->s);
	for (uint32_t i = 0; i < l; ++i) fprintf(out, "\t%u", row[i]);
	if (t->left) XT_dump(t->left, row, l, out);
	if (t->right) XT_dump(t->right, row, l, out);
}

static inline void writeNodes(T *t, char **f) {
	f[t->i] = t->s;
	if (t->left) writeNodes(t->left, f);
	if (t->right) writeNodes(t->right, f);
}

void main(int argc, char *argv[]) {
	if (argc < 3) {puts("Usage: embalmulate in.b6 out.tsv [outTax.tsv] ['GGtrim']"); exit(1);}
	FILE *in = fopen(argv[1],"rb"), *out = fopen(argv[2],"wb");
	FILE *tax = argc > 3 ? fopen(argv[3],"wb") : 0;
	if (!in || !out || (argc > 3 && !tax)) {puts("Can't open file(s)"); exit(1);}
	int ggtrim = 0; if (argc >= 4 && !strcmp(argv[argc-1],"GGtrim")) --argc, ggtrim = 1;
	uint64_t i = 0, ns = 0, nt = 0, nr = 0, ix;
	T *SampT = calloc(1,sizeof(*SampT)), *node = 0;
	XT *RefT = calloc(1,sizeof(*RefT)), *TaxT = calloc(1,sizeof(*TaxT));
	char *line = malloc((2<<16)+1), *lineO = line, *samp;
	while (line = fgets(line,2<<16,in)) {
		char *start = line, *end = start;
		for (; *end && *end != '_' && *end != '\t'; ++end);
		if (!*end) {printf("End of file reached [1]: %llu\n",i); break;}
		if (*end == '_') {
			*end = 0,
			node = T_add(SampT,start,end-start,&ns),
			samp = node->s, ix = node->i;
			for (++end; *end && *end != '\t'; ++end);
			if (!*end) {printf("End of file reached [1b]: %llu\n",i); break;}
		}
		else samp = argv[1], ix = 0;
		start = end+1;
		if (!*start) {printf("End of file reached [2]: %llu\n",i); break;}
		for (end = start; *end && *end != '\t'; ++end);
		if (!*end) {printf("End of file reached [2b]: %llu\n",i); break;}
		*end = 0;
		nr += XT_add(RefT,start,end-start,ix);
		if (tax) {
			end = strchr(end+1,'\0'), *--end = 0;
			for (start = end-1; *start && *start != '\t'; --start);
			char *taxon = ++start; // search taxon for danglies; either expand or contract them.
			if (GGtrim && end > start) {
				while (*(end-1) == '_') {
					do --end; while (end > start && *end != ';');
					*end = 0;
				}
			}
			nt += XT_add(TaxT,taxon,end-taxon,ix);
		}
		++i;
	}
	free(lineO);
	printf("Parsed %llu reads [%llu samples, %llu taxa, %llu refs]. Collating...\n",i,ns,nt,nr);
	char **Samps = malloc(ns*sizeof(*Samps));
	writeNodes(SampT,Samps);
	fputs("#OTU ID",out); if (tax) fputs("#OTU ID",tax);
	for (uint32_t i = 0; i < ns; ++i) fprintf(out, "\t%s", Samps[i]);
	if (tax) for (uint32_t i = 0; i < ns; ++i) fprintf(tax, "\t%s", Samps[i]);
	free(Samps);
	uint32_t *row = malloc(ns*sizeof(*row));
	XT_dump(RefT, row, ns, out);
	if (tax) XT_dump(TaxT, row, ns, tax);
}