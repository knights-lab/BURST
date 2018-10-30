#define VER "v0.92"
#define USAGE "bcov in.alignments.b6 in.table.txt OUT_PREFIX [<VAR>] [PAD <X>] [SPLIT]"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <math.h>
#include <sys/resource.h>

typedef struct {char *S; uint32_t L, *U,*C;} Record;
typedef struct T T;
struct T { char *s; uint64_t i; T *left, *right; uint32_t **Uc, **Cc; };
static inline T* T_add(T *t, char *s, uint64_t l, uint64_t *i) {
	if (!t->s) {
		*t = (T){memcpy(calloc(l+1,1),s,l),(*i)++,0,0,0,0};
		return t;
	}
	int cmp = 0, j = 0;
	for (char *ts = t->s; j < l; ++j) if (s[j]!=ts[j]) {
		cmp = ts[j]-s[j]; break;}
	do {
		if (cmp > 0) {
			if (!t->right) {
				t->right = malloc(sizeof(T));
				*t->right = (T){memcpy(calloc(l+1,1),s,l),(*i)++,0,0,0,0};
				return t->right;
			} t = t->right;
		} else if (cmp < 0) {
			if (!t->left) {
				t->left = malloc(sizeof(T)); 
				*t->left = (T){memcpy(calloc(l+1,1),s,l),(*i)++,0,0,0,0};
				return t->left;
			} t = t->left;
		}
		cmp = 0, j = 0;
		for (char *ts = t->s; j < l; ++j) if (s[j]!=ts[j]) {
			cmp = ts[j]-s[j]; break;}
	} while (cmp); 
	return t;
}
static inline void writeNodes(T *t, char **f, uint32_t ***GU, uint32_t ***GC) {
	f[t->i] = t->s, GU[t->i] = t->Uc, GC[t->i] = t->Cc;
	if (t->left) writeNodes(t->left, f, GU, GC);
	if (t->right) writeNodes(t->right, f, GU, GC);
}
static inline int cmp(const void *a, const void *b) {
	Record *A = (Record *)a, *B = (Record *)b;
	return strcmp(A->S, B->S);
}
static Record * crBST_XAT(char *key, uint32_t li, Record *D) {
	Record *p = D;
	while (li) {
		uint32_t w = li >> 1; 
		char *ref_s = p[w+1].S, *key_s = key;
		while (*key_s == *ref_s++)
			if (*ref_s == '\t' || !*key_s++) return p+w+1;
		if (*(ref_s-1) < *key_s) p += w+1, li -= w+1; 
		else li = w;
	}
	char *ref_s = p->S, *key_s = key;
	while (*key_s == *ref_s++) 
		if (*ref_s == '\t' || !*key_s++) return p;
	return 0;
}
char FILEN[4096] = {0};
int main(int argc, char *argv[]) {
	puts("This is BURSTcoverage (bcov) "VER);
	if (argc < 4) {puts("Usage: " USAGE); exit(1);}
	int inf, inm = -1, DO_SAMP = 0; 
	long PAD = 0;
	FILE *outs, *outu, *outbs, *outbu; 
	char *mm_inf, *mm_inm; 
	double vf = 1.0;
	struct stat sb; uint64_t sz_inf, sz_inm = 0; 

	if (argc > 4 && !strcmp(argv[argc-1],"SPLIT")) --argc, DO_SAMP = 1;
	if (argc > 5 && !strcmp(argv[argc-2],"PAD")) 
		PAD = atoi(argv[argc-1]), argc-=2, printf("Padding %ld bp\n",PAD);
	if (argc > 4) vf = atof(argv[argc-1]);
	inf = open(argv[1],O_RDONLY);
	inm = open(argv[2],O_RDONLY);
	sprintf(FILEN,"%sshared.txt",argv[3]);
	outs = fopen(FILEN,"wb");
	sprintf(FILEN,"%sunique.txt",argv[3]);
	outu = fopen(FILEN,"wb");
	sprintf(FILEN,"%sshared_binary.txt",argv[3]);
	outbs = fopen(FILEN,"wb");
	sprintf(FILEN,"%sunique_binary.txt",argv[3]);
	outbu = fopen(FILEN,"wb");
	
	if (inf==-1 || inm==-1 || !(outs && outu && outbs && outbu)) {
		puts("I/O ERROR. Check filenames."); exit(2);}
	fstat(inm,&sb); sz_inm = sb.st_size; 
	mm_inm = mmap(0,sz_inm+2,PROT_READ|PROT_WRITE,MAP_PRIVATE,inm,0);
	if (mm_inm[sz_inm-1] == '\n') mm_inm[sz_inm-1] = 0;
	madvise(mm_inm, sz_inm+2, MADV_WILLNEED); 
	
	
	fstat(inf,&sb); sz_inf = sb.st_size; 
	mm_inf = mmap(0,sz_inf+4096,PROT_READ|PROT_WRITE,MAP_PRIVATE,inf,0);
	if (mm_inf[sz_inf-1] == '\n') mm_inf[sz_inf-1] = 0;
	madvise(mm_inf,sz_inf+4096,MADV_SEQUENTIAL);
	if (!sz_inf) {puts("B6 file is empty. Exiting."); exit(1);}
	
	setvbuf(outs,0,_IOFBF,1<<20);
	setvbuf(outu,0,_IOFBF,1<<20);
	setvbuf(outbs,0,_IOFBF,1<<20);
	setvbuf(outbu,0,_IOFBF,1<<20);

	// Begin parsing block
	uint32_t rLines = 0, rSz = 1000, nRecs = 0;
	Record *Recs = malloc(sizeof(*Recs)*rSz);
	char *ptr = mm_inm-1;
	do {
		Recs[nRecs].S = ++ptr;
		ptr = strchr(ptr,'\t');
		if (!ptr) {printf("ERROR: map [%u]\n",nRecs+1); exit(1);}
		if (nRecs >= rSz) Recs=realloc(Recs,sizeof(*Recs)*(rSz*=2));
		Recs[nRecs++].L = atol(ptr+1);
	} while (ptr = strchr(ptr+1,'\n'));
	printf("Parsed %u records in map\n",nRecs);

	qsort(Recs,nRecs,sizeof(*Recs),cmp);
	
	// Prepare global structure
	uint64_t whole = 2 * 128 * nRecs + 256;
	for (uint32_t i = 0; i < nRecs; ++i) whole += Recs[i].L << 1;
	uint32_t *HugeBlk = calloc(whole,sizeof(*HugeBlk)), *HgPtr = HugeBlk;
	for (uint32_t i = 0; i < nRecs; ++i)
		Recs[i].U = HgPtr + 128, Recs[i].C = Recs[i].U + Recs[i].L + 128,
		HgPtr = Recs[i].C + Recs[i].L;
	// Parse b6 and increment 
	ptr = mm_inf;
	char *next = 0, *semi, uprv = 1;
	uint64_t six = 0, N = 0;
	--nRecs;
	T *SampT = calloc(1,sizeof(*SampT)), *node = 0;
	uint32_t cnt = 0;
	do {
		if (!(cnt++ & 16383)) printf("Processed %u...\r",cnt);
		int tab = 0; 
		char *q_begin = ptr, *q_end, *r_begin, *r_end;
		while (tab < 1) tab += *ptr++ == '\t';
		q_end = ptr - 1, r_begin = ptr;
		while (tab < 2) tab += *ptr++ == '\t';
		r_end = ptr - 1; 
		while (tab < 8) tab += *ptr++ == '\t';
		long rs = atoi(ptr), re;
		while (tab < 9) tab += *ptr++ == '\t';
		re = atoi(ptr);
		Record *match = crBST_XAT(r_begin, nRecs, Recs);
		next = strchr(ptr,'\n');
		if (!match) {
			fprintf(stderr,"WARNING: couldn't find ref: ");
			fwrite(r_begin,1,r_end-r_begin,stderr);
			fputc('\n',stderr);
			continue;
		}
		char unex; 
		if (!next) unex = 1; 
		else {
			char *q = q_begin, *n = next+1;
			while (q <= q_end) if (*q++ != *n++) break;
			unex = q <= q_end;
		}
		--rs, --re;
		rs -= PAD; re += PAD;
		rs = rs < 0 ? 0 : rs; 
		re = re >= match->L ? match->L : re;
		uint32_t *C = match->C, *U = match->U;
		for (uint32_t s = rs; s < re; ++s) ++C[s];
		if (uprv && unex) for (uint32_t s = rs; s < re; ++s) ++U[s];
		
		if (DO_SAMP && (semi = strchr(q_begin,'_'))) {
			node = T_add(SampT,q_begin,semi-q_begin,&six);
			if (!node->Cc) {
				node->Uc = malloc(sizeof(*node->Uc)*(nRecs+1));
				node->Cc = malloc(sizeof(*node->Cc)*(nRecs+1));
				uint32_t *HugeBlk = calloc(whole,sizeof(*HugeBlk)), *hp = HugeBlk;
				for (uint32_t i = 0; i < nRecs; ++i)
					node->Uc[i] = hp + 128, 
					node->Cc[i] = node->Uc[i] + Recs[i].L + 128,
					hp = node->Cc[i] + Recs[i].L;
			}
			uint32_t *C = node->Cc[match-Recs], *U = node->Uc[match-Recs];
			for (uint32_t s = rs; s < re; ++s) ++C[s];
			if (uprv && unex) for (uint32_t s = rs; s < re; ++s) ++U[s];
		}
		uprv = unex;
	} while ((ptr=next) && ++ptr); 
	
	uint32_t ***GU = malloc(sizeof(*GU)*six), ***GC = malloc(sizeof(*GC)*six);
	char **SampS = malloc(sizeof(*SampS)*six);
	if (DO_SAMP) writeNodes(SampT,SampS,GU,GC), 
		printf("%lu samples found.\n",six);
	
	// Write output files; uniq and shared coverage tables
	printf("Consolidating %u records...\n",cnt);
	fprintf(outs,"#Coverage\tDataset"); fprintf(outu,"#Coverage\tDataset");
	fprintf(outbs,"#Coverage\tDataset"); fprintf(outbu,"#Coverage\tDataset");
	for (uint32_t i = 0; i < six; ++i) fprintf(outs,"\t%s",SampS[i]),
		fprintf(outu,"\t%s",SampS[i]), fprintf(outbs,"\t%s",SampS[i]),
		fprintf(outbu,"\t%s",SampS[i]);
	fputc('\n',outs), fputc('\n',outu), fputc('\n',outbs), fputc('\n',outbu);
	for (uint32_t i = 0; i <= nRecs; ++i) {
		// Shared
		uint32_t *C = Recs[i].C, *U = Recs[i].U, L = Recs[i].L;
		uint64_t tot = 0; uint32_t btot = 0;
		#pragma omp parallel for simd reduction(+:tot,btot)
		for (uint32_t k = 0; k < L; ++k) tot+=C[k], btot+=C[k] != 0;
		if (!tot) continue; // skip reference
		for (char *r = Recs[i].S; *r != '\t'; ++r) // write ref name
			fputc(*r,outs), fputc(*r,outu), fputc(*r,outbs), fputc(*r,outbu);
		double mean = (double)tot/L, bmean = (double)btot/L;
		fprintf(outbs,"\t%.4f",bmean);
		double ssd = 0, d;
		#pragma omp parallel for simd reduction(+:ssd)
		for (uint32_t k = 0; k < L; ++k) d=(double)C[k]-mean, ssd+=d*d; 
		fprintf(outs,"\t%.4f",mean > sqrt(vf*ssd/(L-1)) ? mean : -mean);
		
		// Unique
		tot = 0, btot = 0, ssd = 0;
		#pragma omp parallel for simd reduction(+:tot,btot)
		for (uint32_t k = 0; k < L; ++k) tot+=U[k], btot+=U[k] != 0;
		mean = (double)tot/L, bmean = (double)btot/L;
		fprintf(outbu,"\t%.4f",bmean);
		#pragma omp parallel for simd reduction(+:ssd)
		for (uint32_t k = 0; k < L; ++k) d = (double)U[k]-mean, ssd+=d*d;
		fprintf(outu,"\t%.4f",mean > sqrt(vf*ssd/(L-1)) ? mean : -mean);
		//fprintf(outu,"\t%.4f",sqrt(ssd/(L-1)));
		
		// Loop over all samples as well
		for (uint32_t j = 0; j < six; ++j) {
			uint32_t *C = GC[j][i], *U = GU[j][i];
			tot = 0, btot = 0;
			// Shared
			for (uint32_t k = 0; k < L; ++k) tot+=C[k], btot+=C[k] != 0;
			mean = (double)tot/L, bmean = (double)btot/L, ssd = 0;
			fprintf(outbs,"\t%.4f",bmean);
			for (uint32_t k = 0; k < L; ++k) d = (double)C[k]-mean, ssd+=d*d;
			fprintf(outs,"\t%.4f",mean > sqrt(vf*ssd/(L-1)) ? mean : -mean);
			// Unique
			tot = 0, btot = 0, ssd = 0;
			for (uint32_t k = 0; k < L; ++k) tot+=U[k], btot+=U[k] != 0;
			mean = (double)tot/L, bmean = (double)btot/L;
			fprintf(outbu,"\t%.4f",bmean);
			for (uint32_t k = 0; k < L; ++k) d = (double)U[k]-mean, ssd+=d*d;
			//fprintf(outu,"\t%.4f",sqrt(ssd/(L-1)));
			fprintf(outu,"\t%.4f",mean > sqrt(vf*ssd/(L-1)) ? mean : -mean);
		}
		fputc('\n',outu), fputc('\n',outs), fputc('\n',outbu), fputc('\n',outbs);
	}
}