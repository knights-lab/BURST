#define HOWTO "Usage: ripINSD in.xml out.fasta out.tax [minlev] [cutlast] [noSp.]"
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

void main(int argc, char *argv[]) {
	if (argc < 4 || argc > 7) {puts(HOWTO); exit(1);}
	FILE *in = fopen(argv[1],"rb"), *out = fopen(argv[2],"wb"),
		*tax = fopen(argv[3],"wb");
	int minlev = argc > 4 ? atoi(argv[4]) : 0, 
		cutlast = argc > 5 ? atoi(argv[5]) : 0,
		noSp = argc > 6 ? atoi(argv[6]) : 0;
	if (!in || !out || !tax) {fputs("I/O error\n",stderr); exit(2);}
	char *line = malloc(INT32_MAX), *field, *eol, 
		*acc = malloc(UINT16_MAX), *n = malloc(UINT16_MAX);
	size_t ns = 0; while (++ns, line = fgets(line,INT32_MAX,in)) {
		while (!(field = strstr(line, "<GBSeq_primary-accession>"))) 
			if (!(line = fgets(line,INT32_MAX,in))) {
				printf("Exiting. Parsed %llu records.\n",ns-1);
				exit(0);
			}
		eol = field += 25;
		while (*eol != '<') ++eol;
		*eol = 0;
		strncpy(acc, field, MIN(UINT16_MAX,eol-field+1));

		while (!(field = strstr(line, "<GBSeq_organism>"))) 
			if (!(line = fgets(line,INT32_MAX,in))) {
				fprintf(stderr,"ERROR: incomplete record %llu.\n",ns);
				exit(2);
			}
		eol = field += 16;
		while (*eol != '<') ++eol;
		*eol = 0;
		strncpy(n, field, MIN(UINT16_MAX, eol-field+1));
		char *candidate = strstr(n,"Candidatus ");
		int spaces = 0, semis = 0, s = eol-field; 
		for (char *sp = candidate ? candidate + 11 : n; *sp; ++sp) 
			if (*sp==' ') if (++spaces == 2) *sp = 0, s = sp - n;

		while (!(field = strstr(line, "<GBSeq_taxonomy>"))) 
			if (!(line = fgets(line,INT32_MAX,in))) {
				fprintf(stderr,"ERROR: incomplete record %llu.\n",ns);
				exit(2);
			}
		eol = field += 16;
		while (*eol != '<') semis+=*eol++==';';
		*eol = 0;
		if (semis < minlev || (noSp && ((n[s-1] == '.' && n[s-2] == 'p' && 
				n[s-3]=='s') || !strcmp(n+s-10," bacterium"))))
			{printf("WARNING: '%s' omitted (%s)\n", acc, n); --ns; continue;}
		else if (spaces >=2 && !cutlast) 
			fprintf(tax,"%s\t%s; %s; %s\n", acc, field, n, n+s+1);
		else fprintf(tax,"%s\t%s; %s\n", acc, field, n);

		while (!(field = strstr(line, "<GBSeq_sequence>"))) 
			if (!(line = fgets(line,INT32_MAX,in))) {
				fprintf(stderr,"ERROR: incomplete record %llu.\n",ns);
				exit(2);
			}
		eol = field += 16;
		while (*eol != '<') *eol = *eol > 96 ? *eol - 32 : *eol, ++eol;
		*eol = 0;
		fprintf(out,">%s\n%s\n",acc,field);
	}
}