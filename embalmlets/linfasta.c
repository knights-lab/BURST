#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define LINELEN 2000000000
void main(int argc, char *argv[]) {
	if (argc != 3) {puts("Usage: infile.bad.fasta outfile.good.fasta"); exit(1);}
	char *filename = argv[1], *outfile = argv[2];
	FILE *in = fopen(filename,"rb"), *out = fopen(outfile,"wb");
	char *line = malloc(LINELEN);
	if (!line) {fputs("Out of memory\n",stderr); exit(3);}
	fputs(fgets(line,LINELEN,in),out);
	while (line = fgets(line,LINELEN,in)) {
		if (*line == '>') fputc('\n',out), fputs(line,out);
		else *strchr(line,'\n') = 0, fputs(line,out);
	}
	fputs("\n",out);
}