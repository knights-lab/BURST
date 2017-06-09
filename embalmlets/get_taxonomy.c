#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

// gettax indir outfile.txt
void main(int argc, char *argv[]) {
	if (argc != 3) {puts("Usage: gettax indir outfile.txt"); exit(1);}
	DIR *dir = opendir(argv[1]);
	FILE *out = fopen(argv[2],"wb");
	if (!dir || !out) {fputs("I/O error\n",stderr); exit(2);}
	struct dirent *file;
	char fullPath[UINT16_MAX] = {0};
	int len = strlen(argv[1]);
	strncpy(fullPath, argv[1], len);
	fullPath[len++] = '/';
	size_t flexSz = 16000, numFiles = 0, numRec = 0; //, sz=flexSz;
	char *dump = malloc(flexSz), *init, *accession;
	while (file = readdir(dir)) if (strstr(file->d_name, ".gbff\0")) {
		strncpy(fullPath+len,file->d_name,256);
		printf("Considering '%s' ",file->d_name);
		FILE *in = fopen(fullPath,"rb");
		
		fseek(in,0,SEEK_END);
		size_t sz = (size_t)ftell(in);
		rewind(in);
		if (sz >= flexSz) 
			flexSz = sz+1,
			dump = realloc(dump,flexSz);
		
		fread(dump,1,sz,in);
		dump[sz] = 0;
		fclose(in);

		uint32_t UAs = 0; init = dump;
		while ((accession = strstr(init, "VERSION     ") + 12) - 12) {
			if (!(accession-12)) {fputs("Accession not found.\n",stderr); exit(2);}
			char *offset = strstr(accession, "  ORGANISM  ") + 12;
			if (!(offset-12)) {fputs("Organism not found.\n",stderr); exit(2);}
			char *ending = strstr(offset, "\nREFERENCE");
			if (!ending) ending = strstr(offset, "\nCOMMENT");
			if (!ending) {fputs("Reference/comment not found.\n",stderr); exit(2);}
			char *eoa = strchr(accession, '.'); 
			char *eol = strchr(offset,'\n'); 
			if (!eoa || !eol) {fputs("Cannot find line terminator\n",stderr); exit(2);}
			*eoa = *eol = 0;
			char *next = eol, *curLn;
			*eol = 0, *(ending-1) = ';';
			for (curLn = next + 1; *curLn == ' '; ++curLn);
			fprintf(out,"%s\t", accession);
			// check if n-line organism name
			next = strchr(curLn,'\n');
			*next = 0;
			while (!strchr(curLn,';')) {
				*eol = ' ';
				strncpy(eol + 1, curLn, next - curLn + 1);
				eol += 1 + next - curLn;
				curLn = next + 1;
				next = strchr(curLn,'\n');
				*next = 0;
				for (; *curLn == ' '; ++curLn);
			}
			
			while (curLn < ending) {
				fprintf(out,"%s",curLn);
				for (curLn = next + 1; *curLn == ' '; ++curLn);
				curLn -= curLn < ending;
				next = strchr(curLn, '\n');
				*next = 0;
			}
			// process offset (split after second space)
			char *candidate = strstr(offset,"Candidatus ");
			int spaces = 0, sep; for (char *sp = offset + (candidate? 11: 0); *sp; ++sp) 
				if (*sp==' ') if (++spaces == 2) *sp = 0, sep = sp - offset;
			if (spaces >=2) fprintf(out," %s; %s\n", offset, offset + sep + 1);
			else fprintf(out," %s\n",offset);
			++UAs;
			init = next + 1;
		}
		printf(" [%u records]%s\n", UAs, UAs > 1 ? " NOTE: MULTIPLE RECORDS DETECTED" : "");
		++numFiles, numRec += UAs;
	}
	printf("Considered all .gbff files in directory [%u files, %u records]\n",
		numFiles, numRec);
}