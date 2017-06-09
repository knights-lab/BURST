#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void main(int argc, char *argv[]) {
	char *line = malloc(2000000000);
	FILE *file = fopen(argv[1],"r");
	if (!file || !line) {puts("Invalid file."); exit(1);}
	unsigned rep = argc > 2 ? atoi(argv[2]) : 0, 
		maxlen = 0, i = 0, mi=0, t;
	while (line = fgets(line,2000000000,file)) 
		if (++i, (t=strlen(line)) > maxlen) mi = i, maxlen = t;
	if (rep) printf("%u %u\n",i/rep,maxlen - (mi < i));
	else printf("Length of longest line: %u at line %u\n",maxlen - (mi < i),mi);
}