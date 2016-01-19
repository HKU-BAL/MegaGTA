#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "definitions.h"
#include "utils.h"

int build_lib(int argc, char** argv);
int build_graph(int argc, char** argv);
int read_stat(int argc, char** argv);
int find_start(int argc, char **argv);
int search(int argc, char **argv);
int filter_by_len(int argc, char **argv);
int translate(int argc, char **argv);
int main_assemble(int argc, char **argv);


void show_help(const char *program_name) {
	fprintf(stderr, "Usage: %s <sub_program> [sub options]\n"
					"    sub-programs:\n"
                    "       buildlib              build read library\n"
	                "       buildgraph            build the SdBG\n"
                    "       denovo                de novo assemble contigs from SDBG\n"
	                "       findstart             find starting kmers\n"
	                "       search                A* search\n"
                    "       dumpversion           dump MEGAHIT-GT version\n"
                    "       readstat              get sequence stat from fastq/a files\n"
                    "       filterbylen           filter contigs by length\n"
                    "       translate             translate DNA to Protein\n",
	                program_name);
}

int main(int argc, char **argv) {
	if (argc < 2) {
		show_help(argv[0]);
		exit(1);
	}
    
    if (strcmp(argv[1], "buildlib") == 0) {
        return build_lib(argc - 1, argv + 1); 
    } else if (strcmp(argv[1], "buildgraph") == 0) {
		return build_graph(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "readstat") == 0) {
        return read_stat(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "filterbylen") == 0) {
        return filter_by_len(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "findstart") == 0) {
		return find_start(argc - 1 , argv + 1);
	} else if (strcmp(argv[1], "search") == 0) {
		return search(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "translate") == 0) {
		return translate(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "dumpversion") == 0) {
        printf("%s\n", PACKAGE_VERSION);
        return 0;
    } else if (strcmp(argv[1], "denovo") == 0) {
        return main_assemble(argc - 1, argv + 1);
    } else {
		show_help(argv[0]);
		exit(1);
	}

	return 0;
}
