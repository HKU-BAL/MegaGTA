#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

int build_graph(int argc, char** argv);
int find_start(int argc, char **argv);
int search(int argc, char **argv);

void show_help(const char *program_name) {
	fprintf(stderr, "Usage: %s <sub_program> [sub options]\n"
					"    sub-programs:\n"
	                "       build_graph           build the SdBG\n"
	                "       find_start            find starting kmers\n"
	                "       search                A* search\n",
	                program_name);
}

AutoMaxRssRecorder recorder;

int main(int argc, char **argv) {
	if (argc < 2) {
		show_help(argv[0]);
		exit(1);
	}

	if (strcmp(argv[1], "build") == 0) {
		return build_graph(argc - 1, argv + 1);
	} else if (strcmp(argv[1], "find") == 0) {
		return find_start(argc - 1 , argv + 1);
	} else if (strcmp(argv[1], "search") == 0) {
		return search(argc - 1, argv + 1);
	} else {
		show_help(argv[0]);
		exit(1);
	}

	return 0;
}