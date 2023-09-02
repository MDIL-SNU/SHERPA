#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "files.h"


int main(int argc, char *argv[])
{
    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        printf("Usage: EXTRACTOR {OUTPUT} {COUNT}\n");
        return 0;
    }
    extract_files(argv[1], atoi(argv[2]));
    return 0;
}
