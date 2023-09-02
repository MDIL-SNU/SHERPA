#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "files.h"


int copy_files(char *destination, char *source)
{
    FILE *rp = fopen(source, "r");
    FILE *wp = fopen(destination, "w");
    if (rp == NULL || wp == NULL) {
        printf("Check filenames.");
        return 1;
    }
    char line[65536];
    while (fgets(line, 65536, rp) != NULL) {
        fputs(line, wp);
    }
    fclose(rp);
    fclose(wp);
    return 0;
}


int append_files(char *destination, char *source)
{
    FILE *rp = fopen(source, "r");
    FILE *wp = fopen(destination, "a");
    if (rp == NULL || wp == NULL) {
        printf("Check filenames.");
        return 1;
    }
    char line[65536];
    while (fgets(line, 65536, rp) != NULL) {
        fputs(line, wp);
    }
    fclose(rp);
    fclose(wp);
    return 0;
}


int extract_files(char *filename, int count)
{
    FILE *rp = fopen(filename, "r");
    if (rp == NULL) {
        printf("Check output filename.");
        return 1;
    }
    char tmp_filename[1024], line[1024], tmp_line[1024];
    char *ptr = strrchr(filename, '/');
    char *format1, *format2;
    if (ptr == NULL) {
        format1 = strtok(filename, ".");
    } else {
        format1 = strtok(ptr + 1, ".");
    }
    format2 = strtok(NULL, "\n");
    sprintf(tmp_filename, "%s_%d.%s", format1, count, format2);
    FILE *wp = fopen(tmp_filename, "w");
    int flag = 0;
    while (fgets(tmp_line, 1024, rp) != NULL) {
        strcpy(line, tmp_line);
        if ((strchr(line, '_')) != NULL) {
            if (count == atoi(strtok(tmp_line, "_"))) {
                flag = 1;
            } else {
                if (flag > 0) {
                    break;
                }
            }
        }
        if (flag > 0) {
            fputs(line, wp);
        }
    }
    fclose(rp);
    fclose(wp);
    return 0;
}
