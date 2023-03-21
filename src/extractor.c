#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{
    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        printf("Usage: EXTRACTOR {OUTPUT} {COUNT}\n");
        return 0;
    }
    FILE *rp = fopen(argv[1], "r");
    if (rp == NULL) {
        printf("Check output filename.");
        return 1;
    }
    char filename[1024], line[1024], tmp_line[1024];
    char *ptr = strrchr(argv[1], '/');
    char *format1, *format2;
    if (ptr == NULL) {
        format1 = strtok(argv[1], ".");
    } else {
        format1 = strtok(ptr + 1, ".");
    }
    format2 = strtok(NULL, "\n");
    int count = atoi(argv[2]);
    sprintf(filename, "%s_%d.%s", format1, count, format2);
    FILE *wp = fopen(filename, "w");
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
