#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{
    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        printf("Usage: EXTRACTOR {FILENAME} {COUNT}\n");
        return 0;
    }
    FILE *rp = fopen(argv[1], "r");
    if (rp == NULL) {
        printf("Check input filename.");
        return 1;
    }
    char filename[1024], line[1024], tmp_line[1024], *ptr;
    char *format1 = strtok(argv[1],".");
    char *format2 = strtok(NULL, "\n");
    int count = atoi(argv[2]);
    sprintf(filename, "%s_%d.%s", format1, count, format2);
    FILE *wp = fopen(filename, "w");
    int flag = 0;
    while (fgets(tmp_line, 1024, rp) != NULL) {
        memcpy(line, tmp_line, sizeof(char) * 1024);
        if ((strchr(line, '_')) != NULL) {
            int tmp_count = atoi(strtok(tmp_line, "_"));
            if (count == tmp_count) {
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
