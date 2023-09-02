#ifndef __EXTRACTOR_H__
#define __EXTRACTOR_H__

int copy_files(char *destination, char *source);
int append_files(char *destination, char *source);
int extract_files(char *filename, int count);
#endif
