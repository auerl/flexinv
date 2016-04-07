#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#define PERMS 0666

int openresfile_()
{
char filenm[132];
int fd;

  printf("Open waveform derived C residual file\n");
  printf("File name:\n");
  scanf("%s", filenm);
  printf("File=%s\n",filenm);
  fd=open(filenm,O_RDONLY); 
  return fd;
}

int opennewbin_()
{
char filenm[132];
int fd;

  printf("Open output file\n");
  printf("File name:\n");
  scanf("%s", filenm);
  printf("File=%s",filenm);
  fd =creat(filenm,0666);
  return fd;
}
