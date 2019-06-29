#include <stdio.h>
#include <stdlib.h>

#define N 512
#define PROCX 8
#define PROCY 4
#define PROCZ 4

#define NUMFILES 51
#define INIFILES 500
#define STEPFILES 500

#define NUMVAR 2


void subcopy(double *dest, double *src, int i, int j, int k, int ni,
	     int nj, int nk)
{
  int li, lj, lk, indexSrc = 0;

  for (lk = 0; lk < nk; lk++) {
    for (lj = 0; lj < nj; lj++) {
      for (li = 0; li < ni; li++) {
	int indexDest = (lk + k) * N * N + (lj + j) * N + li + i;
	dest[indexDest] = src[indexSrc];
	indexSrc++;
      }
    }
  }
}

int main(int argc, char **argv)
{
  int readSz;
  double *buf[NUMVAR] = {NULL};
  double *glob[NUMVAR] = {NULL};
  char fname[1000];
  int i, c, l;

  const int lnx = N / PROCX;
  const int lny = N / PROCY;
  const int lnz = N / PROCZ;


  for (l = 0; l < NUMVAR; l++) {
    glob[l] = (double *) calloc(N * N * N, sizeof(double));
    buf[l] = (double *) calloc(lnx * lny * lnz, sizeof(double));
    if (buf[l] == NULL || glob[l] == NULL) {
      fprintf(stderr, "Mem alloc failed var:%d %x %x \n", l, glob, buf);
      exit(4);
    }
  }

  for (i = 0; i < NUMFILES; i++) {
    int num = (1 + i) * STEPFILES + INIFILES;
    int ic, jc, kc;

    c = 0;
    for (kc = 0; kc < PROCZ; kc++)
      for (jc = 0; jc < PROCY; jc++)
	for (ic = 0; ic < PROCX; ic++, c++) {

	  sprintf(fname, "dumpHVAR.%06d.%06d.dat", num, c);
	  printf("Converting:%s at:%d,%d,%d \n", fname, ic, jc, kc);

	  FILE *fp = fopen(fname, "rb");
	  if (!fp) {
	    perror("fopen");
	    exit(2);
	  }

	  for (l = 0; l < NUMVAR; l++) {
	    fread(&readSz, sizeof(readSz), 1, fp);
	    if (readSz != lnx * lny * lnz * sizeof(double)) {
	      fprintf(stderr, "Found %d vs %d \n", readSz, sizeof(double) * lnx * lny * lnz);
	      exit(4);
	    }
	    readSz = fread(buf[l], sizeof(double), lnx * lny * lnz, fp);
	    if (readSz != lnx * lny * lnz) {
	      fprintf(stderr, "Read %d vs %d \n", readSz, sizeof(double) * lnx * lny * lnz);
	      exit(5);
	    }
	    fseek(fp, 4, SEEK_CUR);

	    subcopy(glob[l], buf[l], ic * lnx, jc * lny, kc * lnz, lnx, lny, lnz);
	  }

	  fclose(fp);
	}

    for (l = 0; l < NUMVAR; l++) {
      sprintf(fname, "out-var%d-%06d.bin", l, num);
      printf("Into      :%s \n", fname);
      FILE *fpOut = fopen(fname, "wb");
      if (!fpOut) {
	perror("fopen out file");
	exit(2);
      }

      fwrite(glob[l], N * N * N, sizeof(double), fpOut);
      fclose(fpOut);
    }
  }

  return 0;
}
