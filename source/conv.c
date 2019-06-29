#include <stdio.h>

#define NUMFILES 1
#define INIFILES 0
#define STEPFILES 500

#define NUMVAR 2


int main(int argc, char **argv) {
	int sz[3] = { 64,128,128};
	size_t szPlane = sz[0] * sz[1];
	double buf[szPlane];
	char fname[1000];
	int i, c, l;


	for(i=0; i<NUMFILES; i++) {
	   int num = (1+i)*STEPFILES + INIFILES;

	   sprintf(fname, "dumpHVAR.000%03d.000000.dat", num);
	   printf("Converting:%s \n", fname);

	   FILE *fp=fopen(fname,"rb");
	   sprintf(fname, "out-%03d.bin", num);
	   printf("Into      :%s \n", fname);
	   FILE *fpOut=fopen(fname,"wb");
	   if (!fp) return 2;
	   if (!fpOut) return 3;

	   for(l=0; l<NUMVAR; l++) {
	      fseek(fp, 4, SEEK_CUR);
	      for(c=sz[2]; c--; c>0) {
	         fread (buf, sizeof(buf),1, fp);
	         fwrite(buf, sizeof(buf),1, fpOut);
	      }
	      fseek(fp, 4, SEEK_CUR);
	   }

	   fclose(fp);
	   fclose(fpOut);
	}
	return 0;
}
