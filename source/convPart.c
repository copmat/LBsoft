#include <stdio.h>

#define NUMFILES 51
#define INIFILES 500
#define STEPFILES 500


int main(int argc, char **argv) {
	int sz = 3800;
	double buf[6][sz];
	char fname[1000];
	int i, c, l;


	for(i=0; i<NUMFILES; i++) {
 	   int num = (1+i)*STEPFILES + INIFILES;
	   sprintf(fname, "dumpParticles.%06d.dat", num);
	   printf("Converting:%s \n", fname);

	   FILE *fp=fopen(fname,"rb");
	   sprintf(fname, "part-%06d.csv", num);
	   FILE *fpOut=fopen(fname,"w");
	   if (!fp) return 2;
	   if (!fpOut) return 3;

	   for(l=0; l<6; l++) {
	      fseek(fp, 4, SEEK_CUR);
	      fread (buf+l, sizeof(buf[l]),1, fp);
	      fseek(fp, 4, SEEK_CUR);
	   }
	   for(l=0; l<sz; l++) {
		fprintf(fpOut,"%20.15g,%20.15g,%20.15g,%20.15g,%20.15g,%20.15g \n", buf[0][l], buf[1][l], buf[2][l],
							    buf[3][l], buf[4][l], buf[5][l]);
	   }

	   fclose(fp);
	   fclose(fpOut);
	}
	return 0;
}
