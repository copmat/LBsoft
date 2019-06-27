#include <stdio.h>
#include <stdlib.h>

void subcopy(double *dest, double *src, int i,int j,int k) {
  int li,lj,lk, indexSrc=0;

  for(li=0; li<64; li++) {
    for(lj=0; lj<128; lj++) {
      for(lk=0; lk<128; lk++) {
	 int indexDest = (li + i) * 512*512 + (lj+j)*512 + lk+k;
	 dest[indexDest] = src[indexSrc];
         indexSrc++;
      }
    }
  }
}

int main(int argc, char **argv) {
	int sz[3] = {64,128,128}, readSz;
	double *buf = NULL;
	double *glob = NULL;
	char fname[1000];
	int i, c, l;


	for(l=0; l<1; l++) {
	   glob = (double*) calloc( 512*512*512,sizeof(double) );
	   buf = (double*) calloc( 64*128*128,sizeof(double) );
	   if (buf==NULL || glob==NULL) {
	      fprintf(stderr,"Mem alloc failed %x %x \n", glob, buf);
	      exit(4);
	   }
	}

	for(i=0; i<32; i++) {
	   int num = (1+i)*500 + 10000;
	   int ic,jc,kc;

	   c = 0;
           for(kc=0; kc<4; kc++)
             for(jc=0; jc<4; jc++)
              for(ic=0; ic<8; ic++) {

	      sprintf(fname, "dumpHVAR.0%03d.000%03d.dat", num, c);
	      printf("Converting:%s at:%d,%d,%d \n", fname, ic,jc,kc);

	      FILE *fp=fopen(fname,"rb");
	      if (!fp) return 2;

	      for(l=0; l<1; l++) {
	         fread (&readSz, sizeof(readSz),1, fp);
	         if (readSz != 64*128*128*sizeof(double)) {
	   	    fprintf(stderr,"Found %d vs %d \n", readSz, sizeof(double)*64*128*128);
	   	    exit(4);
	         }
	         fread (buf, 64*128*128,sizeof(double), fp);
	         fseek(fp, 4, SEEK_CUR);

	         subcopy(glob, buf, ic*64,jc*128,kc*128);
	      }

	      fclose(fp);
	      c++;
	   }

	   for(l=0; l<1; l++) {
	       sprintf(fname, "out-var%d-%03d.bin", l,num);
	       printf("Into      :%s \n", fname);
	       FILE *fpOut=fopen(fname,"wb");
	       if (!fpOut) return 3;
	       fwrite (glob, 512*512*512,sizeof(double), fpOut);
	       fclose(fpOut);
	   }
	}

	return 0;
}
