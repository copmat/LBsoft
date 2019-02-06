END {
	nx = 3
	ny = 3
	nz = 3
	eps = 0.0000000001
	OFMT="%15.8f"
	for (i=0; i<nx; i++) {
	 x = 5.8 + i*9.5 + eps
	 for (j=0; j<ny; j++) {
	  y = 6.5 + j*9.5 + eps
	  for (k=0; k<nz; k++) {
		z = 7.1 + k*9.5 + eps
		print "C       ", x,y,z, eps,eps,eps
	  }
	 }
	}
}
	  	
