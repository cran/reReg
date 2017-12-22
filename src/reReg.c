#include <R.h>
#include <Rmath.h>
#include <math.h>

void plLambda(double *sl, double *tij, double *yi, double *weights,
	      int *n, int *N, // n is the length of sl
	      // output
	      double *res) {
  int i, j;
  double dl = 0, Rl = 0;
  for (j = 0; j < *n; j++) {
    res[j] = 1;
    dl = 0;
    Rl = 0;
    for (i = 0; i < *N; i++) {
      if (sl[j] >= tij[i] && sl[j] <= yi[i]) {
	Rl = Rl + weights[i];
      }
      if (sl[j] == tij[i] && tij[i] != yi[i]) {
	dl = dl + weights[i];
      }
    }
    if (Rl > 0) {
      if (j == 0)
	res[j] = (dl / Rl);
      if (j > 0)
	res[j] = res[j - 1] + (dl / Rl);
    }
  }
  res;
}


/* void plLambda(double *sl, double *tij, double *yi, double *weights, */
/* 	      int *n, int *N, // n is the length of sl */
/* 	      // output */
/* 	      double *res) { */
/*   int i, j; */
/*   double dl = 0, Rl = 0; */
/*   for (j = 0; j < *n; j++) { */
/*     res[j] = 1; */
/*     dl = 0; */
/*     Rl = 0; */
/*     for (i = 0; i < *N; i++) { */
/*       if (sl[j] >= tij[i] && sl[j] <= yi[i]) { */
/* 	Rl = Rl + weights[i]; */
/*       } */
/*       if (sl[j] == tij[i] && tij[i] != yi[i]) { */
/* 	dl = dl + weights[i]; */
/*       } */
/*     } */
/*     if (Rl > 0) { */
/*       if (j == 0) */
/* 	res[j] = (1 - dl / Rl); */
/*       if (j > 0) */
/* 	res[j] = res[j - 1] * (1 - dl / Rl); */
/*     } */
/*   } */
/*   res; */
/* } */

void sarm1(double *X, double *Lambda, double *weights, double *gamma, 
	   int *mt, int *n, int *p, int *B, 
	   double *res) {
  int i, r, b, k; 
  double xr; 
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      xr = 0;
      for (k = 0; k < *p; k++) {
	xr += X[i + k * *n] * gamma[k];
      }
      for (r = 0; r < *p; r++) {
	if (Lambda[i] > 0) 
	  res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * (mt[i]/Lambda[i] - exp(xr));
	if (Lambda[i] == 0)
	  res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * (-1 * exp(xr));
      } // end r
    } // end i
  } // end b
  res;
}

// SARM 3
/* void sarm2(double *X, double *T, double *Y, double *weights, double *lambda,  */
/* 	   double *mt, int *n, int *p, int *B, */
/* 	   double *res) { */
/*   double *nu = Calloc(*p, double); */
/*   int i, j, r, b; */
/*   double de; */
/*   for (b = 0; b < *B; b++) { */
/*     for (i = 0; i < *n; i++) { */
/*       for (r = 0; r < *p; r++) { */
/* 	nu[r] = 0; */
/*       } */
/*       de = 0.0; */
/*       for (j = 0; j < *n; j++) { */
/* 	if (T[i] <= Y[j] & lambda[j] > 0) { */
/* 	  for (r = 0; r < *p; r++) { */
/* 	    // nu[r] += mt[j] * X[j + r * *n] / lambda[j]; */
/* 	    nu[r] += X[j + r * *n] / lambda[j]; */
/* 	  } */
/* 	  // de += mt[j] / lambda[j]; */
/* 	  de += 1 / lambda[j]; */
/* 	} */
/*       } */
/*       for (r = 0; r < *p; r ++) { */
/*       	res[r] += X[i + r * *n] - nu[r] / de; */
/*       	// res[r] += (X[i + r * *n] - X[j + r * *n]) * de[r]; */
/*       } */
/*     } */
/*     } */
/*   Free(nu); */
/*   res; */
/* } */

// SARM 2

void sarm2(double *X, double *T, double *Y, double *weights, 
	   int *n, int *p, int *B,
	   double *res) {
  double *nu = Calloc(*p, double);
  int i, j, r, b;
  double de;
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      for (r = 0; r < *p; r++) {
	nu[r] = 0;
      }
      de = 0.0;
      for (j = 0; j < *n; j++) {
	if (T[i] <= Y[j] & T[i] >= T[j]) {
	  for (r = 0; r < *p; r++) {
	    nu[r] += X[j + r * *n];
	  }
	  de += 1;
	}
      }
      for (r = 0; r < *p; r ++) {
      	res[r] += X[i + r * *n] - nu[r] / de;
      	// res[r] += (X[i + r * *n] - X[j + r * *n]) * de[r];
      }
    }
    }
  Free(nu);
  res;
}

/* void alphaEq1(double *X, double *Lambda, double *weights,  */
/* 	      int *mt, int *n, int *p, int *B,  */
/* 	      double *res) { */
/*   int i, j, r, b, iId = 0, jId = 0;  */
/*   for (b = 0; b < *B; b++) { */
/*     for (i = 0; i < *n; i++) { */
/*       for (j = 0; j < *n; j++) { */
/* 	for (r = 0; r < *p; r++) { */
/* 	  if (Lambda[i] != 0 && Lambda[j] != 0) */
/* 	    res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * weights[j + b * *n] * (mt[i] / Lambda[i] -  mt[j] / Lambda[j]); */
/* 	  if (Lambda[i] != 0 && Lambda[j] == 0) */
/* 	    res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * weights[j + b * *n] * (mt[i] / Lambda[i]); */
/* 	  if (Lambda[i] == 0 && Lambda[j] != 0) */
/* 	    res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * weights[j + b * *n] * (0 - mt[j] / Lambda[j]); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   res; */
/* }		  */


void alphaEq(double *X, double *Lambda, int *mt, int *n, int *p, double *res) {
  int i, j, r; 
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      for (r = 0; r < *p; r++) {
	if (Lambda[i] != 0 && Lambda[j] != 0) {
	  res[r] += X[i + r * *n] * (mt[i] / Lambda[i] -  mt[j] / Lambda[j]);
	}
	if (Lambda[i] != 0 && Lambda[j] == 0) {
	  res[r] += X[i + r * *n] * (mt[i] / Lambda[i]);
	}
	if (Lambda[i] == 0 && Lambda[j] != 0) {
	  res[r] += X[i + r * *n] * (0 - mt[j] / Lambda[j]);
	}
      }
    }
  }
  res;
}		 

void betaEst(double *Y, double *X, double *delta, double *z, double *weights,
	     int *n, int *p, int *B, 
	     // output
	     double *res) {
  int i, j, r, b;
  double *nu = Calloc(*p, double); 
  double de;
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      if (delta[i] != 0) {
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	for (j = 0; j < *n; j++) {
	  if (Y[i] <= Y[j]) {
	    for (r = 0; r < *p; r++) {
	      nu[r] += weights[j + b * *n] * z[j] * X[j + r * *n];
	    }
	    de += weights[j + b * *n] * z[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de == 0) {
	    res[r + b * *p] += weights[i + b * *n] *  X[i + r * *n];
	  }
	  if (de != 0) {
	    res[r + b * *p] += weights[i + b * *n] * (X[i + r * *n] - (nu[r] / de));
	  }
	}
      } // end delta
    }
  }
  Free(nu);
  res;
}

void HWb(double *Y, double *X, double *delta, double *z, double *weights,
	     int *n, int *p, int *B, 
	     // output
	 double *res) {
  int i, j, r, b;
  double *nu = Calloc(*p, double); 
  double de;
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      if (delta[i] != 0) {
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	for (j = 0; j < *n; j++) {
	  if (Y[i] <= Y[j]) {
	    for (r = 0; r < *p; r++) {
	      nu[r] += exp(weights[j + b * *n]) * z[j] * X[j + r * *n];
	    }
	    de += exp(weights[j + b * *n]) * z[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de == 0) {
	    res[r + b * *p] += X[i + r * *n];
	  }
	  if (de != 0) {
	    res[r + b * *p] += X[i + r * *n] - (nu[r] / de);
	  }
	}
      } // end delta
    }
  }
  Free(nu);
  res;
}
