#ifndef _HIC_
#define _HIC_

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <R.h>


extern const gsl_rng_type * T;
extern gsl_rng * r;
extern double ubound,lbound;
extern int NK;
extern gsl_matrix_int *Ktable;
extern double epsilon;
extern int HL;
extern gsl_matrix *efl;
extern gsl_matrix *gc;
extern gsl_matrix *map;
extern gsl_matrix *ZZ;
extern gsl_matrix *iZ;
extern gsl_matrix *IH;
extern gsl_matrix *gamma_hat;
extern gsl_vector *dd;
extern gsl_vector *loglambda;

extern double ml,mgc;


extern double MVNsigmax,MVNsigmay,MVNrho;
extern double nI2;


/* number of fit coefficients */
//#define NCOEFFS  11
/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
//#define NBREAK   (NCOEFFS - 2)

extern int ncoeffs;
extern double prior;
extern double anchor;
extern int nbreak;
extern double CEIL;

extern gsl_bspline_workspace *Bw;
extern gsl_matrix *Bx;
extern gsl_vector *B;
extern gsl_vector *Cx;
extern gsl_matrix *Cov;

extern double MaxKnot;



#endif /* #ifndef _HIC_ */
