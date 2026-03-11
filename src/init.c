// new version:

// src/init.c
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*** .C routines from nz.c ***/
void nz (double*, double*, int*, int*, int*, double*, double*, double*,
         int*, double*, int*, int*, double*);
void nz2(double*,          int*, int*, int*, double*, double*, double*,
         int*, double*, int*, int*, double*);
void bn (double*, double*, int*, int*, int*, double*, double*, double*,
         int*, double*, int*, int*, double*);
void bn2(double*,          int*, int*, int*, double*, double*, double*,
         int*, double*, int*, int*, double*);

static const R_CMethodDef CEntries[] = {
  {"nz",  (DL_FUNC) &nz,  13},
  {"nz2", (DL_FUNC) &nz2, 12},
  {"bn",  (DL_FUNC) &bn,  13},
  {"bn2", (DL_FUNC) &bn2, 12},
  {NULL, NULL, 0}
};

/*** Rcpp-generated .Call wrappers from RcppExports.cpp ***/
/* forward-declare exactly what appears in RcppExports.cpp */
SEXP _tRexCAP_rcpp_gen_angle(void);
SEXP _tRexCAP_Rx(SEXP);
SEXP _tRexCAP_Ry(SEXP);
SEXP _tRexCAP_Rz(SEXP);
SEXP _tRexCAP_dist_trans(SEXP, SEXP, SEXP);
SEXP _tRexCAP_beta_func(SEXP, SEXP, SEXP, SEXP);
SEXP _tRexCAP_project_func(SEXP, SEXP);
SEXP _tRexCAP_nll(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _tRexCAP_grad(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _tRexCAP_APG(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _tRexCAP_minimizer(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_tRexCAP_rcpp_gen_angle", (DL_FUNC) &_tRexCAP_rcpp_gen_angle, 0},
  {"_tRexCAP_Rx",             (DL_FUNC) &_tRexCAP_Rx,             1},
  {"_tRexCAP_Ry",             (DL_FUNC) &_tRexCAP_Ry,             1},
  {"_tRexCAP_Rz",             (DL_FUNC) &_tRexCAP_Rz,             1},
  {"_tRexCAP_dist_trans",     (DL_FUNC) &_tRexCAP_dist_trans,     3},
  {"_tRexCAP_beta_func",      (DL_FUNC) &_tRexCAP_beta_func,      4},
  {"_tRexCAP_project_func",   (DL_FUNC) &_tRexCAP_project_func,   2},
  {"_tRexCAP_nll",            (DL_FUNC) &_tRexCAP_nll,            5},
  {"_tRexCAP_grad",           (DL_FUNC) &_tRexCAP_grad,           5},
  {"_tRexCAP_APG",            (DL_FUNC) &_tRexCAP_APG,            7},
  {"_tRexCAP_minimizer",      (DL_FUNC) &_tRexCAP_minimizer,      9},
  {NULL, NULL, 0}
};

/*** single init for the package ***/
void R_init_tRexCAP(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}



// old version:
// // src/init.c
// #include <R.h>
// #include <Rinternals.h>
// #include <R_ext/Rdynload.h>

// /* --- prototypes for your .C routines from nz.c --- */
// void nz (double*, double*, int*, int*, int*, double*, double*, double*,
//          int*, double*, int*, int*, double*);
// void nz2(double*,          int*, int*, int*, double*, double*, double*,
//          int*, double*, int*, int*, double*);
// void bn (double*, double*, int*, int*, int*, double*, double*, double*,
//          int*, double*, int*, int*, double*);
// void bn2(double*,          int*, int*, int*, double*, double*, double*,
//          int*, double*, int*, int*, double*);

// /* .C registration table */
// static const R_CMethodDef CEntries[] = {
//   {"nz",  (DL_FUNC) &nz,  13},
//   {"nz2", (DL_FUNC) &nz2, 12},
//   {"bn",  (DL_FUNC) &bn,  13},
//   {"bn2", (DL_FUNC) &bn2, 12},
//   {NULL, NULL, 0}
// };

// /* forward-declare the Rcpp registrar that compileAttributes() generates */
// void R_init_tRexCAP_RcppExports(DllInfo *dll);

// /* single init for the whole package */
// void R_init_tRexCAP(DllInfo *dll) {
//   /* register your .C routines */
//   R_registerRoutines(dll, CEntries, NULL, NULL, NULL);

//   /* let Rcpp register its .Call routines */
//   R_init_tRexCAP_RcppExports(dll);

//   /* no dynamic symbol lookup */
//   R_useDynamicSymbols(dll, FALSE);
// }
