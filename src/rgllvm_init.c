#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(do_allocations)(void *, void *, void *, void *);
extern void F77_NAME(estimation_logistic)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getvv_terms)(void *, void *, void *);
extern void F77_NAME(store_n1m1)(void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"do_allocations",      (DL_FUNC) &F77_NAME(do_allocations),      4},
  {"estimation_logistic", (DL_FUNC) &F77_NAME(estimation_logistic), 7},
  {"getvv_terms",         (DL_FUNC) &F77_NAME(getvv_terms),         3},
  {"store_n1m1",          (DL_FUNC) &F77_NAME(store_n1m1),          2},
  {NULL, NULL, 0}
};

void R_init_rgllvm(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
