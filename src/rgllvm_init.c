#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rgllvm_efficientscorefc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_esteqallbetaic(SEXP, SEXP, SEXP);
extern SEXP _rgllvm_esteqbetaic(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_esteqc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_esteqfc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_exppartin5c(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_fastLRp_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_firderPCLaikc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_firdertest1(SEXP);
extern SEXP _rgllvm_formula5denoc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_formula5forallkc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_formula5numerc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_getwic(SEXP, SEXP);
extern SEXP _rgllvm_getzilrc(SEXP);
extern SEXP _rgllvm_PCLaikc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_PCLallc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_PCLijklalltc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_PCLijklc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_PCLsandwich(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_pic(SEXP, SEXP);
extern SEXP _rgllvm_scoref1c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_scoref2c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_scoreffinalc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_secderPCLaikallc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_secderPCLaikc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_secdertest1(SEXP);
extern SEXP _rgllvm_sfisherinforinv(SEXP, SEXP);
extern SEXP _rgllvm_sinforcell(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_sumofrowbycol(SEXP, SEXP);
extern SEXP _rgllvm_test1(SEXP);
extern SEXP _rgllvm_varofsthetac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rgllvm_xbetamatrixc(SEXP, SEXP);
extern SEXP _rgllvm_yin5c(SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(do_allocations)(void *, void *, void *, void *);
extern void F77_NAME(estimation_logistic)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getvv_terms)(void *, void *, void *);
extern void F77_NAME(store_n1m1)(void *, void *);

static const R_CallMethodDef CallEntries[] = {
  {"_rgllvm_efficientscorefc", (DL_FUNC) &_rgllvm_efficientscorefc,  9},
  {"_rgllvm_esteqallbetaic",   (DL_FUNC) &_rgllvm_esteqallbetaic,    3},
  {"_rgllvm_esteqbetaic",      (DL_FUNC) &_rgllvm_esteqbetaic,       4},
  {"_rgllvm_esteqc",           (DL_FUNC) &_rgllvm_esteqc,            9},
  {"_rgllvm_esteqfc",          (DL_FUNC) &_rgllvm_esteqfc,           9},
  {"_rgllvm_exppartin5c",      (DL_FUNC) &_rgllvm_exppartin5c,       4},
  {"_rgllvm_fastLRp_",         (DL_FUNC) &_rgllvm_fastLRp_,          6},
  {"_rgllvm_firderPCLaikc",    (DL_FUNC) &_rgllvm_firderPCLaikc,     8},
  {"_rgllvm_firdertest1",      (DL_FUNC) &_rgllvm_firdertest1,       1},
  {"_rgllvm_formula5denoc",    (DL_FUNC) &_rgllvm_formula5denoc,     4},
  {"_rgllvm_formula5forallkc", (DL_FUNC) &_rgllvm_formula5forallkc,  9},
  {"_rgllvm_formula5numerc",   (DL_FUNC) &_rgllvm_formula5numerc,   10},
  {"_rgllvm_getwic",           (DL_FUNC) &_rgllvm_getwic,            2},
  {"_rgllvm_getzilrc",         (DL_FUNC) &_rgllvm_getzilrc,          1},
  {"_rgllvm_PCLaikc",          (DL_FUNC) &_rgllvm_PCLaikc,           8},
  {"_rgllvm_PCLallc",          (DL_FUNC) &_rgllvm_PCLallc,           6},
  {"_rgllvm_PCLijklalltc",     (DL_FUNC) &_rgllvm_PCLijklalltc,      6},
  {"_rgllvm_PCLijklc",         (DL_FUNC) &_rgllvm_PCLijklc,          5},
  {"_rgllvm_PCLsandwich",      (DL_FUNC) &_rgllvm_PCLsandwich,       6},
  {"_rgllvm_pic",              (DL_FUNC) &_rgllvm_pic,               2},
  {"_rgllvm_scoref1c",         (DL_FUNC) &_rgllvm_scoref1c,          8},
  {"_rgllvm_scoref2c",         (DL_FUNC) &_rgllvm_scoref2c,          8},
  {"_rgllvm_scoreffinalc",     (DL_FUNC) &_rgllvm_scoreffinalc,      9},
  {"_rgllvm_secderPCLaikallc", (DL_FUNC) &_rgllvm_secderPCLaikallc,  6},
  {"_rgllvm_secderPCLaikc",    (DL_FUNC) &_rgllvm_secderPCLaikc,     8},
  {"_rgllvm_secdertest1",      (DL_FUNC) &_rgllvm_secdertest1,       1},
  {"_rgllvm_sfisherinforinv",  (DL_FUNC) &_rgllvm_sfisherinforinv,   2},
  {"_rgllvm_sinforcell",       (DL_FUNC) &_rgllvm_sinforcell,        5},
  {"_rgllvm_sumofrowbycol",    (DL_FUNC) &_rgllvm_sumofrowbycol,     2},
  {"_rgllvm_test1",            (DL_FUNC) &_rgllvm_test1,             1},
  {"_rgllvm_varofsthetac",     (DL_FUNC) &_rgllvm_varofsthetac,      6},
  {"_rgllvm_xbetamatrixc",     (DL_FUNC) &_rgllvm_xbetamatrixc,      2},
  {"_rgllvm_yin5c",            (DL_FUNC) &_rgllvm_yin5c,             3},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {"do_allocations",      (DL_FUNC) &F77_NAME(do_allocations),      4},
  {"estimation_logistic", (DL_FUNC) &F77_NAME(estimation_logistic), 7},
  {"getvv_terms",         (DL_FUNC) &F77_NAME(getvv_terms),         3},
  {"store_n1m1",          (DL_FUNC) &F77_NAME(store_n1m1),          2},
  {NULL, NULL, 0}
};

void R_init_rgllvm(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
