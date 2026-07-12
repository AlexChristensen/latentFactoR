#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the C functions you want to make available to R here
extern SEXP r_ziggurat(SEXP n, SEXP r_seed);
extern SEXP r_xoshiro_seeds(SEXP n, SEXP r_seed);
extern SEXP r_xoshiro_uniform(SEXP n, SEXP r_seed);
extern SEXP r_xoshiro_shuffle(SEXP r_vector, SEXP r_seed);
extern SEXP r_xoshiro_weighted_shuffle(SEXP r_vector, SEXP r_prob, SEXP r_seed);
extern SEXP r_xoshiro_shuffle_replace(SEXP r_vector, SEXP r_size, SEXP r_seed);


// Register native routine
static const R_CallMethodDef CallEntries[] = {

    {
        "r_ziggurat",
        (DL_FUNC)&r_ziggurat,
         2
    },
    {
        "r_xoshiro_uniform",
        (DL_FUNC)&r_xoshiro_uniform,
         2
    },
    {
        "r_xoshiro_seeds",
        (DL_FUNC)&r_xoshiro_seeds,
         2
    },
    {
        "r_xoshiro_shuffle",
        (DL_FUNC)&r_xoshiro_shuffle,
         2
    },
    {
        "r_xoshiro_weighted_shuffle",
        (DL_FUNC)&r_xoshiro_weighted_shuffle,
         3
    },
    {
        "r_xoshiro_shuffle_replace",
        (DL_FUNC)&r_xoshiro_shuffle_replace,
         3
    },
    {NULL, NULL, 0}

};

// Set up call in package
void R_init_latentFactoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
