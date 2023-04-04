#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the C functions you want to make available to R here
extern SEXP r_polychoric_correlation_matrix(SEXP r_input_matrix);

// Register native routine
static const R_CallMethodDef CallEntries[] = {
    {
        "r_polychoric_correlation_matrix", // Name of function call in R
        (DL_FUNC)&r_polychoric_correlation_matrix, // Name of C function
         1 // Number of arguments
    },
    {NULL, NULL, 0}
};

// Set up call in package
void R_init_latentFactoR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
