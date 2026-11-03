#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the C functions you want to make available to R here
extern SEXP r_polychoric_correlation_matrix(SEXP r_input_matrix, SEXP r_empty_method, SEXP r_empty_value, SEXP r_rows, SEXP r_cols);
extern SEXP atan_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_);
extern SEXP exp_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_);
extern SEXP gumbel_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_);
extern SEXP log_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_);
extern SEXP weibull_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_, SEXP shape_);


// Register native routine
static const R_CallMethodDef CallEntries[] = {

    {
        "r_polychoric_correlation_matrix", // Name of function call in R
        (DL_FUNC)&r_polychoric_correlation_matrix, // Name of C function
         5 // Number of arguments
    },
    {
        "atan_derivative_c",
        (DL_FUNC)&atan_derivative_c,
         3
    },
    {
        "exp_derivative_c",
        (DL_FUNC)&exp_derivative_c,
         3
    },
    {
        "gumbel_derivative_c",
        (DL_FUNC)&gumbel_derivative_c,
         3
    },
    {
        "log_derivative_c",
        (DL_FUNC)&log_derivative_c,
         3
    },
    {
        "weibull_derivative_c",
        (DL_FUNC)&weibull_derivative_c,
         4
    },
    {NULL, NULL, 0}

};

// Set up call in package
void R_init_L0ggm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
