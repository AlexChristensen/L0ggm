/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%% Derivatives in C %%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Updated 13.03.2026 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>

/* atan_derivative
 *
 * R equivalent:
 *   lambda * (gamma * (gamma + 2 / pi)) / (gamma^2 + x^2)
 */
SEXP atan_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_)
{
    int n = length(x_);
    SEXP out = PROTECT(allocVector(REALSXP, n));

    const double *x      = REAL(x_);
    const double lambda  = REAL(lambda_)[0];
    const double gamma   = REAL(gamma_)[0];

    double gamma2        = gamma * gamma;
    double gamma_scaled  = gamma * (gamma + 2.0 / M_PI);

    for (int i = 0; i < n; i++)
        REAL(out)[i] = lambda * gamma_scaled / (gamma2 + x[i] * x[i]);

    UNPROTECT(1);
    return out;
}


/* exp_derivative
 *
 * R equivalent:
 *   lambda * (1 / gamma) * exp(-(abs(x) / gamma))
 */
SEXP exp_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_)
{
    int n = length(x_);
    SEXP out = PROTECT(allocVector(REALSXP, n));

    const double *x     = REAL(x_);
    const double lambda = REAL(lambda_)[0];
    const double gamma  = REAL(gamma_)[0];

    double inv_gamma    = 1.0 / gamma;

    for (int i = 0; i < n; i++)
        REAL(out)[i] = lambda * inv_gamma * exp(-fabs(x[i]) * inv_gamma);

    UNPROTECT(1);
    return out;
}


/* gumbel_derivative
 *
 * R equivalent:
 *   (lambda / (1 - exp(-1))) * (1 / gamma) * exp(-gamma_x - exp(-gamma_x))
 *   where gamma_x = abs(x) / gamma
 */
SEXP gumbel_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_)
{
    int n = length(x_);
    SEXP out = PROTECT(allocVector(REALSXP, n));

    const double *x     = REAL(x_);
    const double lambda = REAL(lambda_)[0];
    const double gamma  = REAL(gamma_)[0];

    /* Pre-compute the scaling constant: lambda / ((1 - exp(-1)) * gamma) */
    double scale        = lambda / ((1.0 - exp(-1.0)) * gamma);

    for (int i = 0; i < n; i++) {
        double gx       = fabs(x[i]) / gamma;
        REAL(out)[i]    = scale * exp(-gx - exp(-gx));
    }

    UNPROTECT(1);
    return out;
}


/* log_derivative
 *
 * R equivalent:
 *   lambda / ((gamma + abs(x)) * log(1 + 1 / gamma))
 */
SEXP log_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_)
{
    int n = length(x_);
    SEXP out = PROTECT(allocVector(REALSXP, n));

    const double *x     = REAL(x_);
    const double lambda = REAL(lambda_)[0];
    const double gamma  = REAL(gamma_)[0];

    /* log(1 + 1/gamma) is constant across observations */
    double log_term     = log1p(1.0 / gamma);

    for (int i = 0; i < n; i++)
        REAL(out)[i] = lambda / ((gamma + fabs(x[i])) * log_term);

    UNPROTECT(1);
    return out;
}


/* weibull_derivative
 *
 * R equivalent:
 *   lambda * (shape / gamma) * (abs_x / gamma)^(shape - 1) * exp(-(abs_x / gamma)^shape)
 *   where abs_x = pmax(abs(x), .Machine$double.eps)
 */
SEXP weibull_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_, SEXP shape_)
{
    int n = length(x_);
    SEXP out = PROTECT(allocVector(REALSXP, n));

    const double *x     = REAL(x_);
    const double lambda = REAL(lambda_)[0];
    const double gamma  = REAL(gamma_)[0];
    const double shape  = REAL(shape_)[0];

    double scale        = lambda * shape / gamma;

    for (int i = 0; i < n; i++) {
        double abs_x    = fmax(fabs(x[i]), DBL_EPSILON);
        double xg       = abs_x / gamma;
        double xg_s     = pow(xg, shape);
        REAL(out)[i]    = scale * pow(xg, shape - 1.0) * exp(-xg_s);
    }

    UNPROTECT(1);
    return out;
}
