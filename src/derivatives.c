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

    const double gamma2 = gamma * gamma;
    const double scale  = lambda * gamma * (gamma + 2.0 / M_PI);

    for (int i = 0; i < n; i++)
        REAL(out)[i] = scale / (gamma2 + x[i] * x[i]);

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

    const double *x = REAL(x_);
    const double lambda = REAL(lambda_)[0];
    const double gamma = REAL(gamma_)[0];

    const double inv_gamma = 1.0 / gamma;
    const double scale = lambda * inv_gamma;

    for (int i = 0; i < n; i++)
        REAL(out)[i] = scale * exp(-fabs(x[i]) * inv_gamma);

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
    const double scale = lambda / ((1.0 - exp(-1.0)) * gamma);

    for (int i = 0; i < n; i++) {
        double gx = fabs(x[i]) / gamma;
        REAL(out)[i] = scale * exp(-gx - exp(-gx));
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
    const double log_term = log1p(1.0 / gamma);

    for (int i = 0; i < n; i++)
        REAL(out)[i] = lambda / ((gamma + fabs(x[i])) * log_term);

    UNPROTECT(1);
    return out;
}


/* weibull_derivative
 *
 * R equivalent:
 *   abs_x = pmax(abs(x), .Machine$double.eps)
 *
 *   if (shape > 1):
 *     peak       = gamma * ((shape - 1) / shape)^(1 / shape)
 *     peak_value = (shape / gamma) * (peak / gamma)^(shape - 1) * exp(-(peak / gamma)^shape)
 *
 *     lambda * peak_value                                                                      if abs_x <= peak
 *     lambda * (shape / gamma) * (abs_x / gamma)^(shape - 1) * exp(-(abs_x / gamma)^shape)   otherwise
 *
 *   if (shape <= 1):
 *     lambda * (shape / gamma) * (abs_x / gamma)^(shape - 1) * exp(-(abs_x / gamma)^shape)
 */
SEXP weibull_derivative_c(SEXP x_, SEXP lambda_, SEXP gamma_, SEXP shape_)
{
    int n = length(x_);
    SEXP out = PROTECT(allocVector(REALSXP, n));

    const double *x = REAL(x_);
    const double lambda = REAL(lambda_)[0];
    const double gamma = REAL(gamma_)[0];
    const double shape = REAL(shape_)[0];
    
    /* Pre-compute constants */
    const double scale = lambda * (shape / gamma);
    const double shape_one = shape - 1.0;
    
    if (shape > 1.0) {
    
    	/* Pre-compute peak parameters for shape > 1 */
        const double peak = gamma * pow(shape_one / shape, 1.0 / shape);
        const double pg = peak / gamma;
        const double pg_s = pow(pg, shape);
        const double peak_value = pow(pg, shape_one) * exp(-pg_s);
        double val;
        
        for (int i = 0; i < n; i++) {
        
        	/*  Compute absolute without zero check
        	 *  since peak will always be greater than
        	 *  zero when shape > 1 
        	 */
        	double abs_x = fabs(x[i]);

		    if (abs_x <= peak) {
		        val = peak_value;
		    } else {
		    	double xg = abs_x / gamma;
		    	double xg_s = pow(xg, shape);
		        val = pow(xg, shape_one) * exp(-xg_s);
		    }

		    REAL(out)[i] = scale * val;
		}
        
    } else {
    
		for (int i = 0; i < n; i++) {
		    double abs_x = fmax(fabs(x[i]), DBL_EPSILON);
		    double xg = abs_x / gamma;
		    double xg_s = pow(xg, shape);

		    REAL(out)[i] = scale * pow(xg, shape_one) * exp(-xg_s);
		}
    
    }

    UNPROTECT(1);
    return out;
    
}
