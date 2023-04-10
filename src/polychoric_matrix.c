// Headers to include
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdlib.h>
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>

// Define constants in commonly used functions (optimization process)

// Constant for hard cut-off for polychoric
#define CUT 11 // similar to {Turbofuns}

// Constants in `bsm_inverse_cdf`
const double CONST_A[6] = {-39.69683028665376, 220.9460984245205, -275.928510446969, 138.357751867269, -30.66479806614716, 2.506628277459239};
const double CONST_B[5] = {-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572};
const double CONST_C[6] = {-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783};
const double CONST_D[4] = {0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416};

// Constants in `error_function`
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429
#define P 0.3275911

// Constants in `drezner_bivariate_normal`
#define INT_NX 5
#define COR_MAX 0.7
#define BV_FAC1 0.13298076
#define BV_FAC2 0.053051647
const double DOUBLE_X[INT_NX] = {0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
const double DOUBLE_W[INT_NX] = {0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042};

// Constants in `optimize`
#define ZEPS 1e-10

// Using the Beasley-Springer-Moro algorithm to perform `qnorm` function
double bsm_inverse_cdf(double probability){

    // Check for zero and one probabilities
    if(probability == 0){
        return(-INFINITY);
    } else if(probability == 1){
        return(INFINITY);
    }

    // Initialize variables once
    double q, r, x;

    // Determine region
    if(probability >= 0.02425 && probability <= 1 - 0.02425){ // Middle

        // Define q
        q = probability - 0.50;

        // Define r
        r = q * q;

        // Define x
        x = ((((((CONST_A[0] * r + CONST_A[1]) * r + CONST_A[2]) * r + CONST_A[3]) * r + CONST_A[4]) * r + CONST_A[5]) * q) / (((((CONST_B[0] * r + CONST_B[1]) * r + CONST_B[2]) * r + CONST_B[3]) * r + CONST_B[4]) * r + 1);

    } else { // Ends

        // Define q
        if(probability < 0.02425){
            q = sqrt(-2 * log(probability));
        } else {
            q = sqrt(-2 * log(1 - probability));
        }

        // Define x
        x = (((((CONST_C[0] * q + CONST_C[1]) * q + CONST_C[2]) * q + CONST_C[3]) * q + CONST_C[4]) * q + CONST_C[5]) / ((((CONST_D[0] * q + CONST_D[1]) * q + CONST_D[2]) * q + CONST_D[3]) * q + 1);

        // Check if the sign needs to be reversed
        if(probability >= 0.02425){
            x = -x;
        }

    }

    // Return x
    return x;

}

// Find starting index
int starting_index(int* frequency){

    // Initialize start
    int start = 0;

    // Find starting index
    for (int i = 0; i < CUT; i++) {
        if (frequency[i] != 0) {
            start = i;
            break;
        }
    }

    // Return start
    return start;
}

// Find ending index
int ending_index(int* frequency){

    // Initialize end
    int end = 0;

    // Find ending index
    for (int i = CUT - 1; i >= 0; i--) {
        if (frequency[i] != 0) {
            end = i;
            break;
        }
    }

    // Return end
    return end;

}

// Define structure for return values
struct ThresholdsResult {
    int** joint_frequency;
    double* threshold_X;
    double* threshold_Y;
    double* probability_X;
    double* probability_Y;
    int num_elements;
};

// Compute thresholds
struct ThresholdsResult thresholds(int* input_data, int rows, int i, int j) {

    // Initialize iterators and temporary variables
    int k, X, Y;

    // Initialize memory space for frequencies
    int* frequency_X = (int*) calloc(CUT, sizeof(int));
    int* frequency_Y = (int*) calloc(CUT, sizeof(int));

    // Allocate memory space for table
    int** joint_frequency = (int**) malloc(CUT * sizeof(int*));
    for (k = 0; k < CUT; k++) {
        joint_frequency[k] = (int*) calloc(CUT, sizeof(int));
    }

    // Populate joint_frequency and calculate frequencies
    for (k = 0; k < rows; k++) {
        X = input_data[k + i * rows];
        Y = input_data[k + j * rows];
        joint_frequency[X][Y]++;
        frequency_X[X]++;
        frequency_Y[Y]++;
    }

    // Obtain index values
    int min_X = starting_index(frequency_X); // starting index for X
    int max_X = ending_index(frequency_X); // ending index for X
    int min_Y = starting_index(frequency_Y); // starting index for Y
    int max_Y = ending_index(frequency_Y); // ending index for Y
    int min_min = fmin(min_X, min_Y); // minimum of the minimum indices
    int max_max = fmax(max_X, max_Y); // maximum of the maximum indices
    int num_elements = max_max - min_min; // number of total possible elements

    // Initialize non-zero elements
    int* non_zero_freq_X = (int*) malloc(num_elements * sizeof(int));
    int* non_zero_freq_Y = (int*) malloc(num_elements * sizeof(int));

    // Copy non-zero elements for X
    for (k = 0; k < num_elements; k++) {
        non_zero_freq_X[k] = frequency_X[min_min + k];
        non_zero_freq_Y[k] = frequency_Y[min_min + k];
    }

    // Free the memory for frequency_X and frequency_Y
    free(frequency_X);
    free(frequency_Y);

    // Initialize memory space for probabilities and thresholds
    double* probability_X = (double*) malloc((num_elements) * sizeof(double));
    double* probability_Y = (double*) malloc((num_elements) * sizeof(double));
    double* threshold_X = (double*) malloc((num_elements) * sizeof(double));
    double* threshold_Y = (double*) malloc((num_elements) * sizeof(double));

    // Compute probabilities
    for(k = 0; k < num_elements; k++) {
        probability_X[k] = (double) non_zero_freq_X[k] / rows;
        probability_Y[k] = (double) non_zero_freq_Y[k] / rows;
    }

    // Compute cumulative sums
    for (k = 0; k < num_elements - 1; k++) {
        probability_X[k + 1] += probability_X[k];
        probability_Y[k + 1] += probability_Y[k];
    }

    // Obtain thresholds
    for (k = 0; k < num_elements; k++) {
        threshold_X[k] = bsm_inverse_cdf(probability_X[k]);
        threshold_Y[k] = bsm_inverse_cdf(probability_Y[k]);
    }

    // Free memory
    free(non_zero_freq_X);
    free(non_zero_freq_Y);

    // Create structure for return values
    struct ThresholdsResult result;
    result.joint_frequency = joint_frequency;
    result.threshold_X = threshold_X;
    result.threshold_Y = threshold_Y;
    result.probability_X = probability_X;
    result.probability_Y = probability_Y;
    result.num_elements = num_elements;

    // Return
    return result;
}

// Error function
double error_function(double x) {

  // Initialize values
  double t_value, y;

  // Determine sign
  int sign_x = x >= 0 ? 1 : -1;

  // Obtain absolute value
  x = fabs(x);

  // Set t-value
  t_value = 1 / (1 + P * x);

  // Set y
  y = 1 - (((((A5 * t_value + A4) * t_value) + A3) * t_value + A2) * t_value + A1) * t_value * exp(-x * x);

  // Add sign
  return sign_x * y;
}

// Univariate normal CDF
double univariate_normal(double x) {

  // This function is streamlined for use in this function
  //
  // With mean = 0 and sd = 1, then the z-score of x is x

  // With the error function, obtain CDF
  double cdf = 0.5 * (1 + error_function(x / sqrt(2)));

  // Return CDF
  return cdf;

}

// Bivariate normal CDF
// Implements Drezner-Wesolowsky's approximation (Drezner & Wesolowsky, 1990)
// Translated from C++ to C: https://github.com/cran/pbv/blob/master/src/pbv_rcpp_bvnorm.cpp
double drezner_bivariate_normal(double h1, double h2, double rho, double p1, double p2) {

  // Check for infinities
  if(h1 == -INFINITY || h2 == -INFINITY){
    return(0.0);
  } else if(h1 == INFINITY){
    return(p2);
  } else if(h2 == INFINITY){
    return(p1);
  }

  // Initialize probability and h3
  double bv = 0.0;
  double h3;

  // Obtain h12 and absolute correlation
  double h12 = (h1 * h1 + h2 * h2) / 2;
  double rho_abs = fabs(rho);

  // Check for correlation lower than maximum
  if(rho_abs <= COR_MAX) {

    // Initialize r1 and rr2
    double r1, rr2;

    // Compute h3
    h3 = h1 * h2;

    // Compute probability
    for (int i = 0; i < INT_NX; i++) {
      r1 = rho * DOUBLE_X[i];
      rr2 = 1 - r1 * r1;
      bv += DOUBLE_W[i] * exp((r1 * h3 - h12) / rr2) / sqrt(rr2);
    }

    // Finalize probability
    bv = p1 * p2 + rho * bv;

  } else { // Greater than maximum correlation (0.70)

    // Initialize r2 and r3
    double r2 = 1 - rho * rho;
    double r3 = sqrt(r2);

    // Reverse sign for negative correlation
    if(rho < 0.0) {
      h2 = -h2;
      p2 = 1 - p2;
    }

    // Compute h3
    h3 = h1 * h2;

    // Compute h7
    double h7 = exp(-h3 / 2.0);

    // Check for correlation less than 1
    if(rho_abs < 1.0) {

      // Set up variables
      double h6 = fabs(h1 - h2);
      double h5 = h6 * h6 / 2.0;
      h6 = h6 / r3;
      double aa = 0.5 - h3 / 8.0;
      double ab = 3.0 - 2.0 * aa * h5;

      // Compute initial probability
      bv = BV_FAC1 * h6 * ab * (1.0 - univariate_normal(h6)) - exp(-h5 / r2) * (ab + aa * r2) * BV_FAC2;

      double r1, rr;

      // Compute probability
      for (int i = 0; i < INT_NX; i++) {
        r1 = r3 * DOUBLE_X[i];
        rr = r1 * r1;
        r2 = sqrt(1.0 - rr);
        bv += -DOUBLE_W[i] * exp(-h5 / rr) * (exp(-h3 / (1.0 + r2)) / r2 / h7 - 1.0 - aa * rr);
      }

    }

    // Obtain minimum of probabilities
    bv = bv * r3 * h7 + fmin(p1, p2);

    // Adjust for negative correlation
    if (rho < 0) {
      bv = p1 - bv;
    }

  }

  // Return probability
  return bv;

}

// Estimate log-likelihood
double polychoric_log_likelihood(
    double rho, int** joint_frequency,
    double* threshold_X, double* threshold_Y,
    double* probability_X, double* probability_Y,
    int num_elements
) {

    // Set upper and lowers (thresholds and probabilities)
    double lower_ti, upper_ti, lower_pi, upper_pi;
    double lower_tj, upper_tj, lower_pj, upper_pj;

    // Initialize variables
    double log_likelihood = 0.0;
    double probability, log_prob;

      // Compute log-likelihood
      for (int i = 0; i <= num_elements; ++i) {
        for (int j = 0; j <= num_elements; ++j) {

          // Set up thresholds and probabilities
          // X
          if(i == 0){
            lower_ti = -INFINITY; lower_pi = 0;
            upper_ti = threshold_X[i]; upper_pi = probability_X[i];
          }else if(i == num_elements){
            lower_ti = threshold_X[i - 1]; lower_pi = probability_X[i - 1];
            upper_ti = INFINITY; upper_pi = 1;
          }else{
            lower_ti = threshold_X[i - 1]; lower_pi = probability_X[i - 1];
            upper_ti = threshold_X[i]; upper_pi = probability_X[i];
          }
          // Y
          if(j == 0){
            lower_tj = -INFINITY; lower_pj = 0;
            upper_tj = threshold_Y[j]; upper_pj = probability_Y[j];
          }else if(j == num_elements){
            lower_tj = threshold_Y[j - 1]; lower_pj = probability_Y[j - 1];
            upper_tj = INFINITY; upper_pj = 1;
          }else{
            lower_tj = threshold_Y[j - 1]; lower_pj = probability_Y[j - 1];
            upper_tj = threshold_Y[j]; upper_pj = probability_Y[j];
          }

          // Compute bivariate normal CDF
          probability = drezner_bivariate_normal(upper_ti, upper_tj, rho, upper_pi, upper_pj) -
                        drezner_bivariate_normal(lower_ti, upper_tj, rho, lower_pi, upper_pj) -
                        drezner_bivariate_normal(upper_ti, lower_tj, rho, upper_pi, lower_pj) +
                        drezner_bivariate_normal(lower_ti, lower_tj, rho, lower_pi, lower_pj);

          // Handle probabilities equal to zero
          if (probability == 0) {
            probability = DBL_MIN;
          }

          // Compute log probability
          log_prob = log(probability);

          // Handle infinite log probabilities
          if (isinf(log_prob)) {
            log_prob = DBL_MIN;
          }

          // Update log-likelihood
          log_likelihood += joint_frequency[i][j] * log_prob;

        }
      }

      // Return negative log-likelihood
      return -log_likelihood;

}

// Brent's method for optimization
double optimize(double (*f)(double, int**, double*, double*, double*, double*, int),
                int** joint_frequency, double* threshold_X, double* threshold_Y,
                double* probability_X, double* probability_Y,
                int num_elements,
                double lower, double upper, double tol, int max_iter
) {

    // Initialize variables for the optimization algorithm
    double a = lower;
    double b = upper;
    double c = a + (b - a) / 2.0;

    double x = c;
    double w = c;
    double v = c;

    double fx = f(x, joint_frequency, threshold_X, threshold_Y, probability_X, probability_Y, num_elements);
    double fw = fx;
    double fv = fx;

    double d = 0.0;
    double e = 0.0;

    // Iterate using Brent's method until the maximum number of iterations
    // is reached, or until the solution is found within the specified tolerance
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate the midpoint of the current interval and the tolerance
        double mid = (a + b) / 2.0;
        double tol1 = tol * fabs(x) + ZEPS;
        double tol2 = 2.0 * tol1;

        // Check if the solution is within the tolerance
        if (fabs(x - mid) <= (tol2 - (b - a) / 2.0)) {
            break;
        }

        // Initialize variables for the current iteration
        double u;
        double fu;
        bool use_parabola = false;

        // Check if the current iteration can use the parabolic fit
        if (x != w && x != v && w != v) {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;

            q = 2.0 * (q - r);

            if (q > 0.0) {
                p = -p;
            } else {
                q = -q;
            }

            double etemp = e;
            e = d;

            // If the parabolic fit is not appropriate, set the use_parabola flag to false
            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                use_parabola = false;
            } else {
                // Otherwise, set the use_parabola flag to true and calculate the parabolic minimum
                use_parabola = true;
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2) {
                    d = (mid - x >= 0.0) ? tol1 : -tol1;
                }
            }
        }

        // If the parabolic fit is not used, set the variables for the golden section step
        if (!use_parabola) {
            e = (x >= mid) ? a - x : b - x;
            d = 0.3819660 * e;
        }

        // Calculate the function value at the trial point
        u = (fabs(d) >= tol1) ? x + d : x + ((d >= 0.0) ? tol1 : -tol1);
        fu = f(u, joint_frequency, threshold_X, threshold_Y, probability_X, probability_Y, num_elements);

        // Update the intervals based on the function values at the trial point
        if (fu <= fx) {
            if (u >= x) {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }

    // Return the optimal solution
    return x;
}

// Compute polychoric correlation
double polychoric(int* input_data, int rows, int i, int j) {

    // Obtain joint frequency table, probability_X, and probability_Y from thresholds function
    struct ThresholdsResult thresholds_result = thresholds(input_data, rows, i, j);

    // Initialize parameters for optimization
    double lower = -1.0;
    double upper = 1.0;
    double tol = 1e-04; // same tolerance as `optimize` in R
    int max_iter = 100;

    // Perform optimization
    double rho_optimum = optimize(
        polychoric_log_likelihood, thresholds_result.joint_frequency,
        thresholds_result.threshold_X, thresholds_result.threshold_Y,
        thresholds_result.probability_X, thresholds_result.probability_Y,
        thresholds_result.num_elements,
        lower, upper, tol, max_iter
    );

    // Free memory
    for (int i = 0; i <= thresholds_result.num_elements; i++) {
        free(thresholds_result.joint_frequency[i]);
    }
    free(thresholds_result.joint_frequency);
    free(thresholds_result.threshold_X);
    free(thresholds_result.threshold_Y);
    free(thresholds_result.probability_X);
    free(thresholds_result.probability_Y);

    // Return
    return rho_optimum;

}

// The updated polychoric_correlation_matrix function
double** polychoric_correlation_matrix(int* input_data, int rows, int cols) {

    // Initialize iterators
    int i, j;

    // Allocate memory for polychoric_matrix
    double** polychoric_matrix = malloc(cols * sizeof(double*));
    for (i = 0; i < cols; i++) {
      polychoric_matrix[i] = malloc(cols * sizeof(double)); // Change rows to cols
    }

    // Perform polychoric correlations over the input_matrix
    for (i = 0; i < cols; i++) {

        // Loop over other variables
        for (j = 0; j < cols; j++) {

            // Fill diagonal
            if (i == j) {

                polychoric_matrix[i][j] = 1;

            }else if (i > j) { // Fill lower triangle

                // Compute correlation
                double correlation = polychoric(
                    input_data, rows, i, j
                );

                // Add to matrix
                polychoric_matrix[i][j] = correlation;

                // Fill opposite of triangle
                polychoric_matrix[j][i] = correlation;

            }
        }
    }

    // return 0

    return polychoric_matrix;
}

SEXP r_polychoric_correlation_matrix(SEXP r_input_matrix) {

    int rows = nrows(r_input_matrix);
    int cols = ncols(r_input_matrix);

    // Call the C function
    double** c_result = polychoric_correlation_matrix(INTEGER(r_input_matrix), rows, cols);

    // Convert the C result to an R matrix
    SEXP r_result = PROTECT(allocMatrix(REALSXP, cols, cols));
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < cols; j++) {
            REAL(r_result)[i + cols * j] = c_result[i][j];
        }
    }

    // Free the C result memory
    for (int i = 0; i < cols; i++) {
        free(c_result[i]);
    }
    free(c_result);

    UNPROTECT(1);
    return r_result;

}

// To compile as a shared library on Linux
//
// R CMD SHLIB polychoric_matrix.c -o polychoric_matrix.so
//
// To compile as a shared library on Max (through Linux)
//
// R CMD SHLIB polychoric_matrix.c -o polychoric_matrix.dylib
//
// To compile as a shared library on Windows
//
// R CMD SHLIB polychoric_matrix.c -o polychoric_matrix.dll

//    // Print matrix
//    for (i = 0; i < cols; i++) {
//        for (j = 0; j < cols; j++) { // Changed from 'rows' to 'cols'
//            printf("%f ", matrix[i][j]);
//        }
//        printf("\n");
//    }
//
//    // Free allocated matrix memory
//    for (i = 0; i < cols; i++) {
//        free(matrix[i]);
//    }
//    free(matrix);
//
//    // Print vector
//    for (j = i; i < cols; i++) {
//        printf("%f ", vector[i]);
//    }
//
//    // Free allocated matrix memory
//    free(vector);
//
//    // Print value
//    printf("%f ", value);
//
