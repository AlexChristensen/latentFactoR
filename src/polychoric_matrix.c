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

// Constants in `bsm_inverse_cdf`
#define CONST_A {-39.69683028665376, 220.9460984245205, -275.928510446969, 138.357751867269, -30.66479806614716, 2.506628277459239}
#define CONST_B {-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572}
#define CONST_C {-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783}
#define CONST_D {0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416}

// Constants in `error_function`
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429
#define P 0.3275911

// Constants in `genz_bivariate_normal`
// #define PI 3.14159265358979323846
// Use M_PI instead

//// Constants in `drezner_bivariate_normal`
//#define INT_NX 5
//#define DOUBLE_X {0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992}
//#define DOUBLE_W {0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042}
//#define COR_MAX 0.7
//#define BV_FAC1 0.13298076
//#define BV_FAC2 0.053051647

// Constants in `optimize`
#define ZEPS 1e-10

// Obtain maximum categories from a specific column in a matrix
int find_max_categories(int* input_data, int rows, int cols, int col_idx) {

    // Initialize variables
    int arr_max = input_data[col_idx * rows];

    // Find maximum categories
    for (int i = 1; i < rows; i++) {
        if (input_data[i + col_idx * rows] > arr_max) {
            arr_max = input_data[i + col_idx * rows];
        }
    }

    // Return maximum value
    return arr_max;

}

// Define cumulative sum function
void cumulative_sum(double* arr, int n) {

    // Initialize result
    arr[0] = 0;

    // Compute cumulative sum
    for (int i = 1; i < n; i++) {
        arr[i] += arr[i-1];
    }

}

// Obtain joint frequency table
int** joint_frequency_table(int *X, int max_X, int *Y, int max_Y, int n){

    // Allocate memory space for table
    int** joint_frequency = (int**) malloc((max_X + 1) * sizeof(int*));
    for (int i = 0; i <= max_X; i++) {
        joint_frequency[i] = (int*) calloc(max_Y + 1, sizeof(int));
    }

    // Populate table
    for (int i = 0; i < n; i++) {
        joint_frequency[X[i]][Y[i]]++;
    }

    // Return joint frequency table
    return joint_frequency;
}

// Using the Beasley-Springer-Moro algorithm to perform `qnorm` function
double bsm_inverse_cdf(double probability){

    // Check for zero and one probabilities
    if(probability == 0){
        return(-INFINITY);
    } else if(probability == 1){
        return(INFINITY);
    }

    // Define constants
    double const_a[6] = CONST_A;
    double const_b[5] = CONST_B;
    double const_c[6] = CONST_C;
    double const_d[4] = CONST_D;

    // Initialize variables once
    double q, r, x;

    // Determine region
    if(probability >= 0.02425 && probability <= 1 - 0.02425){ // Middle

        // Define q
        q = probability - 0.50;

        // Define r
        r = q * q;

        // Define x
        x = ((((((const_a[0] * r + const_a[1]) * r + const_a[2]) * r + const_a[3]) * r + const_a[4]) * r + const_a[5]) * q) / (((((const_b[0] * r + const_b[1]) * r + const_b[2]) * r + const_b[3]) * r + const_b[4]) * r + 1);

    } else { // Ends

        // Define q
        if(probability < 0.02425){
            q = sqrt(-2 * log(probability));
        } else {
            q = sqrt(-2 * log(1 - probability));
        }

        // Define x
        x = (((((const_c[0] * q + const_c[1]) * q + const_c[2]) * q + const_c[3]) * q + const_c[4]) * q + const_c[5]) / ((((const_d[0] * q + const_d[1]) * q + const_d[2]) * q + const_d[3]) * q + 1);

        // Check if the sign needs to be reversed
        if(probability >= 0.02425){
            x = -x;
        }

    }

    // Return x
    return x;

}

// Define structure for return values
struct ThresholdsResult {
    int** joint_frequency;
    double* threshold_X;
    double* threshold_Y;
    double* probability_X;
    double* probability_Y;
};

// Compute thresholds
struct ThresholdsResult thresholds(int *X, int max_X, int *Y, int max_Y, int n){

    // Obtain joint frequency table
    int** joint_frequency = joint_frequency_table(X, max_X, Y, max_Y, n);

    // Initialize memory space for frequencies
    int* frequency_X = (int*) calloc(max_X, sizeof(int));
    int* frequency_Y = (int*) calloc(max_Y, sizeof(int));

    // Obtain frequencies
    for (int i = 0; i <= max_X; i++) {
        for (int j = 0; j <= max_Y; j++) {
            frequency_X[i] += joint_frequency[i][j];
            frequency_Y[j] += joint_frequency[i][j];
        }
    }

    // Initialize memory space for probabilities and thresholds
    double* probability_X = (double*) malloc((max_X + 1) * sizeof(double));
    double* probability_Y = (double*) malloc((max_Y + 1) * sizeof(double));
    double* threshold_X = (double*) malloc((max_X) * sizeof(double));
    double* threshold_Y = (double*) malloc((max_Y) * sizeof(double));

    // Avoid using two for loops (if possible)
    if(max_X == max_Y){

        // Compute probabilities
        for(int i = 0; i <= max_X; i++) {
            probability_X[i] = (double) frequency_X[i] / n;
            probability_Y[i] = (double) frequency_Y[i] / n;
        }

    }else{

        // Compute probabilities
        for(int i = 0; i <= max_X; i++) {
            probability_X[i] = (double) frequency_X[i] / n;
        }
        for(int j = 0; j <= max_Y; j++) {
            probability_Y[j] = (double) frequency_Y[j] / n;
        }

    }

    // Compute cumulative probabilities
    cumulative_sum(probability_X, max_X);
    cumulative_sum(probability_Y, max_Y);

    // Avoid using two for loops (if possible)
    if(max_X == max_Y){

        // Obtain thresholds
        for (int i = 1; i < max_X; i++) {
            threshold_X[i] = bsm_inverse_cdf(probability_X[i]);
            threshold_Y[i] = bsm_inverse_cdf(probability_Y[i]);
        }

    }else{

        // Obtain thresholds
        for (int i = 1; i < max_X; i++) {
            threshold_X[i] = bsm_inverse_cdf(probability_X[i]);
        }
        for (int j = 1; j < max_Y; j++) {
            threshold_Y[j] = bsm_inverse_cdf(probability_Y[j]);
        }

    }

    // Free memory
    free(frequency_X);
    free(frequency_Y);

    // Create structure for return values
    struct ThresholdsResult result;
    result.joint_frequency = joint_frequency;
    result.threshold_X = threshold_X;
    result.threshold_Y = threshold_Y;
    result.probability_X = probability_X;
    result.probability_Y = probability_Y;

    // Return
    return result;
}

// Error function
double error_function(double x) {

  // Set constants
  double a1 = A1;
  double a2 = A2;
  double a3 = A3;
  double a4 = A4;
  double a5 = A5;
  double p = P;

  // Initialize values
  double t_value, y;

  // Determine sign
  int sign_x = x >= 0 ? 1 : -1;

  // Obtain absolute value
  x = fabs(x);

  // Set t-value
  t_value = 1 / (1 + p * x);

  // Set y
  y = 1 - (((((a5 * t_value + a4) * t_value) + a3) * t_value + a2) * t_value + a1) * t_value * exp(-x * x);

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
// Implements Genz's approximation
double genz_bivariate_normal(double t1, double t2, double rho, double p1, double p2) {

  // Translated to C from MATLAB code: https://www.math.wsu.edu/faculty/genz/software/software.html

    // key
    // dh = t1
    // dk = t2
    // phid(-dh) = p1
    // phid(-dk) = p2

    // Initialize bivariate normal probability
    double bvp = 0.0;

    // Initialize size
    int size_x = 0;
    double w[20];
    double x[20];

    // Check for INFINITY
    if(t1 == INFINITY || t2 == INFINITY){
        return 0;
    }else if(t1 == -INFINITY){
        if(t2 == -INFINITY){
            return 1;
        }else{
            return 1 - p2;
        }
    }else if(t2 == -INFINITY){
        return 1 - p1;
    }else if(rho == 0){
        return (1 - p1) * (1 - p2);
    }else{

        // key
        // h = t1
        // k = t2
        // bvn = bvp

        // Initialize variables
        double tp = 2 * PI;
        double t12 = t1 * t2;
        double abs_rho = fabs(rho);

        // Check for smaller correlation
        if(abs_rho < 0.30){ // Gauss Legendre points and weights n = 6

            // Set size
            size_x = 6;

            // Initialize w and x (performs the [w w] and [1-x 1+x] here)
            w[0] = 0.1713244923791705; w[1] = 0.3607615730481384;
            w[2] = 0.4679139345726904; w[3] = 0.1713244923791705;
            w[4] = 0.3607615730481384; w[5] = 0.4679139345726904;

            x[0] = 0.0675304857968478; x[1] = 0.338790613533735;
            x[2] = 0.761380813916803; x[3] = 1.93246951420315;
            x[4] = 1.66120938646626; x[5] = 1.2386191860832;

            // There's probably a better way to do this step without
            // R getting angry...

//            // Initialize w and x (performs the [w w] and [1-x 1+x] here)
//            double w[] = {
//                0.1713244923791705, 0.3607615730481384, 0.4679139345726904,
//                0.1713244923791705, 0.3607615730481384, 0.4679139345726904
//            };
//            double x[] = {
//                0.0675304857968478, 0.338790613533735, 0.761380813916803,
//                1.93246951420315, 1.66120938646626, 1.2386191860832
//            };

        }else if(abs_rho < 0.75){ // Gauss Legendre points and weights n = 12

            // Set size
            size_x = 12;

            // Initialize w and x (performs the [w w] and [1-x 1+x] here)
            w[0] = 0.04717533638651177; w[1] = 0.1069393259953183;
            w[2] = 0.1600783285433464; w[3] = 0.2031674267230659;
            w[4] = 0.2334925365383547; w[5] = 0.2491470458134029;
            w[6] = 0.04717533638651177; w[7] = 0.1069393259953183;
            w[8] = 0.1600783285433464; w[9] = 0.2031674267230659;
            w[10] = 0.2334925365383547; w[11] = 0.2491470458134029;

            x[0] = 0.0184393657532809; x[1] = 0.095882743629525;
            x[2] = 0.230097325805695; x[3] = 0.412682045713383;
            x[4] = 0.63216850100182; x[5] = 0.874766591488531;
            x[6] = 1.98156063424672; x[7] = 1.90411725637047;
            x[8] = 1.7699026741943; x[9] = 1.58731795428662;
            x[10] = 1.36783149899818; x[11] = 1.12523340851147;

            // There's probably a better way to do this step without
            // R getting angry...

//            // Initialize w and x
//            double w[] = {
//                0.04717533638651177, 0.1069393259953183, 0.1600783285433464, 0.2031674267230659, 0.2334925365383547, 0.2491470458134029,
//                0.04717533638651177, 0.1069393259953183, 0.1600783285433464, 0.2031674267230659, 0.2334925365383547, 0.2491470458134029
//            };
//            double x[] = {
//                0.0184393657532809, 0.095882743629525, 0.230097325805695, 0.412682045713383, 0.63216850100182, 0.874766591488531,
//                1.98156063424672, 1.90411725637047, 1.7699026741943, 1.58731795428662, 1.36783149899818, 1.12523340851147
//            };

        }else{ // Gauss Legendre points and weights n = 20

            // Set size
            size_x = 20;

            // Initialize w and x (performs the [w w] and [1-x 1+x] here)
            w[0] = 0.01761400713915212; w[1] = 0.04060142980038694;
            w[2] = 0.06267204833410906; w[3] = 0.08327674157670475;
            w[4] = 0.1019301198172404; w[5] = 0.1181945319615184;
            w[6] = 0.1316886384491766; w[7] = 0.1420961093183821;
            w[8] = 0.1491729864726037; w[9] = 0.1527533871307259;
            w[10] = 0.01761400713915212; w[11] = 0.04060142980038694;
            w[12] = 0.06267204833410906; w[13] = 0.08327674157670475;
            w[14] = 0.1019301198172404; w[15] = 0.1181945319615184;
            w[16] = 0.1316886384491766; w[17] = 0.1420961093183821;
            w[18] = 0.1491729864726037; w[19] = 0.1527533871307259;

            x[0] = 0.00687140081490512; x[1] = 0.0360280727220862;
            x[2] = 0.0877655717486741; x[3] = 0.160883028177781;
            x[4] = 0.253668093539849; x[5] = 0.363946319273485;
            x[6] = 0.489132998049173; x[7] = 0.62629391128458;
            x[8] = 0.772214148858355; x[9] = 0.923473478866503;
            x[10] = 1.99312859918509; x[11] = 1.96397192727791;
            x[12] = 1.91223442825133; x[13] = 1.83911697182222;
            x[14] = 1.74633190646015; x[15] = 1.63605368072652;
            x[16] = 1.51086700195083; x[17] = 1.37370608871542;
            x[18] = 1.22778585114165; x[19] = 1.0765265211335;

            // There's probably a better way to do this step without
            // R getting angry...

//            // Initialize w and x
//            double w[] = {
//                0.01761400713915212, 0.04060142980038694, 0.06267204833410906, 0.08327674157670475, 0.1019301198172404,
//                0.1181945319615184, 0.1316886384491766, 0.1420961093183821, 0.1491729864726037, 0.1527533871307259,
//                0.01761400713915212, 0.04060142980038694, 0.06267204833410906, 0.08327674157670475, 0.1019301198172404,
//                0.1181945319615184, 0.1316886384491766, 0.1420961093183821, 0.1491729864726037, 0.1527533871307259
//            };
//            double x[] = {
//                0.00687140081490512, 0.0360280727220862, 0.0877655717486741, 0.160883028177781, 0.253668093539849,
//                0.363946319273485, 0.489132998049173, 0.62629391128458, 0.772214148858355, 0.923473478866503,
//                1.99312859918509, 1.96397192727791, 1.91223442825133, 1.83911697182222, 1.74633190646015,
//                1.63605368072652, 1.51086700195083, 1.37370608871542, 1.22778585114165, 1.0765265211335
//            };

        }

        if(abs_rho < 0.925){

            // Initialize variables
            double hs = (t1 * t1 + t2 * t2) / 2;
            double asr = asin(rho) / 2;
            double sn;

            // Compute bvp
            for(int i = 0; i < size_x; i++){
                sn = sin(asr * x[i]);
                bvp += exp((sn * t12 - hs) / (1 - pow(sn, 2))) * w[i];
            }

            // Update bvp
            bvp = bvp * asr / tp +  (1 - p1) * (1 - p2);

        }else{

            // Check for negative correlation
            if(rho < 0){ // Reverse signs on second variable
                t2 = -t2;
                p2 = 1 - p2;
                t12 = -t12;
            }

            if(abs_rho < 1){

                // Initialize variables
                double as = 1 - pow(rho, 2);
                double a = sqrt(as);
                double bs = pow((t1 - t2), 2);
                double asr = -(bs/as + t12) / 2;
                double sp;
                double c = (4 - t12) / 8;
                double d = (12 - t12) /80;

                // Check asr
                if(asr > -100){
                    bvp += a * exp(asr) * (1 - c * (bs - as) * (1 - d * bs) / 3 + c * d * pow(as, 2));
                }

                // Check t12
                if(t12 > -100){
                    double b = sqrt(bs);
                    sp = sqrt(tp) * univariate_normal(-b / a);
                    bvp -= exp(-t12/2) * sp * b * (1 - c * bs * (1 - d * bs) / 3);
                }

                // Initialize and update variables
                a = a / 2;
                double xs, rs, ep;
                double inner_a = 0.0;

                // Compute asr
                for(int i = 0; i < size_x; i++){

                    // MATLAB code vectorizes this `for` loop
                    xs = pow((a * x[i]), 2);
                    asr = -(bs / xs + t12) / 2;

                    if(asr > -100){
                        sp = (1 + c * xs * (1 + 5 * d * xs));
                        rs = sqrt(1 - xs);
                        ep = exp(-(t12 / 2) * xs / pow(1 + rs, 2)) / rs;
                        inner_a += (exp(asr) * (sp - ep)) * w[i];
                    }

                }

                // Compute bvp with inner_a
                bvp = (a * inner_a - bvp) / tp;

            }

            if(rho > 0){

                // Update bvp
                if(t1 > t2){
                    bvp += 1 - p1;
                }else{
                    bvp += 1 - p2;
                }

            }else{

                // Initialize L
                double L;

                if(t1 < 0){
                    L = p2 - p1;
                }else{
                    L = (1 - p1) - (1 - p2);
                }

                // Update bvp
                bvp = L - bvp;

            }

        }

    }

    return fmax(0, fmin(1, bvp));

}

// Estimate log-likelihood
double polychoric_log_likelihood(
    double rho, int** joint_frequency,
    double* threshold_X, double* threshold_Y,
    double* probability_X, double* probability_Y,
    int max_X, int max_Y
) {

  // Set upper and lowers (thresholds and probabilities)
    double lower_ti[max_X], upper_ti[max_X], lower_pi[max_X], upper_pi[max_X];
    double lower_tj[max_Y], upper_tj[max_Y], lower_pj[max_Y], upper_pj[max_Y];

    // Avoid need for double loop (if possible)
    if(max_X == max_Y){

        // Set up bounds for X and Y
        for (int i = 0; i < max_X; i++) {
            if(i == 0){
                // X
                lower_ti[i] = -INFINITY;
                upper_ti[i] = threshold_X[i + 1];
                lower_pi[i] = 0;
                upper_pi[i] = probability_X[i + 1];
                // Y
                lower_tj[i] = -INFINITY;
                upper_tj[i] = threshold_Y[i + 1];
                lower_pj[i] = 0;
                upper_pj[i] = probability_Y[i + 1];
            }else if(i == max_X - 1){
                // X
                lower_ti[i] = threshold_X[i];
                upper_ti[i] = INFINITY;
                lower_pi[i] = probability_X[i];
                upper_pi[i] = 1;
                // Y
                lower_tj[i] = threshold_Y[i];
                upper_tj[i] = INFINITY;
                lower_pj[i] = probability_Y[i];
                upper_pj[i] = 1;
            }else{
                // X
                lower_ti[i] = threshold_X[i];
                upper_ti[i] = threshold_X[i + 1];
                lower_pi[i] = probability_X[i];
                upper_pi[i] = probability_X[i + 1];
                // Y
                lower_tj[i] = threshold_Y[i];
                upper_tj[i] = threshold_Y[i + 1];
                lower_pj[i] = probability_Y[i];
                upper_pj[i] = probability_Y[i + 1];
            }
        }


    }else{

        // Set up bounds for X
        for (int i = 0; i < max_X; i++) {
            if(i == 0){
                lower_ti[i] = -INFINITY;
                upper_ti[i] = threshold_X[i + 1];
                lower_pi[i] = 0;
                upper_pi[i] = probability_X[i + 1];
            }else if(i == max_X - 1){
                lower_ti[i] = threshold_X[i];
                upper_ti[i] = INFINITY;
                lower_pi[i] = probability_X[i];
                upper_pi[i] = 1;
            }else{
                lower_ti[i] = threshold_X[i];
                upper_ti[i] = threshold_X[i + 1];
                lower_pi[i] = probability_X[i];
                upper_pi[i] = probability_X[i + 1];
            }
        }

        // Set up bounds for Y
        for (int j = 0; j < max_Y; j++) {
            if(j == 0){
                lower_tj[j] = -INFINITY;
                upper_tj[j] = threshold_Y[j + 1];
                lower_pj[j] = 0;
                upper_pj[j] = probability_Y[j + 1];
            }else if(j == max_Y - 1){
                lower_tj[j] = threshold_Y[j];
                upper_tj[j] = INFINITY;
                lower_pj[j] = probability_Y[j];
                upper_pj[j] = 1;
            }else{
                lower_tj[j] = threshold_Y[j];
                upper_tj[j] = threshold_Y[j + 1];
                lower_pj[j] = probability_Y[j];
                upper_pj[j] = probability_Y[j + 1];
            }
        }

    }

  // Initialize variables
  double log_likelihood = 0.0;
  double probability, log_prob;

  // Compute log-likelihood
  for (int i = 0; i < max_X; ++i) {
    for (int j = 0; j < max_Y; ++j) {

      // Following MATLAB code: bvnu(xl,yl,r) - bvnu(xu,yl,r) - bvnu(xl,yu,r) + bvnu(xu,yu,r);

      // Compute bivariate normal CDF
      probability = genz_bivariate_normal(lower_ti[i], lower_tj[j], rho, lower_pi[i], lower_pj[j]) -
                    genz_bivariate_normal(upper_ti[i], lower_tj[j], rho, upper_pi[i], lower_pj[j]) -
                    genz_bivariate_normal(lower_ti[i], upper_tj[j], rho, lower_pi[i], upper_pj[j]) +
                    genz_bivariate_normal(upper_ti[i], upper_tj[j], rho, upper_pi[i], upper_pj[j]);

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
      log_likelihood += joint_frequency[i + 1][j + 1] * log_prob;

    }
  }

  // Return negative log-likelihood
  return -log_likelihood;

}

//// Bivariate normal CDF
//// Implements Drezner's approximation (Drezner, 1990)
//double drezner_bivariate_normal(double h1, double h2, double rho, double p1, double p2) {
//
//  // Check for infinities
//  if(h1 == -INFINITY || h2 == -INFINITY){
//    return(0);
//  } else if(h1 == INFINITY){
//    return(p2);
//  } else if(h2 == INFINITY){
//    return(p1);
//  }
//
//  // Set up constants
//  int NX = INT_NX;
//  double X[INT_NX] = DOUBLE_X;
//  double W[INT_NX] = DOUBLE_W;
//  double cor_max = COR_MAX;
//  double bv_fac1 = BV_FAC1;
//  double bv_fac2 = BV_FAC2;
//  double bv = 0;
//
//  // Obtain h12 and absolute correlation
//  double h12 = ((h1 * h1) + (h2 * h2)) / 2;
//  double rho_abs = fabs(rho);
//
//  // Initialize variables on either side of the if statement
//  double h3 = h1 * h2;
//
//  // Check for correlation lower than maximum
//  if(rho_abs <= cor_max) {
//
//    double r1, rr2;
//
//    // Compute probability
//    for (int i = 0; i < NX; i++) {
//      r1 = rho * X[i];
//      rr2 = 1 - r1 * r1;
//      bv += W[i] * exp((r1 * h3 - h12) / rr2) / sqrt(rr2);
//    }
//
//    // Finalize probability
//    bv = p1 * p2 + rho * bv;
//
//  } else {
//
//    double r2 = 1 - rho * rho;
//    double r3 = sqrt(r2);
//
//    // Reverse sign for negative correlation
//    if(rho < 0) {
//      h2 = -h2;
//    }
//
//    // Compute h7
//    double h7 = exp(-h3 / 2.0);
//
//    // Check for correlation less than 1
//    if(rho_abs < 1) {
//      double h6 = fabs(h1 - h2);
//      double h5 = h6 * h6 / 2.0;
//      h6 = h6 / r3;
//      double aa = 0.5 - h3 / 8.0;
//      double ab = 3.0 - 2.0 * aa * h5;
//      bv = bv_fac1 * h6 * ab * (1 - univariate_normal(h6)) - exp(-h5 / r2) * (ab + aa * r2) * bv_fac2;
//
//      double r1, rr;
//
//      // Compute probability
//      for (int i = 0; i < NX; i++) {
//        r1 = r3 * X[i];
//        rr = r1 * r1;
//        r2 = sqrt(1 - rr);
//        bv -= W[i] * exp(-h5 / rr) * (exp(-h3 / (1 + r2)) / r2 / h7 - 1 - aa * rr);
//      }
//
//    }
//
//    // Obtain minimum of probabilities
//    bv = bv * r3 * h7 + fmin(p1, p2);
//
//    // Adjust for negative correlation
//    if (rho < 0) {
//      bv = p1 - bv;
//    }
//
//  }
//
//  // Return probability
//  return bv;
//
//}
//
//// Estimate log-likelihood
//double polychoric_log_likelihood(
//    double rho, int** joint_frequency,
//    double* threshold_X, double* threshold_Y,
//    double* probability_X, double* probability_Y,
//    int max_X, int max_Y
//) {
//
//  // Set upper and lowers (thresholds and probabilities)
//    double lower_ti[max_X], upper_ti[max_X], lower_pi[max_X], upper_pi[max_X];
//    double lower_tj[max_Y], upper_tj[max_Y], lower_pj[max_Y], upper_pj[max_Y];
//
//    // Avoid need for double loop (if possible)
//    if(max_X == max_Y){
//
//        // Set up bounds for X and Y
//        for (int i = 0; i < max_X; i++) {
//            if(i == 0){
//                // X
//                lower_ti[i] = -INFINITY;
//                upper_ti[i] = threshold_X[i + 1];
//                lower_pi[i] = 0;
//                upper_pi[i] = probability_X[i + 1];
//                // Y
//                lower_tj[i] = -INFINITY;
//                upper_tj[i] = threshold_Y[i + 1];
//                lower_pj[i] = 0;
//                upper_pj[i] = probability_Y[i + 1];
//            }else if(i == max_X - 1){
//                // X
//                lower_ti[i] = threshold_X[i];
//                upper_ti[i] = INFINITY;
//                lower_pi[i] = probability_X[i];
//                upper_pi[i] = 1;
//                // Y
//                lower_tj[i] = threshold_Y[i];
//                upper_tj[i] = INFINITY;
//                lower_pj[i] = probability_Y[i];
//                upper_pj[i] = 1;
//            }else{
//                // X
//                lower_ti[i] = threshold_X[i];
//                upper_ti[i] = threshold_X[i + 1];
//                lower_pi[i] = probability_X[i];
//                upper_pi[i] = probability_X[i + 1];
//                // Y
//                lower_tj[i] = threshold_Y[i];
//                upper_tj[i] = threshold_Y[i + 1];
//                lower_pj[i] = probability_Y[i];
//                upper_pj[i] = probability_Y[i + 1];
//            }
//        }
//
//
//    }else{
//
//        // Set up bounds for X
//        for (int i = 0; i < max_X; i++) {
//            if(i == 0){
//                lower_ti[i] = -INFINITY;
//                upper_ti[i] = threshold_X[i + 1];
//                lower_pi[i] = 0;
//                upper_pi[i] = probability_X[i + 1];
//            }else if(i == max_X - 1){
//                lower_ti[i] = threshold_X[i];
//                upper_ti[i] = INFINITY;
//                lower_pi[i] = probability_X[i];
//                upper_pi[i] = 1;
//            }else{
//                lower_ti[i] = threshold_X[i];
//                upper_ti[i] = threshold_X[i + 1];
//                lower_pi[i] = probability_X[i];
//                upper_pi[i] = probability_X[i + 1];
//            }
//        }
//
//        // Set up bounds for Y
//        for (int j = 0; j < max_Y; j++) {
//            if(j == 0){
//                lower_tj[j] = -INFINITY;
//                upper_tj[j] = threshold_Y[j + 1];
//                lower_pj[j] = 0;
//                upper_pj[j] = probability_Y[j + 1];
//            }else if(j == max_Y - 1){
//                lower_tj[j] = threshold_Y[j];
//                upper_tj[j] = INFINITY;
//                lower_pj[j] = probability_Y[j];
//                upper_pj[j] = 1;
//            }else{
//                lower_tj[j] = threshold_Y[j];
//                upper_tj[j] = threshold_Y[j + 1];
//                lower_pj[j] = probability_Y[j];
//                upper_pj[j] = probability_Y[j + 1];
//            }
//        }
//
//    }
//
//  // Initialize variables
//  double log_likelihood = 0.0;
//  double probability, log_prob;
//
//  // Compute log-likelihood
//  for (int i = 0; i < max_X; ++i) {
//    for (int j = 0; j < max_Y; ++j) {
//
//      // Compute bivariate normal CDF
//      probability = drezner_bivariate_normal(upper_ti[i], upper_tj[j], rho, upper_pi[i], upper_pj[j]) -
//                    drezner_bivariate_normal(lower_ti[i], upper_tj[j], rho, lower_pi[i], upper_pj[j]) -
//                    drezner_bivariate_normal(upper_ti[i], lower_tj[j], rho, upper_pi[i], lower_pj[j]) +
//                    drezner_bivariate_normal(lower_ti[i], lower_tj[j], rho, lower_pi[i], lower_pj[j]);
//
//      // Handle probabilities equal to zero
//      if (probability == 0) {
//        probability = DBL_MIN;
//      }
//
//      // Compute log probability
//      log_prob = log(probability);
//
//      // Handle infinite log probabilities
//      if (isinf(log_prob)) {
//        log_prob = DBL_MIN;
//      }
//
//      // Update log-likelihood
//      log_likelihood += joint_frequency[i + 1][j + 1] * log_prob;
//
//    }
//  }
//
//  // Return negative log-likelihood
//  return -log_likelihood;
//
//}

// Brent's method for optimization
double optimize(double (*f)(double, int**, double*, double*, double*, double*, int, int), int** joint_frequency, double* threshold_X, double* threshold_Y, double* probability_X, double* probability_Y, int max_X, int max_Y, double lower, double upper, double tol, int max_iter) {

    // Initialize variables for the optimization algorithm
    double a = lower;
    double b = upper;
    double c = a + (b - a) / 2.0;

    double x = c;
    double w = c;
    double v = c;

    double fx = f(x, joint_frequency, threshold_X, threshold_Y, probability_X, probability_Y, max_X, max_Y);
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
        fu = f(u, joint_frequency, threshold_X, threshold_Y, probability_X, probability_Y, max_X, max_Y);

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
double polychoric(int* X, int max_X, int* Y, int max_Y, int n){

    // Obtain joint frequency table, probability_X, and probability_Y from thresholds function
    struct ThresholdsResult thresholds_result = thresholds(X, max_X, Y, max_Y, n);

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
        max_X, max_Y, lower, upper, tol, max_iter
    );

    // Free memory
    for (int i = 0; i < max_X; ++i) {
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
    int i, j, k;

    // Initialize X and Y
    int X[rows], Y[rows];

    // Initialize maximum categories vector
    int max_categories[cols];

    // Allocate memory for polychoric_matrix
    double** polychoric_matrix = malloc(cols * sizeof(double*));
    for (i = 0; i < cols; i++) {
      polychoric_matrix[i] = malloc(cols * sizeof(double)); // Change rows to cols
      max_categories[i] = find_max_categories(input_data, rows, cols, i); // Saves a loop
    }

    // Perform polychoric correlations over the input_matrix
    for (i = 0; i < cols; i++) {

        // Loop over other variables
        for (j = 0; j < cols; j++) {

            // Fill diagonal
            if (i == j) {

                polychoric_matrix[i][j] = 1;

            }else if (i > j) { // Fill lower triangle

                // Obtain X and Y
                for (k = 0; k < rows; k++) {
                    X[k] = input_data[k + i * rows];
                    Y[k] = input_data[k + j * rows];
                }

                // Compute correlation
                double correlation = polychoric(
                    X, max_categories[i], Y, max_categories[j], rows
                );

                // Fill matrix
                polychoric_matrix[i][j] = correlation;
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
