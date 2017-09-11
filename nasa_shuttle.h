#ifndef NASA_HEADER_FILE

#define NASA_HEADER_FILE
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_blas.h>
#include <string.h>
#include <math.h>
#define COLUMNS 10
#define CLASSES 7
gsl_matrix* read_data();
void read_new_set();
void calculate_probability_matrices(gsl_matrix*);
long double calculate_probability ( double mean , double stdev, double x, int characteristic );

#endif


