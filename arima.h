/*
 *   Copyright (c) 2017 Federico Reghenzani
 *
 *   This is a free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This software is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this software; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Federico Reghenzani, Last modified: Jun 20, 2017 */


#ifndef ARIMA_H
#define ARIMA_H

#include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * The struct containing the model computed coefficient
 */
typedef struct arima_model_s {
    double AR_coeff[ARIMA_DIM][ARIMA_AR_POLES * ARIMA_DIM];
    double MA_coeff[ARIMA_DIM][ARIMA_MA_POLES * ARIMA_DIM];
} arima_model_t;

/**
 * The statistical data of the computed model
 */
typedef struct arima_stats_s {
    int nr_iterations;
    double precision;

    double avpm;
    double loglikelihood;
    double aic;
} arima_stats_t;

/**
 * The result ofthe ARIMA estimation. If the field exit_code is not zero, other fields may have
 * invalid data.
 */
typedef struct arima_result_s {

    enum {
        ARIMA_OK = 0,
        ARIMA_SINGULAR_MATRIX,
        ARIMA_GENERIC_ERROR
    } exit_code;

    arima_model_t model;
    arima_stats_t stats;
} arima_result_t;

/**
 * The variable containing the training dataset. You must fill only the first ARIMA_DIM columns
 * and don't touch the second half of the matrix.
 */
extern double series[2*ARIMA_DIM][ARIMA_NR_SAMPLES];

/**
 * Run the model estimation. After the model estimation you cannot assume that the global variable
 * `series` remains unchanged.
 */
extern void estimate_arima(arima_result_t *result);


/**
 * Print the model in a human-readable form
 */
extern void print_model(const arima_model_t *mod);

/**
 * One-step forecast
 * @param model    The estimated ARIMA model
 * @param max_n The number of elements of next two arrays
 * @param measures A max_n * DIM array containing the y_i(t-1), y_i(t-2), ... measures
 * @param predicted A max_n * DIM array containing the ~y_i(t-1), ~y_i(t-2), ... predicted values
 * @param output a DIM array containing the ~y_1(t), ~y_2(t), ... forecasts
 */
void forecast_arima(const arima_model_t *model, double measures[ARIMA_DIM][ARIMA_AR_POLES + ARIMA_I_POLES],
               double errors[ARIMA_DIM][ARIMA_MA_POLES], double output[ARIMA_DIM]) ;

#ifdef __cplusplus
}
#endif

#endif
