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


#include "arima.h"
#include "global.h"

#include <chrono>
#include <fstream>
#include <stdio.h>
#define TOT_SAMPLES 89783

#define DATA_OFFSET 0

double serie_temp1[TOT_SAMPLES];
double serie_temp2[TOT_SAMPLES];

int main() {
    int i=0;
    arima_result_t res;

    printf("Starting...\n");
    printf("Loading files (1)...\n");

    std::ifstream temp1("temp1.txt");
    while (temp1.good()) {
        temp1 >> serie_temp1[i++];
    }
    temp1.close();

    i = 0;

    printf("Loading files (2)...\n");

    std::ifstream temp2("temp2.txt");
    while (temp2.good()) {
        temp2 >> serie_temp2[i++];
    }
    temp2.close();


    printf("Copying dataset for training...\n");
    for (i=0; i < ARIMA_NR_SAMPLES; i++) {
        series[0][i] = serie_temp1[i+DATA_OFFSET];
        series[1][i] = serie_temp2[i+DATA_OFFSET];
    }


    printf("Estimating ARIMA...\n");
    auto t1 = std::chrono::high_resolution_clock::now();
    estimate_arima(&res);
    auto t2 = std::chrono::high_resolution_clock::now();

    if (res.exit_code) {
        printf("Impossible to determine a model.\n");
    } else {
        print_model(&res.model);
    }
    printf("Estimation Timing: %ld\n", 
        std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());

    double output[ARIMA_DIM];
    double measures[ARIMA_DIM][ARIMA_AR_POLES + ARIMA_I_POLES];
    double errors  [ARIMA_DIM][ARIMA_MA_POLES];

    // To get the first prediction, we cannot know any error term, so just set them to zero
    for ( int i=0; i<ARIMA_DIM; i++) {
        for ( int j=0; j < ARIMA_MA_POLES ; j++ ) {
            errors[i][j] = 0;
        }
    }

    // Then, we set the x(t-1), x(t-2), ... terms
    for ( int j=0; j < ARIMA_AR_POLES + ARIMA_I_POLES; j++ ) {
        measures[0][j] = serie_temp1[ARIMA_NR_SAMPLES+j-1];
        measures[1][j] = serie_temp2[ARIMA_NR_SAMPLES+j-1];
    }

    t1 = std::chrono::high_resolution_clock::now();
    forecast_arima(&res.model, measures, errors, output);
    t2 = std::chrono::high_resolution_clock::now();

    printf("Running Timing: %ld\n", 
        std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());

    for ( int i=0; i<ARIMA_DIM; i++) {
        printf("%f ", output[i]);
    }
    printf("\n");

    return 0;
}


