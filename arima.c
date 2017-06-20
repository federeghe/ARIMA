/*
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Rainer Hegger, Last modified: Feb 6, 2006 */
/*Author: Federico Reghenzani, Last modified: Jun 20, 2017 */

#include "arima.h"
#include "helpers.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#define sqr(x) ((x)*(x))

#ifndef M_PI
    #define M_PI 3.141592653589793238L
#endif

static unsigned int length, poles, offset, series_offset;
static unsigned int aindex[2][MAX_ARMADIM_POLES];

double mat[MAX_ARMADIM_POLES][MAX_ARMADIM_POLES];
double inv[MAX_ARMADIM_POLES][MAX_ARMADIM_POLES];
double vec[MAX_ARMADIM_POLES];
static double coeff[ARIMA_DIM][MAX_ARMADIM_POLES];
static double oldcoeff[ARIMA_DIM][ARMADIM];
static double diffcoeff[ARIMA_MAX_ITERS];
static double diff[ARIMA_DIM][ARIMA_NR_SAMPLES];
static double resi[ARIMA_DIM];
static double xdiff[ARIMA_MAX_ITERS][ARIMA_DIM];

double series[2*ARIMA_DIM][ARIMA_NR_SAMPLES];

static double my_average[ARIMA_DIM];


static void make_difference(void)
{
    unsigned long i,d;

    for (i=length-1;i>0;i--)
        for (d=0;d<ARIMA_DIM;d++)
            series[d][i]=series[d][i]-series[d][i-1];
}

static void make_ar_index(void)
{
    unsigned long i;

    assert(ARDIM <= MAX_ARMADIM_POLES); 


    for (i=0;i < ARDIM;i++) {
        aindex[0][i]=i/poles;
        aindex[1][i]=i%poles;
    }

}

#if ARIMASET
static void make_arima_index(unsigned int ars,unsigned int mas)
{
    unsigned int armad;
    unsigned long i,i0;

    armad=(ars+mas)*ARIMA_DIM;
    assert(armad <= MAX_ARMADIM_POLES); 


    for (i=0;i<ars*ARIMA_DIM;i++) {
        aindex[0][i]=i/ars;
        aindex[1][i]=i%ars;
    }
    i0=ars*ARIMA_DIM;
    for (i=0;i<mas*ARIMA_DIM;i++) {
        aindex[0][i+i0]=ARIMA_DIM+i/mas;
        aindex[1][i+i0]=i%mas;
    }

}
#endif

static void set_averages_to_zero(void)
{
    long i,j,k;

    for (i=0;i<ARIMA_DIM;i++) {

        my_average[i] = 0;
        for (k=0; k < length; k++) {
            my_average[i] += series[i][k+series_offset];
        }
        my_average[i] /= (double)length;

        for (j=0;j<length;j++)
            series[i][j+series_offset] -= my_average[i];
    }
}

static int build_matrix(unsigned int size)
{
    long n,i,j,is,id,js,jd;
    double norm;

    norm=1./((double)length-1.0-(double)poles-(double)offset);

#if ARIMA_DEBUG_EXTRA
    fprintf(stderr,"Normalizing matrix with norm %f\n", norm);
#endif

    for (i=0;i<size;i++) {
        id=aindex[0][i];
        is=aindex[1][i];
        for (j=i;j<size;j++) {
            jd=aindex[0][j];
            js=aindex[1][j];
            mat[i][j]=0.0;
            for (n=offset+poles-1;n<length-1;n++)
                mat[i][j] += series[id][n-is+series_offset]*series[jd][n-js+series_offset];
            mat[i][j] *= norm;
            mat[j][i]=mat[i][j];
        }
    }

#if ARIMA_DEBUG_EXTRA
    fprintf(stderr,"Trying to invert: \n");
    for (i=0;i<size;i++) {
        for (j=0;j<size;j++) {
            fprintf(stderr,"%10f", mat[i][j]);
        }
        fprintf(stderr,"\n");
    }
#endif

  return invert_matrix(size);
}

static void build_vector(unsigned int size,long comp)
{
    long i,is,id,n;
    double norm;

    norm=1./((double)length-1.0-(double)poles-(double)offset);

    for (i=0;i<size;i++) {
        id=aindex[0][i];
        is=aindex[1][i];
        vec[i]=0.0;
        for (n=offset+poles-1;n<length-1;n++)
            vec[i] += series[comp][n+1+series_offset]*series[id][n-is+series_offset];
        vec[i] *= norm;
    }
}

static void multiply_matrix_vector(unsigned int size, int k)
{
  long i,j;

  for (i=0;i<size;i++) {
    coeff[k][i]=0.0;
    for (j=0;j<size;j++)
      coeff[k][i] += inv[i][j]*vec[j];
  }

}

static void make_residuals(unsigned int size)
{
    long n,n1,d,i,is,id;

    for (i=0;i<ARIMA_DIM;i++)
        resi[i]=0.0;

    for (n=poles-1;n<length-1;n++) {
        n1=n+1;
        for (d=0;d<ARIMA_DIM;d++) {
            diff[d][n1]=series[d][n1+series_offset];
            for (i=0;i<size;i++) {
                id=aindex[0][i];
                is=aindex[1][i];
                diff[d][n1] -= coeff[d][i]*series[id][n-is+series_offset];
          }
          resi[d] += sqr(diff[d][n1]);
        }
    }

    for (i=0;i<ARIMA_DIM;i++)
        resi[i]=sqrt(resi[i]/((double)length-(double)poles));

}

#if ARIMASET
static int learn_arima(void) {
    int i, j, iter, hj, last_err, realiter=0;
    double hdiff, alldiff;

    offset=poles;

    make_arima_index(ARIMA_AR_POLES,ARIMA_MA_POLES);

    for (i=0;i<ARIMA_DIM;i++) {
        for (j=0;j<ARMADIM;j++)
            oldcoeff[i][j]=0.0;
    }

    for (iter=1;iter<=ARIMA_MAX_ITERS;iter++) {

        poles=(ARIMA_AR_POLES > ARIMA_MA_POLES)? ARIMA_AR_POLES : ARIMA_MA_POLES;

        offset += poles;
    last_err = build_matrix(ARMADIM);

        for (i=0;i<ARIMA_DIM;i++) {
            build_vector(ARMADIM,i);
            multiply_matrix_vector(ARMADIM,i);
        }

        make_residuals(ARMADIM);

        for (j=0;j<ARIMA_DIM;j++) {
            hdiff=0.0;
            hj=j+ARIMA_DIM;
            for (i=offset;i<length;i++)
                hdiff += sqr(series[hj][i + series_offset]-diff[j][i]);
            for (i=0;i<length;i++) {
                series[hj][i+ series_offset]=diff[j][i];
            }
            xdiff[iter-1][j]=sqrt(hdiff/(double)(length-offset));
        }


        diffcoeff[iter-1]=0.0;
        for (i=0;i<ARIMA_DIM;i++) {
            for (j=0;j<ARIMA_DIM;j++) {
                diffcoeff[iter-1] += sqr(coeff[i][j]-oldcoeff[i][j]);
                oldcoeff[i][j]=coeff[i][j];
            }
        }

        diffcoeff[iter-1]=sqrt(diffcoeff[iter-1]/(double)ARMADIM);
        alldiff=xdiff[iter-1][0];
        for (i=1;i<ARIMA_DIM;i++)
            if (xdiff[iter-1][i] > alldiff)
                alldiff=xdiff[iter-1][i];
        realiter=iter;
        if (alldiff < ARIMA_PRECISION)
            iter=ARIMA_MAX_ITERS;
    }

    if (last_err < 0) {
        return -1;
    }

    return realiter;
}

#endif

static void calculate_errors(arima_result_t *res) {
    int i;

    res->stats.avpm=resi[0]*resi[0];
    res->stats.loglikelihood= -log(resi[0]);
    for (i=1;i<ARIMA_DIM;i++) {
        res->stats.avpm += resi[i]*resi[i];
        res->stats.loglikelihood -= log(resi[i]);
    }
    res->stats.loglikelihood *= ((double)length);
    res->stats.loglikelihood += -((double)length)*((1.0+log(2.*M_PI))*ARIMA_DIM)/2.0;
    res->stats.avpm=sqrt(res->stats.avpm/ARIMA_DIM);

#if ARIMASET
        res->stats.aic=2.0*(ARIMA_AR_POLES+ARIMA_MA_POLES)-2.0*res->stats.loglikelihood;
#else
        res->stats.aic=2.0*poles-2.0*res->stats.loglikelihood;
#endif

}

#if ARIMA_DEBUG
static void print_debug_iterations(int realiter) {
    int i, j;

    fprintf(stderr,"N. Iters = %d\r\n", realiter);

    #if ARIMA_DEBUG_EXTRA

    for (i=0;i<realiter;i++) {
        fprintf(stderr,"#iteration %d ",i+1);
        for (j=0;j<ARIMA_DIM;j++)
            fprintf(stderr,"%e ",xdiff[i][j]);
        fprintf(stderr,"%e",diffcoeff[i]);
        fprintf(stderr,"\r\n");
    }
    #else
    (void) i;
    (void) j;
    #endif

}


static void print_debug_error_stat(const arima_stats_t *err_stat) {
    int i; 
    fprintf(stderr,"#average forcast error= %e\r\n",err_stat->avpm);
    fprintf(stderr,"#individual forecast errors: ");
    for (i=0;i<ARIMA_DIM;i++)
        fprintf(stderr,"%e ",resi[i]);
    
    fprintf(stderr,"\r\n");

    fprintf(stderr,"#Log-Likelihood= %e\t AIC= %e\r\n",err_stat->loglikelihood,err_stat->aic);
}

#endif

static void set_solutions(int size, arima_result_t* result) {
    unsigned int id, is, i, j;
    for (i=0;i<size;i++) {
        id=aindex[0][i];
        is=aindex[1][i];
        if (id < ARIMA_DIM) {
#if ARIMA_DEBUG
            fprintf(stderr,"#x_%u(n-%u) ",id+1,is);
#else
	(void)is;
#endif // ARIMA_DEBUG
            for (j=0;j<ARIMA_DIM;j++)
                result->model.AR_coeff[j][i] = coeff[j][i];
    }
        else {
#if ARIMA_DEBUG
            fprintf(stderr,"#e_%u(n-%u) ",id+1-ARIMA_DIM,is);
#endif // ARIMA_DEBUG
            for (j=0;j<ARIMA_DIM;j++)
                result->model.MA_coeff[j][i - ARIMA_DIM * ARIMA_AR_POLES] = coeff[j][i];
    }
#if ARIMA_DEBUG
        for (j=0;j<ARIMA_DIM;j++)
            fprintf(stderr,"%e ",coeff[j][i]);
        fprintf(stderr,"\r\n");
#endif // ARIMA_DEBUG
    }

#if ARIMA_DEBUG_EXTRA
    for (i=poles;i<length;i++) {
        for (j=0;j<ARIMA_DIM;j++)
            fprintf(stderr,"%e %e ",series[j][i+ series_offset]+my_average[j],diff[j][i]);

        fprintf(stderr,"\r\n");
    }
#endif /* ARIMA_DEBUG_EXTRA */
} 


void estimate_arima(arima_result_t *result)
{
    int i, size, realiter;

    poles   = ARIMA_NR_POLES;
    length  = ARIMA_NR_SAMPLES;

    offset = 0;

    for (i=0;i<ARIMA_I_POLES;i++)
        make_difference();

    series_offset = ARIMA_I_POLES;
    length -= ARIMA_I_POLES;

    set_averages_to_zero();

    make_ar_index();

    if (build_matrix(ARDIM)) {
        result->exit_code = ARIMA_SINGULAR_MATRIX;
    }


    for (i=0;i<ARIMA_DIM;i++) {
        build_vector(ARDIM,i);
        multiply_matrix_vector(ARDIM, i);
    }


    make_residuals(ARDIM);


    size=ARDIM;
  
#if ARIMASET
        size=ARMADIM;
        realiter = learn_arima();
    if (realiter < 0) {
            result->exit_code = ARIMA_SINGULAR_MATRIX;
            return;
    }
#endif

    // Calculate error statistics
    calculate_errors(result);

    set_solutions(size, result);

#if ARIMA_DEBUG
    #if ARIMASET
        print_debug_iterations(realiter);
    #endif

    print_debug_error_stat(&result->stats);
#else
    (void) realiter;
    (void) size;
#endif

    result->exit_code = ARIMA_OK;
    return;
}


extern void print_model(const arima_model_t *model) {

    for (int i=0; i < ARIMA_DIM; i++) {
        printf("x_%d(t+1) = ", i+1);
#if ARIMA_AR_POLES > 0

        for (int j=0; j < ARIMA_AR_POLES * ARIMA_DIM; j++) {

            if (model->AR_coeff[i][j] > 0) {
                printf("+");
            }
            printf("%f x_%d(t", model->AR_coeff[i][j], j/(ARIMA_AR_POLES)+1);
            if ( j % ARIMA_AR_POLES != 0) {
                printf("-%d", j % ARIMA_AR_POLES);
            }
            printf(") ");
        }
#endif
#if ARIMA_MA_POLES > 0
        for (int j=0; j < ARIMA_MA_POLES * ARIMA_DIM; j++) {

            if (model->MA_coeff[i][j] > 0) {
                printf("+");
            }
            printf("%f e_%d(t", model->MA_coeff[i][j], j/(ARIMA_MA_POLES)+1);
            if ( j % ARIMA_MA_POLES != 0) {
                printf(" - %d", j % ARIMA_MA_POLES);
            }
            printf(") ");
        }
#endif

        printf("\n");
    }

}


extern void forecast_arima(const arima_model_t *model, double measures[ARIMA_DIM][ARIMA_AR_POLES + ARIMA_I_POLES],
               double errors[ARIMA_DIM][ARIMA_MA_POLES], double output[ARIMA_DIM]) {
    int i,j,d;
    int offset = 0;

#if ARIMA_I_POLES > 0
    double I_measures[ARIMA_DIM][ARIMA_AR_POLES + ARIMA_I_POLES];

    for (i=0;i<ARIMA_DIM;i++) {
        for (j=0; j<ARIMA_AR_POLES + ARIMA_I_POLES; j++) {
                I_measures[i][j]=measures[i][j];
        }
    }

    for (i=0; i < ARIMA_I_POLES; i++) {
        for (j=ARIMA_AR_POLES + ARIMA_I_POLES-1;j>0;j--) {
                for (d=0;d<ARIMA_DIM;d++) {
                    I_measures[d][j]=I_measures[d][j]-I_measures[d][j-1];
            }
        }
    }
    offset = ARIMA_I_POLES;
#endif

    for (i=0; i < ARIMA_DIM; i++) {

#if ARIMA_I_POLES > 0
        output[i] = measures[i][offset + ARIMA_AR_POLES - 1];
#else
        output[i] = 0;
#endif
        for (j=0; j < ARIMA_AR_POLES * ARIMA_DIM; j++) {
#if ARIMA_I_POLES > 0
            output[i] += model->AR_coeff[i][j] * I_measures[j / ARIMA_AR_POLES][offset + j % ARIMA_AR_POLES];
#else
            output[i] += model->AR_coeff[i][j] * measures[j / ARIMA_AR_POLES][j % ARIMA_AR_POLES];
#endif
        }

        for (j=0; j < ARIMA_MA_POLES * ARIMA_DIM; j++) {
            output[i] += model->MA_coeff[i][j] * errors[j / ARIMA_MA_POLES][j % ARIMA_MA_POLES];
        }
    }

}


