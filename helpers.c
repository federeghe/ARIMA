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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "helpers.h"
#include "global.h"

static double temp_swap[MAX_ARMADIM_POLES];
static double temp_mat[MAX_ARMADIM_POLES][MAX_ARMADIM_POLES];

extern double mat[MAX_ARMADIM_POLES][MAX_ARMADIM_POLES];
extern double inv[MAX_ARMADIM_POLES][MAX_ARMADIM_POLES];
extern double vec[MAX_ARMADIM_POLES];

static int solvele(unsigned int n)
{
    double vswap,max,h,pivot,q;
    int i,j,k,maxi;

    for (i=0;i<n-1;i++) {
        max=fabs(temp_mat[i][i]);
        maxi=i;
        for (j=i+1;j<n;j++)
            if ((h=fabs(temp_mat[j][i])) > max) {
                max=h;
                maxi=j;
            }
        if (maxi != i) {
            for (k=0; k<n; k++) {
                temp_swap[k] = temp_mat[i][k];
                temp_mat[i][k]=temp_mat[maxi][k];
                temp_mat[maxi][k]=temp_swap[k];
            }
            vswap=vec[i];
            vec[i]=vec[maxi];
            vec[maxi]=vswap;
        }

        pivot=temp_mat[i][i];
        if (fabs(pivot) == 0.0) {
            return -1;
        }
        for (j=i+1;j<n;j++) {
            q= -temp_mat[j][i]/pivot;
            temp_mat[j][i]=0.0;
            for (k=i+1;k<n;k++)
                temp_mat[j][k] += q*temp_mat[i][k];
            vec[j] += q*vec[i];
        }
    }

    vec[n-1] /= temp_mat[n-1][n-1];
    for (i=n-2;i>=0;i--) {
        for (j=n-1;j>i;j--)
            vec[i] -= temp_mat[i][j]*vec[j];
        vec[i] /= temp_mat[i][i];
    }

    return 0;
}

int invert_matrix(unsigned int size)
{
    int i,j,k;


    for (i=0;i<size;i++) {
        for (j=0;j<size;j++) {
            vec[j]=(i==j)?1.0:0.0;
            for (k=0;k<size;k++) {
                temp_mat[j][k]=mat[j][k];
            }
        }

        if (solvele(size)) {
            return -1;  // Singular matrix
        }

        for (j=0;j<size;j++) {
            inv[j][i]=vec[j];
        }
    }

    return 0;
}

