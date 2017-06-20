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

#ifndef CONFIG_H
#define CONFIG_H

/**
 * Number of dimensions of the dataset
 */
#define ARIMA_DIM        2

/**
 * Number of samples (for each dimension) to be considered in the arima_series variable. 
 */
#define ARIMA_NR_SAMPLES 2000

/**
 * Number of poles for preliminary AR estimation
 */
#define ARIMA_NR_POLES   1

/**
 * Order for ARIMA model estimation (AR part)
 */
#define ARIMA_AR_POLES   2

/**
 * Degree of differences for ARIMA model estimation (I part)
 */
#define ARIMA_I_POLES    1

/**
 * Order for ARIMA model estimation (MA part)
 */
#define ARIMA_MA_POLES   1

/**
 * Max number of iterations allowed.
 */
#define ARIMA_MAX_ITERS  10000

/**
 * Target precision to reach. The precision may be not reached if ARIMA_MAX_ITERS limit is reached first.
 */
#define ARIMA_PRECISION  1.0e-3


/**
 * Set 1 this option to have debug messages in stderr
 */
#define ARIMA_DEBUG 1

/**
 * Set 1 this option to have insane verbosity in debug messages in stderr
 */
#define ARIMA_DEBUG_EXTRA 1

#endif
