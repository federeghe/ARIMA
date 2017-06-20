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

#ifndef GLOBAL_H
#define GLOBAL_H

#include "config.h"

/** Definition of some useful values */

#define ARDIM      ARIMA_NR_POLES*ARIMA_DIM
#define ARMADIM    (ARIMA_AR_POLES+ARIMA_MA_POLES)*ARIMA_DIM
#define MAX_ARMADIM_POLES (ARDIM > ARMADIM ? ARDIM : ARMADIM)

#if (ARIMA_AR_POLES+ARIMA_I_POLES+ARIMA_MA_POLES) > 0
    #define ARIMASET 1
#else
    #define ARIMASET 0
#endif

/** Some checks **/

#if ARIMA_AR_POLES > ARIMA_NR_SAMPLES || ARIMA_MA_POLES > ARIMA_NR_SAMPLES
    #error "Non-sense value for AR/MA poles"
#endif
#if ARIMA_NR_POLES > ARIMA_NR_SAMPLES 
    #error "Non-sense value for NR poles"
#endif

#if ARIMA_DEBUG_EXTRA
    #if ! ARIMA_DEBUG
        #error "You have to define ARIMA_DEBUG if you declare ARIMA_DEBUG_EXTRA"
    #endif
#endif

#endif
