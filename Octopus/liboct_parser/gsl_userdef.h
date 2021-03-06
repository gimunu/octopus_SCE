/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id: gsl_userdef.h 10978 2013-07-11 15:28:46Z micael $
*/

#ifndef __GSL_USERDEF_H__
#define __GSL_USERDEF_H__

/* Heaviside step function */
gsl_complex gsl_complex_step_real (double a);

/* Minimum and maximum of two arguments (comparing real parts) */  
gsl_complex gsl_complex_min_real (gsl_complex a, gsl_complex b);
gsl_complex gsl_complex_max_real (gsl_complex a, gsl_complex b);

/* arg, abs, abs2 and logabs with complex return values for consistency */
gsl_complex gsl_complex_carg  (gsl_complex a);
gsl_complex gsl_complex_cabs  (gsl_complex a);
gsl_complex gsl_complex_cabs2 (gsl_complex a);
gsl_complex gsl_complex_clogabs (gsl_complex a);

/* error function */
gsl_complex gsl_complex_erf(gsl_complex a);

/* atan2 function */
gsl_complex gsl_complex_arctan2 (gsl_complex a, gsl_complex b);

#endif /* __GSL_USERDEF_H__ */
