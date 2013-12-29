/*
 * Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Mon Mar 24 02:06:57 EST 2003 */

#include "fftw-int.h"
#include "fftw.h"

/* Generated by: /homee/stevenj/cvs/fftw/gensrc/genfft -magic-alignment-check -magic-twiddle-load-all -magic-variables 4 -magic-loopi -hc2real 7 */

/*
 * This function contains 24 FP additions, 19 FP multiplications,
 * (or, 23 additions, 18 multiplications, 1 fused multiply/add),
 * 13 stack variables, and 14 memory accesses
 */
static const fftw_real K2_000000000 =
FFTW_KONST(+2.000000000000000000000000000000000000000000000);
static const fftw_real K1_801937735 =
FFTW_KONST(+1.801937735804838252472204639014890102331838324);
static const fftw_real K445041867 =
FFTW_KONST(+0.445041867912628808577805128993589518932711138);
static const fftw_real K1_246979603 =
FFTW_KONST(+1.246979603717467061050009768008479621264549462);
static const fftw_real K867767478 =
FFTW_KONST(+0.867767478235116240951536665696717509219981456);
static const fftw_real K1_949855824 =
FFTW_KONST(+1.949855824363647214036263365987862434465571601);
static const fftw_real K1_563662964 =
FFTW_KONST(+1.563662964936059617416889053348115500464669037);

/*
 * Generator Id's : 
 * $Id: exprdag.ml,v 1.43 2003/03/16 23:43:46 stevenj Exp $
 * $Id: fft.ml,v 1.44 2003/03/16 23:43:46 stevenj Exp $
 * $Id: to_c.ml,v 1.26 2003/03/16 23:43:46 stevenj Exp $
 */

void fftw_hc2real_7(const fftw_real *real_input,
		    const fftw_real *imag_input, fftw_real *output,
		    int real_istride, int imag_istride, int ostride)
{
     fftw_real tmp9;
     fftw_real tmp13;
     fftw_real tmp11;
     fftw_real tmp1;
     fftw_real tmp4;
     fftw_real tmp2;
     fftw_real tmp3;
     fftw_real tmp5;
     fftw_real tmp12;
     fftw_real tmp10;
     fftw_real tmp6;
     fftw_real tmp8;
     fftw_real tmp7;
     ASSERT_ALIGNED_DOUBLE;
     tmp6 = imag_input[2 * imag_istride];
     tmp8 = imag_input[imag_istride];
     tmp7 = imag_input[3 * imag_istride];
     tmp9 =
	 (K1_563662964 * tmp6) - (K1_949855824 * tmp7) -
	 (K867767478 * tmp8);
     tmp13 =
	 (K867767478 * tmp6) + (K1_563662964 * tmp7) -
	 (K1_949855824 * tmp8);
     tmp11 =
	 (K1_563662964 * tmp8) + (K1_949855824 * tmp6) +
	 (K867767478 * tmp7);
     tmp1 = real_input[0];
     tmp4 = real_input[3 * real_istride];
     tmp2 = real_input[real_istride];
     tmp3 = real_input[2 * real_istride];
     tmp5 =
	 tmp1 + (K1_246979603 * tmp3) - (K445041867 * tmp4) -
	 (K1_801937735 * tmp2);
     tmp12 =
	 tmp1 + (K1_246979603 * tmp4) - (K1_801937735 * tmp3) -
	 (K445041867 * tmp2);
     tmp10 =
	 tmp1 + (K1_246979603 * tmp2) - (K1_801937735 * tmp4) -
	 (K445041867 * tmp3);
     output[4 * ostride] = tmp5 - tmp9;
     output[3 * ostride] = tmp5 + tmp9;
     output[0] = tmp1 + (K2_000000000 * (tmp2 + tmp3 + tmp4));
     output[2 * ostride] = tmp12 + tmp13;
     output[5 * ostride] = tmp12 - tmp13;
     output[6 * ostride] = tmp10 + tmp11;
     output[ostride] = tmp10 - tmp11;
}

fftw_codelet_desc fftw_hc2real_7_desc = {
     "fftw_hc2real_7",
     (void (*)()) fftw_hc2real_7,
     7,
     FFTW_BACKWARD,
     FFTW_HC2REAL,
     169,
     0,
     (const int *) 0,
};
