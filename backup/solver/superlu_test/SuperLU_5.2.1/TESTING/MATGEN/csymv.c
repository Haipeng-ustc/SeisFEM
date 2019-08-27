#include <string.h>
#include "f2c.h"

/* Subroutine */ int csymv_sluslu(char *uplo, integer *n, complex *alpha, complex *
	a, integer *lda, complex *x, integer *incx, complex *beta, complex *y,
	 integer *incy)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CSYMV  performs the matrix-vector  operation   

       y := alpha*A*x + beta*y,   

    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix.   

    Arguments   
    ==========   

    UPLO   - CHARACTER*1   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX array, dimension ( LDA, N )   
             Before entry, with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry, with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced.   
             Unchanged on exit.   

    LDA    - INTEGER   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, N ).   
             Unchanged on exit.   

    X      - COMPLEX array, dimension at least   
             ( 1 + ( N - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the N-   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - COMPLEX   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - COMPLEX array, dimension at least   
             ( 1 + ( N - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

   ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3, q__4;
    /* Local variables */
    static integer info;
    static complex temp1, temp2;
    static integer i, j;
    static integer ix, iy, jx, jy, kx, ky;
    extern int input_error(char *, int *);

#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if ( strncmp(uplo, "U", 1)!=0 && ! strncmp(uplo, "L", 1)!=0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	input_error("CSYMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0.f && alpha->i == 0.f && (beta->r == 1.f && 
	    beta->i == 0.f)) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   

       First form  y := beta*y. */

    if (beta->r != 1.f || beta->i != 0.f) {
	if (*incy == 1) {
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = i;
		    Y(i).r = 0.f, Y(i).i = 0.f;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = i;
		    i__3 = i;
		    q__1.r = beta->r * Y(i).r - beta->i * Y(i).i, 
			    q__1.i = beta->r * Y(i).i + beta->i * Y(i)
			    .r;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = iy;
		    Y(iy).r = 0.f, Y(iy).i = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = iy;
		    i__3 = iy;
		    q__1.r = beta->r * Y(iy).r - beta->i * Y(iy).i, 
			    q__1.i = beta->r * Y(iy).i + beta->i * Y(iy)
			    .r;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }
    if (strncmp(uplo, "U", 1)==0) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		q__1.r = alpha->r * X(j).r - alpha->i * X(j).i, q__1.i =
			 alpha->r * X(j).i + alpha->i * X(j).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(i).r + q__2.r, q__1.i = Y(i).i + q__2.i;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
		    i__3 = i + j * a_dim1;
		    i__4 = i;
		    q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(i).i, 
			    q__2.i = A(i,j).r * X(i).i + A(i,j).i * X(
			    i).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
/* L50: */
		}
		i__2 = j;
		i__3 = j;
		i__4 = j + j * a_dim1;
		q__3.r = temp1.r * A(j,j).r - temp1.i * A(j,j).i, q__3.i = 
			temp1.r * A(j,j).i + temp1.i * A(j,j).r;
		q__2.r = Y(j).r + q__3.r, q__2.i = Y(j).i + q__3.i;
		q__4.r = alpha->r * temp2.r - alpha->i * temp2.i, q__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		Y(j).r = q__1.r, Y(j).i = q__1.i;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, q__1.i =
			 alpha->r * X(jx).i + alpha->i * X(jx).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__3 = iy;
		    i__4 = iy;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(iy).r + q__2.r, q__1.i = Y(iy).i + q__2.i;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    i__3 = i + j * a_dim1;
		    i__4 = ix;
		    q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(ix).i, 
			    q__2.i = A(i,j).r * X(ix).i + A(i,j).i * X(
			    ix).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		i__2 = jy;
		i__3 = jy;
		i__4 = j + j * a_dim1;
		q__3.r = temp1.r * A(j,j).r - temp1.i * A(j,j).i, q__3.i = 
			temp1.r * A(j,j).i + temp1.i * A(j,j).r;
		q__2.r = Y(jy).r + q__3.r, q__2.i = Y(jy).i + q__3.i;
		q__4.r = alpha->r * temp2.r - alpha->i * temp2.i, q__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		q__1.r = alpha->r * X(j).r - alpha->i * X(j).i, q__1.i =
			 alpha->r * X(j).i + alpha->i * X(j).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		i__2 = j;
		i__3 = j;
		i__4 = j + j * a_dim1;
		q__2.r = temp1.r * A(j,j).r - temp1.i * A(j,j).i, q__2.i = 
			temp1.r * A(j,j).i + temp1.i * A(j,j).r;
		q__1.r = Y(j).r + q__2.r, q__1.i = Y(j).i + q__2.i;
		Y(j).r = q__1.r, Y(j).i = q__1.i;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(i).r + q__2.r, q__1.i = Y(i).i + q__2.i;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
		    i__3 = i + j * a_dim1;
		    i__4 = i;
		    q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(i).i, 
			    q__2.i = A(i,j).r * X(i).i + A(i,j).i * X(
			    i).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
/* L90: */
		}
		i__2 = j;
		i__3 = j;
		q__2.r = alpha->r * temp2.r - alpha->i * temp2.i, q__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = Y(j).r + q__2.r, q__1.i = Y(j).i + q__2.i;
		Y(j).r = q__1.r, Y(j).i = q__1.i;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, q__1.i =
			 alpha->r * X(jx).i + alpha->i * X(jx).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		i__2 = jy;
		i__3 = jy;
		i__4 = j + j * a_dim1;
		q__2.r = temp1.r * A(j,j).r - temp1.i * A(j,j).i, q__2.i = 
			temp1.r * A(j,j).i + temp1.i * A(j,j).r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    ix += *incx;
		    iy += *incy;
		    i__3 = iy;
		    i__4 = iy;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(iy).r + q__2.r, q__1.i = Y(iy).i + q__2.i;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    i__3 = i + j * a_dim1;
		    i__4 = ix;
		    q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(ix).i, 
			    q__2.i = A(i,j).r * X(ix).i + A(i,j).i * X(
			    ix).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
/* L110: */
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp2.r - alpha->i * temp2.i, q__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of CSYMV */

} /* csymv_sluslu */
