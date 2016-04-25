#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vec3d.h>
#include "lapacke.h"

extern void printMatrix( char const* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char const* desc, int n, int* a );

#define N 5
#define NRHS 3
#define LDA N
#define LDB NRHS


/* Main program */
int main() {
        /* Locals */
        int n = 5, nrhs = 1, lda = n, ldb = nrhs, info;
        /* Local arrays */
        int ipiv[N];
        double a[LDA*N] = {
            6.80, -6.05, -0.45,  8.32, -9.67,
           -2.11, -3.30,  2.58,  2.71, -5.14,
            5.66, 5.36, -2.70,  4.35, -7.26,
            5.97, -4.44,  0.27, -7.17, 6.08,
            8.23, 1.08,  9.04,  2.14, -6.87
        };
        double* anew = (double*) malloc(LDA*N*sizeof(double));
        for (int i = 0; i <LDA*N; i++){
            anew[i] = a[i];
        }
            
        double b[LDB*N] = {
            4.02,
            6.19,
           -8.22,
           -7.57,
           -3.03
        };
        double* bnew = (double*) malloc(LDB*N*sizeof(double));
        for (int i = 0; i <LDB*N; i++){
            bnew[i] = b[i];
        }
        
        /* Executable statements */
        printf( "LAPACKE_dgesv (row-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */
        info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv,
                        b, ldb );
        /* Check for the exact singularity */
        if( info > 0 ) {
                printf( "The diagonal element of the triangular factor of A,\n" );
                printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
                printf( "the solution could not be computed.\n" );
                exit( 1 );
        }
        /* Print solution */
        printMatrix( "Solution", n, nrhs, b, ldb );
        /* Print details of LU factorization */
        printMatrix( "Details of LU factorization", n, n, a, lda );
        /* Print pivot indices */
        print_int_vector( "Pivot indices", n, ipiv );
        //free(anew); free(bnew);
        exit( 0 );
} /* End of LAPACKE_dgesv Example */

