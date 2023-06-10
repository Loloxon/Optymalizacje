#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define IDX(i, j, n) (((j)+ (i)*(n)))

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */

static double gtod_ref_time_sec = 0.0;

double dclock() {
    double the_time, norm_sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    if (gtod_ref_time_sec == 0.0)
        gtod_ref_time_sec = (double) tv.tv_sec;
    norm_sec = (double) tv.tv_sec - gtod_ref_time_sec;
    the_time = norm_sec + tv.tv_usec * 1.0e-6;
    return the_time;
}

int LUPDecompose(double **A, int N, double Tol, int *P) {

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //decomposition done
}

void permute(double **A, int SIZE, int *P) {
    int i, j, k;
    double *permutationMatrix_;
    double *prepermutedMatrix_;
    double **permutationMatrix;
    double **prepermutedMatrix;
    permutationMatrix_ = (double *) malloc(SIZE * SIZE * sizeof(double));
    prepermutedMatrix_ = (double *) malloc(SIZE * SIZE * sizeof(double));
    permutationMatrix = (double **) malloc(SIZE * sizeof(double));
    prepermutedMatrix = (double **) malloc(SIZE * sizeof(double));
    for (i = 0; i < SIZE; i++) {
        permutationMatrix[i] = permutationMatrix_ + i * SIZE;
        prepermutedMatrix[i] = prepermutedMatrix_ + i * SIZE;
    }
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            permutationMatrix[i][j] = 0;
            prepermutedMatrix[i][j] = A[i][j];
        }
    }
    for (i = 0; i < SIZE; i++) {
        permutationMatrix[i][P[i]] = 1;
    }

    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            double sum = 0;
            for (k = 0; k < SIZE; k++) {
                sum += permutationMatrix[i][k] * prepermutedMatrix[k][j];
            }
            A[i][j] = sum;
        }
    }

//    printf("-------permutation---------\n");
//    for (i = 0; i < SIZE; i++) {
//        for (j = 0; j < SIZE; j++) {
//            printf("%f, ", permutationMatrix[i][j]);
//        }
//        printf("\n");
//    }
//    printf("-------permuted---------\n");
//    for (i = 0; i < SIZE; i++) {
//        for (j = 0; j < SIZE; j++) {
//            printf("%f, ", A[i][j]);
//        }
//        printf("\n");
//    }
}

void multiply(double **A, int SIZE) {
    int i, j, k;
    double *matrixL_;
    double *matrixU_;
    double **matrixL;
    double **matrixU;
    matrixL_ = (double *) malloc(SIZE * SIZE * sizeof(double));
    matrixU_ = (double *) malloc(SIZE * SIZE * sizeof(double));
    matrixL = (double **) malloc(SIZE * sizeof(double));
    matrixU = (double **) malloc(SIZE * sizeof(double));
    for (i = 0; i < SIZE; i++) {
        matrixL[i] = matrixL_ + i * SIZE;
        matrixU[i] = matrixU_ + i * SIZE;
    }
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            if (i > j) {
                matrixL[i][j] = A[i][j];
            } else {
                if (i == j) {
                    matrixL[i][j] = 1;
                }
                matrixU[i][j] = A[i][j];
            }
        }
    }

//    printf("-------lower---------\n");
//    for (i = 0; i < SIZE; i++) {
//        for (j = 0; j < SIZE; j++) {
//            printf("%f, ", matrixL[i][j]);
//        }
//        printf("\n");
//    }
//
//    printf("---------upper--------\n");
//    for (i = 0; i < SIZE; i++) {
//        for (j = 0; j < SIZE; j++) {
//            printf("%f, ", matrixU[i][j]);
//        }
//        printf("\n");
//    }

    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            double sum = 0;
            for (k = 0; k < SIZE; k++) {
                sum += matrixL[i][k] * matrixU[k][j];
            }
            A[i][j] = sum;
        }
    }

//    printf("---------multipied--------\n");
//    for (i = 0; i < SIZE; i++) {
//        for (j = 0; j < SIZE; j++) {
//            printf("%f, ", A[i][j]);
//        }
//        printf("\n");
//    }
}

int check(double **oldA, double **A, int SIZE, int *P) {
    int i, j;
    permute(oldA, SIZE, P);
    multiply(A, SIZE);
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            if (fabs(oldA[i][j] - A[i][j]) > 0.001) {
                printf("Error at: %d %d: %f %f", i, j, oldA[i][j], A[i][j]);
                return 1;
            }
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    double t1, t2;
    int i, j, SIZE, ret;
    double dtime;

    SIZE = atoi(argv[1]);
    double *matrix_;
    double *matrixO_;
    double **matrix;
    double **matrixO;
    matrix_ = (double *) malloc(SIZE * SIZE * sizeof(double));
    matrixO_ = (double *) malloc(SIZE * SIZE * sizeof(double));
    matrix = (double **) malloc(SIZE * sizeof(double));
    matrixO = (double **) malloc(SIZE * sizeof(double));
    for (i = 0; i < SIZE; i++) {
        matrix[i] = matrix_ + i * SIZE;
        matrixO[i] = matrixO_ + i * SIZE;
    }
    srand(1);
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            matrix[i][j] = rand() % 10;
            matrixO[i][j] = matrix[i][j];
//            printf("%f, ", matrix[i][j]);
        }
//        printf("\n");
    }


    int *P;
    P = (int *) malloc((SIZE + 1) * sizeof(int));

    dtime = dclock();
    if (!LUPDecompose(matrix, SIZE, 0.01, P)) {
        fprintf(stderr, "Matrix is not symmetric or not positive definite\n");
    }
    dtime = dclock() - dtime;
    printf("Time: %le \n", dtime);

    if (check(matrixO, matrix, SIZE, P)) {
        printf("!!!!! ERROR WITH DECOMPOSITION !!!!!");
    } else {
        printf("Decomposition works fine!");
    }

//    printf("--------orginal--------\n");
//    for (i = 0; i < SIZE; i++) {
//        for (j = 0; j < SIZE; j++) {
//            printf("%f, ", matrixO[i][j]);
//        }
//        printf("\n");
//    }

    fflush(stdout);

    free(matrix_);
    free(matrix);
    free(P);
    return 0;
}