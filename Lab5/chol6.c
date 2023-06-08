#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <immintrin.h> // OPT3
/* Include PAPI */
#include <papi.h>

/* PAPI macro helpers definitions */
#define NUM_EVENT 2
#define THRESHOLD 100000
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }


#define IDX(i, j, n) (j+ i*n)
//#define IDX(i,j,n) (i*(i+1)/2+j)
#define BLKSIZE 16


static double gtod_ref_time_sec = 0.0;

/* Adapted from the bl2_clock() routine in the BLIS library */

double dclock(){
    double the_time, norm_sec;
    struct timeval tv;
    gettimeofday( &tv, NULL );
    if ( gtod_ref_time_sec == 0.0 )
        gtod_ref_time_sec = ( double ) tv.tv_sec;
    norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;
    the_time = norm_sec + tv.tv_usec * 1.0e-6;
    return the_time;
}

int max(int a, int b){
    if (a > b)
        return (a);
    return (b);
}

int chol(double *A, unsigned int n){
    register unsigned int i,j,k;
    register double tmp;
    register __m256d tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7; // OPT 3
    //register __m128d tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15;

    for (j = 0; j < n; j++) {
        for (i = j; i < n; i++) {
            tmp = A[IDX(i, j, n)];
            for (k = 0; k < j; ) {
                if (k < max(j - BLKSIZE, 0)) {
                    tmp0 = _mm256_loadu_pd(A+IDX(i, k, n)); // <- OPT 4
                    tmp1 = _mm256_loadu_pd(A+IDX(j, k, n));
                    tmp2 = _mm256_loadu_pd(A+IDX(i, k+4, n));
                    tmp3 = _mm256_loadu_pd(A+IDX(j, k+4, n));
                    tmp4 = _mm256_loadu_pd(A+IDX(i, k+8, n));
                    tmp5 = _mm256_loadu_pd(A+IDX(j, k+8, n));
                    tmp6 = _mm256_loadu_pd(A+IDX(i, k+12, n));
                    tmp7 = _mm256_loadu_pd(A+IDX(j, k+12, n));

                    tmp0 = _mm256_mul_pd(tmp0, tmp1); // <- OPT 5
                    tmp2 = _mm256_mul_pd(tmp2, tmp3);
                    tmp4 = _mm256_mul_pd(tmp4, tmp5);
                    tmp6 = _mm256_mul_pd(tmp6, tmp7);

                    tmp0 = _mm256_add_pd(tmp0, tmp2); // <- OPT 6
                    tmp4 = _mm256_add_pd(tmp4, tmp6);

                    tmp0 = _mm256_add_pd(tmp0, tmp4);

                    tmp -= tmp0[0] + tmp0[1] + tmp0[2] + tmp0[3];
                    k += BLKSIZE;
                } else {
                    tmp -= A[IDX(i, k, n)] * A[IDX(j, k, n)];
                    k++;
                }
            }
            A[IDX(i, j, n)] = tmp;
        }

        if (A[IDX(j, j, n)] < 0.0) {
            return (1);
        }

        tmp = sqrt(A[IDX(j, j, n)]);
        for (i = j + 1; i < n; i++){
            A[IDX(i, j, n)] /= tmp;
        }
    }

    return (0);
}

int main(int argc, char ** argv){
    double* A;
    double t1, t2;
    int i, j, n, ret;
    double dtime;

    int measure = 0; // default value of what to measure

    /* PAPI FLOPS variables */
    float real_time, proc_time,mflops;
    long long flpops;
    float ireal_time, iproc_time, imflops;
    long long iflpops;
    int retval;

    /* PAPI counters variables */
    int tmp;
    int EventSet = PAPI_NULL;
    /*must be initialized to PAPI_NULL before calling PAPI_create_event*/

    int event_codes[NUM_EVENT]={PAPI_LST_INS,PAPI_REF_CYC};
    char errstring[PAPI_MAX_STR_LEN];
    long long values[NUM_EVENT];

/* read matrix size */
    n = atoi(argv[1]);
    A = malloc(n * n * sizeof(double));
    assert(A != NULL);

    /* fill matrix */
    for (i = 0; i < n; i++) {
        A[IDX(i, i, n)] = 1.0;
    }


    if (argc>2){
        measure = atoi(argv[2]);
    }

    /* measurments */
    if (measure == 0){
        dtime = dclock();
    }
    if (measure == 1){
        if((retval=PAPI_flops_rate(PAPI_FP_OPS,&ireal_time,&iproc_time,&iflpops,&imflops)) < PAPI_OK){
            printf("Could not initialise PAPI_flops \n");
            printf("Your platform may not support floating point operation event.\n");
            printf("retval: %d\n", retval);
            exit(1);
        }
    }
    if (measure == 2){
        /* initializing library */
        if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT ) {
            fprintf(stderr, "Error: %s\n", errstring);
            exit(1);
        }
        /* Creating event set   */
        if ((retval=PAPI_create_eventset(&EventSet)) != PAPI_OK){
            ERROR_RETURN(retval)
        }
        /* Add the array of events PAPI_TOT_INS and PAPI_TOT_CYC to the eventset*/
        if ((retval=PAPI_add_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK){
            ERROR_RETURN(retval)
        }
        /* Start counting */
        if ( (retval=PAPI_start(EventSet)) != PAPI_OK){
            ERROR_RETURN(retval)
        }
    }

    if (chol(A, n)) {
        fprintf(stderr,"Matrix is not symmetric or not positive definite\n");
    }


    if (measure == 0){
        dtime = dclock()-dtime;
        printf( "Time: %le \n", dtime);
    }
    if (measure == 1){
        if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK) {
            printf("retval: %d\n", retval);
            exit(1);
        }
        printf("Real_time: %f Proc_time: %f flpops: %lld MFLOPS: %f\n", real_time, proc_time,flpops,mflops);
    }
    if (measure == 2){
        /* Stop counting, this reads from the counter as well as stop it. */
        if ( (retval=PAPI_stop(EventSet,values)) != PAPI_OK){
            ERROR_RETURN(retval);
        }
        printf("\nThe total instructions executed are %lld, total cycles %lld\n", values[0],values[1]);
    }

    fflush( stdout );

    free(A);
    return 0;
}



