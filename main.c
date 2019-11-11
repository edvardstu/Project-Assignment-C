#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include "utilities.h"

#define R 10.0

#define N_PARTICLES = 10
#define N_STEPS = 100
#define U_0 = 1.0
#define D_R = 10.0
#define DT = 0.0002

//Diffusive parameters
#define GAMMA_T = 1
#define GAMMA_R = 1

//Boundary interatction
#define LAMBDA_HAR = 10 //FS
#define KAPPA_HAR = 2  //GS

//Particle particle interaction
#define GAMMA_PP 0.1 //GAM
#define R_C 1
//#define SIGMA_PP pow(1/2, 1/6)

extern int errno;

unsigned long int random_seed();
double randDouble(double min, double max, gsl_rng * r);

int main(int argc, char **argv) {
    //////////////////////Setting up variables ////////////////////////
    FILE *fp;
    ///////////////////////////////////////////////////////////////////

    ///////////////////Opening file for results ///////////////////////
    const char * restrict fileName = "data.txt";

    fp = fopen(fileName, "w");
    if (fp == NULL){
        int errnum = errno;
        fprintf(stderr, "Value of errno %d\n", errno);
        fprintf(stderr, "Error: %s\n", strerror(errno));
        return -1;
    }
    ///////////////////////////////////////////////////////////////////

    //////////////////////Setting up GSL RNG ////////////////////////
    // gsl_rng_uniform (r) returns a random float [0,1)
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, random_seed());
    printf ("RNG type: %s\n", gsl_rng_name (r));
    /////////////////////////////////////////////////////////////////

    ///////////////////Closing file with results //////////////////////
    if (fclose(fp)!=0){
        fprintf(stderr, "Error closing file, %c\n", *fileName);
    }
    ///////////////////////////////////////////////////////////////////



    int i, n = 10;
    double min = -3.14;
    double max = -min;
    for (i = 0; i < n; i++)
      {
        double u = randDouble(min, max, r);
        printf ("%.5f\n", u);
      }



    //Freeing RNG
    gsl_rng_free (r);

    return 0;
}

unsigned long int random_seed(){
 struct timeval tv;
 gettimeofday(&tv,0);
 return (tv.tv_sec + tv.tv_usec);
}

double randDouble(double min, double max, gsl_rng * r){
    double u = gsl_rng_uniform(r);
    return u*(max-min)+min;
}
