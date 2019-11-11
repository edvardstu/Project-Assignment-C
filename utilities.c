#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>

#include "utilities.h"

void testFunc(){
    printf("Hello World");
}

unsigned long int random_seed(){
 struct timeval tv;
 gettimeofday(&tv,0);
 return (tv.tv_sec + tv.tv_usec);
}

void setUpRNG(const gsl_rng_type **T, gsl_rng **r){
    gsl_rng_env_setup();
    *T = gsl_rng_default;
    *r = gsl_rng_alloc (*T);
    gsl_rng_set(*r, random_seed());
    printf ("RNG type: %s\n", gsl_rng_name(*r));
}

double randDouble(double min, double max, gsl_rng ** r){
    //Doing alot of extra work here...
    double u = gsl_rng_uniform(*r);
    return u*(max-min)+min;
}

void openFile(const char * restrict fileName, FILE** fp){
    *fp = fopen(fileName, "w");
    if (fp == NULL){
        int errnum = errno;
        fprintf(stderr, "Value of errno %d\n", errno);
        fprintf(stderr, "Error: %s\n", strerror(errno));
        exit(-1);
    } else {
        printf("Opened the file: %s\n", fileName);
    }
}

void closeFile(const char * restrict fileName, FILE** fp){
    if (fclose(*fp)!=0){
        fprintf(stderr, "Error closing file, %s\n", fileName);
        exit(-1);
    } else {
        printf("File closed and saved to: %s\n", fileName);
    }
}

void sunflower(double x[], double y[], int n_particles, int alpha, double R){
    double n_boundary = floor(alpha*sqrt(n_particles));
    double phi = (sqrt(5)+1)/2;
    double r;
    double theta;
    for (int index = 1; index <=n_particles; index++){
        r = radius(index, n_particles, n_boundary, R);
        theta = 2*M_PI*index/pow(phi, 2);
        x[index-1]= r*cos(theta);
        y[index-1]= r*sin(theta);
    }
}

double radius(int index, int n_particles, int n_boundary, double R){
    if (index > (n_particles-n_boundary)){
        return R;
    } else {
        return R*sqrt((double)index-0.5)/sqrt((double)n_particles-(double)(n_boundary-1)/2);
    }
}

double walltime(){
	static struct timeval t;
	gettimeofday ( &t, NULL );
	return ( t.tv_sec + 1e-6 * t.tv_usec );
}
