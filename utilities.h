#ifndef UTILITIES_H
#define UTILITIES_H

#include <gsl/gsl_rng.h>

extern int errno;

void testFunc();

unsigned long int random_seed();

void setUpRNG(const gsl_rng_type **T, gsl_rng **r);

double randDouble(double min, double max, gsl_rng ** r);

void openFile(const char * restrict fileName, FILE** fp);

void closeFile(const char * restrict fileName, FILE** fp);

void sunflower(double x[], double y[], int n_particles, int alpha, double R);

double radius(int index, int n_particles, int n_boundary, double R);

double walltime();



#endif
