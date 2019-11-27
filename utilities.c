#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>

#include "utilities.h"
#include "str_builder.h"


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

void readFile(const char * restrict fileName, FILE** fp){
    *fp = fopen(fileName, "r");
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

const char* restrict createFileNameBase(const char* restrict fileNameBase, bool overwrite){
    const char * restrict fileName;
    const char * TXT = ".txt";
    bool fileExists = true;
    int fileNumber = 0;

    str_builder_t *sb;
    sb = str_builder_create();
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_int(sb,fileNumber);
    str_builder_add_str(sb, TXT, 0);

    if (!overwrite){
        fileName = str_builder_peek(sb);
        while (fileExists){
            if( access( fileName, F_OK ) != -1 ) {
                if (fileNumber==0){
                    printf("Filename: '%s' already exsisted, and a numeral was added\n", fileName);
                }
                fileNumber++;
                str_builder_clear(sb);
                str_builder_add_str(sb, fileNameBase, 0);
                str_builder_add_int(sb,fileNumber);
                str_builder_add_str(sb, TXT, 0);
                fileName = str_builder_peek(sb);
            } else {
                fileExists = false;
            }

        }
    }

    str_builder_clear(sb);
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_int(sb,fileNumber);

    fileName = str_builder_dump(sb, NULL);
    str_builder_destroy(sb);

    return fileName;
}

const char* restrict createFileName(const char* restrict fileNameBase){
    const char * restrict fileName;
    const char * TXT = ".txt";
    str_builder_t *sb;
    sb = str_builder_create();
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_str(sb, TXT, 0);
    fileName = str_builder_dump(sb, NULL);
    str_builder_destroy(sb);
    return fileName;
}

void writeSimulationParameters(const char* restrict fileNameBase, double r, double r_particle, unsigned int n_particles, double u_0, double D_r, unsigned int n_steps, double dt, double gamma_t, double gamma_r, double lambda_har, double kappa_har, double gamma_pp, double r_cut_off_torque_2, double lambda_pp, double r_cut_off_force, double sigma_pp){
    str_builder_t *sb;
    FILE *fp;
    const char * TXT = ".txt";
    sb = str_builder_create();
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_str(sb, "SimulationParameters", 0);
    str_builder_add_str(sb, TXT, 0);
    const char * fileName = str_builder_dump(sb, NULL);
    openFile(fileName, &fp);

    fprintf(fp, "System and particle parameters\n");
    fprintf(fp, "Radius of system: %.1f\n", r);
    fprintf(fp, "Radius of particles: %.3f\n", r_particle);
    fprintf(fp, "Number of particles: %d\n", n_particles);
    fprintf(fp, "Particle velocity: %.1f\n", u_0);
    fprintf(fp, "D_r: %.3f\n", D_r);

    fprintf(fp, "\nNummerical solver parameters\n");
    fprintf(fp, "Number of steps: %d\n", n_steps);
    fprintf(fp, "Time step: %.5f\n", dt);

    fprintf(fp, "\nDiffusion\n");
    fprintf(fp, "Translational diffusion: %.1f\n", gamma_t);
    fprintf(fp, "Rotational diffusion: %.1f\n", gamma_r);

    fprintf(fp, "\nBoundary interaction\n");
    fprintf(fp, "Lambda harmonic: %.1f\n", lambda_har);
    fprintf(fp, "Kappa harmonic: %.1f\n", kappa_har);

    fprintf(fp, "\nParticle-Particle interaction\n");
    fprintf(fp, "Gamma particle-particle: %.1f\n", gamma_pp);
    fprintf(fp, "Cut off radius for torque squared: %.1f\n", r_cut_off_torque_2);
    fprintf(fp, "Lambda particle-particle: %.1f\n", lambda_pp);
    fprintf(fp, "Cut off radius force: %.1f\n", r_cut_off_force);
    fprintf(fp, "Sigma particle-particle: %.3f\n", sigma_pp);



    closeFile(fileName, &fp);
    str_builder_destroy(sb);
}

void writeFinalState(const char* restrict fileNameBase, int n_particles, double time, double x[], double y[], double theta[], double vx[], double vy[], double d_r){
    str_builder_t *sb;
    FILE *fp;
    const char * TXT = ".txt";
    sb = str_builder_create();
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_str(sb, "FinalState", 0);
    str_builder_add_str(sb, TXT, 0);
    const char * fileName = str_builder_dump(sb, NULL);
    openFile(fileName, &fp);

    for (int i=0;i<n_particles;i++) fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", i, time, x[i], y[i], theta[i], vx[i], vy[i], d_r);

    closeFile(fileName, &fp);
    str_builder_destroy(sb);
}

void swapPointers(double **Y_i, double **Y_i_prev){
    double *temp = *Y_i;
    *Y_i = *Y_i_prev;
    *Y_i_prev = temp;
}

const char* restrict createFileNamePrevious(const char* restrict fileNameBase, bool overwrite){
    const char * restrict fileName;
    const char * TXT = ".txt";
    bool fileExists = true;
    int fileNumber = 0;

    str_builder_t *sb;
    sb = str_builder_create();
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_int(sb,fileNumber);
    str_builder_add_str(sb, TXT, 0);

    if (!overwrite){
        fileName = str_builder_peek(sb);
        while (fileExists){
            if( access( fileName, F_OK ) != -1 ) {
                fileNumber++;
                str_builder_clear(sb);
                str_builder_add_str(sb, fileNameBase, 0);
                str_builder_add_int(sb,fileNumber);
                str_builder_add_str(sb, TXT, 0);
                fileName = str_builder_peek(sb);
            } else {
                fileExists = false;
            }

        }
    }

    str_builder_clear(sb);
    str_builder_add_str(sb, fileNameBase, 0);
    str_builder_add_int(sb, fileNumber-1);
    str_builder_add_str(sb, "FinalState.txt", 0);

    fileName = str_builder_dump(sb, NULL);
    str_builder_destroy(sb);

    return fileName;
}

void readInitialState(const char* restrict fileName, int n_particles, double *time, double x[], double y[], double theta[], double vx[], double vy[], double *d_r){
    FILE *fp;
    readFile(fileName, &fp);
    int dummy;
    int status;
    for (int i=0;i<n_particles;i++){
        status = fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", &dummy, time, &x[i], &y[i], &theta[i], &vx[i], &vy[i], d_r);
        if (status != 8){
            printf("Line was read incorrectly, number of elements was %d, not 8\n", status);
        }
    }
    closeFile(fileName, &fp);

}
