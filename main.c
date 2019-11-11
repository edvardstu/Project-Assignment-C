#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#include "utilities.h"
#include "interactions.h"

#define R 10.0

#define N_PARTICLES 10
#define N_STEPS 1000
#define U_0 20.0
#define D_R 10.0
#define DT 0.002

//Diffusive parameters
#define GAMMA_T 1
#define GAMMA_R 1

//Boundary interatction
#define LAMBDA_HAR 20 //FS
#define KAPPA_HAR 2  //GS

//Particle particle interaction
#define GAMMA_PP 0.1 //GAM
#define R_C 1
//#define SIGMA_PP pow(1/2, 1/6)

const double a = sqrt(3);



int main(int argc, char **argv) {
    //////////////////////Setting up variables ////////////////////////
    FILE *fp;
    double x[N_PARTICLES], y[N_PARTICLES], theta[N_PARTICLES];

    unsigned int t, index_p, index_n, i;
    const char * restrict fileName = "data.txt";
    ///////////////////////////////////////////////////////////////////

    ///////////////////Opening file for results ///////////////////////
    openFile(fileName, &fp);
    //fp = fopen(fileName, "w");
    ///////////////////////////////////////////////////////////////////

    //////////////////////Setting up GSL RNG ////////////////////////
    const gsl_rng_type * T;
    gsl_rng * r;
    setUpRNG(&T, &r);
    /////////////////////////////////////////////////////////////////

    //////////Setting up particle position and angle ////////////////
    sunflower(x, y, N_PARTICLES, 2, R);
    for (i=0;i<N_PARTICLES;i++) theta[i]=randDouble(-M_PI, M_PI, &r);
    for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf\n", i, x[i], y[i], theta[i]);
    /////////////////////////////////////////////////////////////////


    //For each time step
    for (t = 1; t < N_STEPS; t++){
        for (index_p = 0; index_p < N_PARTICLES; index_p++){
            //Find forces and torque from wall
            //Assume circular potentail has a centre in (0,0)
            double r_coord, fx_b, fy_b, torque_b; //Could be moved out of scoope depending on OpenMP
            fx_b = 0;
            fy_b = 0;
            torque_b = 0;
            r_coord = sqrt(x[index_p]*x[index_p]+y[index_p]*y[index_p]);
            if (r_coord > R){
                forceHarmonicCircular(&fx_b, &fy_b, r_coord, x[index_p], y[index_p], R, LAMBDA_HAR);
                torqueHarmonicCircular(&torque_b, r_coord, x[index_p], y[index_p], theta[index_p], LAMBDA_HAR, KAPPA_HAR);
            }







            //Update particle parameters
            x[index_p] = x[index_p] + (U_0*cos(theta[index_p]) + fx_b)*DT;
            y[index_p] = y[index_p] + (U_0*sin(theta[index_p]) + fy_b)*DT;
            theta[index_p] = theta[index_p] + sqrt(2*D_R*DT)*randDouble(-a, a, &r);
        }
        for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf\n", i, x[i], y[i], theta[i]);
    }






    ///////////////////Closing file with results //////////////////////
    closeFile(fileName, &fp);
    ///////////////////////////////////////////////////////////////////
    //Freeing RNG
    gsl_rng_free (r);
    //fclose(fp);

    return 0;
}
