#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>

#include "utilities.h"
#include "interactions.h"

#define R 35.0
#define R_PARTICLE 1.0

#define N_PARTICLES 1000
#define N_STEPS 10000
#define U_0 20.0
#define D_R 10.0
#define DT 0.0002

//Diffusive parameters
#define GAMMA_T 1
#define GAMMA_R 1

//Boundary interatction
#define LAMBDA_HAR 20 //FS
#define KAPPA_HAR 2  //GS

//Particle particle interaction
#define GAMMA_PP 0.1 //GAM
#define R_CUT_OFF_SQU 3
//#define SIGMA_PP pow(1/2, 1/6)

const double a = sqrt(3);



int main(int argc, char **argv) {
    double time_start = walltime();
    //Check the number of particles compared to the size of the system

    int max_p = floor((R*R*0.9069)/(R_PARTICLE*R_PARTICLE));
    if (N_PARTICLES>max_p){
        printf("Too many particles\n");
        printf("Maximum number of particles is %d\n", max_p);
        exit(-1);
    }


    //////////////////////Setting up variables ////////////////////////
    FILE *fp;
    double x[N_PARTICLES], y[N_PARTICLES], theta[N_PARTICLES];
    double fs_scale;
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
    sunflower(x, y, N_PARTICLES, 0, R);
    for (i=0;i<N_PARTICLES;i++) theta[i]=randDouble(-M_PI, M_PI, &r);
    for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf\n", i, x[i], y[i], theta[i]);
    /////////////////////////////////////////////////////////////////


    //For each time step
    for (t = 1; t < N_STEPS; t++){
        if (t % 1000 ==0 ) printf("%d\n",t);

        if (t <= 50000){
          fs_scale=0.01+t*0.99/50000.0;
        }


        double fx_n[N_PARTICLES] = {0};
        double fy_n[N_PARTICLES] = {0};
        double torque_n[N_PARTICLES] = {0};
        double number_n[N_PARTICLES] = {0};
        for (index_p = 0; index_p < N_PARTICLES; index_p++){
            double r_pn_2; //radius between particle and neighbour squared
            double delta_x, delta_y, temp_fx_n, temp_fy_n, temp_torque_n;
            for (index_n = index_p + 1; index_n < N_PARTICLES; index_n++){
                delta_x = x[index_p]-x[index_n];
                delta_y = y[index_p]-y[index_n];
                r_pn_2 = delta_x*delta_x + delta_y*delta_y;

                if (r_pn_2 < R_CUT_OFF_SQU){
                    if (r_pn_2 < 1.0){
                    //Here it is assumed that sigma=1/2^(1/6) which gives a particle radius of 1
                        forceWeeksChandlerAndersen(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                        fx_n[index_p] += temp_fx_n;
                        fy_n[index_p] += temp_fy_n;
                        fx_n[index_n] -= temp_fx_n;
                        fy_n[index_n] -= temp_fy_n;
                    }
                    torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], GAMMA_PP);
                    torque_n[index_p] += temp_torque_n;
                    torque_n[index_n] -= temp_torque_n;
                    number_n[index_n]++;
                    number_n[index_p]++;

                }
            }
        }

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
            x[index_p] = x[index_p] + (U_0*cos(theta[index_p]) + fx_b + fs_scale*fx_n[index_p])*DT;
            y[index_p] = y[index_p] + (U_0*sin(theta[index_p]) + fy_b + fs_scale*fy_n[index_p])*DT;
            if (number_n[index_p] > 0){
                theta[index_p] = theta[index_p] + sqrt(2*D_R*DT)*randDouble(-a, a, &r) + (torque_b + fs_scale*torque_n[index_p]/number_n[index_p])*DT;
            } else {
                theta[index_p] = theta[index_p] + sqrt(2*D_R*DT)*randDouble(-a, a, &r) + torque_b*DT;
            }

        }
        if (t % 50 ==0){
            for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf\n", i, x[i], y[i], theta[i]);
        }
    }






    ///////////////////Closing file with results //////////////////////
    closeFile(fileName, &fp);
    ///////////////////////////////////////////////////////////////////
    //Freeing RNG
    gsl_rng_free (r);
    //fclose(fp);

    double time_end = walltime();
    printf("Simulation time: %7.3f ms\n",(time_end-time_start)*1e3);
    return 0;
}
