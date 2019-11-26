#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include <omp.h>
#include <stdbool.h>

#include "utilities.h"
#include "interactions.h"
#include "str_builder.h"



#define R 17.0
#define R_PARTICLE 0.5
#define N_PARTICLES 1000
#define U_0 10.0
#define D_R_C 0.8

#define N_STEPS 50000//400000//300000
#define DT 0.0005 //0.0006

//Diffusive parameters
#define GAMMA_T 1
#define GAMMA_R 1

//Boundary interatction
#define LAMBDA_HAR 200.0//20.0 //FS
#define KAPPA_HAR 10.0  //GS

//Particle particle interaction
#define GAMMA_PP 10.0//1.0 //GAM
#define R_CUT_OFF_TORQUE_2 4.0 //9.0
#define LAMBDA_PP 40.0
#define R_CUT_OFF_FORCE 2.0  // (R_CUT_OFF_FORCE)^2 > R_CUT_OFF_TORQUE_2
//#define SIGMA_PP pow(1/2, 1/6)
//#define SIGMA_PP 1.0 //(SIGMA_PP*2)^2 > R_CUT_OFF_TORQUE_2
const double SIGMA_PP = 1.5;//sqrt(2.0);//1/sqrt(2);

const double a = sqrt(3);



int main(int argc, char **argv) {
    double time_start = walltime();

    bool continueFromPrev = false;
    for (int k=0; k<2; k++){
    printf("Time frame %d of %d\n", k+1, 6);
    //Check the number of particles compared to the size of the system
    double r_particle = (R_CUT_OFF_FORCE-U_0/LAMBDA_PP)/2;
    printf("The particle radius is %f\n", r_particle);
    int max_p = floor((R*R*0.9069)/(R_PARTICLE*R_PARTICLE));
    if (N_PARTICLES>max_p){
        printf("Too many particles\n");
        printf("Maximum number of particles is %d\n", max_p);
        //exit(-1);
    } else {
        printf("The system occupancy is %f\n", (double)N_PARTICLES/max_p);
    }


    //////////////////////Setting up variables ////////////////////////
    FILE *fp;
    double x[N_PARTICLES], y[N_PARTICLES], theta[N_PARTICLES], vx[N_PARTICLES], vy[N_PARTICLES];
    double fs_scale, time;
    //Loop parameters
    unsigned int t, index_p, index_n, i;
    //Parameters for boundary
    double r_coord, fx_b, fy_b, torque_b;
    //Parameters for particle particle interaction
    double delta_x, delta_y, temp_fx_n, temp_fy_n, temp_torque_n, r_pn_2;

    const char * restrict fileNameBase = "transient/vtest";
    //const char * restrict fileNameBase = "/home/edvardst/Documents/NTNU/Programming/Project_Assignment/C_plots/Results/Integrators/AB";
    const bool overwrite = false;
    const char * restrict fileName;

    //Helping variables for Adams_Bashforth
    double *Y_x = malloc(N_PARTICLES*sizeof(double));
    double *Y_x_prev = malloc(N_PARTICLES*sizeof(double));
    double *Y_y = malloc(N_PARTICLES*sizeof(double));
    double *Y_y_prev = malloc(N_PARTICLES*sizeof(double));
    double *Y_th = malloc(N_PARTICLES*sizeof(double));
    double *Y_th_prev = malloc(N_PARTICLES*sizeof(double));

    //Parameters for loading previos Results
    const char * restrict fileNamePrevious;
    //bool continueFromPrev = true;
    double D_R;
    double D_R_I;
    ///////////////////////////////////////////////////////////////////


    //////////////////////Creating file names /////////////////////////
    fileNamePrevious = createFileNamePrevious(fileNameBase, overwrite);
    fileNameBase = createFileNameBase(fileNameBase, overwrite);
    fileName = createFileName(fileNameBase);
    ///////////////////////////////////////////////////////////////////

    //////////////////////Setting up GSL RNG ////////////////////////
    const gsl_rng_type * T;
    gsl_rng * r;
    setUpRNG(&T, &r);
    /////////////////////////////////////////////////////////////////

    //////////Setting up particle position and angle ////////////////

    if (!continueFromPrev){
        time = 0;
        D_R = D_R_C;
        sunflower(x, y, N_PARTICLES, 0, R);
        for (i=0;i<N_PARTICLES;i++){
            theta[i]=randDouble(-M_PI, M_PI, &r);
            vx[i] = cos(theta[i]);
            vy[i] = sin(theta[i]);
        }
    } else {
        readInitialState(fileNamePrevious, N_PARTICLES, &time, x, y, theta, vx, vy, &D_R_I);
        D_R = D_R_I;
        if (D_R_C != D_R_I){
            printf("D_r was initilized to %.3f, but change to %.3f\n", D_R_C, D_R_I);
        }
    }

    /////////////////////////////////////////////////////////////////

    ///////////////////Opening file for results ///////////////////////
    writeSimulationParameters(fileNameBase, R, R_PARTICLE, N_PARTICLES, U_0, D_R, N_STEPS, DT, GAMMA_T, GAMMA_R, LAMBDA_HAR, KAPPA_HAR, GAMMA_PP, R_CUT_OFF_TORQUE_2, LAMBDA_PP, R_CUT_OFF_FORCE, SIGMA_PP);
    openFile(fileName, &fp);
    for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", i, time, x[i], y[i], theta[i], vx[i], vy[i], D_R);
    /////////////////////////////////////////////////////////////////


    fs_scale = 0.0;
    //For each time step
    for (t = 1; t <= N_STEPS; t++){
        if (t % 10000 ==0 ) printf("%d\n",t);

        if (t <= 50000){
            if (!continueFromPrev){
                fs_scale=0.01+t*0.99/50000.0;
            } else {
                fs_scale=1.0;
                D_R = D_R_I + 0.2*(t/50000.0);
            }
        }

        double fx_n[N_PARTICLES] = {0};
        double fy_n[N_PARTICLES] = {0};
        double torque_n[N_PARTICLES] = {0};
        double number_n[N_PARTICLES] = {0};

        for (index_p = 0; index_p < N_PARTICLES; index_p++){
            //Find forces and torque from wall
            //Assume circular potentail has a centre in (0,0)
            fx_b = 0;
            fy_b = 0;
            torque_b = 0;
            r_coord = sqrt(x[index_p]*x[index_p]+y[index_p]*y[index_p]);
            if (r_coord > R){
                forceHarmonicCircular(&fx_b, &fy_b, r_coord, x[index_p], y[index_p], R, LAMBDA_HAR);
                torqueHarmonicCircular(&torque_b, r_coord, x[index_p], y[index_p], theta[index_p], LAMBDA_HAR, KAPPA_HAR);
            }

            for (index_n = index_p + 1; index_n < N_PARTICLES; index_n++){
                delta_x = x[index_p]-x[index_n];
                delta_y = y[index_p]-y[index_n];
                r_pn_2 = delta_x*delta_x + delta_y*delta_y;

                /* For WeeksChandlerAndersen potential
                if (r_pn_2 < R_CUT_OFF_TORQUE_2){
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

                } */
                if (r_pn_2 < R_CUT_OFF_TORQUE_2){
                    //if (r_pn_2 < R_CUT_OFF_FORCE*R_CUT_OFF_FORCE){
                    //if (r_pn_2 < 1.0){
                    if (r_pn_2 < 2*SIGMA_PP*SIGMA_PP) {

                        //forceHarmonicPP(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, R_CUT_OFF_FORCE, LAMBDA_PP);
                        //forceOneOverRQuad(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                        forceOneOverRQuadSig(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, SIGMA_PP);
                        fx_n[index_p] += temp_fx_n;
                        fy_n[index_p] += temp_fy_n;
                        fx_n[index_n] -= temp_fx_n;
                        fy_n[index_n] -= temp_fy_n;
                    }
                    torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], GAMMA_PP, r_pn_2);
                    torque_n[index_p] += temp_torque_n;
                    torque_n[index_n] -= temp_torque_n;
                    number_n[index_n]++;
                    number_n[index_p]++;
                }

            }

            //Update Adams-Bashforth helping parameters
            Y_x[index_p] = U_0*cos(theta[index_p]) + (fx_b + fx_n[index_p])*fs_scale;
            Y_y[index_p] = U_0*sin(theta[index_p]) + (fy_b + fy_n[index_p])*fs_scale;
            if (number_n[index_p] > 0){
                Y_th[index_p] = (torque_b + torque_n[index_p]/number_n[index_p])*fs_scale;
            } else {
                Y_th[index_p] = torque_b*fs_scale;
            }

            if (t==1){
                Y_x_prev[index_p] = Y_x[index_p];
                Y_y_prev[index_p] = Y_y[index_p];
                Y_th_prev[index_p] = Y_th[index_p];
            }

            //Update particle parameters
            x[index_p] = x[index_p] + (3/2*Y_x[index_p] - 1/2*Y_x_prev[index_p])*DT;
            y[index_p] = y[index_p] + (3/2*Y_y[index_p] - 1/2*Y_y_prev[index_p])*DT;
            theta[index_p] = theta[index_p] + (3/2*Y_th[index_p] - 1/2*Y_th_prev[index_p])*DT + sqrt(2*D_R*DT)*randDouble(-a, a, &r);

            vx[i] = 3/2*Y_x[index_p] - 1/2*Y_x_prev[index_p];
            vy[i] = 3/2*Y_y[index_p] - 1/2*Y_y_prev[index_p];


        }
        time += DT;
        if (t % 500 ==0){
        //if (t % 100 ==0){
            for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", i, time, x[i], y[i], theta[i], vx[i], vy[i], D_R);
        }
        swapPointers(Y_x, Y_x_prev);
        swapPointers(Y_y, Y_y_prev);
        swapPointers(Y_th, Y_th_prev);
    }
    writeFinalState(fileNameBase, N_PARTICLES, time, x, y, theta, vx, vy, D_R);



    ///////////////////Closing file with results //////////////////////
    closeFile(fileName, &fp);
    ///////////////////////////////////////////////////////////////////
    //Freeing RNG
    gsl_rng_free (r);
    //fclose(fp);

    double time_end = walltime();
    printf("Simulation time: %7.3f s\n",(time_end-time_start));
    continueFromPrev=true;
}
return 0;
}
