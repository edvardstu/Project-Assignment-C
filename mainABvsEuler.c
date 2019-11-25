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



#define R 6.0//17.0
#define R_PARTICLE 0.5
#define N_PARTICLES 100//1000
#define U_0 10.0
#define D_R 0.0//0.2

#define FACTOR 100
#define N_STEPS 400*FACTOR//400000//300000
#define DT 0.01/FACTOR //0.0006

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
    double x[N_PARTICLES], y[N_PARTICLES], theta[N_PARTICLES];
    double fs_scale;
    //Loop parameters
    unsigned int t, index_p, index_n, i;
    //Parameters for boundary
    double r_coord, fx_b, fy_b, torque_b;
    //Parameters for particle particle interaction
    double delta_x, delta_y, temp_fx_n, temp_fy_n, temp_torque_n, r_pn_2;

    //const char * restrict fileNameBase = "results/testAB";
    const char * restrict fileNameBase = "/home/edvardst/Documents/NTNU/Programming/Project_Assignment/C_plots/Results/Integrators/AB";
    const bool overwrite = true;
    const char * restrict fileName;

    //Helping variables for Adams_Bashforth
    double *Y_x = malloc(N_PARTICLES*sizeof(double));
    double *Y_x_prev = malloc(N_PARTICLES*sizeof(double));
    double *Y_y = malloc(N_PARTICLES*sizeof(double));
    double *Y_y_prev = malloc(N_PARTICLES*sizeof(double));
    double *Y_th = malloc(N_PARTICLES*sizeof(double));
    double *Y_th_prev = malloc(N_PARTICLES*sizeof(double));
    //double Y_x[N_PARTICLES], Y_y[N_PARTICLES], Y_th[N_PARTICLES];
    //double Y_x_prev[N_PARTICLES], Y_y_prev[N_PARTICLES], Y_th_prev[N_PARTICLES];
    ///////////////////////////////////////////////////////////////////


    ///////////////////Opening file for results ///////////////////////
    fileNameBase = createFileNameBase(fileNameBase, overwrite);
    fileName = createFileName(fileNameBase);
    writeSimulationParameters(fileNameBase, R, R_PARTICLE, N_PARTICLES, U_0, D_R, N_STEPS, DT, GAMMA_T, GAMMA_R, LAMBDA_HAR, KAPPA_HAR, GAMMA_PP, R_CUT_OFF_TORQUE_2, LAMBDA_PP, R_CUT_OFF_FORCE, SIGMA_PP);



    openFile(fileName, &fp);
    //fp = fopen(fileName, "w");
    ///////////////////////////////////////////////////////////////////

    //////////////////////Setting up GSL RNG ////////////////////////
    const gsl_rng_type * T;
    gsl_rng * r;
    setUpRNG(&T, &r);
    /////////////////////////////////////////////////////////////////

    //////////Setting up particle position and angle ////////////////
    //sunflower(x, y, N_PARTICLES, 0, R);
    for (i=0;i<N_PARTICLES;i++) theta[i]=randDouble(-M_PI, M_PI, &r);
    for (i=0;i<N_PARTICLES;i++){
        double radius = randDouble(0, R, &r);
        x[i] = radius*cos(theta[i]);
        x[i] = radius*cos(theta[i]);
    }
    //for (i=0;i<N_PARTICLES;i++) theta[i]=atan2(y[i], x[i]) - M_PI/2 + randDouble(-M_PI/10, M_PI/10, &r);
    for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf\n", i, x[i], y[i], theta[i]);

    /*for (i=0;i<N_PARTICLES;i++){
        theta[i] = i*M_PI/N_PARTICLES;
        x[i] = 10.0 - 5.0*cos(theta[i]);
        y[i] = 10.0 - 5.0*sin(theta[i]);
    }*/
    /////////////////////////////////////////////////////////////////


    fs_scale = 0.0;
    //For each time step
    for (t = 1; t <= N_STEPS; t++){
        if (t % 1000 ==0 ) printf("%d\n",t);

        /*if (t <= 50000){
          fs_scale=0.01+t*0.99/50000.0;
        }*/
        fs_scale=1.0;

        double fx_n[N_PARTICLES] = {0};
        double fy_n[N_PARTICLES] = {0};
        double torque_n[N_PARTICLES] = {0};
        double number_n[N_PARTICLES] = {0};
        /*for (index_p = 0; index_p < N_PARTICLES; index_p++){
            fx_n[N_PARTICLES] = 0;// = {0};
            fy_n[N_PARTICLES] = 0;// = {0};
            torque_n[N_PARTICLES] = 0;// = {0};
            number_n[N_PARTICLES] = 0;// = {0};
        }*/

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

            //if (abs((Y_th[index_p] - Y_th_prev[index_p])<0.000001) && (t>1)) printf("%d, %d, AB wrong\n", t, index_p);

            //Update particle parameters
            x[index_p] = x[index_p] + (3/2*Y_x[index_p] - 1/2*Y_x_prev[index_p])*DT;
            y[index_p] = y[index_p] + (3/2*Y_y[index_p] - 1/2*Y_y_prev[index_p])*DT;
            theta[index_p] = theta[index_p] + (3/2*Y_th[index_p] - 1/2*Y_th_prev[index_p])*DT + sqrt(2*D_R*DT)*randDouble(-a, a, &r);



        }
        if (t % FACTOR ==0){
        //if (t % 100 ==0){
            for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf\n", i, x[i], y[i], theta[i]);
        }
//        printf("Y_x\t\t");
//        for (i=0;i<N_PARTICLES;i++) printf("%f\t", Y_x[i]);
//        printf("\n");
        swapPointers(Y_x, Y_x_prev);
        swapPointers(Y_y, Y_y_prev);
        swapPointers(Y_th, Y_th_prev);
//        printf("Y_x swapped\t");
//        for (i=0;i<N_PARTICLES;i++) printf("%f\t", Y_x[i]);
//        printf("\n");
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
