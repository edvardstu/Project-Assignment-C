#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "interactions.h"

//In all the following functions it is assumed that the particles are all ready within range of the potentials
//This is to improve computational time

void forceHarmonicCircular(double *fx_b, double *fy_b, double r_coord, double x, double y, double r_boundary, double lambda_harmonic){
    //Assume the circle is placed in (x,y)=(0,0)
    double beta = atan2(y, x);
    double f = lambda_harmonic*(r_coord-r_boundary);

    *fx_b = -f*x/r_coord;
    *fy_b = -f*y/r_coord;
}

void torqueHarmonicCircular(double *torque_b, double r_coord, double x, double y, double theta, double lambda_harmonic, double kappa_harmonic){
    double beta = atan2(y, x);
    *torque_b = lambda_harmonic*kappa_harmonic*sin(2*(theta-beta));
}

void forceWeeksChandlerAndersen(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y){
    //Here it is assumed that sigma=1/2^(1/6) which gives a particle radius of 1
    double r_pn_6 = r_pn_2*r_pn_2*r_pn_2;
    double f = 12*(1.0/r_pn_6-1.0)/(r_pn_6*r_pn_2);
    *fx_n = f*delta_x;
    *fy_n = f*delta_y;
}

void forceHarmonicPP(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double r_cut_off_force, double lambda_pp){
    double f = -lambda_pp*(sqrt(r_pn_2)-r_cut_off_force);
    *fx_n = f*delta_x;
    *fy_n = f*delta_y;
}

void torqueWeeksChandlerAndersen(double *torque_n, double theta_p, double theta_n, double gamma_pp){
    *torque_n = gamma_pp*sin(theta_n-theta_p);
}
