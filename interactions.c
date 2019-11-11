#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void forceHarmonicCircular(double *fx_b, double *fy_b, double r_coord, double x, double y, double r_boundary, double lambda_harmonic){
    //Assume the circle is placed in (x,y)=(0,0)
    double beta = atan2(y, x);
    double f = lambda_harmonic*(r_coord-r_boundary);

    *fx_b += -f*x/r_coord; //Should be changed to only =
    *fy_b += -f*y/r_coord;
}

void torqueHarmonicCircular(double *torque_b, double r_coord, double x, double y, double theta, double lambda_harmonic, double kappa_harmonic){
    double beta = atan2(y, x);
    *torque_b = lambda_harmonic*kappa_harmonic*sin(2*(theta-beta));
}
