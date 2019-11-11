#ifndef INTERACTIONS_H
#define INTERACTIONS_H

void forceHarmonicCircular(double *fx_p, double *fy_p, double r_coord, double x, double y, double r_boundary, double lambda_harmonic);

void torqueHarmonicCircular(double *torque_b, double r_coord, double x, double y, double theta, double lambda_harmonic, double kappa_harmonic);

#endif
