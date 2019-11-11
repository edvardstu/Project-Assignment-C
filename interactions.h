#ifndef INTERACTIONS_H
#define INTERACTIONS_H

void forceHarmonicCircular(double *fx_p, double *fy_p, double r_coord, double x, double y, double r_boundary, double lambda_harmonic);

void torqueHarmonicCircular(double *torque_b, double r_coord, double x, double y, double theta, double lambda_harmonic, double kappa_harmonic);

void forceWeeksChandlerAndersen(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y);

void torqueWeeksChandlerAndersen(double *torque_n, double theta_p, double theta_n, double gamma_pp);

#endif
