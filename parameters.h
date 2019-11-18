#ifndef PARAMETERS_H
#define PARAMETERS_H

const double R = 17.0;
extern const double R_PARTICLE = 0.5;

extern const unsigned int N_PARTICLES = 1000;
extern const unsigned int N_STEPS = 1000; //500000
extern const double U_0 = 10.0;
extern const double D_R = 0.2;
extern const double DT = 0.0006;

//Diffusive parameters
extern const double GAMMA_T = 1.0;
extern const double GAMMA_R = 1.0;

//Boundary interatction
extern const double LAMBDA_HAR = 200.0; //20.0 //FS
extern const double KAPPA_HAR = 10.0;   //GS

//Particle particle interaction
extern const double GAMMA_PP = 10.0; //1.0 //GAM
extern const double R_CUT_OFF_TORQUE_2 = 4.0; //9.0
extern const double LAMBDA_PP = 40.0;
extern const double R_CUT_OFF_FORCE = 2.0;  // (R_CUT_OFF_FORCE)^2 > R_CUT_OFF_TORQUE_2
//#define SIGMA_PP pow(1/2, 1/6)
extern const double SIGMA_PP = 1.0; //(SIGMA_PP*2)^2 > R_CUT_OFF_TORQUE_2

extern const double a = sqrt(3);

#endif
