#include "ac2dr_aux.h"
#include "ac2dr_config.h"

//use limit less than this angle
//const double limit_ang = 1e-3;
const double limit_ang = 1e-15;
//min angle in the theta direction
const double min_ang = 0.0;

//enforce boundary: 0=no, 1=yes
const int BC_ENFORCE = 1;

//ghost layers: 3 layers = 3
const int ghostL = 3;

//maxp and maxpabs image
const int pmaximgs = 0;

//reduced image size
const int image10 = 0;

#if (TMARCH==1)
//const double dcoef = 3.7e3;
//const double cfl = 1.0;
const double sg_dc = 0.07;
#endif

#if (TMARCH==2)
//const double dcoef = 3.7e3;
const double cfl = 0.1;
const double sg_dc = 0.01;
#endif

//amplitude damping
const double damp_t0 = 100.0;  //start time for damping
const double damp_rise = 10.0; //rising time of damping coefficients
const double damp_vel = 100.0; //assumed propagation speed below which will be damped
const double dcoef = 0.1; //damping coefficient

//wind taper length
//stage 1: no-wind model
//stage 2: wind model with zeros near the axis

const double stg1_t0 = 19.0; //stage 1 starting time
//const int wzero = 3;
//const int wtap = 6;
const int wzero = 20;
const int wtap = 200;
const double wtap2 = 100000;
//const double stg2_t0 = wtap2/damp_vel; //stage 2 starting time
const double stg2_t0 = 1000; //stage 2 starting time





