#ifndef __AC2DR_CONFIG__
#define __AC2DR_CONFIG__

extern const double limit_ang;
extern const double min_ang;
extern const int ghostL, pmaximgs, image10;
extern const double cfl, sg_dc;
extern const int wtap, wzero;
extern const double stg1_t0, stg2_t0;
extern const double damp_t0, damp_rise, damp_vel, dcoef;
extern const double wtap2;

//enforce boundary: 0=no, 1=yes
extern const int BC_ENFORCE;

//using manufactured solution. 0: no, 1: solution 1, 2: solution 2
#define MANUFACTURED_SOL 0

//the order of spatial derivatives
//6=sbp6, 4=4th-order staggered, 2=central finite difference (2nd)
#define FD_ORDER 6

//time marching: 1=RK4 2=FD4
#define TMARCH 1

#endif
