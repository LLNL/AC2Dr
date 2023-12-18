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

//wind taper length
const int wzero = 5;
const int wtap = 50;

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



