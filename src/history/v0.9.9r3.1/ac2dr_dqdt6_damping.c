#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void dqdt6_damping ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage ) {
    double *phix, *phiy, dampux=0.0, dampvx=0.0, dampwx=0.0, dampuy=0.0, dampvy=0.0, dampwy=0.0;
    unsigned long i, j;
    unsigned long nx, ny;
    int gi;
    double r, th, dr, dth, thmin, R0; // thlim;
    double dur, dvr, duth, dvth, dwth; //dwr
    double dwhr, dwhth;
    double c12, c13, c23, c32, c33;
    double phixi, phiyj;
    double a[7] = {-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0};
    double a2[7] = {1.0/90.0, -3.0/20.0, 3.0/2.0, -49.0/18.0, 3.0/2.0, -3.0/20.0, 1.0/90.0};
    double d2vth, d2wth, d2vr, d2wr; //d2uth, d2ur
    double *at, ax;
    double rs, ths, RR, RT;
    int ti;

//    double my_pi;

//    my_pi = 3.141592653589793238462643383279502884;

    nx = minfo->nth;
    ny = minfo->nr;
    dr = minfo->dr;
    dth = minfo->dth;
    thmin = minfo->thmin_global;
    R0 = minfo->R;

    phix = minfo->phix;
    phiy = minfo->phiy;

    ax = minfo->ax;
    at = minfo->at;

    for(j=ghostL; j<(ny-ghostL); j++) {
        for(i=ghostL; i<(nx-ghostL); i++) {
            gi = minfo->xi_start + i;
            phixi=phix[i];
            phiyj=phiy[j];

            dur = (a[0]*u[i][j-3] + a[1]*u[i][j-2] + a[2]*u[i][j-1] + \
                   a[4]*u[i][j+1] + a[5]*u[i][j+2] + a[6]*u[i][j+3])/dr * phiyj;

            dvr = (a[0]*v[i][j-3] + a[1]*v[i][j-2] + a[2]*v[i][j-1] + \
                   a[4]*v[i][j+1] + a[5]*v[i][j+2] + a[6]*v[i][j+3])/dr * phiyj;

//            dwr = (a[0]*w[i][j-3] + a[1]*w[i][j-2] + a[2]*w[i][j-1] + 
//                   a[4]*w[i][j+1] + a[5]*w[i][j+2] + a[6]*w[i][j+3])/dr * phiyj;

            dwhr = (a[0]*what[i][j-3] + a[1]*what[i][j-2] + a[2]*what[i][j-1] + \
                    a[4]*what[i][j+1] + a[5]*what[i][j+2] + a[6]*what[i][j+3])/dr * phiyj;

            duth = (a[0]*u[i-3][j] + a[1]*u[i-2][j] + a[2]*u[i-1][j] + \
                    a[4]*u[i+1][j] + a[5]*u[i+2][j] + a[6]*u[i+3][j])/dth *phixi;

            dwth = (a[0]*w[i-3][j] + a[1]*w[i-2][j] + a[2]*w[i-1][j] + \
                    a[4]*w[i+1][j] + a[5]*w[i+2][j] + a[6]*w[i+3][j])/dth *phixi;

            dvth = (a[0]*v[i-3][j] + a[1]*v[i-2][j] + a[2]*v[i-1][j] + \
                    a[4]*v[i+1][j] + a[5]*v[i+2][j] + a[6]*v[i+3][j])/dth *phixi;

            dwhth = (a[0]*what[i-3][j] + a[1]*what[i-2][j] + a[2]*what[i-1][j] + \
                     a[4]*what[i+1][j] + a[5]*what[i+2][j] + a[6]*what[i+3][j])/dth *phixi;

//            d2uth = (a2[0]*u[i-3][j] + a2[1]*u[i-2][j] + a2[2]*u[i-1][j] + a2[3]*u[i][j] + 
//                     a2[4]*u[i+1][j] + a2[5]*u[i+2][j] + a2[6]*u[i+3][j])/(dth*dth);

            d2vth = (a2[0]*v[i-3][j] + a2[1]*v[i-2][j] + a2[2]*v[i-1][j] + a2[3]*v[i][j] + \
                     a2[4]*v[i+1][j] + a2[5]*v[i+2][j] + a2[6]*v[i+3][j])/(dth*dth);

            d2wth = (a2[0]*w[i-3][j] + a2[1]*w[i-2][j] + a2[2]*w[i-1][j] + a2[3]*w[i][j] + \
                     a2[4]*w[i+1][j] + a2[5]*w[i+2][j] + a2[6]*w[i+3][j])/(dth*dth);

//            d2ur = (a2[0]*u[i][j-3] + a2[1]*u[i][j-2] + a2[2]*u[i][j-1] + a2[3]*u[i][j] + 
//                    a2[4]*u[i][j+1] + a2[5]*u[i][j+2] + a2[6]*u[i][j+3])/(dr*dr);

            d2vr = (a2[0]*v[i][j-3] + a2[1]*v[i][j-2] + a2[2]*v[i][j-1] + a2[3]*v[i][j] + \
                    a2[4]*v[i][j+1] + a2[5]*v[i][j+2] + a2[6]*v[i][j+3])/(dr*dr);

            d2wr = (a2[0]*w[i][j-3] + a2[1]*w[i][j-2] + a2[2]*w[i][j-1] + a2[3]*w[i][j] + \
                    a2[4]*w[i][j+1] + a2[5]*w[i][j+2] + a2[6]*w[i][j+3])/(dr*dr);


            r= R0 + dr*(j-ghostL);
            th = thmin + dth*(gi-ghostL);

            rs = R0+minfo->srcs[0].y;
            ths = minfo->srcs[0].x;

            RR = sqrt( pow(r*sin(th) -rs*sin(ths),2) + pow(r*cos(th) - rs*cos(ths), 2) );
            RT = RR/damp_vel;
            if(RT<damp_t0) {RT=damp_t0;}
            ti = tk-floor(RT/minfo->dt+0.5);
            if(ti < 0) { ti = 0;}
            ti = 0;
            
            //du update
            c12 = 2.0*rhohat[i][j]*chat[i][j]*chat[i][j]/r;

            dampux = dampux4( i, j, u, minfo );
            dampuy = dampuy4( i, j, u, minfo );

            if(gi==ghostL) {

            c13 = rhohat[i][j]*chat[i][j]*chat[i][j]/r;
            u1[i][j] = -rhohat[i][j]*chat[i][j]*chat[i][j]*dvr  \
                       -1.0/r*what[i][j]*duth \
                       -1.0/r*rhohat[i][j]*chat[i][j]*chat[i][j]*dwth \
                       -(c12*v[i][j] + c13*dwth) \
                       + dampux + dampuy;
            //           + at[ti]*ax*pow(dth, 2)*(r*r)*(2.0/(r*r)*d2uth) 
            //           + at[ti]*ax*pow(dr, 2)*(d2ur + 2.0/r * dur);
            } else {
            c13 = rhohat[i][j]*chat[i][j]*chat[i][j]/r*cos(th)/sin(th);
            u1[i][j] = -rhohat[i][j]*chat[i][j]*chat[i][j]*dvr  \
                       -1.0/r*what[i][j]*duth \
                       -1.0/r*rhohat[i][j]*chat[i][j]*chat[i][j]*dwth \
                       -(c12*v[i][j] + c13*w[i][j]) \
                       + dampux + dampuy;
            //           + at[ti]*ax*pow(dth,2)*(r*r)*(1.0/(r*r)*d2uth + cos(th)/sin(th)/(r*r)*duth) 
            //           + at[ti]*ax*pow(dr,2)*(d2ur + 2.0/r *dur);
            }

            //dv update
            c23 = -2.0*what[i][j]/r;

            dampvx = dampux4( i, j, v, minfo );
            dampvy = dampuy4( i, j, v, minfo );

            v1[i][j] = -1.0/rhohat[i][j]*dur \
                       -1.0/r*what[i][j]*dvth \
                       -c23*w[i][j] \
                       + dampvx + dampvy ;

            //dw update

            c32 = what[i][j]/r + dwhr;
            c33 = 1.0/r*dwhth;

            dampwx = dampux4( i, j, w, minfo );
            dampwy = dampuy4( i, j, w, minfo );

            w1[i][j] = -1.0/r/rhohat[i][j]*duth \
                       -1.0/r*what[i][j]*dwth \
                       -c32*v[i][j] \
                       -c33*w[i][j] \
                       + dampwx + dampwy;

            if(gi==ghostL) {
            v1[i][j] += at[ti]*ax*pow(dth,2)*(r*r)*(1.0/(r*r)*d2vth) \
                      + at[ti]*ax*pow(dr,2)*(d2vr);
            w1[i][j] += at[ti]*ax*pow(dth,2)*(r*r)*(1.0/(r*r)*d2wth) \
                      + at[ti]*ax*pow(dr,2)*(d2wr);

            } else {
            v1[i][j] += at[ti]*ax*pow(dth,2)*(r*r)*(1.0/(r*r)*d2vth) \
                      + at[ti]*ax*pow(dr,2)*(d2vr);
            w1[i][j] += at[ti]*ax*pow(dth,2)*(r*r)*(1.0/(r*r)*d2wth) \
                      + at[ti]*ax*pow(dr,2)*(d2wr);
            }
       }
   }

   //source update
   /*
   if(RK4stage == 2 || RK4stage == 3) {
   for(k=0; k<minfo->src_num; k++) {
      i=minfo->srcs[k].xi;
      j=minfo->srcs[k].yi;
      r= R0 + dr*(j-ghostL);
      u1[i][j] += minfo->srcs[k].midq->elements[tk]*rhohat[i][j]*chat[i][j]*chat[i][j] \
                /(M_PI*dr*pow(dth*r/2.0,2))/rhohat[i][j];
      //printf("***********ac2dr_dqdt, source test, %d, %d, %d, %e\n", tk, i, j, minfo->srcs[k].midq->elements[tk]);
   }
   }

   if(RK4stage == 1 || RK4stage == 4) {
   for(k=0; k<minfo->src_num; k++) {
      i=minfo->srcs[k].xi;
      j=minfo->srcs[k].yi;
      r= R0 + dr*(j-ghostL);
      u1[i][j] += minfo->srcs[k].q->elements[tk]*rhohat[i][j]*chat[i][j]*chat[i][j] \
                /(M_PI*dr*pow(dth*r/2.0,2))/rhohat[i][j];
      //printf("***********ac2dr_dqdt, source test, %d\n", tk);
   }
   }
   */
//    printf("***********ac2dr_dqdt.c is working\n");
    return;
}


