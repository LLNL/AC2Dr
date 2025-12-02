#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void dqdt12 ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage ) {
    double *phix, *phiy, dampux=0.0, dampvx=0.0, dampwx=0.0, dampuy=0.0, dampvy=0.0, dampwy=0.0;
    unsigned long i, j;
    unsigned long nx, ny;
    int gi;
    double r, th, dr, dth, thmin, R0; // thlim;
    double dur, dvr, duth, dvth, dwth;
    double dwhr, dwhth;
    double c12, c13, c23, c32, c33;
    double phixi, phiyj;
    double a[13] = {1.0/5544.0, -3.0/1155.0, 1.0/56.0, -5.0/63.0, 15.0/56.0, -6.0/7.0 ,0.0, \
                    6.0/7.0, -15.0/56.0, 5.0/63.0, -1.0/56.0, 3.0/1155.0, -1.0/5544.0};
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

    for(j=6; j<(ny-6); j++) {
        for(i=6; i<(nx-6); i++) {
            gi = minfo->xi_start + i;
            phixi=phix[i];
            phiyj=phiy[j];

            dur = (a[0]*u[i][j-6] + a[1]*u[i][j-5] + a[2]*u[i][j-4] + a[3]*u[i][j-3] + a[4]*u[i][j-2] + a[5]*u[i][j-1] + \
                   a[7]*u[i][j+1] + a[8]*u[i][j+2] + a[9]*u[i][j+3] + a[10]*u[i][j+4] + a[11]*u[i][j+5] + a[12]*u[i][j+6])/dr * phiyj;

            dvr = (a[0]*v[i][j-6] + a[1]*v[i][j-5] + a[2]*v[i][j-4] + a[3]*v[i][j-3] + a[4]*v[i][j-2] + a[5]*v[i][j-1] + \
                   a[7]*v[i][j+1] + a[8]*v[i][j+2] + a[9]*v[i][j+3] + a[10]*v[i][j+4] + a[11]*v[i][j+5] + a[12]*v[i][j+6])/dr * phiyj;

            dwhr = (a[0]*what[i][j-6] + a[1]*what[i][j-5] + a[2]*what[i][j-4] + a[3]*what[i][j-3] + a[4]*what[i][j-2] + a[5]*what[i][j-1] + \
                   a[7]*what[i][j+1] + a[8]*what[i][j+2] + a[9]*what[i][j+3] + a[10]*what[i][j+4] + a[11]*what[i][j+5] + a[12]*what[i][j+6])/dr * phiyj;

            duth = (a[0]*u[i-6][j] + a[1]*u[i-5][j] + a[2]*u[i-4][j] + a[3]*u[i-3][j] + a[4]*u[i-2][j] + a[5]*u[i-1][j] + \
                    a[7]*u[i+1][j] + a[8]*u[i+2][j] + a[9]*u[i+3][j] + a[10]*u[i+4][j] + a[11]*u[i+5][j] + a[12]*u[i+6][j])/dth *phixi;

            dwth = (a[0]*w[i-6][j] + a[1]*w[i-5][j] + a[2]*w[i-4][j] + a[3]*w[i-3][j] + a[4]*w[i-2][j] + a[5]*w[i-1][j] + \
                    a[7]*w[i+1][j] + a[8]*w[i+2][j] + a[9]*w[i+3][j] + a[10]*w[i+4][j] + a[11]*w[i+5][j] + a[12]*w[i+6][j])/dth *phixi;

            dvth = (a[0]*v[i-6][j] + a[1]*v[i-5][j] + a[2]*v[i-4][j] + a[3]*v[i-3][j] + a[4]*v[i-2][j] + a[5]*v[i-1][j] + \
                    a[7]*v[i+1][j] + a[8]*v[i+2][j] + a[9]*v[i+3][j] + a[10]*v[i+4][j] + a[11]*v[i+5][j] + a[12]*v[i+6][j])/dth *phixi;

            dwhth = (a[0]*what[i-6][j] + a[1]*what[i-5][j] + a[2]*what[i-4][j] + a[3]*what[i-3][j] + a[4]*what[i-2][j] + a[5]*what[i-1][j] + \
                    a[7]*what[i+1][j] + a[8]*what[i+2][j] + a[9]*what[i+3][j] + a[10]*what[i+4][j] + a[11]*what[i+5][j] + a[12]*what[i+6][j])/dth *phixi;

            r= R0 + dr*(j-ghostL);
            th = thmin + dth*(gi-ghostL);
            
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

            } else {
            c13 = rhohat[i][j]*chat[i][j]*chat[i][j]/r*cos(th)/sin(th);
            u1[i][j] = -rhohat[i][j]*chat[i][j]*chat[i][j]*dvr  \
                       -1.0/r*what[i][j]*duth \
                       -1.0/r*rhohat[i][j]*chat[i][j]*chat[i][j]*dwth \
                       -(c12*v[i][j] + c13*w[i][j]) \
                       + dampux + dampuy;
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


