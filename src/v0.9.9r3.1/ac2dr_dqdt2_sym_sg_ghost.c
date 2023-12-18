#include <omp.h>
#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void dqdt2_sym_ghost_sg ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage ) {
    double *phix, *phiy, dampux=0.0, dampvx=0.0, dampwx=0.0, dampuy=0.0, dampvy=0.0, dampwy=0.0;
    unsigned long i, j, k;
    unsigned long nx, ny;
    int gi;
    double r, th, dr, dth, thmin, R0; // thlim;
    double dur, dvr, duth, dvth, dwth;
    double dwhr, dwhth;
    double drhocr, drhocth, c11, c12, c13, c21, c23, c31, c32, c33;
    double phixi, phiyj;
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

    for(j=3; j<(ny-3); j++) {
        for(i=3; i<(nx-3); i++) {
            gi = minfo->xi_start + i;
            //du update
            phixi=phix[i];
            phiyj=phiy[j];
            dvr = ((v[i][j] - v[i][j-1])/dr) * phiyj;
            duth = ((u[i+1][j] - u[i-1][j])/(2*dth)) * phixi;
            dwth = ((w[i][j] - w[i-1][j])/(dth)) * phixi;

            drhocth = (((rhohat[i+1][j]*chat[i+1][j]) - (rhohat[i-1][j]*chat[i-1][j]))/(2.0*dth))*phixi;

            r= R0 + dr*(j-ghostL);
            th = thmin + dth*(gi-ghostL);
            
            c11 = 1.0/r*(what[i-1][j]+what[i][j])/2.0/(rhohat[i][j]*chat[i][j])*drhocth;
            c12 = 2.0*chat[i][j]/r;

            dampux = dampux4( i, j, u, minfo );
            dampuy = dampuy4( i, j, u, minfo );

            if(gi==3) {
            c13 = chat[i][j]/r;

            u1[i][j] = -chat[i][j]*dvr  \
                       -1.0/r*(-what[i][j]+what[i][j])/2.0*duth \
                       -1.0/r*chat[i][j]*dwth \
                       -(c11*u[i][j] + c12*(v[i][j]+v[i][j-1])/2.0) \
                       - c13*dwth \
                       + dampux + dampuy;
            } else {
            c13 = chat[i][j]/r*cos(th)/sin(th);
            u1[i][j] = -chat[i][j]*dvr  \
                       -1.0/r*(what[i-1][j]+what[i][j])/2.0*duth \
                       -1.0/r*chat[i][j]*dwth \
                       -(c11*u[i][j] + c12*(v[i][j]+v[i][j-1])/2.0 + c13*(w[i][j]+w[i-1][j])/2.0) \
                       + dampux + dampuy;
            }

            //dv update
            phixi=phix[i];
            phiyj=(phiy[j]+phiy[j+1])/2.0;

            dur  = ((u[i][j+1] - u[i][j])/dr) * phiyj;
            dvth = ((v[i+1][j] - v[i-1][j])/(2.0*dth)) * phixi;
            drhocr = (((rhohat[i][j+1]*chat[i][j+1]) - (rhohat[i][j]*chat[i][j]))/dr)*phiyj;
         
            c21 = 1.0/(0.5*rhohat[i][j]+0.5*rhohat[i][j+1])*drhocr;
            c23 = -2.0*(what[i][j] + what[i-1][j] + what[i][j+1] + what[i-1][j+1])/4.0/(r+dr*0.5);

            dampvx = dampvx4( i, j, v, minfo );
            dampvy = dampvy4( i, j, v, minfo );

            v1[i][j] = -(chat[i][j]+chat[i][j+1])/2.0*dur \
                       -1.0/(r+dr*0.5)*(what[i][j] + what[i-1][j] + what[i][j+1] + what[i-1][j+1])/4.0*dvth \
                       -c21*(u[i][j+1]+u[i][j])/2.0 \
                       -c23*(w[i][j] + w[i-1][j] + w[i][j+1] + w[i-1][j+1])/4.0 \
                       + dampvx + dampvy ;

            //dw update
            phixi=(phix[i]+phix[i+1])/2.0;
            phiyj=phiy[j];

            duth = ((u[i+1][j]-u[i][j])/dth)*phixi;
            dwth = ((w[i+1][j]-w[i-1][j])/(2.0*dth))*phixi;
            drhocth = ((rhohat[i+1][j]*chat[i+1][j] - rhohat[i][j]*chat[i][j])/dth)*phixi;
            dwhr = ((what[i][j+1] - what[i][j-1])/(2.0*dr))*phiyj;
            dwhth = ((what[i+1][j] - what[i-1][j])/(2.0*dth))*phixi;

            c31 = 1.0/r/(0.5*rhohat[i][j]+0.5*rhohat[i+1][j])*drhocth;
            c32 = what[i][j]/r + dwhr;
            c33 = 1.0/r*dwhth;

            dampwx = dampwx4( i, j, w, minfo );
            dampwy = dampwy4( i, j, w, minfo );
           
            w1[i][j] = -1.0/r*(0.5*chat[i][j]+0.5*chat[i+1][j])*duth \
                       -1.0/r*what[i][j]*dwth \
                       -c31*(0.5*u[i+1][j]+0.5*u[i][j]) \
                       -c32*(v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1])/4.0 \
                       -c33*w[i][j] \
                       + dampwx + dampwy;
       }
   }

   //source update
   if(RK4stage == 2 || RK4stage == 3) {
   for(k=0; k<minfo->src_num; k++) {
      i=minfo->srcs[k].xi;
      j=minfo->srcs[k].yi;
      u1[i][j] += minfo->srcs[k].midq->elements[tk]*chat[i][j];
      //printf("***********ac2dr_dqdt, source test, %d, %d, %d, %e\n", tk, i, j, minfo->srcs[k].midq->elements[tk]);
   }
   }

   if(RK4stage == 1 || RK4stage == 4) {
   for(k=0; k<minfo->src_num; k++) {
      i=minfo->srcs[k].xi;
      j=minfo->srcs[k].yi;
      u1[i][j] += minfo->srcs[k].q->elements[tk]*chat[i][j];
      //printf("***********ac2dr_dqdt, source test, %d\n", tk);
   }
   }

//    printf("***********ac2dr_dqdt.c is working\n");
    return;
}


