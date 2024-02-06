#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void setup_sg (fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
    int i, gi;
    int nx, ny, nx_g, nt;
    double xi;
    int zstep;

    nx = minfo->nx;    
    ny = minfo->ny;    
    nt = minfo->nt;
    nx_g = minfo->nx_global;
    minfo->sgrid_length = 40;
    //minfo->sgrid_dc=0.07;
    //minfo->sgrid_dc=0.01;
    minfo->sgrid_dc=sg_dc;
    minfo->sgrid_EL = 1.0/pow(10.0,4);
    minfo->alpha = 1./3.;
    minfo->psix = (double *) malloc( nx * sizeof(double) );
    minfo->psiy = (double *) malloc( ny * sizeof(double) );
    minfo->phix = (double *) malloc( nx * sizeof(double) );
    minfo->phiy = (double *) malloc( ny * sizeof(double) );
    minfo->sigx = (double *) malloc( nx * sizeof(double) );
    minfo->sigy = (double *) malloc( ny * sizeof(double) );
    minfo->etax = (double *) malloc( nx * sizeof(double) );
    minfo->etay = (double *) malloc( ny * sizeof(double) );

    minfo->at = (double *) malloc( nt * sizeof(double) );
    minfo->ax = dcoef/minfo->dt;

    zstep = floor(damp_rise/minfo->dt);

    for(i=0;i<nt;i++) {
      minfo->at[i] = 1.0;
      if((i>=0) & (i<=damp_rise)) {
         minfo->at[i] = pow(sin((double)(i)/zstep*M_PI/2.0), 2);
      }
      minfo->at[i] = 1.0;
    }

    for(i=0;i<nx;i++) {
        gi=i+minfo->xi_start;

        if(gi < (nx_g-minfo->sgrid_length-ghostL)) {
            minfo->psix[i] = 0;
            minfo->etax[i] = 1;
        }
 
        if((gi >= (nx_g-minfo->sgrid_length-ghostL)) & (gi<(nx_g-ghostL)) ) {
            xi = (double)(gi-(nx_g-minfo->sgrid_length-ghostL))/minfo->sgrid_length;
            minfo->psix[i] = pow(xi, 6)*(462.-1980.*xi+3465.*pow(xi,2)-3080.*pow(xi,3) + \
                         1386.*pow(xi,4) - 252.*pow(xi,5));
            minfo->etax[i] = minfo->alpha + (1.-minfo->alpha) * (1-xi);
        }

        if((gi>=nx_g-ghostL) & (gi<nx_g)) {
            minfo->psix[i] = 1;
            minfo->etax[i] = minfo->alpha;
        }
    }


    for(i=0;i<(ny-minfo->sgrid_length-ghostL);i++) {
        minfo->psiy[i] = 0;
        minfo->etay[i] = 1;
//      if(mrank==0) {
//       printf("i1=%d\n",i);
//      }
    }
    for(i=(ny-minfo->sgrid_length-ghostL);i<(ny-ghostL);i++) {
        xi = (double)(i-(ny-minfo->sgrid_length-ghostL))/minfo->sgrid_length;
        minfo->psiy[i] = pow(xi, 6)*(462.-1980.*xi+3465.*pow(xi,2)-3080.*pow(xi,3) + \
                         1386.*pow(xi,4) - 252.*pow(xi,5));
        minfo->etay[i] = minfo->alpha + (1.-minfo->alpha) * (1-xi);
//      if(mrank==0) {
//       printf("i2=%d, psiy=%f\n",i,minfo->psiy[i]);
//      }

    }
    for(i=ny-ghostL;i<ny;i++) {
        minfo->psiy[i] = 1;
        minfo->etay[i] = minfo->alpha;
//      if(mrank==0) {
//       printf("i3=%d\n",i);
//      }
    }

    for(i=0;i<nx;i++) {
        minfo->phix[i] = 1. - (1.-minfo->sgrid_EL)*minfo->psix[i];
        minfo->sigx[i] = minfo->psix[i]/minfo->phix[i];
        //minfo->phix[i] = 1.0;
    }

    for(i=0;i<ny;i++) {
        minfo->phiy[i] = 1. - (1.-minfo->sgrid_EL)*minfo->psiy[i];
        minfo->sigy[i] = minfo->psiy[i]/minfo->phiy[i];
        //minfo->phiy[i] = 1.0;
    }

//   printf("Process %d/%d: sgrid_length = %d, i=%d, ghostL=%d \n", mrank, nRank, minfo->sgrid_length,i,ghostL);

   return;
}



double dampx4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dth,4) / minfo->dt;
    phi = minfo->phix;
    sigma = minfo->sigx;
    val = -eps * pow(-1, 4) * phi[i] * minfo->etay[j] * \
          1./pow(minfo->dth,4) * (sigma[i+1]*u[i+2][j] - 2.*sigma[i+1]*u[i+1][j] + sigma[i+1]*u[i][j] \
                          -2.*sigma[i]*u[i+1][j] + 4.*sigma[i]*u[i][j] - 2.*sigma[i]*u[i-1][j] \
                          +sigma[i-1]*u[i][j] - 2.*sigma[i-1]*u[i-1][j] + sigma[i-1]*u[i-2][j]);
    return val;
}

double dampy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dr,4) / minfo->dt;
    phi = minfo->phiy;
    sigma = minfo->sigy;
    val = -eps * pow(-1, 4) * phi[j] * minfo->etax[i] * \
          1./pow(minfo->dr,4) * (sigma[j+1]*u[i][j+2] - 2.*sigma[j+1]*u[i][j+1] + sigma[j+1]*u[i][j] \
                          -2.*sigma[j]*u[i][j+1] + 4.*sigma[j]*u[i][j] - 2.*sigma[j]*u[i][j-1] \
                          +sigma[j-1]*u[i][j] - 2.*sigma[j-1]*u[i][j-1] + sigma[j-1]*u[i][j-2]);
    return val;
}

double dampux4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dth,4) / minfo->dt;
    phi = minfo->phix;
    sigma = minfo->sigx;
    val = -eps * pow(-1, 4) * phi[i] * minfo->etay[j] * \
          1./pow(minfo->dth,4) * (sigma[i+1]*u[i+2][j] - 2.*sigma[i+1]*u[i+1][j] + sigma[i+1]*u[i][j] \
                          -2.*sigma[i]*u[i+1][j] + 4.*sigma[i]*u[i][j] - 2.*sigma[i]*u[i-1][j] \
                          +sigma[i-1]*u[i][j] - 2.*sigma[i-1]*u[i-1][j] + sigma[i-1]*u[i-2][j]);
    return val;
}

double dampuy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dr,4) / minfo->dt;
    phi = minfo->phiy;
    sigma = minfo->sigy;
    val = -eps * pow(-1, 4) * phi[j] * minfo->etax[i] * \
          1./pow(minfo->dr,4) * (sigma[j+1]*u[i][j+2] - 2.*sigma[j+1]*u[i][j+1] + sigma[j+1]*u[i][j] \
                          -2.*sigma[j]*u[i][j+1] + 4.*sigma[j]*u[i][j] - 2.*sigma[j]*u[i][j-1] \
                          +sigma[j-1]*u[i][j] - 2.*sigma[j-1]*u[i][j-1] + sigma[j-1]*u[i][j-2]);
    return val;
}

double dampvx4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dth,4) / minfo->dt;
    phi = minfo->phix;
    sigma = minfo->sigx;
    val = -eps * pow(-1, 4) * phi[i] * (minfo->etay[j]+minfo->etay[j+1])/2.0 * \
          1./pow(minfo->dth,4) * (sigma[i+1]*u[i+2][j] - 2.*sigma[i+1]*u[i+1][j] + sigma[i+1]*u[i][j] \
                          -2.*sigma[i]*u[i+1][j] + 4.*sigma[i]*u[i][j] - 2.*sigma[i]*u[i-1][j] \
                          +sigma[i-1]*u[i][j] - 2.*sigma[i-1]*u[i-1][j] + sigma[i-1]*u[i-2][j]);
    return val;
}

double dampvy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps, sigp1, sig0, sigm1,phiyj;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dr,4) / minfo->dt;
    phi = minfo->phiy;
    phiyj=(phi[j]+phi[j+1])/2.0;
    sigma = minfo->sigy;
    sigp1 = (sigma[j+1]+sigma[j+2])/2.0;
    sig0 = (sigma[j]+sigma[j+1])/2.0;
    sigm1 = (sigma[j-1]+sigma[j])/2.0;
    val = -eps * pow(-1, 4) * phiyj * minfo->etax[i] * \
          1./pow(minfo->dr,4) * (sigp1*u[i][j+2] - 2.*sigp1*u[i][j+1] + sigp1*u[i][j] \
                          -2.*sig0*u[i][j+1] + 4.*sig0*u[i][j] - 2.*sig0*u[i][j-1] \
                          +sigm1*u[i][j] - 2.*sigm1*u[i][j-1] + sigm1*u[i][j-2]);
    return val;
}

double dampwx4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps, sigp1, sig0, sigm1, phixi;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dth,4) / minfo->dt;
    phi = minfo->phix;
    phixi=(phi[i]+phi[i+1])/2.0;
    sigma = minfo->sigx;
    sigp1 = (sigma[i+1]+sigma[i+2])/2.0;
    sig0 = (sigma[i]+sigma[i+1])/2.0;
    sigm1 = (sigma[i-1]+sigma[i])/2.0;
    val = -eps * pow(-1, 4) * phixi * minfo->etay[j] * \
          1./pow(minfo->dth,4) * (sigp1*u[i+2][j] - 2.*sigp1*u[i+1][j] + sigp1*u[i][j] \
                          -2.*sig0*u[i+1][j] + 4.*sig0*u[i][j] - 2.*sig0*u[i-1][j] \
                          +sigm1*u[i][j] - 2.*sigm1*u[i-1][j] + sigm1*u[i-2][j]);
    return val;
}

double dampwy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo ) {
    double val, eps;
    double *sigma, *phi;
    eps = minfo->sgrid_dc * pow(minfo->dr,4) / minfo->dt;
    phi = minfo->phiy;
    sigma = minfo->sigy;
    val = -eps * pow(-1, 4) * phi[j] * (minfo->etax[i]+minfo->etax[i+1])/2.0 * \
          1./pow(minfo->dr,4) * (sigma[j+1]*u[i][j+2] - 2.*sigma[j+1]*u[i][j+1] + sigma[j+1]*u[i][j] \
                          -2.*sigma[j]*u[i][j+1] + 4.*sigma[j]*u[i][j] - 2.*sigma[j]*u[i][j-1] \
                          +sigma[j-1]*u[i][j] - 2.*sigma[j-1]*u[i][j-1] + sigma[j-1]*u[i][j-2]);
    return val;
}


