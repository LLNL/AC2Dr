#include "ac2dr_aux.h"
#include "ac2dr_config.h"

/*
void initialize_u( double **s, double **u, double **w, double **r, fdmesh *minfo ) {
    int i, j, nx, ny;
    double dh;
    double x,y;

    nx = minfo->nx;
    ny = minfo->ny;
    dh = minfo->dh;

    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            x=i * dh;
            y=j * dh;
            s[i][j] = func_s( x, y, 0.0 );
            u[i][j] = func_u( x, y, 0.0 );
            w[i][j] = func_w( x, y, 0.0 );
            r[i][j] = func_r( x, y, 0.0 );
        }
    }
}
*/

void initialize_u( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
    int i, j, nx, ny;
    double dr, dth;
    double th, r, ths, rs;
    double Amp, w0, R0, RR;
    double chat0;
    int thsi, rsj;
//    double c0;
//  if a process does not have a source, minfo->srcs[0] is not defined.  
    Amp = minfo->srcs[0].p0;
    w0 = minfo->srcs[0].freq;

//    MPI_Barrier(comm2);
//    MPI_Bcast(&Amp, 1, MPI_DOUBLE, 0, comm2);
//    MPI_Bcast(&w0, 1, MPI_DOUBLE, 0, comm2);
//    MPI_Barrier(comm2);

    nx = minfo->nth;
    ny = minfo->nr;
    dr = minfo->dr;
    dth = minfo->dth;
    R0 = minfo->R;

    //ths = 0.25*M_PI/180;
    ths = (minfo->thmin_global+minfo->srcs[0].x);
    thsi = minfo->srcs[0].xi_global;
    //rs = R0+minfo->srcs[0].y;
    //ths = 0.5*M_PI/180;
    rs = R0+minfo->srcs[0].y;
    rsj = minfo->srcs[0].yi;

//    MPI_Barrier(comm2);
//    MPI_Bcast(&ths, 1, MPI_DOUBLE, 0, comm2);
//    MPI_Bcast(&rs, 1, MPI_DOUBLE, 0, comm2);
//    MPI_Bcast(&rsj, 1, MPI_INT, 0, comm2);
//    MPI_Barrier(comm2);

    //rs = 10000.0 + R0;
    //printf("w0 = %f\n", w0);
    //printf("Amp = %f\n", Amp);
//    c0 = minfo->chat[ghostL][rsj];
    chat0 = minfo->chat[thsi][rsj];
//    MPI_Barrier(comm2);
//    MPI_Bcast(&chat0, 1, MPI_DOUBLE, 0, comm2);
//    MPI_Barrier(comm2);
    printf("Process %d/%d: ths = %f, r = %f\n", mrank, nRank, ths, rs);
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            //th=(i+minfo->xi_start) * dth + minfo->thmin_global;
            //r=j * dr + R0;
            //if(ghostL!=0) {
            th=(i+minfo->xi_start-ghostL) * dth + minfo->thmin_global;
            r=(j-ghostL) * dr + R0;
            //}
            RR = sqrt( pow(r*sin(th) -rs*sin(ths),2) + pow(r*cos(th) - rs*cos(ths), 2) ); 
            minfo->u[i][j] = point_src( Amp, RR, w0, chat0 );
            //minfo->u[i][j] = point_src( Amp, RR, w0, minfo->chat[ghostL][rsj] ) / minfo->rhohat[i][j] / minfo->chat[i][j];
            //u[i][j] = cos(10.0*th/(dth*(nx-1))*2.0*M_PI);
            //minfo->u[i][j] = 0.0;
        }
    }
}


double point_src( double A, double R, double f, double c ) {
   double ps;
   double t0, t;
   double w;

//C6smoothbump
  //printf("test A=%f\n", A);
  w = f /3.0/ c;
  t0=-1.0/w/2.0;
  t = R;
  if(t>=t0 && t<=(t0+1.0/w)) {
  ps = A*51480.0*pow(w, 7)*pow(t-t0, 7)*pow(1.0-w*(t-t0), 7)/3.14209;
  } else {
   ps =0.0;
  }

//Hedlin's function
//   if(R<=(4.0*w)) {
//      ps = A*exp(-pow(R/w,2))*cos(R*M_PI/(8.0*w));
//   } else {
//      ps = 0.0;
//   }


/*cos function
   if(R<=(6.0*w)) {
      //ps = A*exp(-pow(R/w,2))*cos(R*M_PI/(8.0*w));
      //ps = A*exp(-pow(R/(2.0*w),2))*cos(R*M_PI/(8.0*w));
      ps = A*cos(R/w*2*M_PI)*cos(R*M_PI/(6.0*w));
   } else {
      ps = 0.0;
   }
*/
   return(ps);
}



void initialize_v( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
    int i, j, nx, ny;

    nx = minfo->nx;
    ny = minfo->ny;
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            minfo->v[i][j] = 0.0;
            //u[i][j] = cos(10.0*th/(dth*(nx-1))*2.0*M_PI);
        }
    }

    return;
}

void initialize_u_zero( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
    int i, j, nx, ny;

    nx = minfo->nx;
    ny = minfo->ny;
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            minfo->u[i][j] = 0.0;
        }
    }

    if(pmaximgs==1) {
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            //minfo->pmax[i][j] = 0.0;
            minfo->pmaxabs[i][j] = 0.0;
        }
    }
    }
    return;
}


void initialize_w( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
    int i, j, nx, ny;

    nx = minfo->nx;
    ny = minfo->ny;
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            minfo->w[i][j] = 0.0;
            //u[i][j] = cos(10.0*th/(dth*(nx-1))*2.0*M_PI);
        }
    }

    return;
}

void init_src( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
   int i, j;

   for(i=0; i<minfo->src_num; i++) {
      for(j=0; j<minfo->nt; j++) {
         minfo->srcs[i].q->elements[j] = \
            minfo->srcs[i].p0*c6smoothbump(minfo->ts->elements[j], minfo->srcs[i].t0, minfo->srcs[i].freq/3.0);
         minfo->srcs[i].midq->elements[j] = \
            minfo->srcs[i].p0*c6smoothbump(minfo->ts->elements[j]+0.5*minfo->dt, minfo->srcs[i].t0, minfo->srcs[i].freq/3.0);
      }
   }
}

