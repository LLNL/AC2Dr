#include <stdio.h>
#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void decompose_nproc (fdmesh *minfo, innout *finfo, int mrank, int totalp, int *nRank, \
                      int *xi_start, int *xi_end, int *nx_local)
{

   int i;
   int nx_tmp;
   int ghost0; //ghost layer for the leftmost block
//   char buf[200];
//   int xi_start[100], xi_end[100], nx_local[100];
//   double *xx, *yy, *c_prof, *rho_prof, *w_prof;
//   char buf[500];

//   sprintf(buf, "totalp = %d", totalp);
//   mpi_printf(comm, buf, mrank, totalp);
   ghost0 = ghostL;
   nx_tmp = ceil((double) (minfo->nx_global-2*ghost0) / (double) totalp);
   if(nx_tmp<(2*ghostL)) nx_tmp=2*ghostL;

//   printf("Process %d/%d: nx_tmp=%d...\n",mrank, totalp,nx_tmp);
   i=0;

   xi_start[i] = nx_tmp*i;
   xi_end[i] = (nx_tmp*(i+1))-1;
   nx_local[i] = xi_end[i] - xi_start[i] + 1;

//   printf("Process %d/%d: nx_global = %ld...\n",mrank, totalp, minfo->nth_global);

   while(xi_end[i] < ((minfo->nx_global-2*ghost0) -1))
   {
      i=i+1;
      xi_start[i] = nx_tmp*i;
      xi_end[i] = (nx_tmp*(i+1))-1;
      nx_local[i] = xi_end[i] - xi_start[i] + 1;
   }

   xi_end[i] = (minfo->nx_global-2*ghost0) - 1;
   nx_local[i] = xi_end[i] - xi_start[i] + 1;

   if(nx_local[i] >= (ghostL*2)) {
      *nRank = i;
   } else {
      *nRank = i-1;
      xi_end[i-1] = xi_end[i];
      nx_local[i-1] = xi_end[i-1] - xi_start[i-1] + 1;
   }

   return;
}



void decompose_domain (fdmesh *minfo, innout *finfo, int mrank, int totalp, int nRank, \
                       int *xi_start, int *xi_end, int *nx_local, MPI_Comm comm)
{

   int i, j;
   int ghost0; //ghost layer for the leftmost block
   char buf[200];

   minfo->nRank = nRank;
   ghost0 = ghostL;
//   printf("Process %d/%d: nRank = %d\n", mrank, totalp, *nRank);
//   printf("Process %d/%d: nx_local = %d...\n",mrank, totalp, nx_local[mrank]);

   //ghost layer adjustment
   minfo->xi_start = xi_start[mrank] + ghost0;
   minfo->xi_end = xi_end[mrank] + ghost0;

   minfo->xi_start = minfo->xi_start - ghost0;
   minfo->xi_end = minfo->xi_end + ghost0;

   minfo->nth = nx_local[mrank]+ghost0+ghost0;
   minfo->nx = nx_local[mrank]+ghost0+ghost0;


//    printf("Process %d/%d: xi_start = %d, xi_end = %d, nx=%d\n",mrank, *nRank, minfo->xi_start, minfo->xi_end, minfo->nx);
    minfo->xx = (double *) malloc ( minfo->nth * sizeof(double));
    minfo->xx_deg = (double *) malloc ( minfo->nth * sizeof(double));
    minfo->yy = (double *) malloc ( minfo->nr * sizeof(double));
    minfo->c_prof = (double *) malloc ( minfo->nr * sizeof(double));
    minfo->rho_prof = (double *) malloc ( minfo->nr * sizeof(double));
    minfo->w_prof = (double *) malloc ( minfo->nr * sizeof(double));

//   MPI_Barrier(comm);
//   printf("Process %d/%d: decompose test...\n",mrank, totalp);
//   MPI_Barrier(comm);

    for(i=0;i<minfo->nth;i++) {
      minfo->xx[i] = (i+minfo->xi_start)*minfo->dth;
      minfo->xx_deg[i] = (i+minfo->xi_start)*minfo->dth*180.0/M_PI;
    }

    for(j=0;j<minfo->nr;j++) {
      minfo->yy[j] = j*minfo->dr;
    }

    minfo->u = (double **) malloc ( minfo->nth * sizeof(double *) );
    minfo->v = (double **) malloc ( minfo->nth * sizeof(double *) );
    minfo->w = (double **) malloc ( minfo->nth * sizeof(double *) );

    minfo->u[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
    minfo->v[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
    minfo->w[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  

    for(i=1;i<minfo->nth;i++) {
        minfo->u[i] = minfo->u[0] + minfo->nr * i;
        minfo->v[i] = minfo->v[0] + minfo->nr * i;
        minfo->w[i] = minfo->w[0] + minfo->nr * i;
    }

    minfo->rhohat = (double **) malloc ( minfo->nth * sizeof(double *) );
    minfo->chat = (double **) malloc ( minfo->nth * sizeof(double *) );
    minfo->what = (double **) malloc ( minfo->nth * sizeof(double *) );

    minfo->rhohat[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
    minfo->chat[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
    minfo->what[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  

    for(i=1;i<minfo->nth;i++) {
        minfo->rhohat[i] = minfo->rhohat[0] + minfo->nr * i;
        minfo->chat[i] = minfo->chat[0] + minfo->nr * i;
        minfo->what[i] = minfo->what[0] + minfo->nr * i;
    }

    if(pmaximgs==1) {
//    minfo->pmax = (double **) malloc ( minfo->nth * sizeof(double *) );
//    minfo->pmax_t = (double **) malloc ( minfo->nth * sizeof(double *) );
    minfo->pmaxabs = (double **) malloc ( minfo->nth * sizeof(double *) );
    minfo->pmaxabs_t = (double **) malloc ( minfo->nth * sizeof(double *) );
//    minfo->pmax[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
//    minfo->pmax_t[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
    minfo->pmaxabs[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  
    minfo->pmaxabs_t[0] = (double *) malloc ( minfo->nth * minfo->nr * sizeof(double) );  

    for(i=1;i<minfo->nth;i++) {
 //       minfo->pmax[i] = minfo->pmax[0] + minfo->nr * i;
//        minfo->pmax_t[i] = minfo->pmax_t[0] + minfo->nr * i;
        minfo->pmaxabs[i] = minfo->pmaxabs[0] + minfo->nr * i;
        minfo->pmaxabs_t[i] = minfo->pmaxabs_t[0] + minfo->nr * i;
    }

    }
   
   sprintf(buf,"xi_start = %d, xi_end = %d, nx=%d", minfo->xi_start, minfo->xi_end, minfo->nx);
   MPI_Barrier(comm);
   //mpi_printf(comm, buf, 200, mrank, nRank);
   print_message_mpi(minfo->logf, comm, buf, 200, mrank, nRank, 0);
   return;
}



void shift_domain ( double **u, int nx, int ny, int byN, int mrank, int nRank, MPI_Comm comm2) {

    int i, j;

   MPI_Barrier(comm2);
   //shift inner domain
    for(j=(ny-byN-1); j>=0; j--) {
       for(i=(nx-byN-1); i>=0; i--) {
         u[i+byN][j+byN] = u[i][j];
       }
    }

/*    for(j=(ny-4); j>=0; j--) {
       for(i=0; i<(nx-3); i++) {
         u[i+3][j+3] = u[i][j];
       }
    }
*/

   MPI_Barrier(comm2);
   return;

}


//void updateBoundary ( fdmesh *minfo, int byN, int mrank, int nRank, MPI_Comm comm2) {
void updateBoundary ( double **u, int nx, int ny, int byN, int mrank, int nRank, MPI_Comm comm2) {

    int xtag1=1;
    MPI_Status status;

   //update boundary

   if(mrank == 0) {
   //send boundary to left and receive from right
    MPI_Sendrecv( &u[ghostL][0], ny*ghostL, MPI_DOUBLE, MPI_PROC_NULL, xtag1,
                 &u[nx-ghostL][0], ny*ghostL, MPI_DOUBLE, mrank+1, xtag1,
                  comm2, &status);
   //send boundary to right and receive from left
    MPI_Sendrecv( &u[nx-2*ghostL][0], ny*ghostL, MPI_DOUBLE, mrank+1, xtag1,
                 &u[0][0], ny*ghostL, MPI_DOUBLE, MPI_PROC_NULL, xtag1,
                  comm2, &status);
   }

   if((mrank > 0) & (mrank < nRank)) {
   //send boundary to left and receive from right
    MPI_Sendrecv( &u[ghostL][0], ny*ghostL, MPI_DOUBLE, mrank-1, xtag1,
                 &u[nx-ghostL][0], ny*ghostL, MPI_DOUBLE, mrank+1, xtag1,
                  comm2, &status);
   //send boundary to right and receive from left
    MPI_Sendrecv( &u[nx-2*ghostL][0], ny*ghostL, MPI_DOUBLE, mrank+1, xtag1,
                 &u[0][0], ny*ghostL, MPI_DOUBLE, mrank-1, xtag1,
                  comm2, &status);
   }

   if(mrank == nRank) {
   //send boundary to left and receive from right
    MPI_Sendrecv( &u[ghostL][0], ny*ghostL, MPI_DOUBLE, mrank-1, xtag1,
                 &u[nx-ghostL][0], ny*ghostL, MPI_DOUBLE, MPI_PROC_NULL, xtag1,
                  comm2, &status);
   //send boundary to right and receive from left
    MPI_Sendrecv( &u[nx-2*ghostL][0], ny*ghostL, MPI_DOUBLE, MPI_PROC_NULL, xtag1,
                 &u[0][0], ny*ghostL, MPI_DOUBLE, mrank-1, xtag1,
                  comm2, &status);
   }
   return;
}



void rigidbound_material ( double **chat, double **rhohat, double **what, int nx, int ny, \
                           int mrank, int nRank, MPI_Comm comm2) {

   int i, j, k;
   

   MPI_Barrier(comm2);

//top
    for(i=0; i<nx; i++) {
      for(k=0; k<ghostL; k++) { 
         chat[i][ny-1-ghostL+(k+1)] = chat[i][ny-1-ghostL-(k+1)];
         rhohat[i][ny-1-ghostL+(k+1)] = rhohat[i][ny-1-ghostL-(k+1)];
         what[i][ny-1-ghostL+(k+1)] = what[i][ny-1-ghostL-(k+1)];
      }
    }

//bottom
    for(i=0; i<nx; i++) {
      for(k=0; k<ghostL; k++) {
         chat[i][ghostL-(k+1)] = chat[i][ghostL+(k+1)];
         rhohat[i][ghostL-(k+1)] = rhohat[i][ghostL+(k+1)];
         what[i][ghostL-(k+1)] = what[i][ghostL+(k+1)];
      }
    }

//left
if(mrank == 0) {
    for(j=0; j<ny; j++) {
         what[ghostL][j] = 0.0;
      for(k=0; k<ghostL; k++) {
         chat[ghostL-(k+1)][j] = chat[ghostL+(k+1)][j];
         rhohat[ghostL-(k+1)][j] = rhohat[ghostL+(k+1)][j];
         what[ghostL-(k+1)][j] = -what[ghostL+(k+1)][j];
      }
    }
}

//right
if(mrank==nRank) {
   for(j=0; j<ny; j++) {
      for(k=0; k<ghostL; k++) {
         chat[nx-1-ghostL+(k+1)][j] = chat[nx-1-ghostL-(k+1)][j];
         rhohat[nx-1-ghostL+(k+1)][j] = rhohat[nx-1-ghostL-(k+1)][j];
         what[nx-1-ghostL+(k+1)][j] = -what[nx-1-ghostL-(k+1)][j];
      }
   }
}
   MPI_Barrier(comm2);
   return;

}



void rigidbound_var ( double **u, double **v, double **w, int nx, int ny, \
                      int mrank, int nRank, MPI_Comm comm2) {

   int i, j, k;
   
   MPI_Barrier(comm2);

//top
    for(i=0; i<nx; i++) {
      for(k=0; k<ghostL; k++) { 
         u[i][ny-1-ghostL+(k+1)] = u[i][ny-1-ghostL-(k+1)];
         v[i][ny-1-ghostL+(k+1)] = -v[i][ny-1-ghostL-(k+1)];
         w[i][ny-1-ghostL+(k+1)] = w[i][ny-1-ghostL-(k+1)];
      }
    }

//bottom
    for(i=0; i<nx; i++) {
      for(k=0; k<ghostL; k++) {
         u[i][ghostL-(k+1)] = u[i][ghostL+(k+1)];
         v[i][ghostL-(k+1)] = -v[i][ghostL+(k+1)];
         w[i][ghostL-(k+1)] = w[i][ghostL+(k+1)];
      }
    }

//left
if(mrank == 0) {
    for(j=0; j<ny; j++) {
         w[ghostL][j] = 0.0;
      for(k=0; k<ghostL; k++) {
         u[ghostL-(k+1)][j] = u[ghostL+(k+1)][j];
         v[ghostL-(k+1)][j] = v[ghostL+(k+1)][j];
         w[ghostL-(k+1)][j] = -w[ghostL+(k+1)][j];
      }
    }
}

//right
if(mrank==nRank) {
   for(j=0; j<ny; j++) {
      for(k=0; k<ghostL; k++) {
         u[nx-1-ghostL+(k+1)][j] = u[nx-1-ghostL-(k+1)][j];
         v[nx-1-ghostL+(k+1)][j] = v[nx-1-ghostL-(k+1)][j];
         w[nx-1-ghostL+(k+1)][j] = -w[nx-1-ghostL-(k+1)][j];
      }
   }
}

   MPI_Barrier(comm2);
   return;

}

