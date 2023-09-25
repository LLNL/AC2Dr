#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void mpi_printf ( MPI_Comm comm, char *msg, int nc, int mrank, int nRank ) {
   int i;
   char buf0[200];
   MPI_Status status;

   if(comm != MPI_COMM_NULL) {
      if(mrank==0) {
         printf("Process %d/%d: %s\n", mrank, nRank, msg);
      }

   for(i=1; i<=nRank; i++) {
      if(i==mrank) {
         MPI_Send( msg, nc, MPI_CHAR, 0, 1, comm);
      }

      if(mrank==0) {
         MPI_Recv( buf0, nc, MPI_CHAR, i, 1, comm, &status);
         printf("Process %d/%d: %s\n", i, nRank, buf0);
      }
      MPI_Barrier(comm);
   }

   }

   return;
}

void wdbin_mpi_array2_reduced( char *fn, double **m, int nx, int ny, int nx_global, int ny_global, \
     double dx, double dy, int mrank, int nRank, MPI_Comm comm ) {

   int nx1, ny1;
   int nx_reduced=0, ny_reduced=0;
   double dx0, dy0;
   double **m1;
   FILE *ptr_file;
   int i, j, indx, si, cRank;
   int tag1=1, tag2=2, tag3=3;
   MPI_Status status;

   dx0 = dx*10.0;
   dy0 = dy*10.0;
   nx1 = nx; ny1 = ny;
//   printf("Process %d/%d: nx1=%d, ny1=%d Set\n", mrank, nRank, nx1, ny1);   
   for(i=ghostL; i<(nx_global-ghostL); i += 10) {
      nx_reduced += 1;
   }

   for(j=ghostL; j<(ny_global-ghostL); j += 10) {
      ny_reduced += 1;
   }
   
   if(mrank == 0) {
      m1 = (double **) malloc ( 2 * nx * sizeof(double *) );
      m1[0] = (double *) malloc ( 2 * nx * ny * sizeof(double) );
      for(i=1; i<(2*nx); i++) {
         m1[i] = m1[0] + ny*i;
      }

      ptr_file = fopen(fn, "wb");
      fwrite(&nx_reduced, sizeof(int), 1, ptr_file);
      fwrite(&ny_reduced, sizeof(int), 1, ptr_file);
      fwrite(&dx0, sizeof(double), 1, ptr_file);
      fwrite(&dy0, sizeof(double), 1, ptr_file);
//      for(indx=0; indx < ((nx1-3) * ny1); indx++) 
//      {
//         fwrite(&m[0][indx], sizeof(double), 1, ptr_file);
//      }
      for(i=ghostL; i<(nx1-ghostL); i+=10) {
         for(j=ghostL; j<(ny1-ghostL); j+=10) {
            fwrite(&m[i][j], sizeof(double), 1, ptr_file);
         }
      }
   }
   
   MPI_Barrier(comm);

   if((mrank <= nRank) & (mrank > 0)) {
       MPI_Send( &nx1, 1, MPI_INT, 0, tag1, comm);
       MPI_Send( &ny1, 1, MPI_INT, 0, tag2, comm);
//       printf("MPI_test.. nx1=%d, ny1=%d\n",nx1,ny1); 
       MPI_Send( m[0], nx1*ny1, MPI_DOUBLE, 0, tag3, comm);
//       printf("Process %d/%d: nx1=%d, ny1=%d Sent\n", mrank, nRank, nx1, ny1);   
   }

   if(mrank == 0) {
      for(cRank=1; cRank <= nRank; cRank++) {
         MPI_Recv( &nx1, 1, MPI_INT, cRank, tag1, comm, &status);
         MPI_Recv( &ny1, 1, MPI_INT, cRank, tag2, comm, &status);
         MPI_Recv( m1[0], nx1*ny1, MPI_DOUBLE, cRank, tag3, comm, &status);
//         printf("Process %d/%d: nx1=%d, ny1=%d Received from %d\n", mrank, nRank, nx1, ny1, cRank);   
         if(cRank==nRank) {
            //for(indx=(ny1*3); indx < (nx1 * ny1); indx++) 
            //{
            //   fwrite(&m1[0][indx], sizeof(double), 1, ptr_file);
            //}
            indx = cRank * (nx - ghostL) - ghostL;
            if((indx%10==0)) { si = ghostL; } else { si = ghostL + (indx/10+1)*10 - indx; };
            for(i=si; i<(nx1-ghostL); i+=10) {
         //printf("wdbin2: Process %d/%d: nx1=%d, ny1=%d, indx=%d, si=%d, i=%d, Received from %d\n", 
         //      mrank, nRank, nx1, ny1, indx, si, i, cRank);   
               for(j=ghostL; j<(ny1-ghostL); j+=10) {
                  fwrite(&m1[i][j], sizeof(double), 1, ptr_file);
               }
            }
         } else {
            //for(indx=(ny1*3); indx < ((nx1-3) * ny1); indx++) 
            //{
            //   fwrite(&m1[0][indx], sizeof(double), 1, ptr_file);
            //}
            indx = cRank * (nx - ghostL) - ghostL;
            if((indx%10==0)) { si = ghostL; } else { si = ghostL + (indx/10+1)*10 - indx; };
            for(i=si; i<(nx1-ghostL); i+=10) {
               for(j=ghostL; j<(ny1-ghostL); j+=10) {
                  fwrite(&m1[i][j], sizeof(double), 1, ptr_file);
               }
            }
         }
      }
   }

   MPI_Barrier(comm);
   if(mrank == 0) {
      fclose(ptr_file);
      free_array2( m1 );
   }

   return;
}


void wdbin_mpi_array2( char *fn, double **m, int nx, int ny, int nx_global, int ny_global, double dx, double dy, \
   int mrank, int nRank, MPI_Comm comm ) {

   int nx1, ny1;
   double dx0, dy0;
   double **m1;
   FILE *ptr_file;
   int i, indx, cRank;
   int tag1=1, tag2=2, tag3=3;
   MPI_Status status;

   nx1 = nx; ny1 = ny;
   dx0 = dx; dy0 = dy;
//   printf("Process %d/%d: nx1=%d, ny1=%d Set\n", mrank, nRank, nx1, ny1);   
   
   if(mrank == 0) {
      m1 = (double **) malloc ( 2 * nx * sizeof(double *) );
      m1[0] = (double *) malloc ( 2 * nx * ny * sizeof(double) );
      for(i=1; i<(2*nx); i++) {
         m1[i] = m1[0] + ny*i;
      }

      ptr_file = fopen(fn, "wb");
      fwrite(&nx_global, sizeof(int), 1, ptr_file);
      fwrite(&ny_global, sizeof(int), 1, ptr_file);
      fwrite(&dx0, sizeof(double), 1, ptr_file);
      fwrite(&dy0, sizeof(double), 1, ptr_file);
      for(indx=0; indx < ((nx1-ghostL) * ny1); indx++) 
      {
         fwrite(&m[0][indx], sizeof(double), 1, ptr_file);
      }
   }
   
   MPI_Barrier(comm);

   if((mrank <= nRank) & (mrank > 0)) {
       MPI_Send( &nx1, 1, MPI_INT, 0, tag1, comm);
       MPI_Send( &ny1, 1, MPI_INT, 0, tag2, comm);
//       printf("MPI_test.. nx1=%d, ny1=%d\n",nx1,ny1); 
       MPI_Send( m[0], nx1*ny1, MPI_DOUBLE, 0, tag3, comm);
//       printf("Process %d/%d: nx1=%d, ny1=%d Sent\n", mrank, nRank, nx1, ny1);   
   }

   if(mrank == 0) {
      for(cRank=1; cRank <= nRank; cRank++) {
         MPI_Recv( &nx1, 1, MPI_INT, cRank, tag1, comm, &status);
         MPI_Recv( &ny1, 1, MPI_INT, cRank, tag2, comm, &status);
         MPI_Recv( m1[0], nx1*ny1, MPI_DOUBLE, cRank, tag3, comm, &status);
//         printf("Process %d/%d: nx1=%d, ny1=%d Received from %d\n", mrank, nRank, nx1, ny1, cRank);   
         if(cRank==nRank) {
            for(indx=(ny1*ghostL); indx < (nx1 * ny1); indx++) 
            {
               fwrite(&m1[0][indx], sizeof(double), 1, ptr_file);
            }
         } else {
            for(indx=(ny1*ghostL); indx < ((nx1-ghostL) * ny1); indx++) 
            {
               fwrite(&m1[0][indx], sizeof(double), 1, ptr_file);
            }
         }
      }
   }

   MPI_Barrier(comm);
   if(mrank == 0) {
      fclose(ptr_file);
      free_array2( m1 );
   }

   return;
}


void wdbin_mpi_array1( char *fn, double *m, int nx, int nx_global, int nb, \
   int mrank, int nRank, MPI_Comm comm ) {

   int nx1, nb1;
   double *m1;
   FILE *ptr_file;
   int i, indx, cRank;
   int tag1=1, tag2=2, tag3=3;
   MPI_Status status;

   nx1 = nx;
   nb1 = nb;
//   printf("Process %d/%d: nx1=%d, ny1=%d Set\n", mrank, nRank, nx1, ny1);   
   
   if(mrank == 0) {
      m1 = (double *) malloc ( 2 * nx * sizeof(double) );
      for(i=1; i<(2*nx); i++) {
         m1[i] = 0.0;
      }

      ptr_file = fopen(fn, "wb");
      fwrite(&nx_global, sizeof(int), 1, ptr_file);
      for(indx=0; indx < (nx1-nb); indx++) 
      {
         fwrite(&m[indx], sizeof(double), 1, ptr_file);
      }
   }
   
   MPI_Barrier(comm);

   if((mrank <= nRank) & (mrank > 0)) {
       MPI_Send( &nx1, 1, MPI_INT, 0, tag1, comm);
       MPI_Send( &nb1, 1, MPI_INT, 0, tag2, comm);
//       printf("MPI_test.. nx1=%d, ny1=%d\n",nx1,ny1); 
       MPI_Send( m, nx1, MPI_DOUBLE, 0, tag3, comm);
//       printf("Process %d/%d: nx1=%d, ny1=%d Sent\n", mrank, nRank, nx1, ny1);   
   }

   if(mrank == 0) {
      for(cRank=1; cRank <= nRank; cRank++) {
         MPI_Recv( &nx1, 1, MPI_INT, cRank, tag1, comm, &status);
         MPI_Recv( &nb1, 1, MPI_INT, cRank, tag2, comm, &status);
         MPI_Recv( m1, nx1, MPI_DOUBLE, cRank, tag3, comm, &status);
//         printf("Process %d/%d: nx1=%d, ny1=%d Received from %d\n", mrank, nRank, nx1, ny1, cRank);   
         if(cRank==nRank) {
            for(indx=nb; indx < nx1; indx++) 
            {
               fwrite(&m1[indx], sizeof(double), 1, ptr_file);
            }
         } else {
            for(indx=nb; indx < (nx1-nb); indx++) 
            {
               fwrite(&m1[indx], sizeof(double), 1, ptr_file);
            }
         }
      }
   }

   MPI_Barrier(comm);
   if(mrank == 0) {
      fclose(ptr_file);
      free((double *) m1);
   }

   return;
}


