#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void material_prop(fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm)
{
   int i, j;
   long yi;
   double temp;
   int c_type, rho_type, w_type, nx, ny;
   int fformat;
   int nx0, ny0;
   double *xx0, *yy0, **zz0, *vel0;

//homogeneous
   nx = minfo->nx;
   ny = minfo->ny;

   c_type = minfo->c_type;
   rho_type = minfo->rho_type;
   w_type = minfo->w_type;

   fformat = finfo->c_format;

   if(c_type==1) {
      for(i=0;i<nx;i++) {
         for(j=0;j<ny;j++) {
            minfo->chat[i][j]=minfo->c_val;
         }
      }
   }

   if(rho_type==1) {
      for(i=0;i<nx;i++) {
         for(j=0;j<ny;j++) {
            minfo->rhohat[i][j]=minfo->rho_val;
         }
      }
   }

   if(w_type==1) {
      for(i=0;i<nx;i++) {
         for(j=0;j<ny;j++) {
            minfo->what[i][j]=minfo->w_val;
         }
      }
   }


//1-d profile
   if(c_type==2) {
     fformat = finfo->c_format;

     if(fformat==1) {
         read_2col(finfo->c_fn, &yy0, &vel0, &ny0);
     }
     if(fformat==0) {
         readbin_2col(finfo->c_fn, &yy0, &vel0, &ny0);
     }
      for(i=0;i<ny;i++) {
         locate(yy0, ny0, minfo->yy[i], &yi);
         if(yi==-1) {
            polint( &yy0[yi+1], &vel0[yi+1], 2, minfo->yy[i], &minfo->c_prof[i], &temp );
/*
polint(float xa[], float ya[], int n, float x, float *y, float *dy)
Given arrays xa[1..n] and ya[1..n], and given a value x,
returns a value y
*/
         } else if(yi==(ny-1)) {
            polint( &yy0[yi-1], &vel0[yi-1], 2, minfo->yy[i], &minfo->c_prof[i], &temp );
         } else {
            polint( &yy0[yi], &vel0[yi], 2, minfo->yy[i], &minfo->c_prof[i], &temp );
         }
      }
   for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
         minfo->chat[i][j] = minfo->c_prof[j];
      }
   }

      //wdbin_array2("test_c.bin", minfo->chat, nx, ny);
      free( yy0 );
      free( vel0 );
   }


   if(rho_type==2) {
     fformat = finfo->rho_format;

     if(fformat==1) {
         read_2col(finfo->rho_fn, &yy0, &vel0, &ny0);
     }
     if(fformat==0) {
         readbin_2col(finfo->rho_fn, &yy0, &vel0, &ny0);
     }
      for(i=0;i<ny;i++) {
         locate(yy0, ny0, minfo->yy[i], &yi);
         if(yi==-1) {
            polint( &yy0[yi+1], &vel0[yi+1], 2, minfo->yy[i], &minfo->rho_prof[i], &temp );
/*
polint(float xa[], float ya[], int n, float x, float *y, float *dy)
Given arrays xa[1..n] and ya[1..n], and given a value x,
returns a value y
*/
         } else if(yi==(ny-1)) {
            polint( &yy0[yi-1], &vel0[yi-1], 2, minfo->yy[i], &minfo->rho_prof[i], &temp );
         } else {
            polint( &yy0[yi], &vel0[yi], 2, minfo->yy[i], &minfo->rho_prof[i], &temp );
         }
      }
   for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
         minfo->rhohat[i][j] = minfo->rho_prof[j];
      }
   }

      //wdbin_array2("test_rho.bin", minfo->rhohat, nx, ny);
      free( yy0 );
      free( vel0 );
   }

   if(w_type==2) {
     fformat = finfo->w_format;

     if(fformat==1) {
         read_2col(finfo->w_fn, &yy0, &vel0, &ny0);
     }
     if(fformat==0) {
         readbin_2col(finfo->w_fn, &yy0, &vel0, &ny0);
     }
      for(i=0;i<ny;i++) {
         locate(yy0, ny0, minfo->yy[i], &yi);
         if(yi==-1) {
            polint( &yy0[yi+1], &vel0[yi+1], 2, minfo->yy[i], &minfo->w_prof[i], &temp );
/*
polint(float xa[], float ya[], int n, float x, float *y, float *dy)
Given arrays xa[1..n] and ya[1..n], and given a value x,
returns a value y
*/
         } else if(yi==(ny-1)) {
            polint( &yy0[yi-1], &vel0[yi-1], 2, minfo->yy[i], &minfo->w_prof[i], &temp );
         } else {
            polint( &yy0[yi], &vel0[yi], 2, minfo->yy[i], &minfo->w_prof[i], &temp );
         }
      }
   for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
         minfo->what[i][j] = minfo->w_prof[j];
      }
   }

      //wdbin_array2("test_wind.bin", minfo->what, nx, ny);
      free( yy0 );
      free( vel0 );
   }



//2-D weather model
   if(c_type==3) {
      if(fformat==0) {
         //readbin_2dfile(outp, fn, nx, ny, xx, yy);
         readbin_2dfile(finfo->c_fn, &nx0, &ny0, &xx0, &yy0, &zz0);

         //wdbin_array2("test_c.bin", zz0, nx0, ny0);
      } else {
         //read_2dfile(outp, fn, nx, ny, xx, yy);
         read_2dfile(finfo->c_fn, &nx0, &ny0, &xx0, &yy0, &zz0);
      }
      material_intp_2d(nx0, ny0, xx0, yy0, zz0, minfo->nx, minfo->ny, minfo->xx_deg, minfo->yy, minfo->chat);
      //wdbin_array2("test_c.bin", minfo->chat, minfo->nx, minfo->ny);
   }

   if(rho_type==3) {
      if(fformat==0) {
         //readbin_2dfile(outp, fn, nx, ny, xx, yy);
         readbin_2dfile(finfo->rho_fn, &nx0, &ny0, &xx0, &yy0, &zz0);
         //wdbin_array2("test_rho.bin", zz0, nx0, ny0);
      } else {
         read_2dfile(finfo->rho_fn, &nx0, &ny0, &xx0, &yy0, &zz0);
      }
      material_intp_2d(nx0, ny0, xx0, yy0, zz0, minfo->nx, minfo->ny, minfo->xx_deg, minfo->yy, minfo->rhohat);
      //wdbin_array2("test_rho.bin", minfo->rhohat, minfo->nx, minfo->ny);
   }

   if(w_type==3) {
      if(fformat==0) {
         //readbin_2dfile(outp, fn, nx, ny, xx, yy);
         readbin_2dfile(finfo->w_fn, &nx0, &ny0, &xx0, &yy0, &zz0);
         //wdbin_array2("test_wind.bin", zz0, nx0, ny0);
      } else {
         //read_2dfile(outp, fn, nx, ny, xx, yy);
         read_2dfile(finfo->w_fn, &nx0, &ny0, &xx0, &yy0, &zz0);
      }
      material_intp_2d(nx0, ny0, xx0, yy0, zz0, minfo->nx, minfo->ny, minfo->xx_deg, minfo->yy, minfo->what);
      //wdbin_array2("test_wind.bin", minfo->what, minfo->nx, minfo->ny);
   }
}

void what_adjust(fdmesh *minfo, int mrank, int nRank) {

   int i, j, gi;
   int zeroi=wzero, taperi=wtap;
   if(mrank <= nRank) {
     for(j=0; j<minfo->ny; j++) {
         for(i=0; i<minfo->nx; i++) {
            gi = i+minfo->xi_start;
              //chat[i][j] = 340.0;
            //chat[i][j] = 340.0 + 40 - 40.0*sin((j*minfo->dr/minfo->rpmax)*M_PI/2.0+M_PI/2.0);
            //what[i][j] = 40.0*sin((j*minfo->dr/minfo->rpmax)*M_PI);
            if(gi<(ghostL+zeroi)) {
               minfo->what[i][j] = 0.0;
            }
            if((gi>=(ghostL+zeroi)) & (gi<=(ghostL+taperi+zeroi))) {
               //what[i][j]=what[i][j]*sin( (minfo->thmin+(minfo->dth*i)) / (minfo->thmin+(minfo->dth*10.0)) * M_PI/2.0);
               minfo->what[i][j]=minfo->what[i][j]*pow(sin((minfo->thmin_global+(minfo->dth*(gi-ghostL-zeroi))) / \
                                 (minfo->thmin_global+(minfo->dth*taperi)) * M_PI/2.0),   2);
               //minfo->what[i][j]=minfo->what[i][j]*pow( minfo->dth*(gi-ghostL-zeroi) / (minfo->dth*taperi), 4) * 
               //pow(sin((minfo->dth*(gi-ghostL-zeroi)) / (minfo->dth*taperi) * M_PI/2.0), 5);
               //minfo->what[i][j]=minfo->what[i][j]*pow( minfo->dth*(gi-ghostL-zeroi) / (minfo->dth*taperi), 2);
            }
         }
     }
   }
   return;
} 


//Read the 2-d file and interpolate it
void read_2dfile(char *fn, int *nx, int *ny, double **xx, double **yy, double ***val)
{
   FILE *ptr_file;
   double *xx0, *yy0, **val0;
   double dx0, dy0;
   int nx0, ny0;
   int i=0, j;
   int k, l;
//   int ftype;
   char *pch, buf[1000], buf2[1000], buf0[1000];

   ptr_file = fopen( fn, "r" );
   if (!ptr_file) {
      sprintf(buf0, "     --- Error: The input file (%s) does not exist.\n", fn);
      elac_error(buf0, -1);
   }

   fgets( buf, 1000, ptr_file );
   strcpy( buf2, buf );
   pch = strtok(buf2, " \n");
//   ftype = atoi(pch);

   fgets( buf, 1000, ptr_file );
   strcpy( buf2, buf );
   pch = strtok(buf2, " \n");
   nx0 = atoi(pch);

   fgets( buf, 1000, ptr_file );
   strcpy( buf2, buf );
   pch = strtok(buf2, " \n");
   ny0 = atoi(pch);

   fgets( buf, 1000, ptr_file );
   strcpy( buf2, buf );
   pch = strtok(buf2, " \n");
   dx0 = (double) atof(pch);

   fgets( buf, 1000, ptr_file );
   strcpy( buf2, buf );
   pch = strtok(buf2, " \n");
   dy0 = (double) atof(pch);

   xx0 = (double*) malloc(nx0*sizeof(double));
   yy0 = (double*) malloc(ny0*sizeof(double));
   val0 = (double **) malloc(nx0*sizeof(double*));

//   elev2x2=(double**) malloc(2*sizeof(double*));
//   for(j=0;j<2;j++) {
//      elev2x2[j] = (double *) malloc(2*sizeof(double));
//   }

   for(j=0;j<nx0;j++) {
      val0[j] = (double *) malloc(ny0*sizeof(double));
      if(j==0) {
         xx0[j] = 0;
      } else {
         xx0[j] = xx0[j-1]+dx0;
      }
   }

   for(j=0;j<ny0;j++) {
      if(j==0) {
         yy0[j] = 0;
      } else {
         yy0[j] = yy0[j-1]+dy0;
      }
   }

   i=0;
   while ( fgets( buf, 1000, ptr_file )!=NULL ) {
      //l = floor(i/nx0); //angle varies faster
      //k = i-l*nx0;
      k = floor(i/ny0); //elev varies faster
      l = i-k*ny0;
      strcpy( buf2, buf );
      pch = strtok(buf2," \n");
      val0[k][l] = (double) atof(pch);
      //cout << xx0[k] << " " << yy0[l] << " " << sur_elev0[k][l] << endl;
      i++;
   }
   if(i<(nx0*ny0)) elac_error("The 2d file does not have enough entries.\n", -1);
   fclose( ptr_file );
   
   *nx = nx0;
   *ny = ny0;
   *xx = xx0;
   *yy = yy0;
   *val = val0;

   return;
}

//Read the 2-d bin file and interpolate it
void readbin_2dfile(char *fn, int *nx, int *ny, double **xx, double **yy, double ***val)
{
   FILE *ptr_file;
   double *xx0, *yy0, **val0;
   double temp;
   double dx0, dy0;
   int f1, nx0, ny0;
   int i=0, j=0;
   int k, l;
   char buf0[1000];
   int NeedSwap = 0;

   ptr_file = fopen( fn, "rb" );
   if (!ptr_file) {
      sprintf(buf0, "     --- Error: The input file (%s) does not exist.\n", fn);
      elac_error(buf0, -1);
   }

   fread(&f1, sizeof(int), 1, ptr_file);

   if(f1 == 1 || f1 == 2) { 
      NeedSwap = 0;
   } else {
      SwapBytes( &f1, 4 );
      if(f1 == 1 || f1 == 2) { 
         NeedSwap = 1;
      } else {
         sprintf(buf0, "     --- Error: The binary input file (%s) must start with an integer value of 1 or 2.\n", fn);
         elac_error(buf0, -1);
      }
   }
   
   fread(&nx0, sizeof(int), 1, ptr_file);
   if(NeedSwap) SwapBytes( &nx0, 4 );
   fread(&ny0, sizeof(int), 1, ptr_file);
   if(NeedSwap) SwapBytes( &ny0, 4 );
   fread(&dx0, sizeof(double), 1, ptr_file);
   if(NeedSwap) SwapBytes( &dx0, 8 );
   fread(&dy0, sizeof(double), 1, ptr_file);
   if(NeedSwap) SwapBytes( &dy0, 8 );
   xx0 = (double*) malloc(nx0*sizeof(double));
   yy0 = (double*) malloc(ny0*sizeof(double));
   val0 = (double **) malloc(nx0*sizeof(double*));
   val0[0] = (double *) malloc(nx0*ny0*sizeof(double));
   //printf("nx0 = %d, ny0 = %d, dx0 = %e, dy0 = %e\n", nx0, ny0, dx0, dy0);
   xx0[0] = 0;
   yy0[0] = 0;

   for(j=1;j<nx0;j++) {
      val0[j] = val0[0] + ny0*j;
      xx0[j] = xx0[j-1]+dx0;
   }
   for(j=1;j<ny0;j++) {
      yy0[j] = yy0[j-1]+dy0;
   }

   for(i=0;i<nx0;i++) {
      for(j=0;j<ny0;j++) {
         val0[i][j] = 0.0;
      }
//      printf("readbin_2dfile test...i, j=%d, %d\n",i,j);
   }

   i=0;
   while ( fread( &temp, sizeof(double), 1, ptr_file )>0 ) {
      if(NeedSwap) SwapBytes( &temp, 8 );
      //l = floor(i/nx0); //angle varies faster
      //k = i-l*nx0;
      k = floor(i/ny0); //elev varies faster
      l = i-k*ny0;

      //strcpy( buf2, buf );
      //pch = strtok(buf2," \n");
      val0[k][l] = temp;
      //cout << xx0[k] << " " << yy0[l] << " " << sur_elev0[k][l] << endl;
      i++;
   }
//      printf("readbin_2dfile test...k, l=%d, %d\n",k,l);
//   printf("i=%d\n",i);
   if(i<(nx0*ny0)) elac_error("The 2d file does not have enough entries.\n", -1);
   fclose( ptr_file );

   *nx = nx0;
   *ny = ny0;
   *xx = xx0;
   *yy = yy0;
   *val = val0;

   return;
}

void material_intp_2d (int nx0, int ny0, double *xx0, double *yy0, double **zz0, \
                       int nx1, int ny1, double *xx1, double *yy1, double **zz1) {
   long int xi, yi;
   double val, temp;
   double **z0_2x2;
   int i, j;

   z0_2x2=(double**) malloc(2*sizeof(double*));
   for(j=0;j<2;j++) {
      z0_2x2[j] = (double *) malloc(2*sizeof(double));
   }

   for(j=0;j<ny1;j++) {
         locate(yy0, ny0, yy1[j], &yi);
         if(yi==-1) yi++;
         if(yi==(ny0-1)) yi--;
      for(i=0;i<nx1;i++) {
         locate(xx0, nx0, xx1[i], &xi);
         if(xi==-1) xi++;
         if(xi==(nx0-1)) xi--;
            z0_2x2[0][0]=zz0[xi][yi];
            z0_2x2[1][0]=zz0[xi+1][yi];
            z0_2x2[0][1]=zz0[xi][yi+1];
            z0_2x2[1][1]=zz0[xi+1][yi+1];
         //printf("test1++++++++++++++\n");
         polin2(&xx0[xi],&yy0[yi],z0_2x2,2,2,xx1[i],yy1[j],&val,&temp); 
         //printf("test2++++++++++++++\n");
         //ll = j*NX+i;
         //sur_elev[ll] = val;
         zz1[i][j] = val;
         //cout << xx[i] << " " << yy[j] <<" " << xi << " " << yi << " " << val << endl;
      }
   }
   for(j=0;j<2;j++) {
      free((double*)z0_2x2[j]);
   }
   free((double**)z0_2x2);
}

void material_minmax(fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2) {
   int i, j;
   double whmin, whmax, cmin, cmax, rhomin, rhomax;
   char buf[200];

    whmin = minfo->what[ghostL][ghostL];
    whmax = minfo->what[ghostL][ghostL];
    cmin = minfo->chat[ghostL][ghostL];
    cmax = minfo->chat[ghostL][ghostL];
    rhomin = minfo->rhohat[ghostL][ghostL];
    rhomax = minfo->rhohat[ghostL][ghostL];

    for(i=ghostL; i<(minfo->nth-ghostL); i++) {
      for(j=ghostL; j<(minfo->nr-ghostL); j++) {
       if(minfo->what[i][j] < whmin) { whmin = minfo->what[i][j]; }
       if(minfo->what[i][j] > whmax) { whmax = minfo->what[i][j]; }
       if(minfo->rhohat[i][j] < rhomin) { rhomin = minfo->rhohat[i][j]; }
       if(minfo->rhohat[i][j] > rhomax) { rhomax = minfo->rhohat[i][j]; }
       if(minfo->chat[i][j] < cmin) { cmin = minfo->chat[i][j]; }
       if(minfo->chat[i][j] > cmax) { cmax = minfo->chat[i][j]; }
      }
    }

    MPI_Barrier(comm2);
    MPI_Allreduce(&whmin, &minfo->whmin, 1, MPI_DOUBLE, MPI_MIN, comm2);
    MPI_Allreduce(&whmax, &minfo->whmax, 1, MPI_DOUBLE, MPI_MAX, comm2);
    MPI_Allreduce(&cmin, &minfo->cmin, 1, MPI_DOUBLE, MPI_MIN, comm2);
    MPI_Allreduce(&cmax, &minfo->cmax, 1, MPI_DOUBLE, MPI_MAX, comm2);
    MPI_Allreduce(&rhomin, &minfo->rhomin, 1, MPI_DOUBLE, MPI_MIN, comm2);
    MPI_Allreduce(&rhomax, &minfo->rhomax, 1, MPI_DOUBLE, MPI_MAX, comm2);

    MPI_Barrier(comm2);
    sprintf(buf, "cmin = %f, cmax=%f, rhomin=%f, rhomax=%f, wmin=%f, wmax=%f",minfo->cmin, minfo->cmax, minfo->rhomin, minfo->rhomax, minfo->whmin, minfo->whmax);
   MPI_Barrier(comm2);
   mpi_printf(comm2, buf, 200, mrank, nRank);

//    printf("Process %d/%d: i=%d, cmin = %f, cmax=%f, rhomin=%f, rhomax=%f, wmin=%f, wmax=%f \n", mrank, nRank, i, cmin, cmax, rhomin, rhomax, whmin, whmax);
   return;
}


void material_smooth3 (double **u, int nx, int ny, int mrank, int nRank, MPI_Comm comm2) {
   int i, j, k, nn=10;
   double sum0, mean0;
   double **u1;

    u1 = (double **) malloc ( nx * sizeof(double *) );
    u1[0] = (double *) malloc ( nx * ny * sizeof(double) );

    for(i=1;i<nx;i++) {
        u1[i] = u1[0] + ny * i;
    }

   //top
   for(i=0; i<nx; i++) {
      u[i][ny-1] = u[i][ny-4];
      u[i][ny-2] = u[i][ny-4];
      u[i][ny-3] = u[i][ny-4];
   }
   //bottom
   for(i=0; i<nx; i++) {
      u[i][0] = u[i][3];
      u[i][1] = u[i][3];
      u[i][2] = u[i][3];
   }

   //left
   if(mrank == 0) {
      for(j=0; j<ny; j++) {
         u[0][j] = u[3][j];
         u[1][j] = u[3][j];
         u[2][j] = u[3][j];
      }
   }

   //right
   if(mrank == nRank) {
      for(j=0; j<ny; j++) {
         u[nx-1][j] = u[nx-4][j];
         u[nx-2][j] = u[nx-4][j];
         u[nx-3][j] = u[nx-4][j];
      }
   }
   
   //interpolation in y direction
   for(i=0; i<nx; i++) {
      for(j=3; j<(ny-3); j++) {
         sum0 = u[i][j];
         for(k=1; k<=nn; k++) {
            if((j-k)<4) {
               sum0 += u[i][4];
            } else {
               sum0 += u[i][j-k];
            }
         }
         for(k=1; k<=nn; k++) {
            if((j+k)>(ny-4)) {
               sum0 += u[i][ny-4];
            } else {
               sum0 += u[i][j+k];
            }
         }
         mean0 = sum0/(nn*2+1);
//        sum0=u[i-3][j-3] + u[i-2][j-3] + u[i-1][j-3] + u[i][j-3] + u[i+1][j-3] + u[i+2][j-3] + u[i+3][j-3] + 
//             u[i-3][j-2] + u[i-2][j-2] + u[i-1][j-2] + u[i][j-2] + u[i+1][j-2] + u[i+2][j-2] + u[i+3][j-2] + 
//             u[i-3][j-1] + u[i-2][j-1] + u[i-1][j-1] + u[i][j-1] + u[i+1][j-1] + u[i+2][j-1] + u[i+3][j-1] + 
//             u[i-3][j] + u[i-2][j] + u[i-1][j] + u[i][j] + u[i+1][j] + u[i+2][j] + u[i+3][j] + 
//             u[i-3][j+1] + u[i-2][j+1] + u[i-1][j+1] + u[i][j+1] + u[i+1][j+1] + u[i+2][j+1] + u[i+3][j+1] + 
//             u[i-3][j+2] + u[i-2][j+2] + u[i-1][j+2] + u[i][j+2] + u[i+1][j+2] + u[i+2][j+2] + u[i+3][j+2] + 
//             u[i-3][j+3] + u[i-2][j+3] + u[i-1][j+3] + u[i][j+3] + u[i+1][j+3] + u[i+2][j+3] + u[i+3][j+3];
//        mean0=sum0/49;
        u1[i][j] = mean0;
      }
   }

   for(i=3; i<(nx-3); i++) {
      for(j=3; j<(ny-3); j++) {
         sum0 = u1[i-3][j] + u1[i-2][j] + u1[i-1][j] + u1[i][j] + u1[i+1][j] + u1[i+2][j] + u1[i+3][j];
         mean0 = sum0/7.0;
         u[i][j] = mean0;
      }
   }

//   for(i=3; i<(nx-3); i++) {
//      for(j=3; j<(ny-3); j++) {
//         u[i][j] = u1[i][j];
//      }
//   }
         
   free_array2( u1 );
   
   return;
}

void writetxt_2dfile(char *fn, int nx, int ny, double *xx, double *yy, double **zz)
{
   FILE *ptr_file;
   int i, j;
   double dx, dy;
   dx = xx[1] - xx[0];
   dy = yy[1] - yy[0];
   ptr_file = fopen(fn, "w");
   fprintf(ptr_file, "%d\n", 2);
   fprintf(ptr_file, "%d\n", nx);
   fprintf(ptr_file, "%d\n", ny);
   fprintf(ptr_file, "%e\n", dx);
   fprintf(ptr_file, "%e\n", dy);
   for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
         fprintf(ptr_file, "%f\n", zz[i][j]);
      }
   }
   fclose(ptr_file);
}

void writebin_2dfile(char *fn, int nx, int ny, double *xx, double *yy, double **zz)
{
   FILE *ptr_file;
   int i, j;
   double dx, dy;
   int ftype;
   ftype=2;
   dx = xx[1] - xx[0];
   dy = yy[1] - yy[0];
   ptr_file = fopen(fn, "wb");

   fwrite(&ftype, sizeof(int), 1, ptr_file);
   fwrite(&nx, sizeof(int), 1, ptr_file);
   fwrite(&ny, sizeof(int), 1, ptr_file);
   fwrite(&dx, sizeof(double), 1, ptr_file);
   fwrite(&dy, sizeof(double), 1, ptr_file);

   for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
         fwrite(&zz[i][j], sizeof(double), 1, ptr_file);
      }
   }
   fclose(ptr_file);
}

void writebin_2col(char *fn, double *xx, double *yy, int imax) {
   FILE *ptr_file;
   int i=0, ftype;
   ftype=1;
   ptr_file = fopen(fn, "wb");
   fwrite(&ftype, sizeof(int), 1, ptr_file);
   while(i<imax) {
      fwrite(&xx[i], sizeof(double), 1, ptr_file);
      fwrite(&yy[i], sizeof(double), 1, ptr_file);
//      printf("%f ", xx[i]);
//      printf("%f\n", yy[i]);
      i++;
   }
   fclose( ptr_file );
}

