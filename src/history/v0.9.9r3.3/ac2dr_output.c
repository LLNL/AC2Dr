#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void ac2dr_out( int mode, unsigned long k, fdmesh *minfo, innout *finfo, \
                 double **u, double **v, double **w, \
                 double **rhohat, double **chat, double **what, int mrank, int nRank, MPI_Comm comm2) {
//mode=0: write pressure on receivers
//mode=1: write pressure snapshot
    int i, j, m, l;
    char fn[500], buf[500];
    double **p;
    p = alloc_array2( minfo->nx, minfo->ny );

    if(mode == 0) {
        for( m=0; m<minfo->sta_num; m++ ) {
            strcpy( fn, finfo->outdir);
            strcat( fn, "/"); 
            strcat( fn, minfo->stas[m].name );
            if(minfo->stas[m].format==0) {
               strcat( fn, ".bin" );
               wdbin_recv2(fn, minfo->stas[m].pout, minfo->stas[m].x*180/M_PI, minfo->stas[m].y, minfo);
            }
            if(minfo->stas[m].format==1) {
               strcat( fn, ".txt" );
               wd_recv2(fn, minfo->stas[m].pout, minfo->stas[m].x*180/M_PI, minfo->stas[m].y, minfo);
            }

        }
    }
 
    if(mode == 1) {
        for( m=0; m<finfo->img_num; m++ ) {
            l = floor(finfo->imgs[m].imgdt/minfo->dt+0.5);
            if(k%l == 0) {
//            printf("img out TEST: %d/%d\n", mrank, nRank);
//    printf("image mode = %s\n", finfo->imgs[i].mode);
                if(strchr(finfo->imgs[m].mode, 'p')!=NULL)
                {
                  strcpy( fn, finfo->outdir);
                  strcat( fn, "/"); 
                  strcat( fn, finfo->imgs[m].file );
                  sprintf( buf, "_p_%08.3f.dat", minfo->dt*k );
                  strcat( fn, buf );
                  //printf("image output = %s\n", fn);

                  //puts(fn);
                  for(j=0; j<minfo->ny; j++) {
                        for(i=0; i<minfo->nx; i++) {
                           p[i][j] = u[i][j];
                        }
                  }
/*
                  if(finfo->imgs[m].format==0) {
                     //wdbin_array2( fn, p, minfo->nx, minfo->ny );
                     if(image10 == 0) {
                     wdbin_mpi_array2(fn, p, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                                      minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     } else {
                     wdbin_mpi_array2_reduced(fn, p, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                                      minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     }
                  }
                  if(finfo->imgs[m].format==1) {
                     //wd_array2( fn, p, minfo->nx, minfo->ny );
                  }
*/
                if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
                     if(finfo->imgs[m].ghostL == 1) {
                           wdbin_mpi_2dimg(fn, p, 2, k*minfo->dt, \
                                           minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                                           minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     } else {
                           wdbin_mpi_2dimg_ngl(fn, p, 2, k*minfo->dt, \
                                          minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                                          minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     }
                }

            //printf("img out TEST2: %d/%d\n", mrank, nRank);
                }

                if(strchr(finfo->imgs[m].mode, 'v')!=NULL)
                {
                  strcpy( fn, finfo->outdir);
                  strcat( fn, "/"); 
                  strcat( fn, finfo->imgs[m].file );
                  sprintf( buf, "_v_%08.3f.dat", minfo->dt*k );
                  strcat( fn, buf );
                  //puts(fn);
                  for(j=0; j<minfo->ny; j++) {
                        for(i=0; i<minfo->nx; i++) {
                           p[i][j] = v[i][j];
                        }
                  }
                  if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
                     if(finfo->imgs[m].ghostL == 1) {
                           wdbin_mpi_2dimg(fn, p, 2, k*minfo->dt, \
                                           minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                                           minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     } else {
                           wdbin_mpi_2dimg_ngl(fn, p, 2, k*minfo->dt, \
                                          minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                                          minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     }
                  }
                }
                if(strchr(finfo->imgs[m].mode, 'w')!=NULL)
                {
                  strcpy( fn, finfo->outdir);
                  strcat( fn, "/"); 
                  strcat( fn, finfo->imgs[m].file );
                  sprintf( buf, "_w_%08.3f.dat", minfo->dt*k );
                  strcat( fn, buf );
                  //puts(fn);
                  for(j=0; j<minfo->ny; j++) {
                        for(i=0; i<minfo->nx; i++) {
                           p[i][j] = w[i][j];
                        }
                  }
                  if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
                     if(finfo->imgs[m].ghostL == 1) {
                           wdbin_mpi_2dimg(fn, p, 2, k*minfo->dt, \
                                           minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                                           minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     } else {
                           wdbin_mpi_2dimg_ngl(fn, p, 2, k*minfo->dt, \
                                          minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                                          minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
                     }
                  }

                }
            }
        }
    }
    free_array2(p);
}


void ac2dr_initout( fdmesh *minfo, innout *finfo, \
                 double **u, double **v, double **w, \
                 double **rhohat, double **chat, double **what, int mrank, int nRank, MPI_Comm comm2 ) {
//mode=0: write pressure on receivers
//mode=1: write pressure snapshot
    int i, j, m;
    char fn[500], buf[500];
    double **p;
    p = alloc_array2( minfo->nx, minfo->ny );

    m=0;

//    printf("inout_image.format=%d: %d/%d\n", finfo->imgs[m].format, mrank, nRank);
    strcpy( fn, finfo->outdir);
    strcat( fn, "/"); 
    strcat( fn, finfo->imgs[m].file );
    sprintf( buf, "_p_%08.3f.dat", 0.0);
    strcat( fn, buf );
    for(j=0; j<minfo->ny; j++) {
      for(i=0; i<minfo->nx; i++) {
         p[i][j] = u[i][j];
      }
    }
    if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
       //wd_array2( fn, p, minfo->nx, minfo->ny );
       if(finfo->imgs[m].ghostL == 1) {
       //wdbin_mpi_array2(fn, p, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
       //                 minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
         wdbin_mpi_2dimg(fn, p, 2, 0.0, \
                         minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                         minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       } else {
         //printf("initout test1: %d/%d\n", mrank, nRank);
         wdbin_mpi_2dimg_ngl(fn, p, 2, 0.0, \
                             minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                             minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
         //printf("initout test2: %d/%d\n", mrank, nRank);
       }
    }
//    if(finfo->imgs[m].format == 0) {
//       //wdbin_mpi_array2( fn, p, minfo->nx, minfo->ny );
//       wdbin_mpi_array2(fn, p, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
//                        minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
//    }

/*
    strcpy( fn, finfo->outdir);
    strcat( fn, "/"); 
    strcat( fn, finfo->imgs[k].file );
    sprintf( buf, "_v_%08.3f.dat", 0.0 );
    strcat( fn, buf );
    for(j=0; j<minfo->ny; j++) {
      for(i=0; i<minfo->nx; i++) {
         p[i][j] = v[i][j];
      }
    }
    if(finfo->imgs[m].format == 1) {
       wd_array2( fn, p, minfo->nx, minfo->ny );
    }
    if(finfo->imgs[m].format == 0) {
       wdbin_array2( fn, p, minfo->nx, minfo->ny );
    }


    strcpy( fn, finfo->outdir);
    strcat( fn, "/"); 
    strcat( fn, finfo->imgs[k].file );
    sprintf( buf, "_w_%08.3f.dat", 0.0 );
    strcat( fn, buf );
    for(j=0; j<minfo->ny; j++) {
      for(i=0; i<minfo->nx; i++) {
         p[i][j] = w[i][j];
      }
    }
    if(finfo->imgs[k].format == 1) {
       wd_array2( fn, p, minfo->nx, minfo->ny );
    }
    if(finfo->imgs[k].format == 0) {
       wdbin_array2( fn, p, minfo->nx, minfo->ny );
    }
*/

    strcpy( fn, finfo->outdir);
    strcat( fn, "/"); 
    strcat( fn, finfo->imgs[m].file );
    sprintf( buf, "_rho_%08.3f.dat", 0.0 );
    strcat( fn, buf );
    for(j=0; j<minfo->ny; j++) {
      for(i=0; i<minfo->nx; i++) {
         p[i][j] = rhohat[i][j];
      }
    }
/*
    if(finfo->imgs[k].format == 1) {
       //wd_array2( fn, p, minfo->nx, minfo->ny );
       wdbin_mpi_array2(fn, p, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                        minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
    }
    if(finfo->imgs[k].format == 0) {
       //wdbin_array2( fn, p, minfo->nx, minfo->ny );
       wdbin_mpi_array2(fn, p, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                        minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
    }
*/
    if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
       if(finfo->imgs[m].ghostL == 1) {
         wdbin_mpi_2dmat(fn, p, 2, \
                         minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                         minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       } else {
         wdbin_mpi_2dmat_ngl(fn, p, 2, \
                             minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                             minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       }
    }

    strcpy( fn, finfo->outdir);
    strcat( fn, "/"); 
    strcat( fn, finfo->imgs[m].file );
    sprintf( buf, "_c_%08.3f.dat", 0.0 );
    strcat( fn, buf );
    for(j=0; j<minfo->ny; j++) {
      for(i=0; i<minfo->nx; i++) {
         p[i][j] = chat[i][j];
      }
    }
    if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
       if(finfo->imgs[m].ghostL == 1) {
         wdbin_mpi_2dmat(fn, p, 2, \
                         minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                         minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       } else {
         wdbin_mpi_2dmat_ngl(fn, p, 2, \
                             minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                             minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       }
    }


    strcpy( fn, finfo->outdir);
    strcat( fn, "/"); 
    strcat( fn, finfo->imgs[m].file );
    sprintf( buf, "_wind_%08.3f.dat", 0.0 );
    strcat( fn, buf );
    for(j=0; j<minfo->ny; j++) {
      for(i=0; i<minfo->nx; i++) {
         p[i][j] = what[i][j];
      }
    }
    if(finfo->imgs[m].format == 0 || finfo->imgs[m].format == 1) {
       if(finfo->imgs[m].ghostL == 1) {
         wdbin_mpi_2dmat(fn, p, 2, \
                         minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                         minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       } else {
         wdbin_mpi_2dmat_ngl(fn, p, 2, \
                             minfo->nx, minfo->ny, minfo->nx_global, minfo->ny,
                             minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
       }
    }

//    printf("inout_TEST2: %d/%d\n", mrank, nRank);
    free_array2(p);
}
