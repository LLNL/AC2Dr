#include "ac2dr_aux.h"
#include "ac2dr_config.h"

#if (FD_ORDER==2)
#define spatial_derivatives dqdt2_sym_ghost_sg
#define spatial_derivatives180 dqdt2_sym_ghost_sg
#endif 

#if (FD_ORDER==6)
#define spatial_derivatives dqdt6_damping
#define spatial_derivatives180 dqdt4_sym_ghost_sg
#endif 


void ac2dr_RK4_damping ( fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm2 ) {

    unsigned long i, j, k, l;
    unsigned long nt, nx, ny;
    double dt;
//    double dr,ds;
//    int xxi, yxi; 
//    double f1, f2;
    //FILE *ptr_file;
    char buf[1000];

    double **u0, **u1, **u2, **u3, **u4;
    double **v0, **v1, **v2, **v3, **v4;
    double **w0, **w1, **w2, **w3, **w4;

    double **u = minfo->u;
    double **v = minfo->v;
    double **w = minfo->w;

    double **chat = minfo->chat;
    double **rhohat = minfo->rhohat;
    double **what = minfo->what;
    double **wh0, **wh1, **wh2;

    nt = minfo->nt;
    dt = minfo->dt;
    nx = minfo->nth;
    ny = minfo->nr;

    u0 = alloc_array2( nx, ny );
    u1 = alloc_array2( nx, ny );
    u2 = alloc_array2( nx, ny );
    u3 = alloc_array2( nx, ny );
    u4 = alloc_array2( nx, ny );

    v0 = alloc_array2( nx, ny );
    v1 = alloc_array2( nx, ny );
    v2 = alloc_array2( nx, ny );
    v3 = alloc_array2( nx, ny );
    v4 = alloc_array2( nx, ny );

    w0 = alloc_array2( nx, ny );
    w1 = alloc_array2( nx, ny );
    w2 = alloc_array2( nx, ny );
    w3 = alloc_array2( nx, ny );
    w4 = alloc_array2( nx, ny );

    wh0 = alloc_array2( nx, ny );
    wh1 = alloc_array2( nx, ny );
    wh2 = alloc_array2( nx, ny );

    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
         wh0[i][j] = 0.0;
         wh1[i][j] = what[i][j];
         wh2[i][j] = what[i][j];
      }
    }


//    printf("TEST: %d/%d\n", mrank, nRank);
    ac2dr_initout(minfo, finfo, u, v, w, rhohat, chat, what, mrank, nRank, comm2 );
//    printf("TEST2: %d/%d\n", mrank, nRank);

   int gi, zeroi=floor(wtap2*0.1/minfo->dr), taperi=floor(wtap2*0.9/minfo->dr);
   if(mrank <= nRank) {
     for(j=0; j<minfo->ny; j++) {
         for(i=0; i<minfo->nx; i++) {
            gi = i+minfo->xi_start;
            if(gi<(ghostL+zeroi)) {
               wh2[i][j] = 0.0;
            }
            if((gi>=(ghostL+zeroi)) & (gi<=(ghostL+taperi+zeroi))) {
               wh2[i][j]=wh2[i][j]*pow(sin((minfo->thmin_global+(minfo->dth*(gi-ghostL-zeroi))) / \
                         (minfo->thmin_global+(minfo->dth*taperi)) * M_PI/2.0),   4);
            }
         }
     }
   }

   what = wh0;

    MPI_Barrier(comm2);
    k=0;
    if(mrank==0) {
    snprintf(buf,998,"k = %lu/%lu", k, nt); //print_message(ptr_file, buf);
    puts(buf);
    }
//    ac2dr_out( 1, k, minfo, finfo, u, v, w, rhohat, chat, what );

    for(l=0; l<minfo->sta_num; l++) {
        i = minfo->stas[l].xi;
        j = minfo->stas[l].yi;
        minfo->stas[l].pout->elements[k] = u[i][j];
    }

    MPI_Barrier(comm2);
    for(l=0; l<finfo->line_num; l++) {
       j = k/floor(finfo->lines[l].linedt / minfo->dt + 0.5);
       for(i=0; i<minfo->nx; i++) {
         finfo->lines[l].lineimg[i][j] = u[i][finfo->lines[l].yi];
       }
    }

    MPI_Barrier(comm2);
    for(k=1; k<nt; k++) {
    if((k*minfo->dt)>stg1_t0) {what = wh1;}
    if((k*minfo->dt)>stg2_t0) {what = wh2;}

    //printf("k = %lu\n",k);
    if( (k%100==0) & (mrank==0) ) {
    snprintf(buf,998,"k = %lu/%lu", k, nt); //print_message(ptr_file, buf);
    puts(buf);
    }

//       dqdt6_sym_sg_test( u, v, w, u0, v0, w0, rhohat, what, chat, minfo);
//not necessary  updatesrc( k, s1, r1, rhohat, chat, jac, minfo, 1);
//    if(k%1==0) {
//       sprintf(buf,"output/test_dU1_snapshot%06d.dat",k);
//       wd_array2( buf, u0, minfo->nx, minfo->ny );
//       sprintf(buf,"output/test_dV1_snapshot%06d.dat",k);
//       wd_array2(buf, v0, minfo->nx, minfo->ny );
//       sprintf(buf,"output/test_dW1_snapshot%06d.dat",k);
//       wd_array2(buf, w0, minfo->nx, minfo->ny );
//    }

//    printf("RK4 mpi test...\n");
        MPI_Barrier(comm2);
        spatial_derivatives( u, v, w, u1, v1, w1, rhohat, what, chat, minfo, k, 1 );

/*
    if(k%1==0) {
   char fn[500];
   strcpy(fn, finfo->outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_u_%06ld.dat", k);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, u1, nx, ny, minfo->nx_global, minfo->ny, mrank, nRank, comm2);
   }
*/
        MPI_Barrier(comm2);
        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt/2.0*u1[i][j];
                v0[i][j] = v[i][j] + dt/2.0*v1[i][j];
                w0[i][j] = w[i][j] + dt/2.0*w1[i][j];
            }
        }

        MPI_Barrier(comm2);
        updateBoundary(u0, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(v0, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(w0, nx, ny, ghostL, mrank, nRank, comm2);

//MPI communication

if(BC_ENFORCE==1) {
        MPI_Barrier(comm2);
   rigidbound_var (u0, v0, w0, minfo->nx, minfo->ny, mrank, nRank, comm2);
}

        MPI_Barrier(comm2);
        spatial_derivatives( u0, v0, w0, u2, v2, w2, rhohat, what, chat, minfo, k, 2 );

        MPI_Barrier(comm2);
        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt/2.0*u2[i][j];
                v0[i][j] = v[i][j] + dt/2.0*v2[i][j];
                w0[i][j] = w[i][j] + dt/2.0*w2[i][j];
            }
        }

        MPI_Barrier(comm2);
        updateBoundary(u0, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(v0, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(w0, nx, ny, ghostL, mrank, nRank, comm2);


if(BC_ENFORCE==1) {
        MPI_Barrier(comm2);
   rigidbound_var (u0, v0, w0, minfo->nx, minfo->ny, mrank, nRank, comm2);
}
 
        MPI_Barrier(comm2);
        spatial_derivatives( u0, v0, w0, u3, v3, w3, rhohat, what, chat, minfo, k, 3 );

        MPI_Barrier(comm2);
        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt*u3[i][j];
                v0[i][j] = v[i][j] + dt*v3[i][j];
                w0[i][j] = w[i][j] + dt*w3[i][j];
            }
        }

        MPI_Barrier(comm2);
        updateBoundary(u0, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(v0, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(w0, nx, ny, ghostL, mrank, nRank, comm2);


if(BC_ENFORCE==1) {
        MPI_Barrier(comm2);
   rigidbound_var (u0, v0, w0, minfo->nx, minfo->ny, mrank, nRank, comm2);
}
 
        MPI_Barrier(comm2);
        spatial_derivatives( u0, v0, w0, u4, v4, w4, rhohat, what, chat, minfo, k, 4 );

        MPI_Barrier(comm2);
        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u[i][j] = u[i][j] + dt/6.0*(u1[i][j] + 2.0*u2[i][j] + 2.0*u3[i][j] + u4[i][j]);
                v[i][j] = v[i][j] + dt/6.0*(v1[i][j] + 2.0*v2[i][j] + 2.0*v3[i][j] + v4[i][j]);
                w[i][j] = w[i][j] + dt/6.0*(w1[i][j] + 2.0*w2[i][j] + 2.0*w3[i][j] + w4[i][j]);
            }
        }

//point src
/*
        for(l=0; l<minfo->src_num; l++) {
            i=minfo->srcs[l].xi;
            j=minfo->srcs[l].yi;
            r= minfo->R + dr*(j-ghostL);
            u[i][j] += dt* minfo->srcs[l].q->elements[k]*rhohat[i][j]*chat[i][j]*chat[i][j] \
                     /(M_PI*dr*pow(dth*r/2.0,2))/rhohat[i][j];
      //printf("***********ac2dr_dqdt, source test, %d\n", tk);
        }
*/
        MPI_Barrier(comm2);
        updateBoundary(u, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(v, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(w, nx, ny, ghostL, mrank, nRank, comm2);

if(BC_ENFORCE==1) {
        MPI_Barrier(comm2);
   rigidbound_var (u, v, w, minfo->nx, minfo->ny, mrank, nRank, comm2);
}       

//Filter
        MPI_Barrier(comm2);
        filter6(u, u1, minfo, 0.2);
        filter6(v, v1, minfo, 0.2);
        filter6(w, w1, minfo, 0.2);

        MPI_Barrier(comm2);
        updateBoundary(u1, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(v1, nx, ny, ghostL, mrank, nRank, comm2);
        MPI_Barrier(comm2);
        updateBoundary(w1, nx, ny, ghostL, mrank, nRank, comm2);
if(BC_ENFORCE==1) {
        MPI_Barrier(comm2);
   rigidbound_var (u1, v1, w1, minfo->nx, minfo->ny, mrank, nRank, comm2);
}
        MPI_Barrier(comm2);
        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u[i][j] = u1[i][j];
                v[i][j] = v1[i][j];
                w[i][j] = w1[i][j];
            }
        }
//Filter

        MPI_Barrier(comm2);
 
        for(l=0; l<minfo->sta_num; l++) {
           i = minfo->stas[l].xi;
           j = minfo->stas[l].yi;
            minfo->stas[l].pout->elements[k] = u[i][j];
        }


    MPI_Barrier(comm2);
    for(l=0; l<finfo->line_num; l++) {
       j = k/floor(finfo->lines[l].linedt / minfo->dt + 0.5);
       for(i=0; i<minfo->nx; i++) {
         finfo->lines[l].lineimg[i][j] = u[i][finfo->lines[l].yi];
       }
    }


//    printf("RK4 test\n");
/*
    if(k%1==0) {
   char fn[500];
   strcpy(fn, finfo->outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_u_%06ld.dat", k);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, u, nx, ny, minfo->nx_global, minfo->ny, mrank, nRank, comm2);
   }
*/
//       sprintf(buf,"output/test_U_snapshot%06d.dat",k);
//       wd_array2( buf, u, minfo->nx, minfo->ny );
//       sprintf(buf,"output/test_V_snapshot%06d.dat",k);
//       wd_array2(buf, v, minfo->nx, minfo->ny );
//       sprintf(buf,"output/test_W_snapshot%06d.dat",k);
//       wd_array2(buf, w, minfo->nx, minfo->ny );
//    }
/*    if(k==(nt-1)) {
       wd_array2( "output/test_U_snapshot100.dat", u, minfo->nx, minfo->ny );
       wd_array2( "output/test_V_snapshot100.dat", v, minfo->nx, minfo->ny );
       wd_array2( "output/test_W_snapshot100.dat", w, minfo->nx, minfo->ny );
    }
*/
        MPI_Barrier(comm2);
        ac2dr_out( 1, k, minfo, finfo, u, v, w, rhohat, chat, what, mrank, nRank, comm2 );

        MPI_Barrier(comm2);
        if(pmaximgs == 1) {
            for(i=0; i<nx; i++) {
               for(j=0; j<ny; j++) {
 //                 if((u[i][j]*rhohat[i][j]*chat[i][j]) > minfo->pmax[i][j]) {
//                     minfo->pmax[i][j] = u[i][j]*rhohat[i][j]*chat[i][j];
//                     minfo->pmax_t[i][j] = k*dt;
//                  }
                  if(fabs(u[i][j]) > minfo->pmaxabs[i][j]) {
                     minfo->pmaxabs[i][j] = fabs(u[i][j]);
                     minfo->pmaxabs_t[i][j] = k*dt;
                  }
               }
            }
        } //pmaximgs
    } //for(k)

    MPI_Barrier(comm2);
//    printf("***********ac2dr_RK4.c is working\n");
    ac2dr_out( 0, k, minfo, finfo, u, v, w, rhohat, chat, what, mrank, nRank, comm2 );
    //fclose(ptr_file);
//    ac2dr_out( 0, k, minfo, finfo, u, v, w, rhohat, chat, what );
//    printf("source number = %d\n", minfo->src_num);
//    wd_vec("./output/srcf.txt",minfo->srcs[0].q);
//    wd_vec("./output/srcf2.txt",minfo->srcs[0].midq);

    if(pmaximgs==1) {
      char fn[500];
/*
      MPI_Barrier(comm2);
      strcpy( fn, finfo->outdir);
      strcat( fn, "/"); 
      strcat( fn, "pmaximg.bin");
      wdbin_mpi_array2(fn, minfo->pmax, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, mrank, nRank, comm2);

      MPI_Barrier(comm2);
      strcpy( fn, finfo->outdir);
      strcat( fn, "/"); 
      strcat( fn, "pmaximg_t.bin");
      wdbin_mpi_array2(fn, minfo->pmax_t, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, mrank, nRank, comm2);
*/
      MPI_Barrier(comm2);
      strcpy( fn, finfo->outdir);
      strcat( fn, "/"); 
      strcat( fn, "pmaxabsimg.bin");
      if(image10 == 0) {
      wdbin_mpi_array2(fn, minfo->pmaxabs, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                       minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
      } else {
      wdbin_mpi_array2_reduced(fn, minfo->pmaxabs, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                       minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
      }

      MPI_Barrier(comm2);
      strcpy( fn, finfo->outdir);
      strcat( fn, "/"); 
      strcat( fn, "pmaxabsimg_t.bin");
      if(image10 == 0) {
      wdbin_mpi_array2(fn, minfo->pmaxabs_t, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                       minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
      } else {
      wdbin_mpi_array2_reduced(fn, minfo->pmaxabs_t, minfo->nx, minfo->ny, minfo->nx_global, minfo->ny, \
                       minfo->dth*180.0/M_PI, minfo->dr, mrank, nRank, comm2);
      }
    }

    MPI_Barrier(comm2);
    if(finfo->line_num > 0) {
      char fn[500];
      for(l=0; l<finfo->line_num; l++) {
         strcpy( fn, finfo->outdir);
         strcat( fn, "/");
         strcat( fn, finfo->lines[l].file );
         sprintf( buf, "_p.dat");
         strcat( fn, buf );
      //printf("image output = %s\n", fn);

      //puts(fn);
         if(finfo->lines[l].format==0) {
      //wdbin_array2( fn, p, minfo->nx, minfo->ny );
            wdbin_mpi_array2(fn, finfo->lines[l].lineimg, minfo->nx, finfo->lines[l].tl, \
                             minfo->nx_global, finfo->lines[l].tl, \
                             minfo->dth*180.0/M_PI, finfo->lines[l].linedt, mrank, nRank, comm2);
         }
      }
    }
    MPI_Barrier(comm2);
    return;
}












