#include "ac2dr_aux.h"
#include "ac2dr_config.h"

//const double gam = 1.4;
//const double R0=8.3143;
//const double M0=0.029;
//const double Rd=8.3143/0.029;
const char ver0[] = "0.9.9r2";
const char rdate[] = "September 20th, 2023";
const char author[] = "Keehoon Kim (LLNL, kim84@llnl.gov)";

#if (TMARCH==1)
#define time_marching ac2dr_RK4
#define time_marching180 ac2dr_RK4180
#endif

#if (TMARCH==2)
#define time_marching ac2dr_FD4
#define time_marching180 ac2dr_RK4180
#endif


int main(int argc, char *argv[]) {
    int mrank, totalp, nRank;
    //finite-difference mesh
    fdmesh minfo;
    //input & output information
    innout finfo;
    time_t t0, t1;
    int *proc_ranks, proc;
//    int i, j, nx, ny;
    char buf[200];

    MPI_Group group1, group2;
    MPI_Comm comm2;

    t0 = time(0);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mrank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalp);

    if(mrank==0) {
      printf("======================================================================\n");
      printf("=      AC2Dr version v%s                     \n", ver0);
      printf("=                                            \n");
      printf("=                    Release date: %s        \n", rdate);
      printf("=                    Contact: %s             \n", author);
      printf("======================================================================\n");
    }

    if(mrank==0) {
      printf("=============================================\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    strcpy(buf, "MPI initialized.");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    MPI_Barrier(MPI_COMM_WORLD);
    if(mrank==0) {
      printf("=============================================\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    strcpy(finfo.logf, "ac2dr.out");

//    ac2dr_input( argv[1], &minfo, &finfo, &rhohat, &chat, &what, mrank); 
    get_path( argv[1], &finfo, mrank );
    strcpy(buf, "Path reading.... Done.");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    if(mrank==0) {
      printf("=============================================\n");
    }

    get_grid( argv[1], &minfo, mrank );
    strcpy(buf, "Grid reading.... Done.");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    if(mrank==0) {
      printf("=============================================\n");
    }

    get_vel( argv[1], &finfo, &minfo, mrank );
    strcpy(buf, "Velocity reading.... Done.");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    if(mrank==0) {
      printf("=============================================\n");
    }

    get_dens( argv[1], &finfo, &minfo, mrank );
    strcpy(buf, "Density reading.... Done.");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    if(mrank==0) {
      printf("=============================================\n");
    }

    get_wind( argv[1], &finfo, &minfo, mrank );
    strcpy(buf, "Wind reading.... Done.");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    if(mrank==0) {
      printf("=============================================\n");
    }

    //domain decomposition

    MPI_Barrier(MPI_COMM_WORLD);
    decompose_domain (&minfo, &finfo, mrank, totalp, &nRank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    proc_ranks = (int *) malloc (totalp * sizeof(int));
    for(proc = 0; proc <= nRank; proc++) {
      proc_ranks[proc] = proc;
    }

   MPI_Comm_group(MPI_COMM_WORLD, &group1);
   MPI_Group_incl(group1, nRank+1, proc_ranks, &group2);
   MPI_Comm_create(MPI_COMM_WORLD, group2, &comm2);

   if(comm2 != MPI_COMM_NULL) {

   strcpy(buf, "Domain decomposition... Done.");
   mpi_printf(comm2, buf, 200, mrank, nRank);
    if(mrank==0) {
      printf("=============================================\n");
    }

   MPI_Barrier(comm2);
    //material specification
   material_prop(&minfo, &finfo, mrank, nRank, comm2);

//   MPI_Barrier(comm2);
//   nx = minfo.nx;
//   ny = minfo.ny;
//   for(j=0; j<ny; j++) {
//      for(i=0; i<nx; i++) {
//         minfo.chat[i][j] = 330.0 + 20.0*sin((j*minfo.dr/3000)*M_PI/2.0+M_PI/2.0);
//      }
//   }
 
   MPI_Barrier(comm2);
   strcpy(buf, "Material model created... Done.");
   mpi_printf(comm2, buf, 200, mrank, nRank);
    if(mrank==0) {
      printf("=============================================\n");
    }

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   strcpy(buf, "Material boundary updated... Done.");
   mpi_printf(comm2, buf, 200, mrank, nRank);
    if(mrank==0) {
      printf("=============================================\n");
    }

//   strcpy(fn, finfo.outdir);
//   strcat(fn, "/");
//   strcat(fn, "snapshot");
//   sprintf(buf, "_w0_%06d.dat", 0);
//   strcat(fn, buf);
//   wdbin_mpi_array2(fn, minfo.what, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank, nRank, comm2);

//   mpi_printf(comm2, "MPI_test.", mrank, nRank);
   //wind field adjustment
   MPI_Barrier(comm2);
   what_adjust(&minfo, mrank, nRank);

   //check min, max, global min/max
   MPI_Barrier(comm2);
   material_minmax(&minfo, mrank, nRank, comm2); 
    if(mrank==0) {
      printf("=============================================\n");
    }

   //matrix shift
   //updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   shift_domain (minfo.chat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   shift_domain (minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   shift_domain (minfo.what, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   rigidbound_material (minfo.chat, minfo.rhohat, minfo.what, minfo.nx, minfo.ny, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );

   //smoothing material models
   MPI_Barrier(comm2);
   material_smooth3(minfo.chat, minfo.nx, minfo.ny, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.rhohat, minfo.nx, minfo.ny, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.what, minfo.nx, minfo.ny, mrank, nRank, comm2);
    
   MPI_Barrier(comm2);
   rigidbound_material (minfo.chat, minfo.rhohat, minfo.what, minfo.nx, minfo.ny, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );

   //smoothing material models (twice)
   MPI_Barrier(comm2);
   material_smooth3(minfo.chat, minfo.nx, minfo.ny, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.rhohat, minfo.nx, minfo.ny, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.what, minfo.nx, minfo.ny, mrank, nRank, comm2);
    
   MPI_Barrier(comm2);
   rigidbound_material (minfo.chat, minfo.rhohat, minfo.what, minfo.nx, minfo.ny, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );

/*
   MPI_Barrier(comm2);
   char fn[500];
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_c_%06d.dat", mrank);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.chat, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.chat, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_rho_%06d.dat", mrank);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.rhohat, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.rhohat, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_wind_%06d.dat", mrank);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.what, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.what, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank, nRank, comm2);
*/
  
   MPI_Barrier(comm2);
   get_time(argv[1], &minfo, mrank, nRank, comm2);
    if(mrank==0) {
      printf("=============================================\n");
    }

   //source initiation
   MPI_Barrier(comm2);
   get_src(argv[1], &minfo, mrank, nRank, comm2);
    if(mrank==0) {
      printf("=============================================\n");
    }


   MPI_Barrier(comm2);
   init_src(&minfo, mrank, nRank, comm2);

   //rec allocation
   MPI_Barrier(comm2);
   get_rec(argv[1], &minfo, mrank, nRank, comm2);
    if(mrank==0) {
      printf("=============================================\n");
    }

   MPI_Barrier(comm2);
   get_img(argv[1], &finfo, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   get_line(argv[1], &finfo, mrank, nRank, comm2);

   MPI_Barrier(comm2);
   set_lineimg(&minfo, &finfo, mrank, nRank, comm2);

   //intialization
   MPI_Barrier(comm2);
   //initialize_u_zero(&minfo, mrank, nRank, comm2);
   initialize_u(&minfo, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   initialize_v(&minfo, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   initialize_w(&minfo, mrank, nRank, comm2);
   MPI_Barrier(comm2);

//   printf("***********MAIN TEST\n");
   rigidbound_var (minfo.u, minfo.v, minfo.w, minfo.nx, minfo.ny, mrank, nRank, comm2);
   MPI_Barrier(comm2);
   updateBoundary(minfo.u, minfo.nx, minfo.ny, ghostL, mrank, nRank, comm2 );

/*
   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_u_%06d.dat", 0);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, minfo.u, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, \
                    mrank, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_v_%06d.dat", 0);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, minfo.v, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, \
                    mrank, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_w_%06d.dat", 0);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, minfo.w, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, \
                    mrank, nRank, comm2);
*/

   MPI_Barrier(comm2);

   setup_sg(&minfo, mrank, nRank, comm2);
   MPI_Barrier(comm2);

//   printf("***********MAIN TEST\n");
/*
   if(mrank == 0) {
      wd_array1("etay0.txt", minfo.etay, minfo.ny);
      wd_array1("psiy0.txt", minfo.psiy, minfo.ny);
      wd_array1("phiy0.txt", minfo.phiy, minfo.ny);
      wd_array1("sigy0.txt", minfo.sigy, minfo.ny);
   }
   if(mrank == 1) {
      wd_array1("etay1.txt", minfo.etay, minfo.ny);
      wd_array1("psiy1.txt", minfo.psiy, minfo.ny);
      wd_array1("phiy1.txt", minfo.phiy, minfo.ny);
      wd_array1("sigy1.txt", minfo.sigy, minfo.ny);
   }
*/
/*
   wdbin_mpi_array1("etax.bin", minfo.etax, minfo.nx, minfo.nx_global, ghostL, mrank, nRank, comm2);
   wdbin_mpi_array1("psix.bin", minfo.psix, minfo.nx, minfo.nx_global, ghostL, mrank, nRank, comm2);
   wdbin_mpi_array1("phix.bin", minfo.phix, minfo.nx, minfo.nx_global, ghostL, mrank, nRank, comm2);
   wdbin_mpi_array1("sigx.bin", minfo.sigx, minfo.nx, minfo.nx_global, ghostL, mrank, nRank, comm2);

if(mrank == 0) {
   wdbin_array1("etay.bin", minfo.etay, minfo.ny);
   wdbin_array1("psiy.bin", minfo.psiy, minfo.ny);
   wdbin_array1("phiy.bin", minfo.phiy, minfo.ny);
   wdbin_array1("sigy.bin", minfo.sigy, minfo.ny);
}
*/
    //time integration
//    if(minfo.ang == 180) {
//    time_marching180( &minfo, &finfo, u, v, w, rhohat, chat, what );
//    } else {
    MPI_Barrier(comm2);
    time_marching( &minfo, &finfo, mrank, nRank, comm2 );
//    }
   } //if(MPI_COMM_NULL != comm2)
    //wd_vec("verify.txt",minfo.srcs[0].q);
//       ac2dr_out( 0, 1, &minfo, &finfo, u, v, w, rhohat, chat, what );

    
   MPI_Barrier(MPI_COMM_WORLD);
    if(mrank == 0) {
      printf("***********Simulation Done\n");
      t1 = time(0);
      printf("Elapsed time (s) : %f\n", difftime(t1, t0)*1.0);
    }

    MPI_Finalize();
//    printf("***********ac2dr_main.c is working\n");
//    free_ac2dr( &minfo, &finfo, rhohat, chat, what );
//    free_array2( u );
//    free_array2( v );
//    free_array2( w );
    return 0;
}
