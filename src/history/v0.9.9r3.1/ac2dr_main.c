#include "ac2dr_aux.h"
#include "ac2dr_config.h"

//const double gam = 1.4;
//const double R0=8.3143;
//const double M0=0.029;
//const double Rd=8.3143/0.029;
const char ver0[] = "0.9.9r3.1";
const char rdate[] = "December 18th, 2023";
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
    int mrank, mrank2, totalp, totalp2, nRank;
    //finite-difference mesh
    fdmesh minfo;
    //input & output information
    innout finfo;
    time_t t0, t1;
    int *proc_ranks, proc;
//    int i, j, nx, ny;
    char buf[200];
    char msg_buf[3000];
    int *xi_start, *xi_end, *nx_local;
    //FILE *ptr_file;
    int indx;
    char FBIN[]="binary", FASC[]="ascii";

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
      sprintf(msg_buf,"=============================================\n");
      printf(msg_buf);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    //strcpy(buf, "MPI initialized.");
    sprintf(buf, "MPI_initialized");
    mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    //print_message_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);

    MPI_Barrier(MPI_COMM_WORLD);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      printf(msg_buf);
    }

    MPI_Barrier(MPI_COMM_WORLD);

//    ac2dr_input( argv[1], &minfo, &finfo, &rhohat, &chat, &what, mrank); 
    get_path( argv[1], &finfo, mrank );
    strcpy(finfo.logf, finfo.outdir);
    strcat(finfo.logf, "/");
    strcat(finfo.logf, "AC2DR.LOG");
    strcpy(minfo.logf, finfo.logf);

    if(mrank==0) {
      //ptr_file = fopen( finfo.logf, "w" );
      sprintf(msg_buf, "======================================================================\n");
      print_messagef(finfo.logf, msg_buf, 1);
      sprintf(msg_buf, "=      AC2Dr version v%s                     \n", ver0);
      print_messagef(finfo.logf, msg_buf, 0);
      sprintf(msg_buf, "=                                            \n");
      print_messagef(finfo.logf, msg_buf, 0);
      sprintf(msg_buf, "=                    Release date: %s        \n", rdate);
      print_messagef(finfo.logf, msg_buf, 0);
      sprintf(msg_buf, "=                    Contact: %s             \n", author);
      print_messagef(finfo.logf, msg_buf, 0);
      sprintf(msg_buf, "======================================================================\n");
      print_messagef(finfo.logf, msg_buf, 0);
    }

    if(mrank==0) {
      sprintf(msg_buf,"=============================================\n");
      print_messagef(finfo.logf, msg_buf, 0);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    //strcpy(buf, "MPI initialized.");
    sprintf(buf, "MPI_initialized");
    //mpi_printf(MPI_COMM_WORLD, buf, 200, mrank, totalp-1);
    print_messagef_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);

    MPI_Barrier(MPI_COMM_WORLD);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      print_messagef(finfo.logf, msg_buf, 0);
    }

    strcpy(buf, "Path reading.... Done.");
    print_message_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

    get_grid( argv[1], &minfo, mrank );
    strcpy(buf, "Grid reading.... Done.");
    print_message_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

    get_vel( argv[1], &finfo, &minfo, mrank );
    strcpy(buf, "Velocity reading.... Done.");
    print_message_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

    get_dens( argv[1], &finfo, &minfo, mrank );
    strcpy(buf, "Density reading.... Done.");
    print_message_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

    get_wind( argv[1], &finfo, &minfo, mrank );
    strcpy(buf, "Wind reading.... Done.");
    print_message_mpi(finfo.logf, MPI_COMM_WORLD, buf, 200, mrank, totalp-1, 0);
    if(mrank==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

    //domain decomposition

    MPI_Barrier(MPI_COMM_WORLD);
    xi_start=(int *) malloc (totalp * sizeof(int));
    xi_end=(int *) malloc (totalp * sizeof(int));
    nx_local=(int *) malloc (totalp * sizeof(int));

    decompose_nproc (&minfo, &finfo, mrank, totalp, &nRank, xi_start, xi_end, nx_local);

    MPI_Barrier(MPI_COMM_WORLD);

    proc_ranks = (int *) malloc (totalp * sizeof(int));
    for(proc = 0; proc <= nRank; proc++) {
      proc_ranks[proc] = proc;
    }

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Comm_group(MPI_COMM_WORLD, &group1);
   MPI_Group_incl(group1, nRank+1, proc_ranks, &group2);
   MPI_Comm_create(MPI_COMM_WORLD, group2, &comm2);

   if(comm2 != MPI_COMM_NULL) {
   MPI_Comm_rank(comm2, &mrank2);
   MPI_Comm_size(comm2, &totalp2);
   //printf("comm1 rank = %d/%d, comm2 rank = %d/%d\n", mrank, totalp, mrank2, totalp2);

   decompose_domain (&minfo, &finfo, mrank2, totalp, nRank, xi_start, xi_end, nx_local, comm2);

    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }
   strcpy(buf, "Domain decomposition... Done.");
   print_message_mpi(finfo.logf, comm2, buf, 200, mrank2, nRank, 0);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

   MPI_Barrier(comm2);
    //material specification
   material_prop(&minfo, &finfo, mrank2, nRank, comm2);

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
   print_message_mpi(finfo.logf, comm2, buf, 200, mrank2, nRank, 0);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );

   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   strcpy(buf, "Material boundary updated... Done.");
   print_message_mpi(finfo.logf, comm2, buf, 200, mrank2, nRank, 0);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

//   strcpy(fn, finfo.outdir);
//   strcat(fn, "/");
//   strcat(fn, "snapshot");
//   sprintf(buf, "_w0_%06d.dat", 0);
//   strcat(fn, buf);
//   wdbin_mpi_array2(fn, minfo.what, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);

//   mpi_printf(comm2, "MPI_test.", mrank2, nRank);
   //wind field adjustment
   MPI_Barrier(comm2);
   what_adjust(&minfo, mrank2, nRank);

   //check min, max, global min/max
   MPI_Barrier(comm2);
   material_minmax(&minfo, mrank2, nRank, comm2); 
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

   //matrix shift
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   shift_domain (minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   shift_domain (minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   shift_domain (minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );

   MPI_Barrier(comm2);
   rigidbound_material (minfo.chat, minfo.rhohat, minfo.what, minfo.nx, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );

   //smoothing material models
   /*
   MPI_Barrier(comm2);
   material_smooth3(minfo.chat, minfo.nx, minfo.ny, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.rhohat, minfo.nx, minfo.ny, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.what, minfo.nx, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   
   MPI_Barrier(comm2);
   rigidbound_material (minfo.chat, minfo.rhohat, minfo.what, minfo.nx, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );

   //smoothing material models (twice)
   MPI_Barrier(comm2);
   material_smooth3(minfo.chat, minfo.nx, minfo.ny, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.rhohat, minfo.nx, minfo.ny, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   material_smooth3(minfo.what, minfo.nx, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
    
   MPI_Barrier(comm2);
   rigidbound_material (minfo.chat, minfo.rhohat, minfo.what, minfo.nx, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   updateBoundary(minfo.chat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.rhohat, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
   MPI_Barrier(comm2);
   updateBoundary(minfo.what, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );
*/
/*
   MPI_Barrier(comm2);
   char fn[500];
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_c_%06d.dat", mrank2);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.chat, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.chat, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_rho_%06d.dat", mrank2);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.rhohat, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.rhohat, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_wind_%06d.dat", mrank2);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.what, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.what, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);
*/  
   MPI_Barrier(comm2);
   get_time(argv[1], &minfo, mrank2, nRank, comm2);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

   //source initiation
   MPI_Barrier(comm2);
   get_src(argv[1], &minfo, mrank2, nRank, comm2);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }


   MPI_Barrier(comm2);
   init_src(&minfo, mrank2, nRank, comm2);

   //rec allocation
   MPI_Barrier(comm2);
   get_rec(argv[1], &minfo, mrank2, nRank, comm2);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

   MPI_Barrier(comm2);
   get_img(argv[1], &finfo, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   get_line(argv[1], &finfo, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   set_lineimg(&minfo, &finfo, mrank2, nRank, comm2);

   //intialization
   MPI_Barrier(comm2);
   //initialize_u_zero(&minfo, mrank2, nRank, comm2);
   initialize_u(&minfo, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   initialize_v(&minfo, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   initialize_w(&minfo, mrank2, nRank, comm2);
   MPI_Barrier(comm2);

//   printf("***********MAIN TEST\n");
   rigidbound_var (minfo.u, minfo.v, minfo.w, minfo.nx, minfo.ny, mrank2, nRank, comm2);
   MPI_Barrier(comm2);
   updateBoundary(minfo.u, minfo.nx, minfo.ny, ghostL, mrank2, nRank, comm2 );

/*
   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_u_%06d.dat", 0);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, minfo.u, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, \
                    mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_v_%06d.dat", 0);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, minfo.v, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, \
                    mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_w_%06d.dat", 0);
   strcat(fn, buf);
   wdbin_mpi_array2(fn, minfo.w, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, \
                    mrank2, nRank, comm2);
*/

   MPI_Barrier(comm2);

   setup_sg(&minfo, mrank2, nRank, comm2);
   MPI_Barrier(comm2);

//   printf("***********MAIN TEST\n");
/*
   if(mrank2 == 0) {
      wd_array1("etay0.txt", minfo.etay, minfo.ny);
      wd_array1("psiy0.txt", minfo.psiy, minfo.ny);
      wd_array1("phiy0.txt", minfo.phiy, minfo.ny);
      wd_array1("sigy0.txt", minfo.sigy, minfo.ny);
   }
   if(mrank2 == 1) {
      wd_array1("etay1.txt", minfo.etay, minfo.ny);
      wd_array1("psiy1.txt", minfo.psiy, minfo.ny);
      wd_array1("phiy1.txt", minfo.phiy, minfo.ny);
      wd_array1("sigy1.txt", minfo.sigy, minfo.ny);
   }
*/
/*
   wdbin_mpi_array1("etax.bin", minfo.etax, minfo.nx, minfo.nx_global, ghostL, mrank2, nRank, comm2);
   wdbin_mpi_array1("psix.bin", minfo.psix, minfo.nx, minfo.nx_global, ghostL, mrank2, nRank, comm2);
   wdbin_mpi_array1("phix.bin", minfo.phix, minfo.nx, minfo.nx_global, ghostL, mrank2, nRank, comm2);
   wdbin_mpi_array1("sigx.bin", minfo.sigx, minfo.nx, minfo.nx_global, ghostL, mrank2, nRank, comm2);

if(mrank2 == 0) {
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

    if(mrank2==0) {
    //ptr_file = fopen( finfo.logf, "w" );
    print_message(finfo.logf, "=============================================\n", 0);
    print_message(finfo.logf, "=============================================\n", 0);
    print_message(finfo.logf, "Input Parameters ... \n", 0);
    print_message(finfo.logf, "=============================================\n", 0);
    print_message(finfo.logf, "=============================================\n", 0);

    sprintf(msg_buf,"path input=%s output=%s \n", finfo.indir, finfo.outdir);
    print_message(finfo.logf, msg_buf, 0);
    sprintf(msg_buf,"grid elevMax=%f angleMax=%f h=%f radius=%f nx=%ld ny=%d\n", \
                     minfo.elev_max, minfo.ang_global, minfo.dr, minfo.R/1000, minfo.nx_global, minfo.ny);
    print_message(finfo.logf, msg_buf, 0);
    sprintf(msg_buf,"time t=%f cfl=%f\n", minfo.T, minfo.cfl); 
    print_message(finfo.logf, msg_buf, 0);

    if(minfo.c_type==1) {
      sprintf(msg_buf,"mspeed value=%f \n", minfo.c_val); print_message(finfo.logf, msg_buf, 0);
    }
    if(minfo.c_type==2) {
      sprintf(msg_buf,"mspeed profile=%s \n", finfo.c_fn); print_message(finfo.logf, msg_buf, 0);
    }
    if(minfo.c_type==3) {
      sprintf(msg_buf,"mspeed 2dfile=%s \n", finfo.c_fn); print_message(finfo.logf, msg_buf, 0);
    }

    if(minfo.rho_type==1) {
      sprintf(msg_buf,"mdensity value=%f \n", minfo.rho_val); print_message(finfo.logf, msg_buf, 0);
    }
    if(minfo.rho_type==2) {
      sprintf(msg_buf,"mdensity profile=%s \n", finfo.rho_fn); print_message(finfo.logf, msg_buf, 0);
    }
    if(minfo.rho_type==3) {
      sprintf(msg_buf,"mdensity 2dfile=%s \n", finfo.rho_fn); print_message(finfo.logf, msg_buf, 0);
    }

    if(minfo.w_type==1) {
      sprintf(msg_buf,"wind value=%f \n", minfo.w_val); print_message(finfo.logf, msg_buf, 0);
    }
    if(minfo.w_type==2) {
      sprintf(msg_buf,"wind profile=%s \n", finfo.w_fn); print_message(finfo.logf, msg_buf, 0);
    }
    if(minfo.w_type==3) {
      sprintf(msg_buf,"wind 2dfile=%s \n", finfo.w_fn); print_message(finfo.logf, msg_buf, 0);
    }

    for(indx=0;indx<minfo.src_num;indx++) {
        sprintf(msg_buf,"asource elev=%f angle=%f p0=%f freq=%f type=%s xi=%ld yi=%ld\n", \
                minfo.srcs[indx].y, minfo.srcs[indx].x*180/M_PI, minfo.srcs[indx].p0, \
                minfo.srcs[indx].freq, minfo.srcs[indx].type, minfo.srcs[indx].xi_global, minfo.srcs[indx].yi);
        print_message(finfo.logf, msg_buf, 0);
    }

    for(indx=0;indx<minfo.sta_num_global;indx++) {
        if(minfo.stas0[indx].format==0) {
        sprintf(msg_buf,"rec name=%s elev=%g angle=%g mode=%s format=%s xi=%ld yi=%ld\n", \
                minfo.stas0[indx].name, minfo.stas0[indx].y, minfo.stas0[indx].x*180/M_PI, minfo.stas0[indx].mode,
                FBIN, minfo.stas0[indx].xi_global, minfo.stas0[indx].yi);
        } else {
        sprintf(msg_buf,"rec name=%s elev=%g angle=%g mode=%s format=%s xi=%ld yi=%ld\n", \
                minfo.stas0[indx].name, minfo.stas0[indx].y, minfo.stas0[indx].x*180/M_PI, minfo.stas0[indx].mode,
                FASC, minfo.stas0[indx].xi_global, minfo.stas0[indx].yi);
        }
        print_message(finfo.logf, msg_buf, 0);
    }

    for(indx=0;indx<finfo.img_num;indx++) {
        sprintf(msg_buf,"image timeInterval=%f mode=%s file=%s ghostL=%d\n", \
                finfo.imgs[indx].imgdt, finfo.imgs[indx].mode, finfo.imgs[indx].file, finfo.imgs[indx].ghostL);
        print_message(finfo.logf, msg_buf, 0);
    }
    } //mrank2==0




    MPI_Barrier(comm2);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }
    sprintf(buf, "Time Stepping Starts...");
    print_message_mpi(finfo.logf, comm2, buf, 200, mrank2, nRank, 0);

    time_marching( &minfo, &finfo, mrank2, nRank, comm2 );

    MPI_Barrier(comm2);
    sprintf(buf, "Time Stepping Ends");
    print_message_mpi(finfo.logf, comm2, buf, 200, mrank2, nRank, 0);
    if(mrank2==0) {
      sprintf(msg_buf, "=============================================\n");
      print_message(finfo.logf, msg_buf, 0);
    }

/*
   MPI_Barrier(comm2);
   char fn[500];
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_c_%06d.dat", mrank2);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.chat, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.chat, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_rho_%06d.dat", mrank2);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.rhohat, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.rhohat, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);

   MPI_Barrier(comm2);
   strcpy(fn, finfo.outdir);
   strcat(fn, "/");
   strcat(fn, "snapshot");
   sprintf(buf, "_wind_%06d.dat", mrank2);
   strcat(fn, buf);
   wdbin_array2(fn, minfo.what, minfo.nx, minfo.ny );
   //wdbin_mpi_array2(fn, minfo.what, minfo.nx, minfo.ny, minfo.nx_global, minfo.ny, mrank2, nRank, comm2);
*/
//    }
   } //if(MPI_COMM_NULL != comm2)
    //wd_vec("verify.txt",minfo.srcs[0].q);
//       ac2dr_out( 0, 1, &minfo, &finfo, u, v, w, rhohat, chat, what );

    
   MPI_Barrier(MPI_COMM_WORLD);
    if(mrank == 0) {
      sprintf(msg_buf, "***********Simulation Done\n");
      print_message(finfo.logf, msg_buf, 0);
      t1 = time(0);
      sprintf(msg_buf, "Elapsed time (s) : %f\n", difftime(t1, t0)*1.0);
      print_message(finfo.logf, msg_buf, 0);
    }

    MPI_Finalize();
//    printf("***********ac2dr_main.c is working\n");
//    free_ac2dr( &minfo, &finfo, rhohat, chat, what );
//    free_array2( u );
//    free_array2( v );
//    free_array2( w );
    return 0;
}
