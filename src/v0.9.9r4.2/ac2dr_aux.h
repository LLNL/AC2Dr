#ifndef _AC2DR_
#define _AC2DR_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

extern const double gam, R0, M0, Rd;

//Declaration of structures
typedef struct double_vector {
    unsigned long length;
    double  *elements;
} dvec;

typedef struct double_matrix {
    unsigned long nrows, ncols;
    double  *elements;
} dmat;

typedef struct msource {
    double  x,y,freq,t0,p0;
    unsigned long xi, xi_global, yi;
    char    type[1000];
    dvec    *q, *midq;
} msource;

typedef struct fsource {
    char    name[100];
    double  x,y,freq,t0,p0;
    unsigned long xi, yi;
    char    type[1000];
    dvec    *fx, *fy;
} fsource;


typedef struct station {
    char    name[100];
    double  x,y;
    unsigned long xi, yi, xi_global;
    dvec    *pout, *rhout, *uout, *wout;
    char mode[10];
    int format;
} station;

typedef struct image {
    char    file[1000];
    char    mode[100];
    double  imgdt;
    double  imgt;
    int     format, ghostL;
} image;

typedef struct lineout {
    char    file[1000];
    char    mode[100];
    double  elev;
    int  yi, tl;
    double  linedt;
    double  **lineimg;
    int     format;
} lineout;


typedef struct innout {
    char indir[1000], outdir[1000];
    char inputf[1000];
    char logf[1000];
    //char presf[1000];
    //char tempf[1000];
    //char windxf[1000];
    //char windyf[1000];
    char c_fn[1000];
    char rho_fn[1000];
    char w_fn[1000];
    int  img_num, line_num;
    double *imgdt;
    int c_format, rho_format, w_format;
    image *imgs;
    lineout *lines;
} innout;

typedef struct fdmesh {
    //global
    unsigned long nth_global, nx_global;
    //local 
    unsigned int nr, ny, nth, nx, nt, nRank;
    unsigned int xi_start, xi_end;
    //axis
    double *xx, *yy, *xx_deg;
    //local variables
    double **u, **v, **w;
    double **pmax, **pmaxabs, **pmax_t, **pmaxabs_t;
    double **rhohat, **chat, **what;
    double R, elev_max, rpmin, rpmax, ang, th, dth, dr, angmin, angmax;
    double ang_global, th_global, angmin_global, angmax_global, thmin_global, thmax_global;
    //double R=6371*1000, elev_max, rpmin, rpmax, smin, smax, th, ds, dr;
    double T, dt, cfl;
    dvec   *ts;
    int src_num, src_num_global, sta_num, sta_num_global;
    //double phmin, phmax, thmin, thmax, rhomin, rhomax, cmin, cmax, uhmin, uhmax, whmin, whmax;
    double rhomin, rhomax, cmin, cmax, whmin, whmax;
    //iteration of moving average
    int c_smooth, rho_smooth, w_smooth;
    int sgrid_length;
    double sgrid_dc, sgrid_EL, alpha;
    double *psix, *psiy, *phix, *phiy, *sigx, *sigy, *etax, *etay;
    double *psix_local, *phix_local, *sigx_local, *etax_local;
    double *at, ax;
    double c_val, rho_val, w_val;
    double *c_prof, *rho_prof, *w_prof;
//    double *c_prof0, *rho_prof0, *w_prof0;
//    double **chat0, **what0, **rhohat0;
//    double *xx_c0, *xx_w0, **xx_rho0;
//    double *yy_c0, *yy_w0, **yy_rho0;
    int c_type, rho_type, w_type, effc;
    char logf[1000];
    msource *srcs;
    station *stas;
    station *stas0;
    lineout *louts;
} fdmesh;

//ac2dr_aux.c
double** alloc_array2( unsigned long, unsigned long );
void free_array2( double ** );
double **rd_array2( char * );
void wd_array1( char *fn, double *m, unsigned long n);
void wdbin_array1( char *fn, double *m, unsigned long n);
void wd_array2( char *, double **, unsigned long, unsigned long );
void wdbin_array2( char *fn, double **m, unsigned long nrows, unsigned long ncols );
dvec* alloc_dvec( unsigned long );
void free_dvec( dvec * );
double getval_dvec( dvec *, unsigned long );
dmat* alloc_dmat( unsigned long, unsigned long );
void free_dmat( dmat * );
double getval_dmat( dmat *, unsigned long, unsigned long );
void copy_dvec( dvec *, dvec * );
void copy_mat( dmat *, dmat * );
dmat* rd_mat( char * );
dvec* rd_vec( char * );
void wd_rec( char *, dvec *, fdmesh *);
void wd_recv2( char *fn, dvec *d_vec, double rx, double ry, fdmesh *minfo );
void wdbin_rec( char *fn, dvec *d_vec, fdmesh *minfo );
void wdbin_recv2( char *fn, dvec *d_vec, double rx, double ry, fdmesh *minfo );
void wd_vec( char *, dvec * );
void wd_mat( char *, dmat * );
void elac_error( char *, int );
void free_ac2dr( fdmesh *, innout *, double **, double **, double **);
void free_fdmesh( fdmesh * );
double c6smoothbump( double, double, double );
double c6smoothbump3d( double, double, double );

//ac2dr_input.c
//void ac2dr_input( char *, fdmesh *, innout *, double ***, double ***, double ***, int);
void get_path( char *fn, innout *finfo, int );
void get_grid( char *, fdmesh *, int );
void get_vel( char *, innout *, fdmesh *, int );
void get_dens( char *, innout *, fdmesh *, int );
void get_wind( char *, innout *, fdmesh *, int );

void get_time( char *, fdmesh *, int, int, MPI_Comm );
void get_src( char *, fdmesh *, int, int, MPI_Comm );
void get_rec( char *, fdmesh *, int, int, MPI_Comm );
void get_img( char *, innout *, int, int, MPI_Comm);
void get_line( char *, innout *, int, int, MPI_Comm);
void set_lineimg( fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm2 );
void print_message( char *, char *, int );
void print_messagef( char *, char *, int );
void print_message_mpi ( char *, MPI_Comm comm, char *msg, int nc, int mrank, int nRank, int );
void print_messagef_mpi ( char *, MPI_Comm comm, char *msg, int nc, int mrank, int nRank, int );


//ac2dr_RK4

void ac2dr_RK4 ( fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm2 );
void ac2dr_RK4_damping ( fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm2 );
void ac2dr_FD4 ( fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm2 );
void ac2dr_RK4180 ( fdmesh *, innout *, double **, double **, double **, double **, \
                  double **, double **);
void dqdt2_sym_ghost_sg ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage );
void dqdt4_sym_ghost_sg ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage );
void dqdt4 ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage );
void dqdt12 ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage );
void dqdt6 ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage );
void dqdt6_damping ( double **u, double **v, double **w, double **u1, double **v1, double **w1, \
            double **rhohat, double **what, double **chat, fdmesh *minfo, int tk, int RK4stage );
void getSBP6_in( double * );
void getSBP6_left( double [][9] );
void getSBP6_right( double [][9] );

void getsymmetric_left( int [][6] );
void getsymmetric_right( int [][6], unsigned long );
void getantisymmetric_left( int iLB[][6] );
void getantisymmetric_right( int iRB[][6] );
void bcenforce_rigid_bottom( double **v, fdmesh *minfo );
void bcenforce_rigid_top( double **v, fdmesh *minfo );
void bcenforce_rigid_right( double **v, fdmesh *minfo );
void bcenforce_rigid_left( double **v, fdmesh *minfo );
void bcenforce_rigid_bc1( double **u, double **v, double **w, fdmesh *minfo );
void bcenforce_rigid_bc3( double **u, double **v, double **w, fdmesh *minfo );
void bcenforce_sg_bc3( double **u, double **v, double **w, fdmesh *minfo );
void bcenforce_rigid_material3( double **chat, double **rhohat, double **what, fdmesh *minfo );

void initialize_u( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 );
void initialize_u_zero( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 );
void initialize_v( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 );
void initialize_w( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 );
void init_src( fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 );
double point_src( double A, double R, double w, double c );


//ac2dr_src.c
double get_psi( double );

//ac2dr_output
void ac2dr_out( int , unsigned long , fdmesh *, innout *, \
                 double **, double **, double **, double **, \
                 double **, double **, int, int, MPI_Comm );
void ac2dr_initout(fdmesh *, innout *, \
                 double **, double **, double **, double **, \
                 double **, double **, int, int, MPI_Comm );


//supergrid
double dampx4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampux4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampuy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampvx4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampvy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampwx4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
double dampwy4( unsigned long i, unsigned long j, double **u, fdmesh *minfo );
void setup_sg (fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2);


void bcenforce_symmetric_left3( double **w, fdmesh *minfo );
void bcenforce_symmetric_right3( double **w, fdmesh *minfo );
void bcenforce_symmetric_bottom3( double **w, fdmesh *minfo );
void bcenforce_symmetric_top3( double **w, fdmesh *minfo );
void bcenforce_antisymmetric_left3( double **w, fdmesh *minfo );
void bcenforce_antisymmetric_right3( double **w, fdmesh *minfo );
void bcenforce_antisymmetric_bottom3( double **w, fdmesh *minfo );
void bcenforce_antisymmetric_top3( double **w, fdmesh *minfo );

//material
void material_prop(fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm);
void material_effc(fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm);
void what_adjust(fdmesh *minfo, int mrank, int nRank);
void material_intp_2d (int nx0, int ny0, double *xx0, double *yy0, double **zz0, \
                       int nx1, int ny1, double *xx1, double *yy1, double **zz1);

void write_2col(char *fn, double *xx, double *yy, int n);
void writebin_2col(char *fn, double *xx, double *yy, int imax);
void read_2col(char *fn, double **xx, double **yy, int *n);
void readbin_2col(char *fn, double **xx, double **yy, int *n);

void writetxt_2dfile(char *fn, int nx, int ny, double *xx, double *yy, double **zz);
void writebin_2dfile(char *fn, int nx, int ny, double *xx, double *yy, double **zz);
void read_2dfile(char *fn, int *nx, int *ny, double **xx, double **yy, double ***val);
void readbin_2dfile(char *fn, int *nx, int *ny, double **xx, double **yy, double ***val);

void polint(double *xa, double *ya, int n, double x, double *y, double *dy );
void polin2(double *x1a, double *x2a, double **ya, int m, int n, double x1, double x2, \
            double *y, double *dy);
void locate(double *xx, unsigned long n, double x, long *j);


//mpi

void decompose_nproc (fdmesh *minfo, innout *finfo, int mrank, int totalp, int *nRank, \
                      int *xi_start, int *xi_end, int *nx_local);
void decompose_domain (fdmesh *minfo, innout *finfo, int mrank, int totalp, int nRank, \
                       int *xi_start, int *xi_end, int *nx_local, MPI_Comm comm);
void material_minmax(fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2);
void updateBoundary ( double **u, int nx, int ny, int byN, int mrank, int nRank, MPI_Comm comm2) ;
void shift_domain ( double **u, int nx, int ny, int byN, int mrank, int nRank, MPI_Comm comm2);
void rigidbound_material ( double **chat, double **rhohat, double **what, int nx, int ny, \
                           int mrank, int nRank, MPI_Comm comm2);
void rigidbound_var ( double **u, double **v, double **w, int nx, int ny, int mrank, int nRank, MPI_Comm comm2);
void material_smooth3 (double **u, int nx, int ny, int mrank, int nRank, MPI_Comm comm2);


//mpi_aux
void mpi_printf ( MPI_Comm comm, char *msg, int nc, int mrank, int nRank );
void wdbin_mpi_array2( char *fn, double **m, int nx, int ny, int nx_global, int ny_global, \
   double dx, double dy, int mrank, int nRank, MPI_Comm comm );
void wdbin_mpi_array2_reduced( char *fn, double **m, int nx, int ny, int nx_global, int ny_global, \
   double dx, double dy, int mrank, int nRank, MPI_Comm comm );
void wdbin_mpi_array1( char *fn, double *m, int nx, int nx_global, int nb, \
   int mrank, int nRank, MPI_Comm comm );

void wdbin_mpi_2dimg( char *fn, double **m, int f1, double imgt, \
   int nx, int ny, int nx_global, int ny_global, double dx, double dy, \
   int mrank, int nRank, MPI_Comm comm );

void wdbin_mpi_2dimg_ngl( char *fn, double **m, int f1, double imgt, \
   int nx, int ny, int nx_global, int ny_global, double dx, double dy, \
   int mrank, int nRank, MPI_Comm comm );

void wdbin_mpi_2dmat( char *fn, double **m, int f1, \
   int nx, int ny, int nx_global, int ny_global, double dx, double dy, \
   int mrank, int nRank, MPI_Comm comm );

void wdbin_mpi_2dmat_ngl( char *fn, double **m, int f1, \
   int nx, int ny, int nx_global, int ny_global, double dx, double dy, \
   int mrank, int nRank, MPI_Comm comm );


//filter
void filter6 ( double **u, double **u1, fdmesh *minfo, double d0 );
void filter12 ( double **u, double **u1, fdmesh *minfo, double d0 );

//Auxiliary

void SwapBytes(void *pv, size_t n);







#endif
