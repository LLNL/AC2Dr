#include "ac2dr_aux.h"

dvec* alloc_dvec( unsigned long length ) {
    dvec *m;
    m = (dvec *) malloc( sizeof(dvec) );
    m->elements = (double *) malloc( length * sizeof(double) );
    m->length = length;
    return m;
}

void free_dvec( dvec  *m )
{
    free((double *) m->elements);
    free((dvec *) m);
    return;
}

double getval_dvec( dvec* m, unsigned long xi) {
    return(m->elements[xi]);
}

dmat* alloc_dmat( unsigned long nrows, unsigned long ncols) {
    dmat *m;
    m = (dmat *) malloc( sizeof(dmat) );
    m->elements = (double *) malloc( nrows*ncols*sizeof(double) );
    m->nrows = nrows;
    m->ncols = ncols;
    return m;
}

void free_dmat( dmat* m ) {
    free((double *) m->elements);
    free((dmat *) m);
    return;
}

double getval_dmat( dmat* m, unsigned long xi, unsigned long yi ) {
    unsigned long idx;
    idx = xi + m->nrows*yi;
    return(m->elements[idx]);
}

//
double** alloc_array2( unsigned long nrows, unsigned long ncols) {
    unsigned long i;
    double **m;
    m = (double **) malloc ( nrows * sizeof(double *) );
    m[0] = (double *) malloc ( nrows * ncols * sizeof(double) );
    for(i=1;i<nrows;i++) {
        m[i] = m[0] + ncols * i;
    }
    return m;
}

void free_array2( double **m ) {
    free((double *) m[0]);
    free((double **) m);
    return;
}


//copy from n to m
void copy_dvec( dvec* m, dvec* n ) {
    unsigned long i;
    for(i=0;i<m->length;i++) {
       m->elements[i] = n->elements[i]; 
    }
    return;
}

void copy_mat( dmat* m, dmat* n ) {
    unsigned long i;
    for(i=0;i<(m->nrows*m->ncols);i++) {
       m->elements[i] = n->elements[i]; 
    }
    return;
}

double **rd_array2( char *fn )
{
    FILE *ptr_file;
    char buf[1000];
    long num_rows, num_cols, indx;
    double **BB;

    ptr_file =fopen(fn,"r");

    fgets(buf,1000,ptr_file);
    num_rows = atol( buf );
    fgets(buf,1000,ptr_file);
    num_cols = atol( buf );

    BB = alloc_array2( num_rows, num_cols );
    indx=0;
    while (fgets(buf,1000, ptr_file)!=NULL)
    {
        //printf("%s",buf);
        BB[0][indx] = atof( buf );
        indx++;
    }
    fclose(ptr_file);
    return BB;
}



dmat *rd_mat( char *fn )
{
    FILE *ptr_file;
    char buf[1000];
    long num_rows, num_cols, indx;
    dmat *BB;

    ptr_file =fopen(fn,"r");

    fgets(buf,1000,ptr_file);
    num_rows = atol( buf );
    fgets(buf,1000,ptr_file);
    num_cols = atol( buf );

    BB = (dmat *) alloc_dmat( num_rows, num_cols );
    indx=0;
    while (fgets(buf,1000, ptr_file)!=NULL)
    {
        //printf("%s",buf);
        BB->elements[indx] = atof( buf );
        indx++;
    }
    fclose(ptr_file);
    return BB;
}


dvec *rd_vec( char *fn )
{
    FILE *ptr_file;
    char buf[1000];
    long indx;
    dvec *d_vec;

    ptr_file =fopen(fn,"r");

    indx=0;
    while (fgets(buf,1000, ptr_file)!=NULL)
    {
        indx++;
    }

    d_vec = (dvec *) alloc_dvec( indx );
    d_vec->length = indx;

    indx = 0;
    rewind( ptr_file );

    while (fgets(buf,1000, ptr_file)!=NULL)
    {
        d_vec->elements[indx]=atof( buf );
        indx++;
    }

    fclose(ptr_file);
    return d_vec;
}

void wd_rec( char *fn, dvec *d_vec, fdmesh *minfo )
{
    FILE *ptr_file;
    long indx;

    ptr_file =fopen(fn,"w");

    fprintf(ptr_file,"%.17g\n",minfo->dt);
    for(indx=0 ; indx < d_vec->length ; indx++)
    {
        fprintf(ptr_file,"%.17g\n",d_vec->elements[indx]);
    }

    fclose(ptr_file);
    return;
}

void wd_recv2( char *fn, dvec *d_vec, double rx, double ry, fdmesh *minfo )
{
    FILE *ptr_file;
    long indx;

    ptr_file =fopen(fn,"w");

    fprintf(ptr_file,"%.17g\n",rx);
    fprintf(ptr_file,"%.17g\n",ry);
    fprintf(ptr_file,"%.17g\n",minfo->dt);
    fprintf(ptr_file,"%ld\n",d_vec->length);
    for(indx=0 ; indx < d_vec->length ; indx++)
    {
        fprintf(ptr_file,"%.17g\n",d_vec->elements[indx]);
    }

    fclose(ptr_file);
    return;
}

void wdbin_rec( char *fn, dvec *d_vec, fdmesh *minfo )
{
    FILE *ptr_file;
    long indx;

    ptr_file =fopen(fn,"wb");

    //fprintf(ptr_file,"%.17g\n",minfo->dt);
    fwrite(&minfo->dt, sizeof(double), 1, ptr_file);
    for(indx=0 ; indx < d_vec->length ; indx++)
    {
        //fprintf(ptr_file,"%.17g\n",d_vec->elements[indx]);
        fwrite(&d_vec->elements[indx], sizeof(double), 1, ptr_file);
    }

    fclose(ptr_file);
    return;
}

void wdbin_recv2( char *fn, dvec *d_vec, double rx, double ry, fdmesh *minfo )
{
    FILE *ptr_file;
    long indx;
    int f1=1, tlen;

    tlen = d_vec->length;

    ptr_file =fopen(fn,"wb");

    //fprintf(ptr_file,"%.17g\n",minfo->dt);
    fwrite(&f1, sizeof(int), 1, ptr_file);
    fwrite(&rx, sizeof(double), 1, ptr_file);
    fwrite(&ry, sizeof(double), 1, ptr_file);
    fwrite(&minfo->dt, sizeof(double), 1, ptr_file);
    fwrite(&tlen, sizeof(int), 1, ptr_file);
    for(indx=0 ; indx < d_vec->length ; indx++)
    {
        //fprintf(ptr_file,"%.17g\n",d_vec->elements[indx]);
        fwrite(&d_vec->elements[indx], sizeof(double), 1, ptr_file);
    }

    fclose(ptr_file);
    return;
}


void wd_vec( char *fn, dvec *d_vec )
{
    FILE *ptr_file;
    long indx;

    ptr_file =fopen(fn,"w");

    for(indx=0 ; indx < d_vec->length ; indx++)
    {
        fprintf(ptr_file,"%g\n",d_vec->elements[indx]);
    }

    fclose(ptr_file);
    return;
}

void wd_mat( char *fn, dmat *d_mat )
{
    FILE *ptr_file;
    long indx;

    ptr_file =fopen(fn,"w");

    fprintf(ptr_file,"%lu\n",d_mat->nrows);
    fprintf(ptr_file,"%lu\n",d_mat->ncols);
    for(indx=0 ; indx < d_mat->nrows * d_mat->ncols ; indx++)
    {
        fprintf(ptr_file,"%g\n",d_mat->elements[indx]);
    }

    fclose(ptr_file);
    return;
}

void wd_array2( char *fn, double **m, unsigned long nrows, unsigned long ncols )
{
    FILE *ptr_file;
    unsigned long indx;

    ptr_file =fopen(fn,"w");

    fprintf(ptr_file,"%lu\n",nrows);
    fprintf(ptr_file,"%lu\n",ncols);
    for(indx=0 ; indx < nrows * ncols ; indx++)
    {
        fprintf(ptr_file,"%.24g\n",m[0][indx]);
    }

    fclose(ptr_file);
    return;
}

void wdbin_array2( char *fn, double **m, unsigned long nrows, unsigned long ncols )
{
    FILE *ptr_file;
    unsigned long indx;

    ptr_file =fopen(fn,"wb");

//    fprintf(ptr_file,"%lu\n",nrows);
//    fprintf(ptr_file,"%lu\n",ncols);
    fwrite(&nrows, sizeof(int), 1, ptr_file);
    fwrite(&ncols, sizeof(int), 1, ptr_file);
    for(indx=0 ; indx < nrows * ncols ; indx++)
    {
        //fprintf(ptr_file,"%.24g\n",m[0][indx]);
        fwrite(&m[0][indx], sizeof(double), 1, ptr_file);
    }

    fclose(ptr_file);
    return;
}


void wd_array1( char *fn, double *m, unsigned long n)
{
    FILE *ptr_file;
    unsigned long indx;

    ptr_file =fopen(fn,"w");

    for(indx=0 ; indx < n ; indx++)
    {
        fprintf(ptr_file,"%.17g\n",m[indx]);
    }

    fclose(ptr_file);
    return;
}

void wdbin_array1( char *fn, double *m, unsigned long n)
{
    FILE *ptr_file;
    unsigned long indx;

    ptr_file =fopen(fn,"wb");

    fwrite(&n, sizeof(int), 1, ptr_file);
    for(indx=0 ; indx < n ; indx++)
    {
        //fprintf(ptr_file,"%.24g\n",m[0][indx]);
        fwrite(&m[indx], sizeof(double), 1, ptr_file);
    }

    fclose(ptr_file);
    return;
}


void elac_error( char  *msg, int   code  )
{
    fprintf(stderr, "%s\n", msg);
    exit(code);
}

void free_fdmesh( fdmesh *minfo ) {
    int i;
    free_dvec( minfo->ts ); 
    for(i=0;i<minfo->src_num;i++) {
        free_dvec( minfo->srcs[i].q );
    }
    for(i=0;i<minfo->sta_num;i++) {
        free_dvec( minfo->stas[i].pout );
        free_dvec( minfo->stas[i].rhout );
        free_dvec( minfo->stas[i].uout );
        free_dvec( minfo->stas[i].wout );
    }
    return;
}


void free_ac2dr( fdmesh *minfo, innout *finfo, double **rhohat, double **chat, \
                  double **what) {
    free_array2( rhohat );
    free_array2( chat );
    free_array2( what );
    free_fdmesh ( minfo );
    return;
}


double c6smoothbump( double t, double t0, double w ) {
    double val;

    if(t<t0) { val = 0; }
    if((t>=t0) & (t<=(t0+1/w))) {
        val = 51480*pow(w,7)*pow(t-t0,7)*pow(1-w*(t-t0),7);
    }
    if(t>(t0+1/w)) { val = 0; }
    return val;
}

double c6smoothbump3d( double t, double t0, double w ) {
    double val;
    if(t<t0) { val = 0; }
    if((t>=t0) & (t<=(t0+1/w))) {
        val = 51480*7*pow(w,7)*(30*pow(t-t0,4)*pow(1-w*(t-t0),7) +
                -3*42*w*pow(t-t0,5)*pow(1-w*(t-t0),6) +
                +3*42*pow(w,2)*pow(t-t0,6)*pow(1-w*(t-t0),5) +
                -30*pow(w,3)*pow(t-t0,7)*pow(1-w*(t-t0),4));
    }
    if(t>(t0+1/w)) { val = 0; }

    return val;
}


void read_2col(char *fn, double **xx, double **yy, int *n) {
   FILE *ptr_file;
   unsigned long i=0;
   double *aa, *bb;
   char *pch, buf[1000], buf2[1000], buf0[1000];

   ptr_file=fopen(fn, "r");
   if (!ptr_file) {
      sprintf(buf0, "     --- Error: The input file (%s) does not exist.\n", fn);
      elac_error(buf0, -1);
   }

   while ( fgets( buf, 1000, ptr_file )!=NULL ) {
      i++;
   }

   rewind(ptr_file);

   *n = i;
   aa = (double *) malloc(i*sizeof(double));
   bb = (double *) malloc(i*sizeof(double));
   i=0;

   while ( fgets( buf, 1000, ptr_file )!=NULL ) {
      strcpy( buf2, buf );
      pch = strtok(buf2, " \n");
      aa[i] = atof(pch); 
      //printf("%f ", aa[i]);
      pch = strtok(NULL, " \n");
      bb[i] = atof(pch);
      //printf("%f\n", bb[i]);
      pch = strtok(NULL, " \n");
      i++;
   }
   fclose( ptr_file );

   *xx=aa;
   *yy=bb;
}

void readbin_2col(char *fn, double **xx, double **yy, int *n) {
   FILE *ptr_file;
   unsigned long i=0;
   double *aa, *bb, val, N;
   char buf0[1000];
   int f1;
   int NeedSwap = 0;

   ptr_file=fopen(fn, "rb");
   if (!ptr_file) {
      sprintf(buf0, "     --- Error: The input file (%s) does not exist.\n", fn);
      elac_error(buf0, -1);
   }

   fread(&f1, sizeof(int), 1, ptr_file);
   while ( fread(&val, sizeof(double), 1, ptr_file )>0 ) {
      i++;
   }
   fclose( ptr_file );

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

   //printf("fread N = %ld\n", i);
   N = i/2;
   *n = i/2;
   aa = (double *) malloc(N*sizeof(double));
   bb = (double *) malloc(N*sizeof(double));
   i=0;

   ptr_file=fopen(fn, "rb");
   fread(&f1, sizeof(int), 1, ptr_file);
   for(i=0;i<N;i++) {
      fread(&val, sizeof(double), 1, ptr_file );
      if(NeedSwap) SwapBytes( &val, 8 );
      aa[i] = val; 
      fread(&val, sizeof(double), 1, ptr_file );
      if(NeedSwap) SwapBytes( &val, 8 );
      bb[i] = val;
      //printf("%f %f\n", aa[i], bb[i]);
   }
   fclose( ptr_file );

   *xx=aa;
   *yy=bb;
}


void write_2col(char *fn, double *xx, double *yy, int imax)
{
   FILE *ptr_file;
   int i = 0;
   ptr_file = fopen(fn, "w");
   while(i<imax)
   {
      fprintf(ptr_file, "%f %f\n", xx[i], yy[i]);
      i++;
   }
   fclose( ptr_file );
}

void locate(double *xx, unsigned long n, double x, long *j) {
   unsigned long ju, jm, jl;
   int ascnd;

   jl=0;
   ju=n+1;
   ascnd=(xx[n-1] >= xx[1-1]);
   while(ju-jl>1) {
      jm=(ju+jl) >> 1;
      if((x>=xx[jm-1]) == ascnd)
         jl=jm;
      else
         ju=jm;
   }
   if (x==xx[1-1]) *j=1-1;
   else if(x==xx[n-1]) *j=n-1-1;
   else *j=jl-1;
}

void polint( double *xa, double *ya, int n, double x, double *y, double *dy ) {
   int i, m, ns=1;
   double den, dif, dift, ho, hp, w;
   double *c, *d;

   dif=fabs(x-xa[0]);
   c=(double *) malloc(n*sizeof(double));
   d=(double *) malloc(n*sizeof(double));
   for(i=1;i<=n;i++) {
      if((dift=fabs(x-xa[i-1]))<dif) {
         ns=i-1;
         dif=dift;
      }
      c[i-1]=ya[i-1];
      d[i-1]=ya[i-1];
   }
   *y=ya[-1+ns--];
   for(m=1;m<n;m++) {
      for(i=1;i<=n-m;i++) {
         ho=xa[i-1]-x;
         hp=xa[i+m-1]-x;
         w=c[i+1-1]-d[i-1];
         if((den=ho-hp)==0.0) elac_error("Error in routine polint",-1);
         den=w/den;
         d[i-1]=hp*den;
         c[i-1]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1-1] : d[-1+ns--]));
   }
    free((double *) c);
    free((double *) d);
}

void polin2(double *x1a, double *x2a, double **ya, int m, int n, double x1, double x2, \
   double *y, double *dy) {
   int j;
   double *ymtmp;

   ymtmp=(double *) malloc(m*sizeof(double));
   for(j=1;j<=m;j++) {
      polint(x2a,&ya[j-1][0],n,x2,&ymtmp[j-1],dy);
   }
   polint(x1a,ymtmp,m,x1,y,dy);
   free((double *) ymtmp);
}

void filter6 ( double **u, double **u1, fdmesh *minfo, double d0 ) {

   int i,j;
   int nx, ny;
   double dux, duy;
   double dj[7] = {-1.0/64.0, 3.0/32.0, -15.0/64.0, 5.0/16.0, -15.0/64.0, 3.0/32.0, -1.0/64.0};
   nx = minfo->nx;
   ny = minfo->ny;

   for(j=3; j<(ny-3); j++) {
      for(i=3; i<(nx-3); i++) {
         dux = dj[0]*u[i-3][j] + dj[1]*u[i-2][j] + dj[2]*u[i-1][j] + dj[3]*u[i][j] + \
              dj[4]*u[i+1][j] + dj[5]*u[i+2][j] + dj[6]*u[i+3][j];
         duy = dj[0]*u[i][j-3] + dj[1]*u[i][j-2] + dj[2]*u[i][j-1] + dj[3]*u[i][j] + \
              dj[4]*u[i][j+1] + dj[5]*u[i][j+2] + dj[6]*u[i][j+3];

         u1[i][j] = u[i][j] - d0*dux - d0*duy;
      }
   }

}

void filter12 ( double **u, double **u1, fdmesh *minfo, double d0 ) {
   int i,j;
   int nx, ny;
   double dux, duy;
   double dj[13] = {1.0/4096.0, -3.0/1024.0, 33.0/2048.0, -55.0/1024.0, 495.0/4096.0, -99.0/512.0,231.0/1024.0, \
                    -99.0/512.0, 495.0/4096.0, -55.0/1024.0, 33.0/2048.0, -3.0/1024.0, 1.0/4096.0};
   nx = minfo->nx;
   ny = minfo->ny;

   for(j=6; j<(ny-6); j++) {
      for(i=6; i<(nx-6); i++) {
         dux = dj[0]*u[i-6][j] + dj[1]*u[i-5][j] + dj[2]*u[i-4][j] + dj[3]*u[i-3][j] + dj[4]*u[i-2][j] + dj[5]*u[i-1][j] + dj[6]*u[i][j] + \
               dj[7]*u[i+1][j] + dj[8]*u[i+2][j] + dj[9]*u[i+3][j] + dj[10]*u[i+4][j] + dj[11]*u[i+5][j] + dj[12]*u[i+6][j];
         duy = dj[0]*u[i][j-6] + dj[1]*u[i][j-5] + dj[2]*u[i][j-4] + dj[3]*u[i][j-3] + dj[4]*u[i][j-2] + dj[5]*u[i][j-1] + dj[6]*u[i][j] + \
               dj[7]*u[i][j+1] + dj[8]*u[i][j+2] + dj[9]*u[i][j+3] + dj[10]*u[i][j+4] + dj[11]*u[i][j+5] + dj[12]*u[i][j+6];
         u1[i][j] = u[i][j] - d0*dux - d0*duy;
      }
   }

}


void SwapBytes(void *pv, size_t n) {
   char *p = pv;
   size_t lo, hi;
   char tmp;
   for(lo=0, hi=n-1; hi>lo; lo++, hi--) {
      tmp = p[lo];
      p[lo]=p[hi];
      p[hi]=tmp;
   }
}
