#include <stdio.h>
#include <dirent.h>
#include "ac2dr_aux.h"
#include "ac2dr_config.h"

void get_path( char *fn, innout *finfo, int mrank)
{
    FILE    *ptr_file;
    char    buf[1000],buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file.\n", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "path" ) == 0 ) id_num++;
        }
    }

    if( id_num != 1 )
        elac_error("path is not defined or defined repeatedly.\n",-1);

    rewind( ptr_file );
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "path" ) == 0 ) {
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"input")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(finfo->indir, pch2);
                    }
                    if(strcmp(pch2,"output")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(finfo->outdir, pch2);
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );

    DIR *dip;
    if((dip=opendir(finfo->indir)) == NULL)
    {
        elac_error("     --- Error: Input directory does not exist.\n", -1);
    } else { closedir(dip); }
    if((dip=opendir(finfo->outdir)) == NULL)
    {
        elac_error("     --- Error: Output directory does not exist.\n", -1);
    } else { closedir(dip); }

    return;
}

//////////////////////////////////////////////////////////


void get_grid( char *fn, fdmesh *minfo, int mrank )
{
    FILE    *ptr_file;
    char    buf[1000],buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file.\n", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "grid" ) == 0 ) id_num++;
        }
    }

    if( id_num != 1 )
        elac_error("GRID is not defined or defined repeatedly.\n",-1);

    rewind( ptr_file );
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "grid" ) == 0 ) {
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"radius")==0) {
                        pch2=strtok(NULL,"=");
                        minfo->R = atof(pch2)*1000.0;
                    }
                    if(strcmp(pch2,"elev")==0) {
                        pch2=strtok(NULL,"=");
                        minfo->elev_max = atof(pch2);
                    }
                    if(strcmp(pch2,"angle")==0) {
                        pch2=strtok(NULL,"=");
                        minfo->ang_global = atof(pch2);
                        minfo->angmin_global = min_ang;
                        minfo->angmax_global = minfo->angmin_global + minfo->ang_global;

                    }
                    if(strcmp(pch2,"h")==0) {
                        pch2=strtok(NULL,"=");
                        minfo->dr = atof(pch2);
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );
   
    minfo->thmax_global=minfo->angmax_global/180.0*M_PI;
    minfo->thmin_global=minfo->angmin_global/180.0*M_PI;
    minfo->th_global=minfo->thmax_global-minfo->thmin_global;

    if(minfo->ang==180) {
    minfo->rpmin=0.0;
    minfo->rpmax=minfo->elev_max;
    minfo->nr=floor(minfo->rpmax/minfo->dr+0.5);
    minfo->dth=minfo->dr/(minfo->rpmax+minfo->R);
    minfo->nth_global=floor(minfo->th_global/minfo->dth);
    minfo->dth=(180.0)/(minfo->nth_global-1)*M_PI/180.0;
    } else {
    minfo->rpmin=0.0;
    minfo->rpmax=minfo->elev_max;
    minfo->dth=minfo->dr/(minfo->rpmax+minfo->R);
    minfo->nth_global=floor(minfo->th_global/minfo->dth+0.5);
    //minfo->dr = minfo->dr / 2.0;
    minfo->nr=floor(minfo->rpmax/minfo->dr+0.5);
    minfo->ny= minfo->nr;
    }

    minfo->nx_global=minfo->nth_global;
    if(mrank == 0) {
      //printf("R = %g\n",minfo->R);
      printf("Modeling Parameters:\n");
      printf("     --- dth degree = %g\n",minfo->dth);
      printf("     --- nth = %ld\n",minfo->nth_global);
      printf("     --- nr = %d\n",minfo->nr);
      printf("================================================\n");
    }
    return;
}

void get_vel( char *fn, innout *finfo, fdmesh *minfo, int mrank)
{
    FILE    *ptr_file;
    char    buf[1000],buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    minfo->c_type = 0;
    minfo->c_val = 0;
    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file.\n", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "mspeed" ) == 0 ) id_num++;
        }
    }

    if( id_num != 1 )
        elac_error("mspeed (sound speed) is not defined or defined repeatedly.\n",-1);

    rewind( ptr_file );
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "mspeed" ) == 0 ) {
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"value")==0) {
                        pch2=strtok(NULL,"=");
                        if(minfo->c_type!=0)
                           elac_error("mspeed: the sound speed is already defined.\n",-1);
                        minfo->c_val = atof(pch2);
                        minfo->c_type = 1;
                    }
                    if(strcmp(pch2,"profile")==0) {
                        pch2=strtok(NULL,"=");
                        if(minfo->c_type!=0)
                           elac_error("mspeed: the sound speed is already defined.\n",-1);
                        strcpy(finfo->c_fn, finfo->indir);
                        strcat(finfo->c_fn, "/");
                        strcat(finfo->c_fn, pch2);
                        minfo->c_type = 2;
                    }
                    if(strcmp(pch2,"2dfile")==0) {
                        pch2=strtok(NULL,"=");
                        //strcpy(*c_fn, pch2);
                        if(minfo->c_type!=0)
                           elac_error("mspeed: the sound speed is already defined.\n",-1);
                        strcpy(finfo->c_fn, finfo->indir);
                        strcat(finfo->c_fn, "/");
                        strcat(finfo->c_fn, pch2);
                        minfo->c_type = 3;
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        finfo->c_format = 0;
                        if(strcmp(pch2,"binary")==0) finfo->c_format = 0;
                        if(strcmp(pch2,"ascii")==0) finfo->c_format = 1;
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );
    return;
}

void get_dens( char *fn, innout *finfo, fdmesh *minfo, int mrank )
{
    FILE    *ptr_file;
    char    buf[1000],buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    minfo->rho_type = 0;
    minfo->rho_val = 0;
    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file.\n", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "mdensity" ) == 0 ) id_num++;
        }
    }

    if( id_num != 1 )
        elac_error("mdensity (material density) is not defined or defined repeatedly.\n",-1);

    rewind( ptr_file );
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "mdensity" ) == 0 ) {
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"value")==0) {
                        pch2=strtok(NULL,"=");
                        if(minfo->rho_type!=0)
                           elac_error("mdensity: the material density is already defined.\n",-1);
                        minfo->rho_val = atof(pch2);
                        minfo->rho_type = 1;
                    }
                    if(strcmp(pch2,"profile")==0) {
                        pch2=strtok(NULL,"=");
                        if(minfo->rho_type!=0)
                           elac_error("mdensity: the material density is already defined.\n",-1);
                        strcpy(finfo->rho_fn, finfo->indir);
                        strcat(finfo->rho_fn, "/");
                        strcat(finfo->rho_fn, pch2);
                        minfo->rho_type = 2;
                    }
                    if(strcmp(pch2,"2dfile")==0) {
                        pch2=strtok(NULL,"=");
                        //strcpy(*rho_fn, pch2);
                        if(minfo->rho_type!=0)
                           elac_error("mdensity: the material density is already defined.\n",-1);
                        strcpy(finfo->rho_fn, finfo->indir);
                        strcat(finfo->rho_fn, "/");
                        strcat(finfo->rho_fn, pch2);
                        minfo->rho_type = 3;
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        finfo->rho_format = 0;
                        if(strcmp(pch2,"binary")==0) finfo->rho_format = 0;
                        if(strcmp(pch2,"ascii")==0) finfo->rho_format = 1;
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );
    return;
}

void get_wind( char *fn, innout *finfo, fdmesh *minfo, int mrank )
{
    FILE    *ptr_file;
    char    buf[1000],buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    minfo->w_type = 0;
    minfo->w_val = 0;
    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file.\n", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "wind" ) == 0 ) id_num++;
        }
    }

    if( id_num != 1 )
        elac_error("wind (background flow) is not defined or defined repeatedly.\n",-1);

    rewind( ptr_file );
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "wind" ) == 0 ) {
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"value")==0) {
                        pch2=strtok(NULL,"=");
                        if(minfo->w_type!=0)
                           elac_error("wind: wind speed is already defined.\n",-1);
                        minfo->w_val = atof(pch2);
                        minfo->w_type = 1;
                    }
                    if(strcmp(pch2,"profile")==0) {
                        pch2=strtok(NULL,"=");
                        if(minfo->w_type!=0)
                           elac_error("wind: wind speed is already defined.\n",-1);
                        strcpy(finfo->w_fn, finfo->indir);
                        strcat(finfo->w_fn, "/");
                        strcat(finfo->w_fn, pch2);
                        minfo->w_type = 2;
                    }
                    if(strcmp(pch2,"2dfile")==0) {
                        pch2=strtok(NULL,"=");
                        //strcpy(*w_fn, pch2);
                        if(minfo->w_type!=0)
                           elac_error("wind: wind speed is already defined.\n",-1);
                        strcpy(finfo->w_fn, finfo->indir);
                        strcat(finfo->w_fn, "/");
                        strcat(finfo->w_fn, pch2);
                        minfo->w_type = 3;
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        finfo->w_format = 0;
                        if(strcmp(pch2,"binary")==0) finfo->w_format = 0;
                        if(strcmp(pch2,"ascii")==0) finfo->w_format = 1;
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );
    return;
}

void get_time( char *fn, fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 )
{
    FILE    *ptr_file;
    char    buf[1000],buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;
    double tt;
    int i;

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file.\n", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "time" ) == 0 ) id_num++;
        }
    }

    if( id_num != 1 )
        elac_error("TIME is not defined or defined repeatedly.\n",-1);

    rewind( ptr_file );
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "time" ) == 0 ) {
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"t")==0) {
                        pch2=strtok(NULL,"=");
                        minfo->T = atof(pch2);
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );

    minfo->dt = cfl * (minfo->R*minfo->dth) / (minfo->cmax);
    tt=0;
    for(i=0; tt<=minfo->T; i++) {
      tt=i*minfo->dt;
    }
    minfo->nt = i;
    minfo->ts = (dvec *) alloc_dvec( minfo->nt);
    for(i=0; i<(minfo->nt); i++) {
      tt = i*minfo->dt;
      minfo->ts->elements[i] = tt;
    }

    sprintf(buf, "dt = %f, nt = %d", minfo->dt, minfo->nt);
    MPI_Barrier(comm2);
    mpi_printf(comm2, buf, 200, mrank, nRank);

    return;
}


void get_src( char *fn, fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2)
{
    FILE    *ptr_file;
    char    buf[1000], buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    int     i, id_num=0, local_id;
    double  x_g=0.0, src_p0=0, src_freq=0, src_y=0.0;
    int     xi_g=0, xi_l;
    char    src_type[100];

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "asource" ) == 0 ) id_num++;
        }
    }

    if( id_num < 1 )
        elac_error("No source is given.\n",-1);

   minfo->src_num_global = id_num;
//    minfo->src_num = id_num;
//    minfo->srcs = (msource *) malloc( id_num*sizeof(msource) );

   rewind( ptr_file );
    id_num = 0;
    local_id = 0;
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "asource" ) == 0 ) {
                id_num++;
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"angle")==0) {
                        pch2=strtok(NULL,"=");
                        x_g = atof(pch2)/180*M_PI;
                        xi_g = floor(x_g/minfo->dth+0.5) + ghostL;
                        if((xi_g>=(minfo->xi_start+ghostL)) &  (xi_g<=(minfo->xi_end-ghostL))) {
                           local_id++;
                        }
                        //minfo->srcs[id_num-1].x=atof(pch2)*M_PI/180;
                        //minfo->srcs[id_num-1].xi=floor(minfo->srcs[id_num-1].x/minfo->dth+0.5)+ghostL+1;
                    }
                    if(strcmp(pch2,"elev")==0) {
                        pch2=strtok(NULL,"=");
                        //minfo->srcs[id_num-1].y=atof(pch2);
                        //minfo->srcs[id_num-1].yi=floor(minfo->srcs[id_num-1].y/minfo->dr+0.5)+ghostL+1;
                    }
                    if(strcmp(pch2,"p0")==0) {
                        pch2=strtok(NULL,"=");
                        //minfo->srcs[id_num-1].p0=atof(pch2);
                    }
                    if(strcmp(pch2,"freq")==0) {
                        pch2=strtok(NULL,"=");
                        //minfo->srcs[id_num-1].freq=atof(pch2);
                        //printf("freq = %s\n",pch2);
                        //printf("freq = %f\n",atof(pch2));
                    }
                    if(strcmp(pch2,"type")==0) {
                        pch2=strtok(NULL,"=");
                        //strcpy(minfo->srcs[id_num-1].type,pch2);
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }

    minfo->src_num = local_id;
    minfo->srcs = (msource *) malloc( local_id * sizeof(msource) );
    for(i=0; i<minfo->src_num; i++) {
      minfo->srcs[i].q = (dvec *) alloc_dvec (minfo->nt);
      minfo->srcs[i].midq = (dvec *) alloc_dvec (minfo->nt);
    }

   rewind( ptr_file );
   id_num = 0;
   local_id = 0;
   while ( fgets( buf, 1000, ptr_file )!=NULL ) {
      strcpy( buf2, buf );
      pch = strtok(buf2," \n");
      itok=1;
      if( pch != NULL && pch[0] != '#') {
         if( strcmp( pch, "asource" ) == 0 ) {
            id_num++;
            pch = strtok( NULL, " \n");
            while ( pch != NULL ) {
               strcpy( opt, pch ); pch2 = strtok( opt, "=");
                  if(strcmp(pch2,"angle")==0) {
                     pch2=strtok(NULL,"=");
                     x_g = atof(pch2)/180*M_PI;
                     xi_g = floor(x_g/minfo->dth+0.5) + ghostL;
                        //minfo->srcs[id_num-1].x=atof(pch2)*M_PI/180;
                        //minfo->srcs[id_num-1].xi=floor(minfo->srcs[id_num-1].x/minfo->dth+0.5)+ghostL+1;
                  }
                    if(strcmp(pch2,"elev")==0) {
                        pch2=strtok(NULL,"=");
                        src_y = atof(pch2);
                        //minfo->srcs[id_num-1].y=atof(pch2);
                        //minfo->srcs[id_num-1].yi=floor(minfo->srcs[id_num-1].y/minfo->dr+0.5)+ghostL+1;
                    }
                    if(strcmp(pch2,"p0")==0) {
                        pch2=strtok(NULL,"=");
                        src_p0 = atof(pch2);
                        //minfo->srcs[id_num-1].p0=atof(pch2);
                    }
                    if(strcmp(pch2,"freq")==0) {
                        pch2=strtok(NULL,"=");
                        src_freq = atof(pch2);
                        //minfo->srcs[id_num-1].freq=atof(pch2);
                        //printf("freq = %s\n",pch2);
                        //printf("freq = %f\n",atof(pch2));
                    }
                    if(strcmp(pch2,"type")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(src_type, pch2);
                        //strcpy(minfo->srcs[id_num-1].type,pch2);
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
                
                if((xi_g>=(minfo->xi_start+ghostL)) & (xi_g<=(minfo->xi_end-ghostL))) {
                     local_id++;
                     xi_l = xi_g - minfo->xi_start;
                     minfo->srcs[local_id-1].x = x_g;
                     minfo->srcs[local_id-1].xi_global = xi_g;
                     minfo->srcs[local_id-1].xi = xi_l; 

                     minfo->srcs[local_id-1].y = src_y;
                     minfo->srcs[local_id-1].yi = floor(src_y/minfo->dr)+ghostL;
 
                     minfo->srcs[local_id-1].p0 = src_p0;
                     minfo->srcs[local_id-1].freq = src_freq;
                     minfo->srcs[local_id-1].t0 = 0.0;
                     
                     strcpy(minfo->srcs[local_id-1].type, src_type);
                }
            }
        }
    }

    fclose( ptr_file );
    sprintf(buf, "global src=%d, local src=%d", minfo->src_num_global, minfo->src_num);
    MPI_Barrier(comm2);
    mpi_printf(comm2, buf, 200, mrank, nRank);

    return;
}


void get_rec( char *fn, fdmesh *minfo, int mrank, int nRank, MPI_Comm comm2 )
{
    FILE    *ptr_file;
    char    buf[1000], buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0, local_id;
    int     i;
    double  x_g=0;
    int     xi_g=0, xi_l;
    char   sta_mode[100], sta_name[100];
    int    sta_format=0;
    double sta_y=0;

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "rec" ) == 0 ) id_num++;
        }
    }

    if( id_num < 1 )
        elac_error("No receiver given.\n",-1);

    minfo->sta_num_global = id_num;

//    minfo->sta_num = id_num;
//    minfo->stas = (station *) malloc( id_num*sizeof(station) );

   rewind( ptr_file );
    id_num = 0;
    local_id = 0;
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "rec" ) == 0 ) {
                id_num++;
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"angle")==0) {
                        pch2=strtok(NULL,"=");
                        x_g=atof(pch2)/180.0*M_PI;
                        xi_g=floor(x_g/minfo->dth + 0.5)+ghostL;
                        if((xi_g>=(minfo->xi_start+ghostL)) &  (xi_g<=(minfo->xi_end-ghostL))) {
                        local_id++;
                        //minfo->stas[id_num-1].x=atof(pch2)/180.0*M_PI;
                        //minfo->stas[id_num-1].xi=floor(minfo->stas[id_num-1].x/minfo->dth+0.5)+ghostL+1;
                        }

                    }
                    if(strcmp(pch2,"elev")==0) {
                        pch2=strtok(NULL,"=");
                        //minfo->stas[id_num-1].y=atof(pch2);
                        //minfo->stas[id_num-1].yi=floor(minfo->stas[id_num-1].y/minfo->dr+0.5)+ghostL+1;

                    }
                    if(strcmp(pch2,"mode")==0) {
                        pch2=strtok(NULL,"=");
                        //strcpy(minfo->stas[id_num-1].mode,pch2);
                    }
                    if(strcmp(pch2,"name")==0) {
                        pch2=strtok(NULL,"=");
                        //strcpy(minfo->stas[id_num-1].name,pch2);
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        //minfo->stas[id_num-1].format=atoi(pch2);
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }

    minfo->sta_num = local_id;
    minfo->stas = (station *) malloc( local_id*sizeof(station) );
    for(i=0; i<minfo->sta_num; i++) {
      minfo->stas[i].pout = (dvec *) alloc_dvec (minfo->nt);
      minfo->stas[i].rhout = (dvec *) alloc_dvec (minfo->nt);
      minfo->stas[i].uout = (dvec *) alloc_dvec (minfo->nt);
      minfo->stas[i].wout = (dvec *) alloc_dvec (minfo->nt);
    }

   rewind( ptr_file );
    id_num = 0;
    local_id = 0;
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "rec" ) == 0 ) {
                id_num++;
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"angle")==0) {
                        pch2=strtok(NULL,"=");
                        x_g=atof(pch2)/180.0*M_PI;
                        xi_g=floor(x_g/minfo->dth + 0.5)+ghostL;
                    }
                    if(strcmp(pch2,"elev")==0) {
                        pch2=strtok(NULL,"=");
                        sta_y=atof(pch2);
                    }
                    if(strcmp(pch2,"mode")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(sta_mode,pch2);
                    }
                    if(strcmp(pch2,"name")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(sta_name,pch2);
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        //sta_format=atoi(pch2);
                        sta_format=0;
                        if(strcmp(pch2,"binary")==0) sta_format = 0;
                        if(strcmp(pch2,"ascii")==0) sta_format = 1;
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
                    if((xi_g>=(minfo->xi_start+ghostL)) &  (xi_g<=(minfo->xi_end-ghostL))) {
                           local_id++;
                           xi_l = xi_g - minfo->xi_start;
                           minfo->stas[local_id-1].x=x_g;
                           minfo->stas[local_id-1].xi_global=xi_g;
                           minfo->stas[local_id-1].xi=xi_l;

                           minfo->stas[local_id-1].y=sta_y;
                           minfo->stas[local_id-1].yi=floor(sta_y/minfo->dr)+ghostL;

                           strcpy(minfo->stas[local_id-1].mode,sta_mode);
                           strcpy(minfo->stas[local_id-1].name,sta_name);
                           minfo->stas[local_id-1].format=sta_format;
                    }
            }
        }
    }

    fclose( ptr_file );

    sprintf(buf,"global rec=%d, local rec=%d", minfo->sta_num_global, minfo->sta_num);
    MPI_Barrier(comm2);
    mpi_printf(comm2, buf, 200, mrank, nRank);

    return;
}

void get_img( char *fn, innout *finfo, int mrank, int nRank, MPI_Comm comm2)
{
    FILE    *ptr_file;
    char    buf[1000], buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "image" ) == 0 ) id_num++;
        }
    }

    finfo->img_num = id_num;
    finfo->imgs = (image *) malloc ( finfo->img_num * sizeof(image) );

    if( id_num < 1 )
        return;
        //elac_error("No image given.\n",-1);

    rewind( ptr_file );
    id_num = 0;
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "image" ) == 0 ) {
                id_num++;
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"timeInterval")==0) {
                        pch2=strtok(NULL,"=");
                        finfo->imgs[id_num-1].imgdt=atof(pch2);
                    }
                    if(strcmp(pch2,"mode")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(finfo->imgs[id_num-1].mode,pch2);
                    }
                    if(strcmp(pch2,"file")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(finfo->imgs[id_num-1].file,pch2);
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        //finfo->imgs[id_num-1].format=atoi(pch2);
                        finfo->imgs[id_num-1].format=0;
                        if(strcmp(pch2,"binary")==0) finfo->imgs[id_num-1].format = 0;
                        if(strcmp(pch2,"ascii")==0) finfo->imgs[id_num-1].format = 1;
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );
    return;
}

void get_line( char *fn, innout *finfo, int mrank, int nRank, MPI_Comm comm2)
{
    FILE    *ptr_file;
    char    buf[1000], buf2[1000];
    long    indx=0,itok=0;
    char    *pch, opt[1000], *pch2;
    long    id_num=0;

    ptr_file = fopen( fn, "r" );
    if (!ptr_file) elac_error("ac2dr: cannot find the input file", -1);

    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        pch = strtok(buf," \n");
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "line" ) == 0 ) id_num++;
        }
    }

    finfo->line_num = id_num;
    finfo->lines = (lineout *) malloc ( finfo->line_num * sizeof(lineout) );

    if( id_num < 1 )
        return;
        //elac_error("No image given.\n",-1);

    rewind( ptr_file );
    id_num = 0;
    while ( fgets( buf, 1000, ptr_file )!=NULL ) {
        strcpy( buf2, buf );
        pch = strtok(buf2," \n");
        itok=1;
        if( pch != NULL && pch[0] != '#') {
            if( strcmp( pch, "line" ) == 0 ) {
                id_num++;
                pch = strtok( NULL, " \n");
                while ( pch != NULL ) {
                    strcpy( opt, pch ); pch2 = strtok( opt, "=");
                    if(strcmp(pch2,"timeInterval")==0) {
                        pch2=strtok(NULL,"=");
                        finfo->lines[id_num-1].linedt=atof(pch2);
                    }
                    if(strcmp(pch2,"elev")==0) {
                        pch2=strtok(NULL,"=");
                        finfo->lines[id_num-1].elev=atof(pch2);
                    }
                    if(strcmp(pch2,"mode")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(finfo->lines[id_num-1].mode,pch2);
                    }
                    if(strcmp(pch2,"file")==0) {
                        pch2=strtok(NULL,"=");
                        strcpy(finfo->lines[id_num-1].file,pch2);
                    }
                    if(strcmp(pch2,"format")==0) {
                        pch2=strtok(NULL,"=");
                        //finfo->lines[id_num-1].format=atoi(pch2);
                        finfo->lines[id_num-1].format=0;
                        if(strcmp(pch2,"binary")==0) finfo->lines[id_num-1].format = 0;
                        if(strcmp(pch2,"ascii")==0) finfo->lines[id_num-1].format = 1;
                    }
                    itok++;
                    strcpy( buf2, buf );
                    pch = strtok( buf2, " \n" );
                    for(indx=0;indx<itok;indx++) {
                        pch = strtok( NULL, " \n");
                    }
                }
            }
        }
    }
    fclose( ptr_file );
    return;
}

void set_lineimg( fdmesh *minfo, innout *finfo, int mrank, int nRank, MPI_Comm comm2 ) {
   int i, k;
   int l, tl;

   for(i=0; i<finfo->line_num; i++) {
      finfo->lines[i].yi = floor(finfo->lines[i].elev/minfo->dr)+ghostL;
      l = floor(finfo->lines[i].linedt / minfo->dt + 0.5);
      tl = 0;
      for(k=0;k<minfo->nt;k++) {
         if(k%l == 0) tl = tl+1;
      } 
      finfo->lines[i].tl = tl; 
      finfo->lines[i].lineimg = alloc_array2( minfo->nx, tl );
   }
   return;
}


void print_message (FILE *ptr_file, char *msg) {
    fprintf(ptr_file,msg);
    fprintf(stdout,msg);

}


