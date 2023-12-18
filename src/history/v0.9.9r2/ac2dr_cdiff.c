#include "ac2dr_aux.h"
#include "ac2dr_config.h"

#if (FD_ORDER==6)
#define spatial_derivatives dqdt6_sym
#define spatial_derivatives180 dqdt6_sym180
#endif 

#if (FD_ORDER==2)
#define spatial_derivatives dqdt2_sym
#define spatial_derivatives180 dqdt6_sym180
#endif 


void ac2dr_cdiff ( fdmesh *minfo, innout *finfo, double **u, double **v, double **w, \
                  double **rhohat, double **chat, double **what ) {

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

    nt = minfo->nt;
    dt = minfo->dt;
    nx = minfo->nth;
    ny = minfo->nr;
//    ds = minfo->ds;
//    dr = minfo->dr;

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

//    ac2dr_out( 1, 0, minfo, finfo, u, v, w, rhohat, chat, what );
    wd_array2( "output/test_U_snapshot000000.dat", u, minfo->nx, minfo->ny );
    wd_array2( "output/test_V_snapshot000000.dat", v, minfo->nx, minfo->ny );
    wd_array2( "output/test_W_snapshot000000.dat", w, minfo->nx, minfo->ny );
    wd_array2( "output/test_chat_snapshot000000.dat", chat, minfo->nx, minfo->ny );
    wd_array2( "output/test_rhohat_snapshot000000.dat", rhohat, minfo->nx, minfo->ny );
    wd_array2( "output/test_what_snapshot000000.dat", what, minfo->nx, minfo->ny );

    k=0;

    for(l=0; l<minfo->sta_num; l++) {
        i = minfo->stas[l].xi;
        j = minfo->stas[l].yi;
        minfo->stas[l].pout->elements[k] = u[i][j];
    }

    for(k=1; k<nt; k++) {
    //printf("k = %lu\n",k);
    snprintf(buf,998,"k = %lu", k); //print_message(ptr_file, buf);
    puts(buf);


        spatial_derivatives( u, v, w, u1, v1, w1, rhohat, what, chat, minfo);
//not necessary  updatesrc( k, s1, r1, rhohat, chat, jac, minfo, 1);

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt/2.0*u1[i][j];
                v0[i][j] = v[i][j] + dt/2.0*v1[i][j];
                w0[i][j] = w[i][j] + dt/2.0*w1[i][j];
            }
        }

if(BC_ENFORCE==1) {
        bcenforce_rigid_bottom( v0, minfo );
        bcenforce_rigid_top( v0, minfo );
        bcenforce_rigid_right( w0, minfo );
//        bcenforce_rigid_left( w0, minfo );
}

        spatial_derivatives( u0, v0, w0, u2, v2, w2, rhohat, what, chat, minfo );

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt/2.0*u2[i][j];
                v0[i][j] = v[i][j] + dt/2.0*v2[i][j];
                w0[i][j] = w[i][j] + dt/2.0*w2[i][j];
            }
        }
if(BC_ENFORCE==1) {
        bcenforce_rigid_bottom( v0, minfo );
        bcenforce_rigid_top( v0, minfo );
        bcenforce_rigid_right( w0, minfo );
//        bcenforce_rigid_left( w0, minfo );
}
        spatial_derivatives( u0, v0, w0, u3, v3, w3, rhohat, what, chat, minfo );

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt*u3[i][j];
                v0[i][j] = v[i][j] + dt*v3[i][j];
                w0[i][j] = w[i][j] + dt*w3[i][j];
            }
        }
if(BC_ENFORCE==1) {
        bcenforce_rigid_bottom( v0, minfo );
        bcenforce_rigid_top( v0, minfo );
        bcenforce_rigid_right( w0, minfo );
//        bcenforce_rigid_left( w0, minfo );
}
        spatial_derivatives( u0, v0, w0, u4, v4, w4, rhohat, what, chat, minfo );

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u[i][j] = u[i][j] + dt/6.0*(u1[i][j] + 2.0*u2[i][j] + 2.0*u3[i][j] + u4[i][j]);
                v[i][j] = v[i][j] + dt/6.0*(v1[i][j] + 2.0*v2[i][j] + 2.0*v3[i][j] + v4[i][j]);
                w[i][j] = w[i][j] + dt/6.0*(w1[i][j] + 2.0*w2[i][j] + 2.0*w3[i][j] + w4[i][j]);
            }
        }
if(BC_ENFORCE==1) {
        bcenforce_rigid_bottom( v, minfo );
        bcenforce_rigid_top( v, minfo );
        bcenforce_rigid_right( w, minfo );
//        bcenforce_rigid_left( w, minfo );
}       
        for(l=0; l<minfo->sta_num; l++) {
            i = minfo->stas[l].xi;
            j = minfo->stas[l].yi;
            minfo->stas[l].pout->elements[k] = u[i][j]*rhohat[i][j]*chat[i][j];
        }
//    printf("test\n");
    if(k%100==0) {
       sprintf(buf,"output/test_U_snapshot%06d.dat",k);
       wd_array2( buf, u, minfo->nx, minfo->ny );
       sprintf(buf,"output/test_V_snapshot%06d.dat",k);
       wd_array2(buf, v, minfo->nx, minfo->ny );
       sprintf(buf,"output/test_W_snapshot%06d.dat",k);
       wd_array2(buf, w, minfo->nx, minfo->ny );
    }
/*    if(k==(nt-1)) {
       wd_array2( "output/test_U_snapshot100.dat", u, minfo->nx, minfo->ny );
       wd_array2( "output/test_V_snapshot100.dat", v, minfo->nx, minfo->ny );
       wd_array2( "output/test_W_snapshot100.dat", w, minfo->nx, minfo->ny );
    }
*/

//        ac2dr_out( 1, k, minfo, finfo, u, v, w, rhohat, chat, what );
    }

//    printf("***********ac2dr_RK4.c is working\n");
    ac2dr_out( 0, k, minfo, finfo, u, v, w, rhohat, chat, what );
    //fclose(ptr_file);
//    ac2dr_out( 0, k, minfo, finfo, u, v, w, rhohat, chat, what );
//    printf("source number = %d\n", minfo->src_num);
//    wd_vec("./output/srcf.txt",minfo->srcs[0].q);
//    wd_vec("./output/srcf2.txt",minfo->srcs[0].midq);
    return;
}












//////////////////////////////////////////////////////////////////////////////////////////

void ac2dr_RK4180 ( fdmesh *minfo, innout *finfo, double **u, double **v, double **w, \
                  double **rhohat, double **chat, double **what ) {

    unsigned long i, j, k;
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

    nt = minfo->nt;
    dt = minfo->dt;
    nx = minfo->nth;
    ny = minfo->nr;
//    ds = minfo->ds;
//    dr = minfo->dr;

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

//    ac2dr_out( 1, 0, minfo, finfo, u, v, w, rhohat, chat, what );
    wd_array2( "output/test_U_snapshot0.dat", u, minfo->nx, minfo->ny );
    wd_array2( "output/test_V_snapshot0.dat", v, minfo->nx, minfo->ny );
    wd_array2( "output/test_W_snapshot0.dat", w, minfo->nx, minfo->ny );


    printf("Running RK4_180 \n");
    k=0;
//    for(l=0; l<minfo->sta_num; l++) {
//        i = minfo->stas[l].xi;
//        j = minfo->stas[l].yi;
//        minfo->stas[l].pout->elements[k] = u[i][j]*rhohat[i][j]*chat[i][j];
//    }

    for(k=1; k<nt; k++) {
    //printf("k = %lu\n",k);
    snprintf(buf,998,"k = %lu", k); //print_message(ptr_file, buf);
    puts(buf);


        spatial_derivatives180( u, v, w, u1, v1, w1, rhohat, what, chat, minfo);
//not necessary  updatesrc( k, s1, r1, rhohat, chat, jac, minfo, 1);

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt/2.0*u1[i][j];
                v0[i][j] = v[i][j] + dt/2.0*v1[i][j];
                w0[i][j] = w[i][j] + dt/2.0*w1[i][j];
            }
        }
        bcenforce_rigid_bottom( v0, minfo );
        bcenforce_rigid_top( v0, minfo );
//        bcenforce_rigid_right( w0, minfo );
//        bcenforce_rigid_left( w0, minfo );

        spatial_derivatives180( u0, v0, w0, u2, v2, w2, rhohat, what, chat, minfo );

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt/2.0*u2[i][j];
                v0[i][j] = v[i][j] + dt/2.0*v2[i][j];
                w0[i][j] = w[i][j] + dt/2.0*w2[i][j];
            }
        }
        bcenforce_rigid_bottom( v0, minfo );
        bcenforce_rigid_top( v0, minfo );
//        bcenforce_rigid_right( w0, minfo );
//        bcenforce_rigid_left( w0, minfo );

        spatial_derivatives180( u0, v0, w0, u3, v3, w3, rhohat, what, chat, minfo );

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u0[i][j] = u[i][j] + dt*u3[i][j];
                v0[i][j] = v[i][j] + dt*v3[i][j];
                w0[i][j] = w[i][j] + dt*w3[i][j];
            }
        }
        bcenforce_rigid_bottom( v0, minfo );
        bcenforce_rigid_top( v0, minfo );
//        bcenforce_rigid_right( w0, minfo );
//        bcenforce_rigid_left( w0, minfo );

        spatial_derivatives180( u0, v0, w0, u4, v4, w4, rhohat, what, chat, minfo );

        for(j=0; j<ny; j++) {
            for(i=0; i<nx; i++) {
                u[i][j] = u[i][j] + dt/6.*(u1[i][j] + 2.*u2[i][j] + 2.*u3[i][j] + u4[i][j]);
                v[i][j] = v[i][j] + dt/6.*(v1[i][j] + 2.*v2[i][j] + 2.*v3[i][j] + v4[i][j]);
                w[i][j] = w[i][j] + dt/6.*(w1[i][j] + 2.*w2[i][j] + 2.*w3[i][j] + w4[i][j]);
            }
        }
        bcenforce_rigid_bottom( v, minfo );
        bcenforce_rigid_top( v, minfo );
//        bcenforce_rigid_right( w, minfo );
//        bcenforce_rigid_left( w, minfo );
       
//        for(l=0; l<minfo->sta_num; l++) {
//            i = minfo->stas[l].xi;
//            j = minfo->stas[l].yi;
//            minfo->stas[l].pout->elements[k] = u[i][j];
//        }
//    if(k%5==0) {
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

//        ac2dr_out( 1, k, minfo, finfo, u, v, w, rhohat, chat, what );
    }
    printf("***********ac2dr_RK4.c is working\n");
    //fclose(ptr_file);
//    ac2dr_out( 0, k, minfo, finfo, u, v, w, rhohat, chat, what );
//    printf("source number = %d\n", minfo->src_num);
//    wd_vec("./output/srcf.txt",minfo->srcs[0].q);
//    wd_vec("./output/srcf2.txt",minfo->srcs[0].midq);
    return;
}





