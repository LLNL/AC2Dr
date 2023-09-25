#include "ac2dr_aux.h"

void getsymmetric_left( int iLB[][6] ) {
    iLB[0][0] = 3;
    iLB[0][1] = 2;
    iLB[0][2] = 1;
    iLB[0][3] = 1;
    iLB[0][4] = 2;
    iLB[0][5] = 3;

    iLB[1][0] = 2;
    iLB[1][1] = 1;
    iLB[1][2] = 0;
    iLB[1][3] = 2;
    iLB[1][4] = 3;
    iLB[1][5] = 4;

    iLB[2][0] = 1;
    iLB[2][1] = 0;
    iLB[2][2] = 1;
    iLB[2][3] = 3;
    iLB[2][4] = 4;
    iLB[2][5] = 5;
    return;
}

void getsymmetric_right( int iRB[][6], unsigned long nx ) {
    iRB[0][0] = nx-6;
    iRB[0][1] = nx-5;
    iRB[0][2] = nx-4;
    iRB[0][3] = nx-2;
    iRB[0][4] = nx-1;
    iRB[0][5] = nx-2;

    iRB[1][0] = nx-5;
    iRB[1][1] = nx-4;
    iRB[1][2] = nx-3;
    iRB[1][3] = nx-1;
    iRB[1][4] = nx-2;
    iRB[1][5] = nx-3;

    iRB[2][0] = nx-4;
    iRB[2][1] = nx-3;
    iRB[2][2] = nx-2;
    iRB[2][3] = nx-2;
    iRB[2][4] = nx-3;
    iRB[2][5] = nx-4;
    return;
}


void getantisymmetric_left( int iLB[][6] ) {
    iLB[0][0] = -1.0;
    iLB[0][1] = -1.0;
    iLB[0][2] = -1.0;
    iLB[0][3] = 1.0;
    iLB[0][4] = 1.0;
    iLB[0][5] = 1.0;

    iLB[1][0] = -1.0;
    iLB[1][1] = -1.0;
    iLB[1][2] = 1.0;
    iLB[1][3] = 1.0;
    iLB[1][4] = 1.0;
    iLB[1][5] = 1.0;

    iLB[2][0] = -1.0;
    iLB[2][1] = 1.0;
    iLB[2][2] = 1.0;
    iLB[2][3] = 1.0;
    iLB[2][4] = 1.0;
    iLB[2][5] = 1.0;
    return;
}

void getantisymmetric_right( int iRB[][6] ) {
    iRB[0][0] = 1.0;
    iRB[0][1] = 1.0;
    iRB[0][2] = 1.0;
    iRB[0][3] = 1.0;
    iRB[0][4] = 1.0;
    iRB[0][5] = -1.0;

    iRB[1][0] = 1.0;
    iRB[1][1] = 1.0;
    iRB[1][2] = 1.0;
    iRB[1][3] = 1.0;
    iRB[1][4] = -1.0;
    iRB[1][5] = -1.0;

    iRB[2][0] = 1.0;
    iRB[2][1] = 1.0;
    iRB[2][2] = 1.0;
    iRB[2][3] = -1.0;
    iRB[2][4] = -1.0;
    iRB[2][5] = -1.0;

    return;
}


void bcenforce_rigid_bottom( double **v, fdmesh *minfo ) {
    int i, j, nx;
    nx = minfo->nx;
    j = 0;
    for(i=0; i<nx; i++) {
         v[i][j] = 0.0;
    }
}

void bcenforce_rigid_top( double **v, fdmesh *minfo ) {
    int i, j, nx;
    nx = minfo->nx;
    j = minfo->ny-1;
    for(i=0; i<nx; i++) {
         v[i][j] = 0.0;
    }
}

void bcenforce_rigid_right( double **w, fdmesh *minfo ) {
    int i, j, ny;
    ny = minfo->ny;
    i = minfo->nx-1;
    for(j=0; j<ny; j++) {
         w[i][j] = 0.0;
    }
}

void bcenforce_rigid_left( double **w, fdmesh *minfo ) {
    int i, j, ny;
    ny = minfo->ny;
    i = 0;
    for(j=0; j<ny; j++) {
         w[i][j] = 0.0;
    }
}

//////////////////////////////////////////////////////
void bcenforce_symmetric_left3( double **w, fdmesh *minfo ) {
    int i, j, ny;
    ny = minfo->ny;
    i = 2;
    for(j=0; j<ny; j++) {
         w[i][j] = w[4][j];
    }
    i = 1;
    for(j=0; j<ny; j++) {
         w[i][j] = w[5][j];
    }
    i = 0;
    for(j=0; j<ny; j++) {
         w[i][j] = w[6][j];
    }
}

void bcenforce_antisymmetric_left3( double **w, fdmesh *minfo ) {
    int i, j, ny;
    ny = minfo->ny;
    i = 2;
    for(j=0; j<ny; j++) {
         w[i][j] = -w[4][j];
    }
    i = 1;
    for(j=0; j<ny; j++) {
         w[i][j] = -w[5][j];
    }
    i = 0;
    for(j=0; j<ny; j++) {
         w[i][j] = -w[6][j];
    }
}


void bcenforce_symmetric_right3( double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;
    i = nx-3;
    for(j=0; j<ny; j++) {
         w[i][j] = w[nx-5][j];
    }
    i = nx-2;
    for(j=0; j<ny; j++) {
         w[i][j] = w[nx-6][j];
    }
    i = nx-1;
    for(j=0; j<ny; j++) {
         w[i][j] = w[nx-7][j];
    }
}

void bcenforce_symmetric_bottom3( double **w, fdmesh *minfo ) {
    int i, j, nx;
    nx = minfo->nx;
    j = 2;
    for(i=0; i<nx; i++) {
         w[i][j] = w[i][4];
    }
    j = 1;
    for(i=0; i<nx; i++) {
         w[i][j] = w[i][5];
    }
    j = 0;
    for(i=0; i<nx; i++) {
         w[i][j] = w[i][6];
    }
}


void bcenforce_symmetric_top3( double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;
    j = ny-3;
    for(i=0; i<nx; i++) {
         w[i][j] = w[i][ny-5];
    }
    j = ny-2;
    for(i=0; i<nx; i++) {
         w[i][j] = w[i][ny-6];
    }
    j = ny-1;
    for(i=0; i<nx; i++) {
         w[i][j] = w[i][ny-7];
    }
}

void bcenforce_antisymmetric_right3( double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;
    i = nx-3;
    for(j=0; j<ny; j++) {
         w[i][j] = -w[nx-5][j];
    }
    i = nx-2;
    for(j=0; j<ny; j++) {
         w[i][j] = -w[nx-6][j];
    }
    i = nx-1;
    for(j=0; j<ny; j++) {
         w[i][j] = -w[nx-7][j];
    }
}

void bcenforce_antisymmetric_bottom3( double **w, fdmesh *minfo ) {
    int i, j, nx;
    nx = minfo->nx;
    j = 2;
    for(i=0; i<nx; i++) {
         w[i][j] = -w[i][4];
    }
    j = 1;
    for(i=0; i<nx; i++) {
         w[i][j] = -w[i][5];
    }
    j = 0;
    for(i=0; i<nx; i++) {
         w[i][j] = -w[i][6];
    }
}


void bcenforce_antisymmetric_top3( double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;
    j = ny-3;
    for(i=0; i<nx; i++) {
         w[i][j] = -w[i][ny-5];
    }
    j = ny-2;
    for(i=0; i<nx; i++) {
         w[i][j] = -w[i][ny-6];
    }
    j = ny-1;
    for(i=0; i<nx; i++) {
         w[i][j] = -w[i][ny-7];
    }
}



//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void bcenforce_rigid_bc1( double **u, double **v, double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;

//left
    for(j=0; j<ny; j++) {
         u[0][j] = u[2][j];
         v[0][j] = v[2][j];
         w[0][j] = -w[1][j];
    }

//bottom
    for(i=0; i<nx; i++) {
         u[i][0] = u[i][1];
         v[i][0] = 0.0;
         w[i][0] = w[i][1];
    }

//right
    for(j=0; j<ny; j++) {
         u[nx-1][j] = u[nx-2][j];
         v[nx-1][j] = v[nx-2][j];
         w[nx-1][j] = -w[nx-3][j];
         w[nx-2][j] = 0.0;
    }

//top
    for(i=0; i<nx; i++) {
         u[i][ny-1] = u[i][ny-2];
         v[i][ny-1] = -v[i][ny-3];
         w[i][ny-1] = w[i][ny-2];
    }
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void bcenforce_rigid_bc3( double **u, double **v, double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;

//left
    for(j=0; j<ny; j++) {
         u[0][j] = u[6][j];
         u[1][j] = u[5][j];
         u[2][j] = u[4][j];

         v[0][j] = v[6][j];
         v[1][j] = v[5][j];
         v[2][j] = v[4][j];

         w[0][j] = -w[5][j];
         w[1][j] = -w[4][j];
         w[2][j] = -w[3][j];
    }

//right
   for(j=0; j<ny; j++) {
         u[nx-1][j] = u[nx-6][j];
         u[nx-2][j] = u[nx-5][j];
         u[nx-3][j] = u[nx-4][j];

         v[nx-1][j] = v[nx-6][j];
         v[nx-2][j] = v[nx-5][j];
         v[nx-3][j] = v[nx-4][j];

         w[nx-1][j] = -w[nx-7][j];
         w[nx-2][j] = -w[nx-6][j];
         w[nx-3][j] = -w[nx-5][j];
   }

//bottom
    for(i=0; i<nx; i++) {
         u[i][0] = u[i][7];
         u[i][1] = u[i][6];
         u[i][2] = u[i][5];
         u[i][3] = u[i][4];

         w[i][0] = w[i][7];
         w[i][1] = w[i][6];
         w[i][2] = w[i][5];
         w[i][3] = w[i][4];

         v[i][0] = -v[i][6];
         v[i][1] = -v[i][5];
         v[i][2] = -v[i][4];;
         v[i][3] = 0.0;
    }

//top
    for(i=0; i<nx; i++) {
         u[i][ny-1] = u[i][ny-6];
         u[i][ny-2] = u[i][ny-5];
         u[i][ny-3] = u[i][ny-4];

         w[i][ny-1] = w[i][ny-6];
         w[i][ny-2] = w[i][ny-5];
         w[i][ny-3] = w[i][ny-4];

         v[i][ny-1] = -v[i][ny-7];
         v[i][ny-2] = -v[i][ny-6];
         v[i][ny-3] = -v[i][ny-5];
         v[i][ny-4] = -0.0;
    }
}



//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void bcenforce_sg_bc3( double **u, double **v, double **w, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;

//left
    for(j=0; j<ny; j++) {
         u[0][j] = u[6][j];
         u[1][j] = u[5][j];
         u[2][j] = u[4][j];

         v[0][j] = v[6][j];
         v[1][j] = v[5][j];
         v[2][j] = v[4][j];

         w[0][j] = -w[5][j];
         w[1][j] = -w[4][j];
         w[2][j] = -w[3][j];
    }

//bottom
    for(i=0; i<nx; i++) {
         u[i][0] = u[i][7];
         u[i][1] = u[i][6];
         u[i][2] = u[i][5];
         u[i][3] = u[i][4];

         w[i][0] = w[i][7];
         w[i][1] = w[i][6];
         w[i][2] = w[i][5];
         w[i][3] = w[i][4];

         v[i][0] = -v[i][6];
         v[i][1] = -v[i][5];
         v[i][2] = -v[i][4];;
         v[i][3] = 0.0;
    }

//right
   for(j=0; j<ny; j++) {
         u[nx-1][j] = 0.0;
         u[nx-2][j] = 0.0;
         u[nx-3][j] = 0.0;

         v[nx-1][j] = 0.0;
         v[nx-2][j] = 0.0;
         v[nx-3][j] = 0.0;

         w[nx-1][j] = 0.0;
         w[nx-2][j] = 0.0;
         w[nx-3][j] = 0.0;
   }

//top
    for(i=0; i<nx; i++) {
         u[i][ny-1] = 0.0;
         u[i][ny-2] = 0.0;
         u[i][ny-3] = 0.0;

         w[i][ny-1] = 0.0;
         w[i][ny-2] = 0.0;
         w[i][ny-3] = 0.0;

         v[i][ny-1] =  0.0;
         v[i][ny-2] =  0.0;
         v[i][ny-3] =  0.0;
         v[i][ny-4] = -0.0;
    }
}


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void bcenforce_rigid_material3( double **chat, double **rhohat, double **what, fdmesh *minfo ) {
    int i, j, nx, ny;
    nx = minfo->nx;
    ny = minfo->ny;

    for(j=(ny-4); j>=0; j--) {
       for(i=0; i<(nx-3); i++) {
         chat[i+3][j+3] = chat[i][j];
         rhohat[i+3][j+3] = rhohat[i][j];
         what[i+3][j+3] = what[i][j];
       }
    }

//left
    for(j=0; j<ny; j++) {
         chat[0][j] = chat[6][j];
         chat[1][j] = chat[5][j];
         chat[2][j] = chat[4][j];

         rhohat[0][j] = rhohat[6][j];
         rhohat[1][j] = rhohat[5][j];
         rhohat[2][j] = rhohat[4][j];

         what[0][j] = -what[5][j];
         what[1][j] = -what[4][j];
         what[2][j] = -what[3][j];
    }

//bottom
    for(i=0; i<nx; i++) {
         chat[i][0] = chat[i][7];
         chat[i][1] = chat[i][6];
         chat[i][2] = chat[i][5];
         chat[i][3] = chat[i][4];

         rhohat[i][0] = rhohat[i][7];
         rhohat[i][1] = rhohat[i][6];
         rhohat[i][2] = rhohat[i][5];
         rhohat[i][3] = rhohat[i][4];

         what[i][0] = what[i][7];
         what[i][1] = what[i][6];
         what[i][2] = what[i][5];
         what[i][3] = what[i][4];
    }


//right
   for(j=0; j<ny; j++) {
         chat[nx-1][j] = chat[nx-4][j];
         chat[nx-2][j] = chat[nx-4][j];
         chat[nx-3][j] = chat[nx-4][j];

         rhohat[nx-1][j] = rhohat[nx-4][j];
         rhohat[nx-2][j] = rhohat[nx-4][j];
         rhohat[nx-3][j] = rhohat[nx-4][j];

         what[nx-1][j] = what[nx-4][j];
         what[nx-2][j] = what[nx-4][j];
         what[nx-3][j] = what[nx-4][j];
   }

//top
    for(i=0; i<nx; i++) {
         chat[i][ny-1] = chat[i][ny-4];
         chat[i][ny-2] = chat[i][ny-4];
         chat[i][ny-3] = chat[i][ny-4];

         rhohat[i][ny-1] = rhohat[i][ny-4];
         rhohat[i][ny-2] = rhohat[i][ny-4];
         rhohat[i][ny-3] = rhohat[i][ny-4];

         what[i][ny-1] = what[i][ny-4];
         what[i][ny-2] = what[i][ny-4];
         what[i][ny-3] = what[i][ny-4];
    }
}



