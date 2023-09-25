#include "ac2dr_aux.h"

double get_psi( double xi0 ) {
    double psi, xi;
    xi = fabs(xi0);
    psi = 0;
    if(xi<2) {
        psi = 1.0/2.0 - 1.0/12.0*xi - 5.0/32.0*pow(xi,2) + \
              5.0/192.0*pow(xi,3) + 1.0/128.0*pow(xi,4) - 1.0/768.0*pow(xi,5);
    }
    if((xi>=2) & (xi<4)) {
        psi = 1.0/2.0 - 13.0/48.0*xi - 5.0/64.0*pow(xi,2) + 25.0/384.0*pow(xi,3) - \
              3.0/256.0*pow(xi,4) + 1.0/1536.0*pow(xi,5);
    }
    if((xi>=4) & (xi<6)) {
        psi = 1.0/2.0 - 137.0/240.0*xi + 15.0/64.0*pow(xi,2) - 17.0/384.0*pow(xi,3) + \
              1.0/256.0*pow(xi,4) - 1.0/7680.0*pow(xi,5);
    }
    return psi;
}
