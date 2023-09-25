-------------------------------------------------------------------------------
AC2Dr: Acoustic Codes in Spherical Coordinates

AC2Dr is a 2-D numerical solver for the acoustic wave equation using the finite
difference method. The acoustic wave equation is split into two first-order
differential equations for pressure and particle velocity, and approximated by
the sixth-order accurate central difference in space and time. The equations
are solved in the axisymmetric spherical coordinates, and hence, are able to
simulate 3-D spherical spreading of wavefields by 2-D. The absorbing boundary,
which prevents outgoing waves from reflecting at the computational domain
boundary, is realized by the super-grid method based on coordinate stretching.
AC2Dr supports the message passing interface (MPI) on multicore CPU,
improving simulation performance dramatically.

-------------------------------------------------------------------------------
REQUIREMENTS

C compilers compatible with the GNU C and environments enabling the Message
Passing Interface (MPI)

-------------------------------------------------------------------------------
INSTALLATION

make all

-------------------------------------------------------------------------------
DOCUMENTATION

See the user guide included in the the doc folder

-------------------------------------------------------------------------------
LICENSE

AC2Dr is distributed under the terms of the MIT license. All new contributions
must be made under this license.

SPDX-License-Identifier: MIT

LLNL-CODE-854516
