#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

#include "nsopde_headers/macros.h"
#include "nsopde_headers/io.h"
#include "nsopde_headers/space_discretization.h"
#include "nsopde_headers/ivp.h"
#include "nsopde_headers/initial_functions.h"
#include "nsopde_headers/lax_friedrichs_1d.h"
#include "nsopde_headers/boundary_conditions_1d.h"

int main(){
    vr16 domainX = {0, 2*3.14159};

    int N1 = 100;
    int N2 = 1000;

    r16 dx1 = (domainX[1]-domainX[0])/N1;
    r16 dx2 = (domainX[1]-domainX[0])/N2;

    r16 speed = 0.5;

    r16 Tf = 10;

    vr16 grid = discretize1d(domainX, N2);
    vr16 U1 = initialize_ivp(grid, &sine_wave);

    writeVector(grid, "grid.txt");

    // animate the solution for the 1d linear advection equation using 1d Lax friedrichs scheme
    lf1order_linear_advection_1d(U1, dx1, speed, Tf, true, periodic);

    return 0;
}