/*
    Header file contains implementstion of the lax friedrichs scheme for 1d for various ivp's
*/
#pragma once

#include<iostream>
#include<vector>
#include<functional>

#include "macros.h"

using namespace std;

void lf1order_linear_advection_1d(vr16& U, r16 dx, r16 speed, r16 Tf, bool animate, function<void(vr16&)> bc){
    int N = U.size()-1;
    
    vr16 U0(N+3), U1(N+3);

    for(int i=1 ; i<=N+1 ; i++){
        U0[i] = U[i-1];
    }

    r16 t=0;
    r16 dt = CFL1*dx/speed;

    while(t < Tf){
        cout << "t=" << t << " | dt=" << dt << endl;
        t += dt;

        bc(U0);

        for(int i=1 ; i<=N+1 ; i++){
            U1[i] = 0.5*( U0[i+1] + U0[i-1] ) - 0.5*speed*(dt/dx)*( U0[i+1] - U0[i-1] );
        }

        for(int i=1 ; i<=N+1 ; i++){
            U0[i] = U1[i];
        }

        if(animate){
            writeVector(U1, "u.txt");
        }
    }

    for(int i=1 ; i<=N+1 ; i++){
        U[i-1] = U1[i];
    }
}
