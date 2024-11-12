/*
    Contains functions which can implement the boundary conditions as required in 1D problems
*/
#pragma once

#include<iostream>
#include<vector>

using namespace std;

#include "macros.h"

void periodic(vr16& U){
    U[0] = U[U.size()-2];
    U[U.size()-1] = U[0];
}

void free(vr16& U){
    U[0] = U[1];
    U[U.size()-1] = U[U.size()-2];
}
