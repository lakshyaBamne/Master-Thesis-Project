/*
    this header file contains functions to discretize the space domain
*/
#pragma once

#include<iostream>
#include<vector>

#include "macros.h"

using namespace std;

// function to discretize 1d domain
vector<real16> discretize1d(vector<real16> x, int N){
    real16 dx = (x[1]-x[0])/N;

    vector<real16> grid;
    for(int i=0 ; i<=N ; i++){
        grid.push_back( x[0] + i*dx );
    }

    return grid;
}

