/*
    Class to represent an initial value problem
*/
#pragma once

#include<iostream>
#include<utility>
#include<vector>
#include<functional>

#include "macros.h"

using namespace std;

vector<real16> initialize_ivp(vector<real16> grid, function<vector<real16>(vector<real16>&)> F){
    return F(grid);
}


