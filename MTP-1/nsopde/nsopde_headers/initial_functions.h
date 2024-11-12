/*
    Header file contains some example functions to initialize the problems
*/
#pragma once

#include<iostream>
#include<vector>
#include<cmath>
#include<random>

using namespace std;

#include "macros.h"

// trigonometric functions which are periodic and always smooth
vr16 sine_wave(vr16& grid){
    vr16 wave;

    for(auto i : grid){
        wave.push_back(sin(i));
    }

    return wave;
}

vr16 cosine_wave(vr16& grid){
    vr16 wave;
    
    for(auto i : grid){
        wave.push_back(cos(i));
    }

    return wave;
}

// periodic function contains discontinuities at integer values
vr16 square_wave(vr16& grid){
    vr16 wave;

    for(auto i : grid){
        if((int)floor(i)%2 == 0){
            wave.push_back(1);
        }
        else{
            wave.push_back(-1);
        }
    }

    return wave;
}

vr16 signum_at_0(vr16& grid){
    vr16 wave;
    
    for(auto i : grid){
        if(i < 0){
            wave.push_back(-1);
        }
        else{
            wave.push_back(1);
        }
    }

    return wave;
}
