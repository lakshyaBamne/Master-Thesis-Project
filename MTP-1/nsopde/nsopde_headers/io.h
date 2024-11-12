/*
    Header file to handle input output operations
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<fstream>

#include "macros.h"

using namespace std;

void writeVector(vector<real16>& x, string file_name){
    ofstream fout(file_name, ios::app);

    stringstream ss;
    for(auto i : x){
        ss << i << " ";
    }
    ss << endl;

    fout << ss.str();

    fout.close();
}

