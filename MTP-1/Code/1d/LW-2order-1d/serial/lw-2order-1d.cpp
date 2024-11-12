/*
    Lax Wendorff Scheme is a 2nd order explicit finite difference scheme 
    -> second order taylor series expansion is used for both spatial and time domains
*/

#include<iostream>
#include<vector>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<string>
#include<utility>
#include<sstream>
#include<chrono>

using namespace std;

#define ld long double
#define vld vector<long double>
#define vvld vector< vector<long double> >
#define pldld pair<long double, long double>
#define pss pair<string,string>

// overloading the output stream operator to log a vector
ostream & operator<<(ostream &out, vld &v){
    out << "-> ";
    for(auto i : v){
        out << i << ", ";
    }

    return out;
}

// constants used in the program
namespace cts{
    int N = 100;

    // change the initial conditions for the problem being run
    ld L = 0l;
    ld R = 1l;
    
    pldld TIME = {0l, 2l}; // MCW
    // pldld TIME = {0l, 0.012l}; // SCW
    // pldld TIME = {0l, 0.038l}; // BLW

    ld DX = (R-L)/N;

    string PROBLEM = "MCW";
    // string PROBLEM = "SCW";
    // string PROBLEM = "BLW";

    string BC = "FREE";
    // string BC = "REFLECTIVE";

    ld GAMMA = 1.4l;
    ld THETA = 1.3l;
    ld CFL = 0.45l; // for 2nd order scheme
}

vld make_grid();
void initialize_variables(vld& grid, vvld& V, vvld& U, vvld& F);
void extend_cells(vvld& V, vvld& U, vvld& F);

void write(vld& v, string filename);


int main(){
    // Initialize the computational grid 
    vld grid = make_grid();

    write(grid, "grid.txt");

    // Initialize the primitive variables based on the problem
    vvld V(3, vld(cts::N+1, 0l)); //  primitive variables
    vvld U(3, vld(cts::N+1, 0l)); // conserved variables
    vvld F(3, vld(cts::N+1, 0l)); // flux variables

    vvld Unew(3, vld(cts::N+1, 0l)); // updating conserved variables

    initialize_variables(grid, V, U, F);

    write(V[0], "density.txt");

    ld T = cts::TIME.first;
    ld DT;
    ld max_speed=0;
    ld a; // characteristic speed

    cout << "Serial execution starting for - 1D 1st order LF Scheme (" << cts::PROBLEM << ")" << endl;

    // time calculation
    auto time_start = chrono::high_resolution_clock::now();

    while(T <= cts::TIME.second){
        cout << "T= " << T << " | DT= " << DT << endl; // log

        // first calculate the time step DT
        max_speed = 0;
        for(int i=1 ; i<=cts::N-1 ; i++){
            max_speed = max(max_speed, abs(V[1][i])+sqrtl(cts::GAMMA*V[2][i]/V[0][i]));
        }
        
        DT = cts::CFL*cts::DX/max_speed;

        // update time
        T += DT;

        /*
            Code for Lax-Wendorff scheme (2nd order taylor series expansion)
        */
        // Using the Primitive variables initialize other variables
        for(int i=1 ; i<=cts::N-1 ; i++){
            // conserved variables
            U[0][i] = V[0][i];
            U[1][i] = V[0][i]*V[1][i];
            U[2][i] = (V[2][i]/(cts::GAMMA-1l)) + (0.5l*V[0][i]*pow(V[1][i],2));

            // flux variables
            F[0][i] = V[0][i]*V[1][i];
            F[1][i] = V[2][i] + V[0][i]*pow(V[1][i],2);
            F[2][i] = V[1][i]*( (cts::GAMMA*V[2][i]/(cts::GAMMA-1l)) + (0.5l*V[0][i]*pow(V[1][i],2)) );
        }

        // extend cells
        extend_cells(V, U, F);

        // update the conserved variables
        for(int u=0 ; u<=2 ; u++){
            for(int i=1 ; i<=cts::N-1 ; i++){
                // character speed calculation

                // Update equations
                
            }
        }

        // update the primitive variables for the next iteration
        for(int i=1 ; i<cts::N-1 ; i++){
            V[0][i] = Unew[0][i];
            V[1][i] = Unew[1][i]/Unew[0][i];
            V[2][i] = (cts::GAMMA-1l)*(Unew[2][i] - 0.5l*(pow(Unew[1][i],2)/Unew[0][i]));
        }

    }

    // time calculation
    auto time_end = chrono::high_resolution_clock::now();

    chrono::duration<long double> elapsed_time = time_end - time_start;
    cout << "Runtime for serial program : " << elapsed_time.count() << " seconds." << endl;

    write(U[0], "density.txt");


    return 0;
}

vld make_grid(){
    vld grid(cts::N+1);

    for(int i=0 ; i<=cts::N ; i++){
        grid[i] = cts::L + i*cts::DX;
    }

    return grid;
}

void initialize_variables(vld& grid, vvld& V, vvld& U, vvld&F){
    // Initialize the primitive variables firs
    if(cts::PROBLEM == "MCW"){
        for(int i=1 ; i<=cts::N-1 ; i++){
            if(grid[i] <= 0.3){
                V[0][i] = 1.4l;
                V[1][i] = 0.1l;
                V[2][i] = 1l;
            }
            else{
                V[0][i] = 1l;
                V[1][i] = 0.1l;
                V[2][i] = 1l;
            }
        }
    }
    else if(cts::PROBLEM == "SCW"){
        for(int i=1 ; i<=cts::N-1 ; i++){
            if(grid[i] <= 0.8){
                V[0][i] = 1l;
                V[1][i] = -19.59745l;
                V[2][i] = 1000l;
            }
            else{
                V[0][i] = 1l;
                V[1][i] = -19.59745l;
                V[2][i] = 0.01l;
            }
        }
    }
    else if(cts::PROBLEM == "BLW"){
        for(int i=0 ; i<grid.size() ; i++){
            V[0][i] = 1.0l;
            V[1][i] = 0.0l;
            
            if( grid[i] < 0.1 ){
                V[2][i] = 1000.0l;
            }
            else if( grid[i]>=0.1 && grid[i]<=0.9 ){
                V[2][i] = 0.01l;
            }
            else{
                V[2][i] = 100.0l;
            }
        }
    }
    else{
        cout << "[ERROR] Enter correct problem" << endl;
    }
}

void extend_cells(vvld& V, vvld& U, vvld& F){
    if(cts::BC == "FREE"){
        for(int i=0 ; i<=2 ; i++){
            V[i][0] = V[i][1];
            V[i][cts::N] = V[i][cts::N-1];

            U[i][0] = U[i][1];
            U[i][cts::N] = U[i][cts::N-1];
        
            F[i][0] = F[i][1];
            F[i][cts::N] = F[i][cts::N-1];
        }
    }
    else if(cts::BC == "REFLECTIVE"){
        // V[0], U[0], F[0]
        V[0][0] = V[0][1];
        V[0][cts::N] = V[0][cts::N-1];
        
        U[0][0] = U[0][1];
        U[0][cts::N] = U[0][cts::N-1];
    
        F[0][0] = -F[0][1];
        F[0][cts::N] = -F[0][cts::N-1];

        // V[1], U[1], F[1]
        V[1][0] = -V[1][1];
        V[1][cts::N] = -V[1][cts::N-1];

        U[1][0] = -U[1][1];
        U[1][cts::N] = -U[1][cts::N-1];
    
        F[1][0] = -F[1][1];
        F[1][cts::N] = -F[1][cts::N-1];
        
        
        // V[2], U[2], F[2]
        V[2][0] = V[2][1];
        V[2][cts::N] = V[2][cts::N-1];

        U[2][0] = U[2][1];
        U[2][cts::N] = U[2][cts::N-1];
    
        F[2][0] = -F[2][1];
        F[2][cts::N] = -F[2][cts::N-1];
    }
    else{
        cout << "[ERROR] Please enter correct boundary conditions" << endl;
    }
}

void write(vld& vec, string filename){
    stringstream ss;
    for(auto v : vec){
        ss << v << " ";
    }

    ofstream fout(filename, ios::app);
    fout << ss.str() << endl;

    fout.close();
}
