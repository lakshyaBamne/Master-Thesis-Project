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
    int N = 1000;

    // change the initial conditions for the problem being run
    ld L = 0l;
    ld R = 1l;
    
    pldld TIME = {0l, 2l}; // MCW
    // pldld TIME = {0l, 0.012l}; // SCW
    // pldld TIME = {0l, 0.038l}; // BLW

    ld DX = (R-L)/(N-1);

    string PROBLEM = "MCW";
    // string PROBLEM = "SCW";
    // string PROBLEM = "BLW";

    string BC = "FREE";
    // string BC = "REFLECTIVE";

    string SCHEME = "1D Lax Friedrichs Scheme (1st order)";

    ld GAMMA = 1.4l;
    ld THETA = 1.3l;
    ld CFL = 0.45l; // for 1st order scheme
}

vld make_grid();
void initialize_variables(vld& grid, vvld& V);
void extend_cells(vvld& V);
void write(vld& v, string filename);

void write_parameters();

int main(){
    vld grid = make_grid();

    write(grid, "grid.txt"); // write the grid values

    // Initialize the primitive variables based on the problem
    vvld V(3, vld(cts::N+2, 0l)); //  primitive variables
    vvld U(3, vld(cts::N+2, 0l)); // conserved variables
    vvld F(3, vld(cts::N+2, 0l)); // flux variables

    vvld Un(3, vld(cts::N+2, 0l)); // updated conserved variables

    initialize_variables(grid, V);

    write(V[0], "density.txt"); // write initial density

    ld T = cts::TIME.first;
    ld DT;
    ld max_speed;

    cout << "Serial execution starting for - 1D 1st order LF Scheme (" << cts::PROBLEM << ")" << endl;

    // time calculation
    auto time_start = chrono::high_resolution_clock::now();

    while(T <= cts::TIME.second){
        cout << "T= " << T << " | DT= " << DT << endl; // log

        // first calculate the time step DT
        max_speed = 0;
        for(int i=1 ; i<=cts::N ; i++){
            max_speed = max(max_speed, abs(V[1][i])+sqrtl(cts::GAMMA*V[2][i]/V[0][i]));
        }
        
        DT = cts::CFL*cts::DX/max_speed;

        // update time
        T += DT;

        // extend cells
        extend_cells(V);

        // Using the Primitive variables initialize other variables
        for(int i=0 ; i<=cts::N+1 ; i++){
            // conserved variables
            U[0][i] = V[0][i];
            U[1][i] = V[0][i]*V[1][i];
            U[2][i] = (V[2][i]/(cts::GAMMA-1l)) + (0.5l*V[0][i]*pow(V[1][i],2));

            // flux variables
            F[0][i] = V[0][i]*V[1][i];
            F[1][i] = V[2][i] + V[0][i]*pow(V[1][i],2);
            F[2][i] = V[1][i]*( (cts::GAMMA*V[2][i]/(cts::GAMMA-1l)) + (0.5l*V[0][i]*pow(V[1][i],2)) );
        }

        // update the conserved variables
        for(int u=0 ; u<=2 ; u++){
            for(int i=1 ; i<=cts::N ; i++){
                Un[u][i] = 0.5l*(U[u][i+1]+U[u][i-1]) - 0.5l*(DT/cts::DX)*(F[u][i+1] - F[u][i-1]);
            }
        }

        // update the primitive variables for the next iteration
        for(int i=1 ; i<cts::N ; i++){
            V[0][i] = Un[0][i];
            V[1][i] = Un[1][i]/Un[0][i];
            V[2][i] = (cts::GAMMA-1l)*(Un[2][i] - 0.5l*(pow(Un[1][i],2)/Un[0][i]));
        }

    }

    // time calculation
    auto time_end = chrono::high_resolution_clock::now();

    chrono::duration<long double> elapsed_time = time_end - time_start;
    cout << "Runtime for serial program : " << elapsed_time.count() << " seconds." << endl;

    write(V[0], "density.txt"); // write final density

    // write the parameters used in the problem
    write_parameters();

    return 0;
}

vld make_grid(){
    vld grid(cts::N+2);

    for(int i=1 ; i<=cts::N ; i++){
        grid[i] = cts::L + (i-1)*cts::DX;
    }

    return grid;
}

void initialize_variables(vld& grid, vvld& V){
    // Initialize the primitive variables firs
    if(cts::PROBLEM == "MCW"){
        for(int i=1 ; i<=cts::N ; i++){
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
        for(int i=1 ; i<=cts::N ; i++){
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
        for(int i=1 ; i<=cts::N ; i++){
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

void extend_cells(vvld& V){
    if(cts::BC == "FREE"){
        for(int i=0 ; i<=2 ; i++){
            V[i][0] = V[i][1];
            V[i][cts::N+1] = V[i][cts::N];
        }
    }
    else if(cts::BC == "REFLECTIVE"){
        // density
        V[0][0] = V[0][1];
        V[0][cts::N+1] = V[0][cts::N];

        // velocity
        V[1][0] = -V[1][1];
        V[1][cts::N+1] = -V[1][cts::N];

        // pressure
        V[2][0] = V[2][1];
        V[2][cts::N+1] = V[2][cts::N];
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

    cout << "[LOG] file write successful in file [" << filename << "]" << endl;
}

void write_parameters(){
    ofstream fout("parameters.txt", ios::app);

    fout << cts::N << endl;
    fout << cts::PROBLEM << endl;
    fout << cts::BC << endl;
    fout << cts::TIME.first << " " << cts::TIME.second << endl;
    fout << cts::SCHEME << endl;

    fout.close();
    cout << "[LOG] file write successful in file [parameters.txt]" << endl;
}

