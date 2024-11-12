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
#define vvvld vector< vector< vector<long double> > >
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
    int N = 500;

    // change the initial conditions for the problem being run
    ld Xl = -0.2l;
    ld Xr = 0.2l;

    ld Yl = 0l;
    ld Yr = 0.8l;
    
    ld DX = (Xr-Xl)/(N-1);
    ld DY = (Yr-Yl)/(N-1);

    pldld TIME = {0l, 0.5l}; // MCW

    string PROBLEM = "MCW";

    string BC = "FREE";

    string SCHEME = "2D Lax Friedrichs Scheme (1st order)";

    ld GAMMA = 1.4l;
    ld THETA = 1.3l;
    ld CFL = 0.45l; // for 1st order scheme
}

vvld make_grids();
void initialize_variables(vld& X, vld& Y, vvvld& V);
void extend_cells(vvvld& V);
void write_matrix(vvld& matrix, string filename);
void write_vector(vld& vec, string filename);
void write_parameters();

int main(){
    vvld grids = make_grids();

    vvvld V(4, vvld(cts::N+2, vld(cts::N+2, 0l)));
    initialize_variables(grids[0], grids[1], V);

    // write output
    write_parameters();
    write_vector(grids[0], "grid.txt");
    write_vector(grids[1], "grid.txt");
    write_matrix(V[0], "initial_density.txt");

    // now we can start producing the solutions
    vvvld U(4, vvld(cts::N+2, vld(cts::N+2, 0l))); // conserved variables
    vvvld F(4, vvld(cts::N+2, vld(cts::N+2, 0l))); // flux variables 1 
    vvvld G(4, vvld(cts::N+2, vld(cts::N+2, 0l))); // flux variables 2

    vvvld Un(4, vvld(cts::N+2, vld(cts::N+2, 0l))); // updated conserved variables

    ld T = cts::TIME.first;
    ld DT;
    ld amax=0, bmax=0, acos_speed;

    cout << "Serial execution starting for - 2D 1st order LF Scheme (" << cts::PROBLEM << ")" << endl;

    // time calculation
    auto time_start = chrono::high_resolution_clock::now();

    while(T <= cts::TIME.second){
        cout << "T= " << T << " | DT= " << DT << endl; // log

        // first calculate the time step DT
        for(int i=1 ; i<=cts::N ; i++){
            for(int j=1 ; j<=cts::N ; j++){
                acos_speed = sqrtl(cts::GAMMA*V[3][i][j]/V[0][i][j]);
                amax = max(amax, abs(V[1][i][j])+acos_speed);
                bmax = max(bmax, abs(V[2][i][j])+acos_speed);
            }
        }

        DT = cts::CFL*min(cts::DX/amax, cts::DY/bmax);

        // update time
        T += DT;

        // extend cells
        extend_cells(V);

        // Using the Primitive variables initialize other variables
        for(int i=0 ; i<=cts::N+1 ; i++){
            for(int j=0 ; j<=cts::N+1 ; j++){
                // conserved variables
                U[0][i][j] = V[0][i][j];
                U[1][i][j] = V[0][i][j]*V[1][i][j];
                U[2][i][j] = V[0][i][j]*V[2][i][j];
                U[3][i][j] = (V[3][i][j]/(cts::GAMMA-1)) + 0.5*V[0][i][j]*(pow(V[1][i][j],2)+pow(V[2][i][j],2));

                // flux 1
                F[0][i][j] = U[1][i][j];
                F[1][i][j] = V[0][i][j]*pow(V[1][i][j],2) + V[3][i][j];
                F[2][i][j] = V[0][i][j]*V[1][i][j]*V[2][i][j];
                F[3][i][j] = V[1][i][j]*(U[3][i][j] + V[3][i][j]);

                // flux 2
                G[0][i][j] = U[2][i][j];
                G[1][i][j] = V[0][i][j]*V[1][i][j]*V[2][i][j];
                G[2][i][j] = V[0][i][j]*pow(V[2][i][j],2) + V[3][i][j];
                G[3][i][j] = V[2][i][j]*(U[3][i][j] + V[3][i][j]);
            }
        }

        // update the conserved variables
        for(int i=1 ; i<=cts::N ; i++){
            for(int j=1 ; j<=cts::N ; j++){
                for(int u=0 ; u<4 ; u++){
                    Un[u][i][j] = 0.25*(U[u][i+1][j]+U[u][i-1][j]+U[u][i][j+1]+U[u][i][j-1]) - (DT/(2*cts::DX))*(F[u][i+1][j]-F[u][i-1][j]) - (DT/(2*cts::DY))*(G[u][i][j+1]-G[u][i][j-1]);
                }
            }
        }
        

        // update the primitive variables for the next iteration
        for(int i=1 ; i<=cts::N ; i++){
            for(int j=1 ; j<=cts::N ; j++){
                V[0][i][j] = Un[0][i][j];
                V[1][i][j] = Un[1][i][j]/Un[0][i][j];
                V[2][i][j] = Un[2][i][j]/Un[0][i][j];
                V[3][i][j] = (cts::GAMMA-1)*(0.5*Un[0][i][j]*(pow(V[1][i][j],2)+pow(V[2][i][j],2)));
            }
        }

    }

    // time calculation
    auto time_end = chrono::high_resolution_clock::now();

    chrono::duration<long double> elapsed_time = time_end - time_start;
    cout << "Runtime for serial program : " << elapsed_time.count() << " seconds." << endl;

    write_matrix(V[0], "final_density.txt"); // write final density

    return 0;
}

vvld make_grids(){
    vld gridX(cts::N+2, 0l);
    vld gridY(cts::N+2, 0l);

    for(int i=1 ; i<=cts::N ; i++){
        gridX[i] = cts::Xl + (i-1)*cts::DX;
        gridY[i] = cts::Yl + (i-1)*cts::DY;
    }

    return {gridX, gridY};
}

void initialize_variables(vld& X, vld& Y, vvvld& V){
    if(cts::PROBLEM == "MCW"){
        for(int i=1 ; i<=cts::N ; i++){
            for(int j=1 ; j<=cts::N ; j++){
                if( (X[i]>-0.1l && X[i]<0.1l && Y[j]>0l && Y[j]<0.02l) || (X[i]>-0.02l && X[i]<0.02l && Y[j]>0.02l && Y[j]<0.1l) || (pow(X[i]+0.02,2) + pow(Y[j]-0.02,2) < pow(0.08,2)) || (pow(X[i]-0.02,2) + pow(Y[j]-0.02,2) < pow(0.08,2)) ){
                    V[0][i][j] = 1.4l;
                    V[1][i][j] = 0l;
                    V[2][i][j] = 0.2l;
                    V[3][i][j] = 1.0l;
                }
                else{
                    V[0][i][j] = 1.0l;
                    V[1][i][j] = 0l;
                    V[2][i][j] = 0.2l;
                    V[3][i][j] = 1.0l;
                }
            }
        }
    }
    else{
        cout << "[ERROR] Enter correct problem" << endl;
    }
}

void extend_cells(vvvld& V){
    if(cts::BC == "FREE"){
        for(int i=1 ; i<=cts::N ; i++){
            for(int u=0 ; u<4 ; u++){
                V[u][0][i] = V[u][1][i];
                V[u][cts::N+1][i] = V[u][cts::N][i];

                V[u][i][0] = V[u][i][1];
                V[u][i][cts::N+1] = V[u][i][cts::N];
            }
        }
    }
    else{
        cout << "[ERROR] Enter correct boundary conditions" << endl;
    }
}

void write_vector(vld& vec, string filename){
    ofstream fout(filename, ios::app);

    for(int i=1 ; i<=cts::N ; i++){
        fout << vec[i] << " ";
    }
    fout << endl;

    fout.close();

    cout << "[LOG] file write successful in file [" << filename << "]" << endl;
}

void write_matrix(vvld& matrix, string filename){
    ofstream fout(filename, ios::app);

    for(int j=cts::N ; j>=1 ; j--){
        for(int i=1 ; i<=cts::N ; i++){
            fout << matrix[i][j] << " ";
        }
        fout << endl;
    }

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
