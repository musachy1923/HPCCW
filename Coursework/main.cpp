#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

void printMatrix(double* A, const int& M, const int& N);
void writeFile(string fileName, int Ny, int Nx, double* x);

int main(int argc, char **argv)
{
    //
    double a = 5;
    double b = 0.06;
    double epsilon = 50;
    double mu1 = 50;
    double mu2 = 0.0;
    //
    
	double dx = 1;
    double dy = 1;
    double h = dx;
    int Nx = 101;
    int Ny = 50;
    double Lx = dx*Nx;
    double Ly = dy*Ny;
    
    
    //meshgrid
    
    double * x = new double [Nx*Ny];
    double * y = new double [Nx*Ny];
    
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            x[Nx*i+j] = dx*j;
            y[Nx*i+j] = dy*i;
        }
    }
    //printMatrix(x, Ny, Nx);
    //printMatrix(y, Ny, Nx);
    
    //boundary conditions from the handout
    double dt = 0.001;
    double t = 0.001;
    int tn = t/dt; //number of time nodes
    
    double * u = new double[Nx*Ny*(tn)];
    double * v = new double[Nx*Ny*(tn)];
    
    
    
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            if (y[Nx*i+j] > Ly/2) {
                u[Nx*i+j] = 1;
            }
            else {
                u[Nx*i+j] = 0;
            }
            
            if (x[Nx*i+j] < Lx/2) {
                v[Nx*i+j] = a/2;
            }
            else {
                v[Nx*i+j] = 0;
            }
            
        }
    }
    
    
    //Neuman boundary conditions
    
    //x = 0
    //Merge with upper code when done
    
    //do for i = 1, come back later for i = 0
    //eg issue with 
    
    for (int i = 1; i < Ny+abs(Nx-Ny); i++) {
        if (i < Ny){
            //Left
            u[i*Nx] = 1/h/h*(u[1+Nx*i]+u[Nx*(i+1)]+u[Nx*(i-1)]-3*u[Nx*i]);
            //v[i*Nx] = 1/h/h*(v[1+Nx*i]+v[Nx*(i+1)]+v[Nx*(i-1)]-3*v[Nx*i]);
            //Right
            u[Nx-1+i*Nx] = 1/h/h*(u[Nx-2+Nx*i]+u[Nx-1+Nx*(i+1)]+u[Nx-1+Nx*(i-1)]-3*u[Nx-1+Nx*i]);
            //v[Nx-1+i*Nx] = 1/h/h*(v[Nx-2+Nx*i]+v[Nx-1+Nx*(i+1)]+v[Nx-1+Nx*(i-1)]-3*v[Nx-1+Nx*i]);
        }
        if (i < Nx) {
            //Lower
            u[i] = 1/h/h*(u[1+i]+u[(i+1)]+u[(i-1)]-3*u[i]);
            //v[i] = 1/h/h*(v[1+i]+v[(i+1)]+v[(i-1)]-3*v[i]);
            //Upper
            u[Nx*(Ny-1)+i] = 1/h/h*(u[Nx*(Ny-2)+i]+u[Nx*(Ny-1)+i+1]+u[Nx*(Ny-1)+i-1]-3*u[Nx*(Ny-1)+i]);
            //v[Nx*(Ny-1)+i] = 1/h/h*(v[Nx*(Ny-2)+i]+v[Nx*(Ny-1)+i+1]+v[Nx*(Ny-1)+i-1]-3*v[Nx*(Ny-1)+i]);
            
        }
        
        
}

        
        
        
        

    
    
    
    
    //printMatrix(u, Ny, Nx);
    //printMatrix(v, Ny, Nx);
    
    writeFile("x.txt", Ny, Nx, x);
    writeFile("y.txt", Ny, Nx, y);
    writeFile("u.txt", Ny, Nx, u);
    writeFile("v.txt", Ny, Nx, v);
    
    delete[] x;
    delete[] y;
    delete[] u;
    delete[] v;
}

 
void printMatrix(double* A, const int& M, const int& N) {
    
    for (int i = 0; i < M; i++)  {
        for (int j = 0; j < N; j++)  {
            cout << setw(10) << A[j*M+i];
        }
        cout << endl;
    }
    cout << endl;
}

void writeFile(string fileName, int Ny, int Nx, double * x) {
    ofstream file;
    file.open(fileName);
    for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx;  j++) {
                    file << x[Nx*i+j]<<" ";
    }
    file << "\n";
}
    file.close();
    
    
    //MPI cartesian topology
}