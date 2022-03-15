#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#define u_n_1(x) for (1/h/h*(u[x]+u[x]+u[x]-3*u[x])+epsilon*u[x]*(1-u[x])*(u[x]-(v[x]+b)/a))*dt+u[x];
#define v_n_1(x) for (1/h/h*(v[x]+v[x]+v[x]-3*v[x])+u[x]*u[x]*u[x]-v[x)*dt+v[x];
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
for (int k = 0; k < tn; k++) {
   ts = Nx*Ny*k;
   t2 = Nx*Ny*(k+1)
        for (int i = 1; i < Ny+abs(Nx-Ny); i++) {
        if (i < Ny){
            //Left
            u[t2+i*Nx] = (1/h/h*(u[ts+1+Nx*i]+u[ts+Nx*(i+1)]+u[ts+Nx*(i-1)]-3*u[ts+Nx*i])+epsilon*u[ts+Nx*i]*(1-u[ts+Nx*i])*(u[ts+Nx*i]-(v[ts+i*Nx]+b)/a))*dt+u[ts+Nx*i];
            v[t2+i*Nx] = (1/h/h*(v[ts+1+Nx*i]+v[ts+Nx*(i+1)]+v[ts+Nx*(i-1)]-3*v[ts+Nx*i])+u[ts+Nx*i]*u[ts+Nx*i]*u[ts+Nx*i]-v[ts+i*Nx])*dt+v[ts+Nx*i];
            //Right
            u[t2+Nx-1+i*Nx] = (1/h/h*(u[ts+Nx-2+Nx*i]+u[ts+Nx-1+Nx*(i+1)]+u[ts+Nx-1+Nx*(i-1)]-3*u[ts+Nx-1+Nx*i])+epsilon*u[ts+Nx-1+Nx*i]*(1-u[ts+Nx-1+Nx*i])*(u[ts+Nx-1+Nx*i]-(v[ts+Nx-1+Nx*i]+b)/a))*dt+u[ts+Nx-1+Nx*i];
            v[t2+Nx-1+i*Nx] = (1/h/h*(v[ts+Nx-2+Nx*i]+v[ts+Nx-1+Nx*(i+1)]+v[ts+Nx-1+Nx*(i-1)]-3*v[ts+Nx-1+Nx*i])-3*v[ts+Nx-1+Nx*i])+u[ts+Nx-1+Nx*i]*u[ts+Nx-1+Nx*i]*u[ts+Nx-1+Nx*i]-v[ts+Nx-1+Nx*i])*dt+v[ts+Nx-1+Nx*i];
        }
        if (i < Nx) {
            //Lower
            u[t2+i] = (1/h/h*(u[ts+1+i]+u[ts+(i+1)]+u[ts+(i-1)]-3*u[ts+i])+epsilon*u[ts+i]*(1-u[ts+i])*(u[ts+i]-(v[ts+i]+b)/a))*dt+u[ts+i];
            v[t2+i] = (1/h/h*(v[ts+1+i]+v[ts+(i+1)]+v[(ts+i-1)]-3*v[ts+i])+u[ts+i]*u[ts+i]*u[ts+i]-v[ts+i])*dt+v[ts+i];
            //Upper
            u[t2+Nx*(Ny-1)+i] = (1/h/h*(u[ts+Nx*(Ny-2)+i]+u[ts+Nx*(Ny-1)+i+1]+u[ts+Nx*(Ny-1)+i-1]-3*u[ts+Nx*(Ny-1)+i])+epsilon*u[ts+Nx*(Ny-1)+i]*(1-u[ts+Nx*(Ny-1)+i])*(u[ts+Nx*(Ny-1)+i]-(v[ts+Nx*(Ny-1)+i]+b)/a))*dt+u[ts+Nx*(Ny-1)+i];
            v[t2+Nx*(Ny-1)+i] = (1/h/h*(v[ts+Nx*(Ny-2)+i]+v[ts+Nx*(Ny-1)+i+1]+v[ts+Nx*(Ny-1)+i-1]-3*v[ts+Nx*(Ny-1)+i])+u[ts+Nx*(Ny-1)+i]*u[ts+Nx*(Ny-1)+i]*u[ts+Nx*(Ny-1)+i]-v[ts+Nx*(Ny-1)+i])*dt+v[ts+Nx*(Ny-1)+i];;
        }
}
u[t2] = (1/h/h*(u[ts+1]+u[ts+Nx]-2*u[ts])+epsilon*u[ts]*(1-u[ts])*(u[ts]-(v[ts]+b)/a))*dt+u[ts];
u[t2+Nx-1] = (1/h/h*(u[ts+Nx-1]+u[ts+Nx-1]-2*u[ts+Nx-1])+epsilon*u[ts+Nx-1]*(1-u[ts+Nx-1])*(u[ts+Nx-1]-(v[ts+Nx-1]+b)/a))*dt+u[ts+Nx-1];
u[t2+Nx*(Ny-1)] = (1/h/h*(u[ts+Nx*(Ny-1)]+u[ts+Nx*(Ny-1)]-2*u[ts+Nx*(Ny-1)])+epsilon*u[ts+Nx*(Ny-1)]*(1-u[ts+Nx*(Ny-1)])*(u[ts+Nx*(Ny-1)]-(v[ts+Nx*(Ny-1)]+b)/a))*dt+u[ts+Nx*(Ny-1)];
u[t2+Nx*Ny] = (1/h/h*(u[ts+Nx*Ny]+u[ts+Nx*Ny]-2*u[ts+Nx*Ny])+epsilon*u[ts+Nx*Ny]*(1-u[ts+Nx*Ny])*(u[ts+Nx*Ny]-(v[ts+Nx*Ny]+b)/a))*dt+u[ts+Nx*Ny];

v[t2] = (1/h/h*(v[ts]+v[ts]-2*v[ts])+u[ts]*u[ts]*u[ts]-v[ts])*dt+v[ts];
v[t2+Nx-1] = (1/h/h*(v[ts+Nx-1]+v[ts+Nx-1]-2*v[ts+Nx-1])+u[ts+Nx-1]*u[ts+Nx-1*u[ts+Nx-1]-v[ts+Nx-1])*dt+v[ts+Nx-1];

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