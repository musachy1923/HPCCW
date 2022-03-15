#include <iostream>
#include <iomanip>
#include <fstream>
//useful macros
#define f1(x) (epsilon*u[x]*(1-u[x])*(u[x]-(v[x]+b)/a))*dt+u[x]
#define f2(x) (u[x]*u[x]*u[x]-v[x])*dt+v[x]
#define neuman(x, a, b, c, d) (1/h/h*(x[a]+x[b]+x[c]-3*x[d]))*dt
#define neuman2(x, a, b, c) (1/h/h*(x[a]+x[b]-2*x[c]))
#define F77NAME(x) x##_
using namespace std;


extern "C" {void F77NAME(dgemv) (const char& trans,
                             const int& M, const int& N,
                             const double& alpha, double * A, const int& lda, 
                             const double * x, const int& incx,
                             const double& beta, double * Y, const int& incy);
}
                 
                             
void printMatrix(double* A, const int& M, const int& N);
void writeFile(string fileName, int Ny, int Nx, int tn, double* x);
void singleUpperBandedOptimised(int n, double* H, double alpha, double beta);
void singleUpperBanded(int nsv, double* H, double alpha, double beta);
void arrayPrint(double * a, int n);

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
    int Nx = 10;
    int Ny = 10;
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
    double t = 0.002;
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
for (int k = 1; k < tn; k++) {
   int ts = Nx*Ny*(k-1);
   int t2 = Nx*Ny*(k);
        for (int i = 1; i < Ny+abs(Nx-Ny); i++) {
        if (i < Ny){
            //Left
            u[t2+i*Nx] = (1/h/h*(u[ts+1+Nx*i]+u[ts+Nx*(i+1)]+u[ts+Nx*(i-1)]-3*u[ts+Nx*i]))*dt+f1(ts+Nx*i);
            v[t2+i*Nx] = (1/h/h*(v[ts+1+Nx*i]+v[ts+Nx*(i+1)]+v[ts+Nx*(i-1)]-3*v[ts+Nx*i]))*dt+f2(ts+Nx*i);
            //Right
            
            u[t2+Nx-1+i*Nx] = (1/h/h*(u[ts+Nx-2+Nx*i]+u[ts+Nx-1+Nx*(i+1)]+u[ts+Nx-1+Nx*(i-1)]-3*u[ts+Nx-1+Nx*i]))*dt+f1(ts+Nx-1+Nx*i);
            v[t2+Nx-1+i*Nx] = (1/h/h*(v[ts+Nx-2+Nx*i]+v[ts+Nx-1+Nx*(i+1)]+v[ts+Nx-1+Nx*(i-1)]-3*v[ts+Nx-1+Nx*i]))*dt+f2(ts+Nx-1+Nx*i);
        
          }
        if (i < Nx) {
            //Lower
            u[t2+i] = (1/h/h*(u[ts+1+i]+u[ts+(i+1)]+u[ts+(i-1)]-3*u[ts+i]))*dt+f1(ts+i);
            v[t2+i] = (1/h/h*(v[ts+1+i]+v[ts+(i+1)]+v[(ts+i-1)]-3*v[ts+i]))*dt+f2(ts+i);
            //Upper
            u[t2+Nx*(Ny-1)+i] = (1/h/h*(u[ts+Nx*(Ny-2)+i]+u[ts+Nx*(Ny-1)+i+1]+u[ts+Nx*(Ny-1)+i-1]-3*u[ts+Nx*(Ny-1)+i]))*dt+f1(ts+Nx*(Ny-1)+i);
            v[t2+Nx*(Ny-1)+i] = (1/h/h*(v[ts+Nx*(Ny-2)+i]+v[ts+Nx*(Ny-1)+i+1]+v[ts+Nx*(Ny-1)+i-1]-3*v[ts+Nx*(Ny-1)+i]))*dt+f2(ts+Nx*(Ny-1)+i);
         }
        
}

u[t2] = (1/h/h*(u[ts+1]+u[ts+Nx]-2*u[ts]))*dt+f1(ts);
u[t2+Nx-1] = (1/h/h*(u[ts+Nx-1]+u[ts+Nx-1]-2*u[ts+Nx-1]))*dt+f1(ts+Nx-1);
u[t2+Nx*(Ny-1)] = (1/h/h*(u[ts+Nx*(Ny-1)]+u[ts+Nx*(Ny-1)]-2*u[ts+Nx*(Ny-1)]))*dt+f1(ts+Nx*(Ny-1));
u[t2+Nx*Ny] = (1/h/h*(u[ts+Nx*Ny]+u[ts+Nx*Ny]-2*u[ts+Nx*Ny]))*dt+f1(ts+Nx*Ny);

v[t2] = (1/h/h*(v[ts]+v[ts]-2*v[ts]))*dt+f2(ts);
v[t2+Nx-1] = (1/h/h*(v[ts+Nx-1]+v[ts+Nx-1]-2*v[ts+Nx-1]))*dt+f2(ts+Nx-1);
v[t2+Nx*(Ny-1)] = (1/h/h*(v[ts+Nx*(Ny-1)]+v[ts+Nx*(Ny-1)]-2*v[ts+Nx*(Ny-1)]))*dt+f2(ts+Nx*(Ny-1));
v[t2+Nx*Ny] = (1/h/h*(v[ts+Nx*Ny]+v[ts+Nx*Ny]-2*v[ts+Nx*Ny]))*dt+f2(ts+Nx*Ny);

//for this timestep the conditions at boundaries are set, just need to evaluate the inner grid


double * A = new double[(Nx-2)*(Ny-2)];
double * y = new double[Nx-2];
double * u_j = new double[(Nx-2)];
double * u_jp1 = new double[(Nx-2)];
double * u_jm1 = new double[(Nx-2)];
double * fu1 = new double[(Nx-2)];
double * v_j = new double[(Nx-2)];
double * v_jp1 = new double[(Nx-2)];
double * v_jm1 = new double[(Nx-2)];
double * fv2 = new double[(Nx-2)];
singleUpperBanded(Ny-2, A, -4, 1);

for (int j = 1; j < Ny-1; j++) {
    for (int i = 1; i < Nx-1; i++) {
        u_j[i-1] = u[i+Nx*j];
        u_jp1[i-1] = u[i+Nx*(j+1)]+u[i+Nx*(j-1)]+f1(i+Nx*j);
        v_j[i-1] = v[i+Nx*j];
        v_jp1[i-1] = v[i+Nx*(j+1)]+v[i+Nx*(j-1)]+f2(i+Nx*j);
    }
    F77NAME(dgemv) ('N', Ny-2, Nx-2, 1, A, Nx-2, u_j, 1, 0, y, 1);
    
    for (int i = 1; i < Nx-1; i++) {
        u[t2+Nx*(j)+i] = y[i-1]+u_jp1[i-1];
        fu1[i-1] = y[i-1]+u_jp1[i-1];
    }
    
}





arrayPrint(fu1, Nx-2);


//printMatrix(A, Ny-2, Nx-2);
} 







    
    
    //printMatrix(u, Ny, Nx);
    //printMatrix(v, Ny, Nx);
    
    writeFile("x.txt", Ny, Nx, tn, x);
    writeFile("y.txt", Ny, Nx, tn, y);
    writeFile("u.txt", Ny, Nx, tn, u);
    writeFile("v.txt", Ny, Nx, tn, v);
    
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

void writeFile(string fileName, int Ny, int Nx, int tn,  double * x) {
    ofstream file;
    file.open(fileName);
    for (int k = 0; k < tn; k++) {
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx;  j++) {
                    
                file << x[Nx*Ny*k+Nx*i+j]<<" ";
    }
    file << "\n";
        
    }
    file << "-----------------------------------------------------------------------------------------------" << endl;
}
    file.close();
    
    
    //MPI cartesian topology
}

void singleUpperBandedOptimised(int n, double* H, double alpha, double beta) {
    const int ldh = 2;      //3 Diagonal and upper diagonal
    H[1] = alpha;
    for (int i = 1; i < n; ++i) {
        H[i*ldh] = beta;
        H[i*ldh + 1] = alpha;
        //H[i*ldh+2]= gamma;
    }
}

void singleUpperBanded(int nsv, double* H, double alpha, double beta) {
    H[0] = alpha;
    for (int i = 1; i < nsv; ++i) {
        H[i*nsv + i - 1] = beta;
        H[i*nsv + i] = alpha;
    }
    //printMatrix(H, nsv, nsv);
}

void arrayPrint(double * a, int n) {
    for (int i = 0; i < n; i++) {
    
        cout << a[i] << endl;    
    
    }
}