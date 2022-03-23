#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
//useful macros
//export OMP NUM_THREADS
#define laplacianu(a, b, c, d) ((u[a]+u[b]+u[c]-3*u[d]))
#define laplacianv(a, b, c, d) ((v[a]+v[b]+v[c]-3*v[d]))
#define laplacianuc(a, b, c) ((u[a]+u[b]-2*u[c]))
#define laplacianvc(a, b, c) ((v[a]+v[b]-2*v[c]))
#define f1(x) (epsilon*u[x]*(1-u[x])*(u[x]-(v[x]+b)*diva))*dt
#define f2(x) (u[x]*u[x]*u[x]-v[x])*dt


#define F77NAME(x) x##_
using namespace std;

void printMatrix(double* A, const int& M, const int& N);
void writeFile(string fileName, int Ny, int Nx, int tn, double* x);
void arrayPrint(double * a, int n);

int main(int argc, char **argv)
{
    //
    double a = 0.75;
    double b = 0.06;
    double epsilon = 50;
    double mu1 = 5;
    double mu2 = 0.0;
    //
	double dx = 1;
    double dy = 1;
    double h = dx;
    int Nx = 50;
    int Ny = 50;
    

    double Lx = dx*(Nx-1);
    double Ly = dy*(Ny-1);
    double dt = 0.001;
    double t =  100;
    int tn = t/dt; //number of time nodes
    
    double hmu1dt = mu1/h/h*dt;
    double hmu2dt = mu2/h/h*dt;
    //meshgrid
    double * x = new double [Nx*Ny];
    double * y = new double [Nx*Ny];
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            x[Nx*i+j] = dx*j;
            y[Nx*i+j] = dy*i;
        }
    }
    writeFile("xr.txt", Ny, Nx, 1, x);
    writeFile("yr.txt", Ny, Nx, 1, y);
    //boundary conditions from the handout
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
    //arrayPrint(u, Nx*Ny);
    //Neuman boundary conditions
    //x = 0
    //Merge with upper code when done
    //do for i = 1, come back later for i = 0
    //eg issue with 
int NxNy = Nx*Ny;
int Nx2 = 2*Nx;
double diva = 1/a;

for (int k = 1; k < tn; k++) {
   int ts = Nx*Ny*(k-1);
   int t2 = Nx*Ny*(k);
   int i = 0;
   int j = 0;
   #pragma omp parallel for private(i)
        for (i = 1; i < Ny+abs(Nx-Ny); i++) {
            int Nxi = i*Nx;

        if (i < Ny){
            //Left
            u[t2+Nxi] = (hmu1dt*laplacianu(ts+Nxi+1, ts+Nxi+Nx, ts+Nxi-Nx, ts+Nxi)+f1(ts+Nxi))+u[ts+Nxi];
            v[t2+Nxi] = (hmu2dt*laplacianv(ts+Nxi+1, ts+Nxi+Nx, ts+Nxi-Nx, ts+Nxi)+f2(ts+Nxi))+v[ts+Nxi];
            //Right
            u[t2+Nxi+Nx-1] = (hmu1dt*laplacianu(ts+Nxi+Nx-2, ts+Nxi+Nx+Nx-1, ts+Nxi-Nx+Nx-1, ts+Nxi+Nx-1)+f1(ts+Nxi+Nx-1))+u[ts+Nxi+Nx-1];
            v[t2+Nxi+Nx-1] = (hmu2dt*laplacianv(ts+Nxi+Nx-2, ts+Nxi+Nx+Nx-1, ts+Nxi-Nx+Nx-1, ts+Nxi+Nx-1)+f2(ts+Nxi+Nx-1))+v[ts+Nxi+Nx-1];
          }
        if (i < Nx) {
            //Lower
            u[t2+i] = (hmu1dt*laplacianu(ts+i+1, ts+i-1, ts+Nx+i, ts+i)+f1(ts+i))+u[ts+i];
            v[t2+i] = (hmu2dt*laplacianv(ts+i+1, ts+i-1, ts+Nx+i, ts+i)+f2(ts+i))+v[ts+i];
            //Upper
            
            u[t2+NxNy-Nx+i] = (hmu1dt*laplacianu(ts+NxNy-Nx+i+1, ts+NxNy-Nx+i-1, ts+NxNy-Nx2+i, ts+NxNy-Nx+i)+f1(ts+NxNy-Nx+i))+u[ts+NxNy-Nx+i];
            v[t2+NxNy-Nx+i] = (hmu2dt*laplacianv(ts+NxNy-Nx+i+1, ts+NxNy-Nx+i-1, ts+NxNy-Nx2+i, ts+NxNy-Nx+i)+f2(ts+NxNy-Nx+i))+v[ts+NxNy-Nx+i];
         }
}

//corners
u[t2] = (hmu1dt*laplacianuc(ts+1, ts+Nx, ts)+f1(ts))+u[ts];
u[t2+NxNy-Nx] = (hmu1dt*laplacianuc(ts+NxNy-Nx, ts+NxNy-Nx2, ts+NxNy-Nx)+f1(ts+NxNy-Nx))+u[ts+NxNy-Nx];
u[t2+Nx-1] = (hmu1dt*laplacianuc(ts+Nx-2, ts+Nx-1+Nx, ts+Nx-1)+f1(ts+Nx-1))+u[ts+Nx-1];
u[t2+NxNy-Nx+Nx-1] = (hmu1dt*laplacianuc(ts+NxNy-Nx+Nx-2, ts+NxNy-Nx2+Nx-1, ts+NxNy-Nx+Nx-1)+f1(ts+NxNy-Nx+Nx-1))+u[ts+NxNy-Nx+Nx-1]; 

v[t2] = (hmu2dt*laplacianvc(ts+1, ts+Nx, ts)+f2(ts))+v[ts];
v[t2+NxNy-Nx] = (hmu2dt*laplacianvc(ts+NxNy-Nx, ts+NxNy-Nx2, ts+NxNy-Nx)+f2(ts+NxNy-Nx))+v[ts+NxNy-Nx];
v[t2+Nx-1] = (hmu2dt*laplacianvc(ts+Nx-2, ts+Nx-1+Nx, ts+Nx-1)+f2(ts+Nx-1))+v[ts+Nx-1];
v[t2+NxNy-Nx+Nx-1] = (hmu2dt*laplacianvc(ts+NxNy-Nx+Nx-2, ts+NxNy-Nx2+Nx-1, ts+NxNy-Nx+Nx-1)+f2(ts+NxNy-Nx+Nx-1))+v[ts+NxNy-Nx+Nx-1]; 

//for this timestep the conditions at boundaries are set, just need to evaluate the inner grid

#pragma omp parallel for 
for (j = 1; j < Ny-1; j++) {
    int Nxj = Nx*j;
    for (i = 1; i < Nx-1; i++) {
        u[t2+Nxj+i] = hmu1dt*((u[ts+i+1+Nxj]-4*u[ts+i+Nxj]+u[ts+i-1+Nxj]+(u[ts+i+Nxj+Nx]+u[ts+i+Nxj-Nx])))+f1(ts+i+Nxj)+u[ts+Nxj+i];
        v[t2+Nxj+i] = hmu2dt*((v[ts+i+1+Nxj]-4*v[ts+i+Nxj]+v[ts+i-1+Nxj]+(v[ts+i+Nxj+Nx]+v[ts+i+Nxj-Nx])))+f2(ts+i+Nxj)+v[ts+Nxj+i];
    }
}

if (k == tn-1) {
    double * u2 = new double[Nx*Ny];
    double * v2 = new double[Nx*Ny];
   for (j = 0; j < Ny; j++) {
    for (i = 0; i < Nx; i++) {
        u2[Nx*(j)+i] = u[t2+Nx*(j)+i];
        v2[Nx*(j)+i] = v[t2+Nx*(j)+i];
            
    }
    
}
    writeFile("ur.txt", Ny, Nx, 1, u2);
    writeFile("vr.txt", Ny, Nx, 1, v2);
    delete[] u2;
    delete[] v2;
}
 

} 
//Write solution to file
//writeFile("x.txt", Ny, Nx, 1, x);
//writeFile("y.txt", Ny, Nx, 1, y);
//writeFile("u.txt", Ny, Nx, tn, u);
//writeFile("v.txt", Ny, Nx, tn, v);

    
    
    delete[] x;
    delete[] y;
    delete[] u;
    delete[] v;

}

 
void printMatrix(double* A, const int& M, const int& N) {
    
    for (int i = 0; i < M; i++)  {
        for (int j = 0; j < N; j++)  {
            cout << setw(2) << A[j*M+i];
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
    //file << "-----------------------------------------------------------------------------------------------" << endl;
}
    file.close();
    
    
    //MPI cartesian topology
}







