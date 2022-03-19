#include <iostream>
#include <iomanip>
#include <fstream>
//useful macros
#define laplacianu(a, b, c, d) (1/h/h*(u[a]+u[b]+u[c]-3*u[d]))
#define laplacianv(a, b, c, d) (1/h/h*(v[a]+v[b]+v[c]-3*v[d]))
#define laplacianuc(a, b, c) (1/h/h*(u[a]+u[b]-2*u[c]))
#define laplacianvc(a, b, c) (1/h/h*(v[a]+v[b]-2*v[c]))
#define f1(x) (epsilon*u[x]*(1-u[x])*(u[x]-(v[x]+b)/a))
#define f2(x) (u[x]*u[x]*u[x]-v[x])


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
void singleBanded(int nsv, double* H, double alpha, double beta);
void arrayPrint(double * a, int n);
void Banded(int n, int m, int NS, double* H, double alpha, double beta);


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
    int Nx = 10;
    int Ny = 10;
    double Lx = dx*Nx;
    double Ly = dy*Ny;
    
    double dt = 0.001;
    double t =  0.002;
    int tn = t/dt; //number of time nodes
    
    //
    double * Au = new double[(Nx-2)*(Ny-2)*(Nx-2)*(Ny-2)];
    double * Av = new double[(Nx-2)*(Ny-2)*(Nx-2)*(Ny-2)];
    Banded(Nx-2, Ny-2, (Nx-2)*(Ny-2), Au, mu1/h/h*(-4*dt)+1, mu1/h/h*dt);
    Banded(Nx-2, Ny-2, (Nx-2)*(Ny-2), Av, mu2/h/h*(-4*dt)+1, mu2/h/h*dt);
    printMatrix(Au, (Nx-2)*(Ny-2), (Nx-2)*(Ny-2));
    //
    double * yu = new double[(Nx-2)*(Ny-2)];
    double * yv = new double[(Nx-2)*(Ny-2)];
    double * u_ = new double[(Nx-2)*(Ny-2)];
    double * v_ = new double[(Nx-2)*(Ny-2)];
    
    
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

for (int k = 1; k < tn; k++) {

   int ts = Nx*Ny*(k-1);
   int t2 = Nx*Ny*(k);
   

        for (int i = 1; i < Ny+abs(Nx-Ny); i++) {
        if (i < Ny){
            //Left
            u[t2+i*Nx] = (mu1*laplacianu(ts+Nx*i+1, ts+Nx*(i+1), ts+Nx*(i-1), ts+Nx*i)+f1(ts+Nx*i))*dt+u[ts+Nx*i];
            v[t2+i*Nx] = (mu2*laplacianv(ts+Nx*i+1, ts+Nx*(i+1), ts+Nx*(i-1), ts+Nx*i)+f2(ts+Nx*i))*dt+v[ts+Nx*i];
            //Right
            u[t2+i*Nx+Nx-1] = (mu1*laplacianu(ts+i*Nx+Nx-2, ts+(i+1)*Nx+Nx-1, ts+(i-1)*Nx+Nx-1, ts+i*Nx+Nx-1)+f1(ts+i*Nx+Nx-1))*dt+u[ts+i*Nx+Nx-1];
            v[t2+i*Nx+Nx-1] = (mu2*laplacianv(ts+i*Nx+Nx-2, ts+(i+1)*Nx+Nx-1, ts+(i-1)*Nx+Nx-1, ts+i*Nx+Nx-1)+f2(ts+i*Nx+Nx-1))*dt+v[ts+i*Nx+Nx-1];
          }
        if (i < Nx) {
            //Lower
            u[t2+i] = (mu1*laplacianu(ts+i+1, ts+i-1, ts+Nx+i, ts+i)+f1(ts+i))*dt+u[ts+i];
            v[t2+i] = (mu2*laplacianv(ts+i+1, ts+i-1, ts+Nx+i, ts+i)+f2(ts+i))*dt+v[ts+i];
            //Upper
            
            u[t2+Nx*(Ny-1)+i] = (mu1*laplacianu(ts+Nx*(Ny-1)+i+1, ts+Nx*(Ny-1)+i-1, ts+Nx*(Ny-2)+i, ts+Nx*(Ny-1)+i)+f1(ts+Nx*(Ny-1)+i))*dt+u[ts+Nx*(Ny-1)+i];
            v[t2+Nx*(Ny-1)+i] = (mu2*laplacianv(ts+Nx*(Ny-1)+i+1, ts+Nx*(Ny-1)+i-1, ts+Nx*(Ny-2)+i, ts+Nx*(Ny-1)+i)+f2(ts+Nx*(Ny-1)+i))*dt+v[ts+Nx*(Ny-1)+i];
         }
}


//corners
u[t2] = (mu1*laplacianuc(ts+1, ts+Nx, ts)+f1(ts))*dt+u[ts];
u[t2+Nx*(Ny-1)] = (mu1*laplacianuc(ts+Nx*(Ny-1), ts+Nx*(Ny-2), ts+Nx*(Ny-1))+f1(ts+Nx*(Ny-1)))*dt+u[ts+Nx*(Ny-1)];
u[t2+Nx-1] = (mu1*laplacianuc(ts+Nx-2, ts+Nx-1+Nx, ts+Nx-1)+f1(ts+Nx-1))*dt+u[ts+Nx-1];
u[t2+Nx*(Ny-1)+Nx-1] = (mu1*laplacianuc(ts+Nx*(Ny-1)+Nx-2, ts+Nx*(Ny-2)+Nx-1, ts+Nx*(Ny-1)+Nx-1)+f1(ts+Nx*(Ny-1)+Nx-1))*dt+u[ts+Nx*(Ny-1)+Nx-1]; 

v[t2] = (mu2*laplacianvc(ts+1, ts+Nx, ts)+f2(ts))*dt+v[ts];
v[t2+Nx*(Ny-1)] = (mu2*laplacianvc(ts+Nx*(Ny-1), ts+Nx*(Ny-2), ts+Nx*(Ny-1))+f2(ts+Nx*(Ny-1)))*dt+v[ts+Nx*(Ny-1)];
v[t2+Nx-1] = (mu2*laplacianvc(ts+Nx-2, ts+Nx-1+Nx, ts+Nx-1)+f2(ts+Nx-1))*dt+v[ts+Nx-1];
v[t2+Nx*(Ny-1)+Nx-1] = (mu2*laplacianvc(ts+Nx*(Ny-1)+Nx-2, ts+Nx*(Ny-2)+Nx-1, ts+Nx*(Ny-1)+Nx-1)+f2(ts+Nx*(Ny-1)+Nx-1))*dt+v[ts+Nx*(Ny-1)+Nx-1]; 

//for this timestep the conditions at boundaries are set, just need to evaluate the inner grid







for (int j = 1; j < Nx-1; j++) {
    for (int i = 1; i < (Nx-1); i++) {
      u_[(j-1)*(Nx-2)+i-1] = u[ts + Nx*j + i];
      v_[(j-1)*(Nx-2)+i-1] = v[ts + Nx*j + i];   
    }
}
//printMatrix(u, Ny, Nx);
//printMatrix(x, Ny, Nx);
//printMatrix(y, Ny, Nx);
//printMatrix(u_, (Ny-2), (Nx-2));
//arrayPrint(u_, (Nx-2)*(Ny-2));



F77NAME(dgemv) ('N', (Nx-2)*(Ny-2), (Nx-2)*(Ny-2), 1, Au, (Nx-2)*(Ny-2), u_, 1, 0, yu, 1);
F77NAME(dgemv) ('N', (Nx-2)*(Ny-2), (Nx-2)*(Ny-2), 1, Av, (Nx-2)*(Ny-2), v_, 1, 0, yv, 1);

for (int j = 1; j < (Ny-1); j++) {
    for (int i = 1; i < (Nx-1); i++) {
        if (i ==1) {
            u[t2+Nx*j+i] = yu[(j-1)*(Nx-2)+i-1]+f1(ts+Nx*j+i)*dt+mu1*dt/h/h*u[ts+Nx*j+i-1];
            v[t2+Nx*j+i]=  yv[(j-1)*(Nx-2)+i-1]+f2(ts+Nx*j+i)*dt+mu2*dt/h/h*v[ts+Nx*j+i-1];
        }
        else if (i==Nx-2) {
            u[t2+Nx*j+i] = yu[(j-1)*(Nx-2)+i-1]+f1(ts+Nx*j+i)*dt+mu1*dt/h/h*u[ts+Nx*j+i+1];
            v[t2+Nx*j+i]=  yv[(j-1)*(Nx-2)+i-1]+f2(ts+Nx*j+i)*dt+mu2*dt/h/h*v[ts+Nx*j+i+1];
            
        }
        else {
            u[t2+Nx*j+i] = yu[(j-1)*(Nx-2)+i-1]+f1(ts+Nx*j+i)*dt;
            v[t2+Nx*j+i]=  yv[(j-1)*(Nx-2)+i-1]+f2(ts+Nx*j+i)*dt;
        }
        
        
        
    // << yu(ts+Nx*j+i)*dt << endl;
    //cout << f2(ts+Nx*j+i)*dt << endl;
}
    
    
}


//arrayPrint(yu, (Nx-2)*(Ny-2));
//printMatrix(yu, Ny-2, Nx-2);
//printMatrix(Au, (Nx-2)*(Ny-2), (Nx-2)*(Ny-2));

//printMatrix(A, (Ny-2)*(Nx-2), (Ny-2)*(Nx-2));

//arrayPrint(yu, (Nx-2)*(Ny-2));
//


/*
double * A = new double[(Nx-2)*(Ny-2)];
double * yu = new double[Nx-2];
double * yv = new double[Nx-2];
double * u_j = new double[(Nx-2)];
double * u_jp1 = new double[(Nx-2)];
double * v_j = new double[(Nx-2)];
double * v_jp1 = new double[(Nx-2)];

singleBanded(Ny-2, A, -4, 1);

for (int j = 1; j < Ny-1; j++) {
    for (int i = 2; i < Nx-2; i++) {
        u_j[i-1] = mu1/h/h*u[ts+i+Nx*j];
        u_jp1[i-1] = mu1/h/h*(u[ts+i+Nx*(j+1)]+u[i+Nx*(j-1)])*dt+f1(ts+i+Nx*j);
        v_j[i-1] = mu2/h/h*v[ts+i+Nx*j];
        v_jp1[i-1] = mu2/h/h*(v[ts+i+Nx*(j+1)]+v[ts+i+Nx*(j-1)])*dt+f2(ts+i+Nx*j);
    }
    F77NAME(dgemv) ('N', Ny-2, Nx-2, 1, A, Nx-2, u_j, 1, 0, yu, 1);
    F77NAME(dgemv) ('N', Ny-2, Nx-2, 1, A, Nx-2, v_j, 1, 0, yv, 1);
    for (int i = 2; i < Nx-2; i++) {
        u[t2+Nx*(j)+i] = yu[i-1]*dt+u_jp1[i-1];
        v[t2+Nx*(j)+i] = yv[i-1]*dt+v_jp1[i-1];
    }
    u[t2+Nx*(j)+1] = yu[0]*dt+u_jp1[0]+mu1/h/h*u[ts+Nx*(j)];
    v[t2+Nx*(j)+Nx-2] = yv[Nx-2-1]*dt+v_jp1[Nx-2-1]+mu2/h/h*v[ts+Nx*(j)];
}
*/
if (k == tn-1) {
    double * u2 = new double[Nx*Ny];
    double * v2 = new double[Nx*Ny];
   for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
        u2[Nx*(j)+i] = u[t2+Nx*(j)+i];
        v2[Nx*(j)+i] = v[t2+Nx*(j)+i];
    }
}
 
    writeFile("xr.txt", Ny, Nx, 1, x);
    writeFile("yr.txt", Ny, Nx, 1, y);
    writeFile("ur.txt", Ny, Nx, 1, u2);
    writeFile("vr.txt", Ny, Nx, 1, v2);
     
 
}
 

} 







//Write solution at last timestep to file

    //printMatrix(u, Ny, Nx);
    //printMatrix(v, Ny, Nx);
writeFile("x.txt", Ny, Nx, 1, x);
writeFile("y.txt", Ny, Nx, 1, y);
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

void singleBanded(int nsv, double* H, double alpha, double beta) {
    H[0] = alpha;
    H[1] = beta;
    H[(nsv-1)*nsv + nsv-1 - 1] = beta;
    H[(nsv-1)*nsv + nsv-1] = alpha;
    for (int i = 1; i < nsv-1; ++i) {
        H[i*nsv + i - 1] = beta;
        H[i*nsv + i] = alpha;
        H[i*nsv + i + 1] = beta;
    }
    
}

void Banded(int n, int m, int NS, double* H, double alpha, double beta) {
    for (int j = 0; j < m; j++){
        //H[j*NS*n+j*n+NS*n] = alpha;
        
        H[j*NS*n+j*n] = alpha;
        H[j*NS*n+j*n+1] = beta;
        H[(n-1)*NS + n-1 - 1+NS*n*j+j*n] = beta;
        H[(n-1)*NS + n-1+NS*n*j+j*n] = alpha;
        for (int i = 1; i < n-1; ++i) {
            H[i*NS + +j*NS*n+n*j+i- 1] = beta;
            H[i*NS + +j*NS*n+n*j+i] = alpha;
            H[i*NS + +j*NS*n+n*j+i+ 1] = beta;
        }
        if (j < m-1) {
            for (int i = 0; i < n; i++) {
            H[j*NS*n+j*n+NS*n+NS*(i)+i] = beta;
            H[j*NS*n+j*n+n+NS*(i)+i] = beta;
        }
        }
    }
}

void singleUpperBanded(int nsv, double* H, double alpha, double beta) {
    H[0] = alpha;
    H[1] = beta;
    for (int i = 1; i < nsv; ++i) {
        H[i*nsv + i - 1] = beta;
        H[i*nsv + i] = alpha;
        H[i*nsv + i + 1] = beta;
    }
    //printMatrix(H, nsv, nsv);
}

void arrayPrint(double * a, int n) {
    for (int i = 0; i < n; i++) {
    
        cout << a[i] << endl;    
    
    }
}