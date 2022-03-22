#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
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
                             
            void F77NAME(dsymv) (const char& uplo,
                             const int& N,
                             const double& alpha, double * A, const int& lda, 
                             const double * x, const int& incx,
                             const double& beta, double * Y, const int& incy);
            void F77NAME(dspmv) (const char& uplo,
                             const int& N,
                             const double& alpha, double * A, 
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
void sym2PackedConverter(int n, double* Hp, double* H);


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
    int Nx = 51;
    int Ny = 51;
    double Lx = dx*(Nx-1);
    double Ly = dy*(Ny-1);
    
    double dt = 0.001;
    double t =  1;
    int tn = t/dt; //number of time nodes

    //
    //int Nx = 10;
    //int Ny = 10;
    int nm   =(Nx-2)*(Ny-2);
    int    rank   = 0;
    int    size   = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double * yur = new double[nm];
    double * yvr = new double[nm];
    double * up = new double [3*(Nx-2)];
    double * vp = new double [3*(Nx-2)];
    double * yup = new double [(Nx-2)];
    double * yvp = new double [(Nx-2)];
    double * Aup = new double [2*(Nx-2)*(Nx-2)];
    double * Aup2 = new double [2*(Nx-2)*(Nx-2)];
    double * Aupr = new double [3*(Nx-2)*(Nx-2)];
    double * Avp = new double [2*(Nx-2)*(Nx-2)];
    double * Avp2 = new double [2*(Nx-2)*(Nx-2)];
    double * Avpr = new double [3*(Nx-2)*(Nx-2)];
    double * yu = new double[nm];
    double * yv = new double[nm];
    double * u_ = new double[nm];
    double * v_ = new double[nm];
    double * u = new double[Nx*Ny*(tn)];
    double * v = new double[Nx*Ny*(tn)];
    double * x = new double [Nx*Ny];
    double * y = new double [Nx*Ny];
    if (rank == 0) {
        for (int i = 0; i < nm; i++) {
            if (i < nm) {
        yur[i] = 0, yvr[i] = 0; 
        yu[i] = 0; yv[i] = 0;
            u_[i] = 0; v_[i] = 0;}
        if (i < 3*(Nx-2)) {
        up[i] = 0, vp[i] = 0;} 
        if (i < Nx-2) {
        yup[i] = 0, yvp[i] = 0;}
        if (i < 2*(Nx-2)*(Nx-2)) {
        Aup[i] = 0, Avp[i] = 0;
        Aup2[i] = 0, Avp2[i] = 0;}
        if (i < 3*(Nx-2)*(Nx-2)) {
        Aupr[i] = 0, Avpr[i] = 0;}
        //cout << yur[i] << endl;
    }
    //
    //double * Au = new double [nm*nm];
    //double * Av = new double [nm*nm];
    //
    //int n = (Nx-2)*(Ny-2);
    //double * Aup = new double[n*(n+1)/2];
    //double * Avp = new double[n*(n+1)/2];
        double * Au = new double [nm*nm];
        double * Av = new double [nm*nm];
        for (int i = 0; i < nm*nm; i++) {
            Au[i] = 0;
            Av[i] = 0;
        }
    Banded(Nx-2, Ny-2, (Nx-2)*(Ny-2), Au, mu1/h/h*(-4*dt)+1, mu1/h/h*dt);
    Banded(Nx-2, Ny-2, (Nx-2)*(Ny-2), Av, mu2/h/h*(-4*dt)+1, mu2/h/h*dt);
    //sym2PackedConverter(n, Aup, Au);
    //sym2PackedConverter(n, Avp, Av);
    //printMatrix(Au, (Nx-2)*(Ny-2), (Nx-2)*(Ny-2));
    //double * yu = new double[(Nx-2)*(Ny-2)];
    //double * yv = new double[(Nx-2)*(Ny-2)];
    //double * u_ = new double[(Nx-2)*(Ny-2)];
    //double * v_ = new double[(Nx-2)*(Ny-2)];
    //meshgrid
    double * x = new double [Nx*Ny];
    double * y = new double [Nx*Ny];
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            x[Nx*i+j] = dx*j;
            y[Nx*i+j] = dy*i;
        }
    }
    //boundary conditions from the handout

    

    
    
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            if (y[Nx*i+j] > Ly/2) {
                u[Nx*i+j] = 1;
                //cout << u[Nx*i+j]  << endl;
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
        for (int j = 0; j < 2*(Nx-2); j++) {
        for (int i = 0; i < Nx-2; i++ ) {
            Aup2[j*(Nx-2)+i] = Au[(Ny-4)*(Nx-2)*nm+nm-Nx+2+j*nm+i];
            Avp2[j*(Nx-2)+i] = Av[(Ny-4)*(Nx-2)*nm+nm-Nx+2+j*nm+i];
        }
    }
    //repeated matrix
    for (int j = 0; j < 3*(Nx-2); j++) {
        for (int i = 0; i < Nx-2; i++ ) {
            Aupr[j*(Nx-2)+i] = Au[Nx-2+j*nm+i];
            Avpr[j*(Nx-2)+i] = Av[Nx-2+j*nm+i];
        }
    }
    
    }//if rank=0
    MPI_Bcast(x, nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //arrayPrint(u, Nx*Ny);
    //Neuman boundary conditions
    //x = 0
    //Merge with upper code when done
    //do for i = 1, come back later for i = 0
    //eg issue with 
//cout << "rank" << rank << endl;
for (int k = 1; k < tn; k++) {

   int ts = Nx*Ny*(k-1);
   int t2 = Nx*Ny*(k);
   
if (rank == 0) {
    for (int i = 0; i < Nx*Ny*(tn); i++) {
        //cout << u[i] << endl;
    }
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
            //cout << u[t2+i]<< endl;
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

    //cout << u[ts + Nx*j + i] << endl;//----------------------------------------------------------------------
    }
}
}//end ifrank == 0

int np = Ny-4;
    
MPI_Comm_size(MPI_COMM_WORLD, &size);

MPI_Bcast(Aupr, 3*(Nx-2)*(Nx-2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(Avpr, 3*(Nx-2)*(Nx-2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(u_, nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(v_, nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(yu, nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(yv, nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//F77NAME(dspmv) ('l', n, 1, Aup, u_, 1, 0, yu, 1);
//F77NAME(dspmv) ('l', n, 1, Avp, v_, 1, 0, yv, 1);

//
    int    r     = np % size;            // remainder
    int    kp     = (np - r) / size;      // minimum size of number of blocks(chunks)
    //cout << kp << endl;//-------------------------------------------------------------------------------------------------
    int    ustart = 0;                   // start index of number of blocks(chunks) n
    int    uend   = 0;                   // end index of number of blocks(chunks)
    int    nob = 0;
    if (rank < (np % size)) {            // for ranks < r, number of block(chunk) is size k + 1
        kp++;
        nob = kp;
        ustart = kp * rank;
        uend   = kp * (rank + 1)-1;
    }
    else {                              // for ranks > r, chunk size is k
        ustart = (kp+1) * r + kp * (rank - r);
        uend   = (kp+1) * r + kp * (rank - r + 1)-1;
        nob = kp;
    }
    
    //allocate u to up in a for loop from 0 to k,
    //each loop returns a yp a up gets multiplied which gets added onto a bigger y vector that contains yp for all blocks in cpu
    //
    for (int j = 0; j < kp; j++) {
         for (int i = 0; i < 3*(Nx-2); i++) {
         up[i] = u_[ustart*(Nx-2)+j*(Nx-2)+i];
         vp[i] = v_[ustart*(Nx-2)+j*(Nx-2)+i];
         //cout << up[i] << endl;------------------------------------------------------------------------------
    }
    F77NAME(dgemv) ('N', Nx-2, 3*(Nx-2), 1, Aupr, Nx-2, up, 1, 0, yup, 1);
    F77NAME(dgemv) ('N', Nx-2, 3*(Nx-2), 1, Avpr, Nx-2, vp, 1, 0, yvp, 1);
         
         for (int i = 0; i < (Nx-2); i++) {
         yu[Nx-2+ustart*(Nx-2)+(Nx-2)*j+i] = yup[i];
         yv[Nx-2+ustart*(Nx-2)+(Nx-2)*j+i] = yvp[i];
    }
    //cout << "------------" << endl;
    }
    
    for (int i = 0; i < nm; i++) {
                 if (rank==0) {
             //cout << yu[i] << endl;//----------------------------------------------------------------------------------
         }
    }
    //Calculate matrix local dgemv
    /*
    cout << "rank " << rank << " has " << nob << " blocks to multiply" <<endl;
    cout << "ustart: " << ustart*(Nx-2) <<endl;
    cout << "uend: " << uend*(Nx-2) <<endl;
     * */
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(yu, yur, nm, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        double * upc = new double[2*(Nx-2)];
        double * vpc = new double[2*(Nx-2)];
        for (int i = 0; i < 2*(Nx-2); i++) {
            upc[i] = 0, vpc[i] = 0;
        }
        for (int i = 0; i < 2*(Nx-2); i++) {
             upc[i] = u_[i];
             vpc[i] = v_[i];
             //cout << upc[i] << endl;
        }
        for (int i = 0; i < 2*(Nx-2)*(Nx-2); i++) {
            //cout << Aup[i] << endl;
        }
        F77NAME(dgemv) ('N', Nx-2, 2*(Nx-2), 1, Aup, Nx-2, upc, 1, 0, yup, 1);
        F77NAME(dgemv) ('N', Nx-2, 2*(Nx-2), 1, Avp, Nx-2, upc, 1, 0, yvp, 1);
        for (int i = 0; i < (Nx-2); i++) {
             yur[i] = yup[i];
             yvr[i] = yvp[i];
             //cout << yup[i] << endl;
        }
        //bottom right corner
        for (int i = 0; i < 2*(Nx-2); i++) {
             upc[i] = u_[nm-2*(Nx-2)+i];
             vpc[i] = v_[nm-2*(Nx-2)+i];
             //cout << u_[nm-2*(Nx-2)+i] << endl;
        }
        F77NAME(dgemv) ('N', Nx-2, 2*(Nx-2), 1, Aup2, Nx-2, upc, 1, 0, yup, 1);
        F77NAME(dgemv) ('N', Nx-2, 2*(Nx-2), 1, Avp2, Nx-2, vpc, 1, 0, yvp, 1);
        for (int i = 0; i < (Nx-2); i++) {
             yur[nm-(Nx-2)+i] = yup[i];
             yvr[nm-(Nx-2)+i] = yvp[i];
             //cout << yup[i] << endl; //----------------------------------------------------------------------------
             //cout << yu[nm-(Nx-2)+i] << endl;
        }
    for (int i = 0; i < nm; i++) {
    //cout << yur[i] << endl;
}    
    
for (int j = 1; j < (Ny-1); j++) {
    for (int i = 1; i < (Nx-1); i++) {
        if (i ==1) {
            u[t2+Nx*j+i] = yur[(j-1)*(Nx-2)+i-1]+f1(ts+Nx*j+i)*dt+mu1*dt/h/h*u[ts+Nx*j+i-1];
            v[t2+Nx*j+i]=  yvr[(j-1)*(Nx-2)+i-1]+f2(ts+Nx*j+i)*dt+mu2*dt/h/h*v[ts+Nx*j+i-1];
        }
        else if (i==Nx-2) {
            u[t2+Nx*j+i] = yur[(j-1)*(Nx-2)+i-1]+f1(ts+Nx*j+i)*dt+mu1*dt/h/h*u[ts+Nx*j+i+1];
            v[t2+Nx*j+i]=  yvr[(j-1)*(Nx-2)+i-1]+f2(ts+Nx*j+i)*dt+mu2*dt/h/h*v[ts+Nx*j+i+1];
            
        }
        else {
            u[t2+Nx*j+i] = yur[(j-1)*(Nx-2)+i-1]+f1(ts+Nx*j+i)*dt;
            v[t2+Nx*j+i]=  yvr[(j-1)*(Nx-2)+i-1]+f2(ts+Nx*j+i)*dt;
        }
    }
}
//

if (k == tn-1) {
    double * u2 = new double[Nx*Ny];
    double * v2 = new double[Nx*Ny];
   for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
        u2[Nx*(j)+i] = u[t2+Nx*(j)+i];
        v2[Nx*(j)+i] = v[t2+Nx*(j)+i];
        //cout << u2[Nx*(j)+i] << endl;//----------------------------------------------------------------------------
    }
}
 
    writeFile("xr.txt", Ny, Nx, 1, x);
    writeFile("yr.txt", Ny, Nx, 1, y);
    writeFile("ur.txt", Ny, Nx, 1, u2);
    writeFile("vr.txt", Ny, Nx, 1, v2);
     
 
}
    }
} 
//Write solution at last timestep to file
    delete[] x;
    delete[] y;
    delete[] u;
    delete[] v;
    MPI_Finalize();
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

void sym2PackedConverter(int n, double* Hp, double* H) {
    //size of Hp must be n*(n+1)/2, where n is dimension of symmetric matrix
    int k = 0;
    for (int j = 0; j < n; j++) {
            for (int i = 0; i < n-j; i++){
               Hp[k] = H[j*n+j+i]; 
               k = k + 1;
            }
    }
}