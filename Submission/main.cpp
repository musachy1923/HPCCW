#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <boost/program_options.hpp>
//useful macros for the equations, just makes readability easier
#define laplacianu(a, b, c, d) ((u[a]+u[b]+u[c]-3*u[d]))
#define laplacianv(a, b, c, d) ((v[a]+v[b]+v[c]-3*v[d]))
#define laplacianuc(a, b, c) ((u[a]+u[b]-2*u[c]))
#define laplacianvc(a, b, c) ((v[a]+v[b]-2*v[c]))
#define f1(a) (epsilon*u[a]*(1-u[a])*(u[a]-(v[a]+b)*diva))*dt
#define f2(a) (u[a]*u[a]*u[a]-v[a])*dt
using namespace std;
namespace po = boost::program_options;

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
    //for (int k = 0; k < tn; k++) {
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx;  j++) {
                file << x[Nx*i+j]<<" ";
    }
    file << "\n";
    }
    file.close();

}

//ReactionDiffusion class, encapsulates both variables u and v. Hold other constants
class ReactionDiffusion {
    
public:
ReactionDiffusion() = default;
~ReactionDiffusion() {};
    //public methods
    void SetParameters(int ac, char* av[]);
    void SetInitialConditions();
    void TimeIntegrate();
private:
    int Nx, Ny, ts, t2, tn;
    double dt, t, a, b, mu1, mu2, epsilon, dx=1, dy=1, h, hmu1dt, hmu2dt;
    double * u = new double[2*Nx*Ny], * v = new double[2*Nx*Ny], * x = new double[Nx*Ny], * y= new double[Nx*Ny];
    double * u2 = new double[Nx*Ny], * v2 = new double[Nx*Ny];
    
    
    
};

//SetParameters is made to accept command line arguements
void ReactionDiffusion::SetParameters(int ac, char* av[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("dt", po::value<double>() ->default_value(0.001), "Time-step to use")
        ("T", po::value<double>()->default_value(100), "Total integration time")
        ("Nx", po::value<int>()->default_value(101), "Number of grid points in x")
        ("Ny", po::value<int>()->default_value(101), "Number of grid points in y")
        ("a", po::value<double>() ->default_value(0.75), "Value of parameter a")
        ("b", po::value<double>()->default_value(0.06), "Value of parameter b")
        ("mu1", po::value<double>()->default_value(5), "Value of parameter mu1")
        ("mu2", po::value<double>()->default_value(0), "Value of parameter mu2")
        ("eps", po::value<double>()->default_value(50), "Value of parameter epsilon");
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).
                options(desc).run(), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        cout << desc << endl;
    
    }
    a = vm["a"].as<double>();
    b = vm["b"].as<double>();
    epsilon = vm["eps"].as<double>();
    mu1 = vm["mu1"].as<double>();
    mu2 = vm["mu2"].as<double>();
    Nx = vm["Nx"].as<int>();
    Ny = vm["Ny"].as<int>();
    dt = vm["dt"].as<double>();
    t =  vm["T"].as<double>();
    tn = t/dt;
    
}

void ReactionDiffusion::SetInitialConditions() {
    //constants
    h = dx;
    double Lx = dx*(Nx-1);
    double Ly = dy*(Ny-1);
    
    u = new double[Nx*Ny*(2)];
    v = new double[Nx*Ny*(2)];
    //meshgrid
    x = new double [Nx*Ny];
    y = new double [Nx*Ny];
    
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            x[Nx*i+j] = dx*j;
            y[Nx*i+j] = dy*i;
        }
    }
    //intial boundary conditions from handout are applied
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
    
    
    
}

void ReactionDiffusion::TimeIntegrate() {
    //to optimise calculations, most multiplications are divisions are initiated here
    hmu1dt=mu1/h/h*dt;
    hmu2dt=mu2/h/h*dt;
    int NxNy = Nx*Ny;
    int Nx2 = 2*Nx;
    double diva = 1/a;
    
    //This time loops imlements the 2D-Diffusuion equation that was given in the handout
    for (int k = 1; k < tn; k++) {
       int ts = 0;
       int t2 = Nx*Ny;
       int j = 0;
       int i = 0;
       //The boundary conditions for the corners using neumann boundary conditions
        u[t2] = (hmu1dt*laplacianuc(ts+1, ts+Nx, ts)+f1(ts))+u[ts];
        u[t2+NxNy-Nx] = (hmu1dt*laplacianuc(ts+NxNy-Nx, ts+NxNy-Nx2, ts+NxNy-Nx)+f1(ts+NxNy-Nx))+u[ts+NxNy-Nx];
        u[t2+Nx-1] = (hmu1dt*laplacianuc(ts+Nx-2, ts+Nx-1+Nx, ts+Nx-1)+f1(ts+Nx-1))+u[ts+Nx-1];
        u[t2+NxNy-Nx+Nx-1] = (hmu1dt*laplacianuc(ts+NxNy-Nx+Nx-2, ts+NxNy-Nx2+Nx-1, ts+NxNy-Nx+Nx-1)+f1(ts+NxNy-Nx+Nx-1))+u[ts+NxNy-Nx+Nx-1]; 
        v[t2] = (hmu2dt*laplacianvc(ts+1, ts+Nx, ts)+f2(ts))+v[ts];
        v[t2+NxNy-Nx] = (hmu2dt*laplacianvc(ts+NxNy-Nx, ts+NxNy-Nx2, ts+NxNy-Nx)+f2(ts+NxNy-Nx))+v[ts+NxNy-Nx];
        v[t2+Nx-1] = (hmu2dt*laplacianvc(ts+Nx-2, ts+Nx-1+Nx, ts+Nx-1)+f2(ts+Nx-1))+v[ts+Nx-1];
        v[t2+NxNy-Nx+Nx-1] = (hmu2dt*laplacianvc(ts+NxNy-Nx+Nx-2, ts+NxNy-Nx2+Nx-1, ts+NxNy-Nx+Nx-1)+f2(ts+NxNy-Nx+Nx-1))+v[ts+NxNy-Nx+Nx-1]; 

    //conditions at boundaries are set, the following evaluates the inner grid. 
    //The neumann boundary conditions for the next time step calculates with the rest of grid and indexing are distinguished using if statements.
#pragma omp parallel for
    for (j = 1; j < Ny; j++) {
        
        int Nxj = Nx*j;
        
        if (j < Ny){
                //Left side of grid neumann boundary conditions
                u[t2+Nxj] = (hmu1dt*laplacianu(ts+Nxj+1, ts+Nxj+Nx, ts+Nxj-Nx, ts+Nxj)+f1(ts+Nxj))+u[ts+Nxj];
                v[t2+Nxj] = (hmu2dt*laplacianv(ts+Nxj+1, ts+Nxj+Nx, ts+Nxj-Nx, ts+Nxj)+f2(ts+Nxj))+v[ts+Nxj];
                //Right side of grid boundary conditions
                u[t2+Nxj+Nx-1] = (hmu1dt*laplacianu(ts+Nxj+Nx-2, ts+Nxj+Nx+Nx-1, ts+Nxj-Nx+Nx-1, ts+Nxj+Nx-1)+f1(ts+Nxj+Nx-1))+u[ts+Nxj+Nx-1];
                v[t2+Nxj+Nx-1] = (hmu2dt*laplacianv(ts+Nxj+Nx-2, ts+Nxj+Nx+Nx-1, ts+Nxj-Nx+Nx-1, ts+Nxj+Nx-1)+f2(ts+Nxj+Nx-1))+v[ts+Nxj+Nx-1];
              }
              
        if (j < Nx) {
                //Lower side of grid neumann boundary conditions
                u[t2+j] = (hmu1dt*laplacianu(ts+j+1, ts+j-1, ts+Nx+j, ts+j)+f1(ts+j))+u[ts+j];
                v[t2+j] = (hmu2dt*laplacianv(ts+j+1, ts+j-1, ts+Nx+j, ts+j)+f2(ts+j))+v[ts+j];
                //Upper side of grid neumann boundary conditions
                u[t2+NxNy-Nx+j] = (hmu1dt*laplacianu(ts+NxNy-Nx+j+1, ts+NxNy-Nx+j-1, ts+NxNy-Nx2+j, ts+NxNy-Nx+j)+f1(ts+NxNy-Nx+j))+u[ts+NxNy-Nx+j];
                v[t2+NxNy-Nx+j] = (hmu2dt*laplacianv(ts+NxNy-Nx+j+1, ts+NxNy-Nx+j-1, ts+NxNy-Nx2+j, ts+NxNy-Nx+j)+f2(ts+NxNy-Nx+j))+v[ts+NxNy-Nx+j];
             }
             
        
        if (j<Ny-1 && j>0) {
        for (i = 1; i < Nx-1; i++) {
            //main for loop that evaluates the u and v values in the next time step.
            u[t2+Nxj+i] = hmu1dt*((u[ts+i+1+Nxj]-4*u[ts+i+Nxj]+u[ts+i-1+Nxj]+(u[ts+i+Nxj+Nx]+u[ts+i+Nxj-Nx])))+f1(ts+i+Nxj)+u[ts+Nxj+i];
            v[t2+Nxj+i] = hmu2dt*((v[ts+i+1+Nxj]-4*v[ts+i+Nxj]+v[ts+i-1+Nxj]+(v[ts+i+Nxj+Nx]+v[ts+i+Nxj-Nx])))+f2(ts+i+Nxj)+v[ts+Nxj+i];
        }
        }
    }
    

    for (i = 0; i < Nx*Ny; i++) {

        u[i] = u[t2+i];
        v[i] = v[t2+i];
        
    }
    //write the solution at the last timestep
    if (k == tn-1) {
        u2 = new double[Nx*Ny];
        v2 = new double[Nx*Ny];
        
       for (j = 0; j < NxNy; j++) {
            u2[j] = u[t2+j];
            v2[j] = v[t2+j];
    }
        writeFile("ur.txt", Ny, Nx, 1, u2);
        writeFile("vr.txt", Ny, Nx, 1, v2);
        delete[] u2;
        delete[] v2;
    }
     

    }
        delete[] x;
        delete[] y;
        delete[] u;
        delete[] v;
        
    }

int main(int ac, char* av[]) {
    
    //calls class methods to solve the 2D Diffusion equation 
    ReactionDiffusion r;
    
    r.SetParameters(ac, av);
    
    r.SetInitialConditions();
    
    r.TimeIntegrate();

    
}