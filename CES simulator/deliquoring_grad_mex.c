#include <math.h>
#include <stdlib.h>
#include <mex.h>


// functions
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))

/* Input Arguments */
#define T prhs[0]
#define X prhs[1]
#define NUMBER_NODES prhs[2]
#define STEP_GRID_DIMLESS prhs[3]
#define PB prhs[4]
#define LAMBDA prhs[5]
#define SCALING prhs[6]
#define PGIN prhs[7]
#define PGOUT prhs[8]
#define PG prhs[9]



double fluxes_HRFVM_vL(double x_minus1, double x,double x_plus1,int j,int number_nodes) {
    double r, r1, phi_r;
    r=(x-x_minus1+1e-8)/(x_plus1-x+1e-8);
    r1=r;
    if (r<0) {
      r1=-r1;
    }

    phi_r=(r+r1)/(1+r1);

    return x+0.5*phi_r*(x_plus1-x);
  }



void deliquoring_pde_model(double* t, double SR[], int number_nodes, double step_grid_dimless,
                  double Pb, double lambda, double scaling, double Pgin, double Pgout,
                  double Pg[], double dSdtheta[])

                  {
                    double ulin, Pl[number_nodes], kl[number_nodes],
                            dPldz[number_nodes], duldz[number_nodes], flux_Pl_in,
                            fluxes_ul[number_nodes], ul[number_nodes], fluxes_Pl[number_nodes],
                            flux_ul_in  ;
                    int j;
                           // First Loop
                           // first node
                           ulin=0;
                           flux_ul_in=ulin;
                           Pl[0]=(Pg[0]-Pb*pow(SR[0],-1/lambda))/Pb/scaling;
                           kl[0]=pow(SR[0],(2+3*lambda)/lambda);
                           Pl[1]=(Pg[1]-Pb*pow(SR[1],-1/lambda))/Pb/scaling;
                           kl[1]=pow(SR[1],(2+3*lambda)/lambda);
                           fluxes_Pl[0]=0.5*(Pl[0]+Pl[1]);
                           flux_Pl_in=Pl[0];
                           dPldz[0]=(fluxes_Pl[0]-flux_Pl_in)*2/step_grid_dimless;
                           ul[0]=-kl[0]*dPldz[0];

                           //second up to number_nodes-1
                           for (j=1; j<number_nodes-1;j++) {
                              Pl[j+1]=(Pg[j+1]-Pb*pow(SR[j+1],-1/lambda))/Pb/scaling;
                              kl[j+1]=pow(SR[j+1],(2+3*lambda)/lambda);
                              fluxes_Pl[j]=fluxes_HRFVM_vL(Pl[j-1], Pl[j], Pl[j+1], j, number_nodes);
                              dPldz[j]=(fluxes_Pl[j]-fluxes_Pl[j-1])/step_grid_dimless;
                              ul[j]=-kl[j]*dPldz[j];

                           }

                           // last node
                           fluxes_Pl[number_nodes-1]=Pl[number_nodes-1];

                           // Second Loop
                           // last node
                           dPldz[number_nodes-1]=(fluxes_Pl[number_nodes-1]-fluxes_Pl[number_nodes-2])*2/step_grid_dimless; //step of dPldz halved for the last volume
                           ul[number_nodes-1]=-kl[number_nodes-1]*dPldz[number_nodes-1];
                           fluxes_ul[number_nodes-1]=ul[number_nodes-1];
                           // first node
                           fluxes_ul[0]=0.5*(ul[0]+ul[1]);
                           duldz[0]=(fluxes_ul[0]-flux_ul_in)/step_grid_dimless;
                           dSdtheta[0] = -duldz[0];
                           

                           for (j=1; j<number_nodes-1;j++) {
                              fluxes_ul[j]=fluxes_HRFVM_vL(ul[j-1], ul[j], ul[j+1], j, number_nodes);
                              duldz[j]=(fluxes_ul[j]-fluxes_ul[j-1])/step_grid_dimless;
                              dSdtheta[j] = -duldz[j];
                           }
                           duldz[number_nodes-1]=(fluxes_ul[number_nodes-1]-fluxes_ul[number_nodes-2])/step_grid_dimless*2; //step for duldz is halved for the last volume
                           dSdtheta[number_nodes-1] = -duldz[number_nodes-1];
                          // mexPrintf("fluxes_ul1=%.2e, fluxes_ul2=%.2e, fluxes_ulN-1=%.2e, fluxes_ulN=%.2e \n", fluxes_ul[0], fluxes_ul[1], fluxes_ul[number_nodes-2],fluxes_ul[number_nodes-1]);
                          // mexPrintf("duldz1=%.2e, duldz2=%.2e, duldzN-1=%.2e, duldzN=%.2e \n", duldz[0], duldz[1], duldz[number_nodes-2],duldz[number_nodes-1]);


                         }


#define DSDTHETA plhs[0]

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{
    // Variables Declarations

    int number_nodes;

    double* dSdtheta;

    double *t, *SR, step_grid_dimless, Pb, lambda, scaling, Pgin, Pgout, *Pg;

    t = mxGetPr(T);
    SR = mxGetPr(X);
    number_nodes = mxGetScalar(NUMBER_NODES);
    step_grid_dimless = mxGetScalar(STEP_GRID_DIMLESS);
    Pb = mxGetScalar(PB);
    lambda = mxGetScalar(LAMBDA);
    scaling = mxGetScalar(SCALING);
    Pgin = mxGetScalar(PGIN);
    Pgout = mxGetScalar(PGOUT);
    Pg = mxGetPr(PG);

    /* Create a matrix for the return argument */
    DSDTHETA = mxCreateDoubleMatrix((mwSize) number_nodes, (mwSize) 1, mxREAL);
    dSdtheta = mxGetPr(DSDTHETA);

    deliquoring_pde_model(t, SR, number_nodes, step_grid_dimless,
                      Pb, lambda, scaling, Pgin, Pgout,
                      Pg, dSdtheta);

    return;
}

/* LocalWords:  yp maxlhs
 */
