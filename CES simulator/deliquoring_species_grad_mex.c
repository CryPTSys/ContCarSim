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
#define NUMBER_COMPONENTS prhs[10]
#define S_INF prhs[11]



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



void deliquoring_pde_model(double* t, double x[], int number_nodes, double step_grid_dimless,
                  double Pb, double lambda, double scaling, double Pgin, double Pgout,
                  double Pg[], double dXdtheta[], int number_components, double S_inf)

                  {
                    double ulin, Pl[number_nodes], kl[number_nodes],
                            dPldz[number_nodes], duldz[number_nodes], flux_Pl_in,
                            fluxes_ul[number_nodes], ul[number_nodes], fluxes_Pl[number_nodes],
                            flux_ul_in, SR[number_nodes], S[number_nodes], wi[number_nodes*number_components] ,
                            dwidz[number_nodes*number_components] ;

                    int i, j;
                           // First Loop
                           // first node
                           ulin=0;
                           flux_ul_in=ulin;
                           // get SR, S and wi for the first two nodes
                           SR[0]=x[0];
                           SR[1]=x[1];
                           S[0]=S_inf+SR[0]*(1-S_inf);
                           S[1]=S_inf+SR[1]*(1-S_inf);
                           for (i=0;i<number_components;i++) {
                             wi[0+number_nodes*i]=x[number_nodes+0+number_nodes*i];
                             dwidz[0+number_nodes*i]=0;
                             wi[1+number_nodes*i]=x[number_nodes+2+number_nodes*i];
                             dwidz[1+number_nodes*i]=(wi[1+number_nodes*i]-wi[0+number_nodes*i])/step_grid_dimless;
                           }

                           // calculations first two nodes
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
                             for (i=0;i<number_components;i++) {
                               wi[j+number_nodes*i]=x[number_nodes+j+number_nodes*i];

                               dwidz[j+number_nodes*i]=(wi[j+number_nodes*i]-wi[j-1+number_nodes*i])/step_grid_dimless;
                              }

                              SR[j+1]=x[j+1];
                              S[j+1]=S_inf+SR[j+1]*(1-S_inf);
                              Pl[j+1]=(Pg[j+1]-Pb*pow(SR[j+1],-1/lambda))/Pb/scaling;
                              kl[j+1]=pow(SR[j+1],(2+3*lambda)/lambda);
                              fluxes_Pl[j]=fluxes_HRFVM_vL(Pl[j-1], Pl[j], Pl[j+1], j, number_nodes);
                              dPldz[j]=(fluxes_Pl[j]-fluxes_Pl[j-1])/step_grid_dimless;
                              ul[j]=-kl[j]*dPldz[j];
                              //mexPrintf("dwidz1=%.2e\n", Pl[j]);
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
                           dXdtheta[0] = -duldz[0];


                           for (j=1; j<number_nodes-1;j++) {
                              fluxes_ul[j]=fluxes_HRFVM_vL(ul[j-1], ul[j], ul[j+1], j, number_nodes);
                              duldz[j]=(fluxes_ul[j]-fluxes_ul[j-1])/step_grid_dimless;
                              dXdtheta[j] = -duldz[j];

                              for (i=0;i<number_components;i++) {
                                  dXdtheta[number_nodes+j+i*number_nodes]=-(1-S_inf)*ul[j]/S[j]*dwidz[j+number_nodes*i];

                                }
                              //  mexPrintf("wj1=%.3e wj2=%.3e wj3=%.3e\n", dXdtheta[j+number_nodes*0], dXdtheta[j+number_nodes*1], dXdtheta[j+number_nodes*2]);
                              //mexPrintf("ul=%.3e S=%.3e dwidz=%.3e\n", ul[j], S[j], S_inf);

                           }
                            duldz[number_nodes-1]=(fluxes_ul[number_nodes-1]-fluxes_ul[number_nodes-2])/step_grid_dimless*2; //step for duldz is halved for the last volume
                            dXdtheta[number_nodes-1] = -duldz[number_nodes-1];
                          // mexPrintf("duldz1=%.2e, duldz2=%.2e, duldzN-1=%.2e, duldzN=%.2e \n", duldz[0], duldz[1], duldz[number_nodes-2],duldz[number_nodes-1]);
                            for (i=0;i<number_components;i++) {
                                wi[number_nodes-1+number_nodes*i]=x[2*number_nodes-1+number_nodes*i];
                                dwidz[number_nodes-1+number_nodes*i]=(wi[number_nodes-1+number_nodes*i]-wi[number_nodes-2+number_nodes*i])/step_grid_dimless;
                                j=0;
                                dXdtheta[number_nodes+j+i*number_nodes]=0;
                                j=number_nodes-1;
                                dXdtheta[number_nodes+j+i*number_nodes]=-(1-S_inf)*ul[j]/S[j]*dwidz[j+number_nodes*i];
                          }
                  /*        mexPrintf("ul0=%.2e, ul1=%.2e, ul2=%.2e, , dul0=%.2e, dul1=%.2e, dul2=%.2e \n", ul[0], ul[1],ul[2], duldz[0], duldz[1],duldz[2]);
                          mexPrintf("ul3=%.2e, ul4=%.2e, ul5=%.2e, , dul3=%.2e, dul4=%.2e, dul5=%.2e \n", ul[3], ul[4],ul[5], duldz[3], duldz[4],duldz[5]);
                          mexPrintf("fluxes0=%.5e, fluxes1=%.5e, fluxes2=%.5e,", fluxes_ul[0], fluxes_ul[1],fluxes_ul[2]);
                          mexPrintf("fluxes3=%.5e, fluxes4=%.5e, fluxes5=%.5e,", fluxes_ul[3], fluxes_ul[4],fluxes_ul[5]); */

                         }


#define DXDTHETA plhs[0]

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{
    // Variables Declarations

    int number_nodes, number_components;

    double* dXdtheta;

    double *t, *SR, step_grid_dimless, Pb, lambda, scaling, Pgin, Pgout, *Pg, S_inf;


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
    number_components = mxGetScalar(NUMBER_COMPONENTS);
    S_inf = mxGetScalar(S_INF);

    /* Create a matrix for the return argument */
    DXDTHETA = mxCreateDoubleMatrix((mwSize) number_nodes*(1+number_components), (mwSize) 1, mxREAL);
    dXdtheta = mxGetPr(DXDTHETA);

    deliquoring_pde_model(t, SR, number_nodes, step_grid_dimless,
                      Pb, lambda, scaling, Pgin, Pgout,
                      Pg, dXdtheta, number_components, S_inf);

    return;
}
