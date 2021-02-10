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
#define NUMBER_VOLATILE_COMPONENTS prhs[3]
#define EPSL_NON_VOL prhs[4]
#define TINLET_DRYING prhs[5]
#define CP_N2 prhs[6]
#define RHO_SOL prhs[7]
#define EE prhs[8]
#define UG0 prhs[9]
#define PPROFILE prhs[10]
#define MW_INERT prhs[11]
#define STEP_GRID_DRYING prhs[12]
#define VL_CRIT prhs[13]
#define VL_EQ prhs[14]
#define COEFF_ANTOINE prhs[15]
#define MW_COMPONENTS prhs[16]
#define H_M prhs[17]
#define A_V prhs[18]
#define RHOLCOMP prhs[19]
#define LATHEAT prhs[20]
#define CPLCOMP prhs[21]
#define CP_S prhs[22]
#define H_T_J prhs[23]


void drying_pde_model(double* t, double x[], int number_nodes,
                  int number_volatile_components, double epsL_non_vol[],
                  double Tinlet_drying, double cp_N2[],
                  double rho_sol, double E, double ug0, double Pprofile[],
                  double MW_inert, double step_grid_drying, double vl_crit[],
                  double vl_eq[], double coeff_antoine[], double MW_components[],
                  double h_M[], double a_V, double RhoLComp[],
                  double LatHeat[], double CpLComp[], double cp_s,
                  double dxdt[], double h_T_j)
                  {
                    double Tg[number_nodes], Ts[number_nodes], epsL_volatile[number_nodes],
                           epsL[number_nodes], cpgTg_in, ug_j[number_nodes], rho_g_j[number_nodes], cp_g_j[number_nodes],
                           eps_g_j, cpgTg[number_nodes], dcpgTgdz[number_nodes],
                           dcpLdt[number_nodes], Psat_volatiles, EB_S[number_nodes],
                           EB_G[number_nodes], dcp_gasdT_j[number_nodes], sum_DR_j,
                           sum_heat_rxn, wig_in, f, drying_driving_force, rho_l[number_nodes], cp_l[number_nodes],
                           wig[number_nodes*number_volatile_components], wil[number_nodes*number_volatile_components],conc_il[number_nodes*number_volatile_components],
                           xv_liq[number_nodes*number_volatile_components],DR_j[number_nodes*number_volatile_components],
                           MB_components_G[number_nodes*number_volatile_components],
                           MB_components_L[number_nodes*number_volatile_components], dwig_dz[number_nodes*number_volatile_components],
                           MWgas[number_nodes], wg_N2[number_nodes], dTgdz[number_nodes],
                           a1, a2, a3, a4, a5, b1, b2;

                    int i, j;

                    (void)t;
     
                    // drying effectiveness parameters regressed from TGA experiments on PCM (See CES paper)
                    a1=-3.56837E7;
                    a2=1.7128E6;
                    a3=-3.0190E4;
                    a4=233.6379;
                    a5=0.0086;
                    b1=3.3735;
                    b2=0.6299;
                    
                    // Initialization
                    for (j=0; j<number_nodes; j++) {
                        Tg[j]=x[j];
                        Ts[j]=Tg[j]; 
                        epsL_volatile[j]=0;
                        MWgas[j]=0;
                        wg_N2[j]=1;
                        dcpLdt[j]=0;
                            for(i=0;i<number_volatile_components;i++) {
                                  wig[j+i*number_nodes]=x[number_nodes*(1+i)+j];
                                  xv_liq[j+i*number_nodes]=x[number_nodes*(1+number_volatile_components+i)+j];

                                  MWgas[j] += wig[j+i*number_nodes]*MW_components[i];
                                  wg_N2[j] -= wig[j+i*number_nodes];
                                  epsL_volatile[j] += xv_liq[j+i*number_nodes];

                             }
                        MWgas[j] += wg_N2[j]*MW_inert;
                        epsL[j]=epsL_volatile[j]+epsL_non_vol[j];
                    }
                    cpgTg_in=Tinlet_drying*(cp_N2[0]+cp_N2[1]*Tinlet_drying+cp_N2[2]*pow(Tinlet_drying,2)+
                          cp_N2[3]*pow(Tinlet_drying,3));
                    wig_in=0;


                    // nodes from 1 to end
                    for (j=0; j< number_nodes; j++) {
                          ug_j[j]=ug0*sqrt(Tg[j]/298.);
                          rho_g_j[j]=Pprofile[j]/(8.314*Tg[j])*MWgas[j];
                          cp_g_j[j]=cp_N2[0]+cp_N2[1]*Tg[j]+cp_N2[2]*pow(Tg[j],2)+cp_N2[3]*pow(Tg[j],3);

                          eps_g_j=E-epsL[j];
                          cpgTg[j]=Tg[j]*cp_g_j[j];
                          
                          // Tg, upwind 1
                          dTgdz[j]=(Tg[j]-Tg[j-1])/(step_grid_drying);
                          dTgdz[0]=(Tg[0]-Tinlet_drying)/(step_grid_drying);                         
                          
                          dcp_gasdT_j[j]=cp_N2[1]+2*cp_N2[2]*Tg[j]+3*cp_N2[3]*pow(Tg[j],2);
                          sum_DR_j=0;
                          sum_heat_rxn=0;

                          // rho and cpL
                          rho_l[j]=0;
                          cp_l[j]=0;

                          for (i=0; i<number_volatile_components; i++) {
                              conc_il[j+i*number_nodes]=xv_liq[j+i*number_nodes]*RhoLComp[i]/epsL[j];
                              rho_l[j] += conc_il[j+i*number_nodes];
                                
                              // 1st order
                              dwig_dz[j+i*number_nodes]=(wig[j+i*number_nodes]-wig[j-1+i*number_nodes])/step_grid_drying;
                              
                              dwig_dz[0+i*number_nodes]=(wig[0+i*number_nodes]-wig_in)/step_grid_drying;

                              if (xv_liq[j+i*number_nodes]<0.016) 
                              {
                                 f=a1*pow(xv_liq[j+i*number_nodes],4)+a2*pow(xv_liq[j+i*number_nodes],3)+
                                  a3*pow(xv_liq[j+i*number_nodes],2)+a4*pow(xv_liq[j+i*number_nodes],1)+
                                  a5;
                                }
                              else if (xv_liq[j+i*number_nodes]<0.11) {
                                 f=b1*xv_liq[j+i*number_nodes]+b2;
                                }
                              
                              f=MIN(f,1);
                              f=MAX(f,0);
                                
                              // Antoine
                              Psat_volatiles = exp(coeff_antoine[0+i*5]+coeff_antoine[1+i*5]/Tg[j]+
                                coeff_antoine[2+i*5]*log(Tg[j])+coeff_antoine[3+i*5]*pow(Tg[j],coeff_antoine[4+i*5]));

                              drying_driving_force=MAX(Psat_volatiles-wig[j+i*number_nodes]/MW_components[i]*MWgas[j]*Pprofile[j],0);

                              DR_j[j+i*number_nodes]=h_M[i]*a_V*drying_driving_force*f;

                              sum_DR_j += DR_j[j+i*number_nodes];
                              sum_heat_rxn += DR_j[j+i*number_nodes]*LatHeat[i];

                              dxdt[j+number_nodes*(1+number_volatile_components+i)]= -DR_j[j+i*number_nodes]/RhoLComp[i]; //MB_components_L

                          for (i=0; i<number_volatile_components; i++) {

                                     wil[j+i*number_nodes]=conc_il[j+i*number_nodes]/rho_l[j];
                                     cp_l[j] += wil[j+i*number_nodes]*CpLComp[i];

                                     dxdt[j+number_nodes*(1+i)]= (-ug_j[j]*rho_g_j[j]*dwig_dz[j+i*number_nodes]+DR_j[j+i*number_nodes]-
                                         wig[j+i*number_nodes]*sum_DR_j)/(rho_g_j[j]*eps_g_j); //MB_components_G

                            }

                          dxdt[j]= (-sum_heat_rxn-ug_j[j]*rho_g_j[j]*cp_g_j[j]*dTgdz[j]) //EB_G
                              /(rho_g_j[j]*eps_g_j*cp_g_j[j]+rho_g_j[j]*
                              eps_g_j*Tg[j]*dcp_gasdT_j[j]+cp_s*rho_sol*(1-E)+cp_l[j]*rho_l[j]*epsL[j]);
                     }
                   }


#define DXDT plhs[0]

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{
    // Variables Declarations

    int number_nodes, number_volatile_components;

    double* dxdt;

    double *t, Tinlet_drying, *x, h_T_j,
      rho_sol, E, ug0, MW_inert, step_grid_drying, a_V, cp_s,
      *vl_crit, *vl_eq, *epsL_non_vol, *Pprofile, *cp_N2, *coeff_antoine,
      *MW_components, *h_M, *RhoLComp, *LatHeat, *CpLComp ;


    size_t m, n;


    t = mxGetPr(T);
    x = mxGetPr(X);
    number_nodes = mxGetScalar(NUMBER_NODES);
    number_volatile_components = mxGetScalar(NUMBER_VOLATILE_COMPONENTS);
    epsL_non_vol = mxGetPr(EPSL_NON_VOL);
    Tinlet_drying = mxGetScalar(TINLET_DRYING);
    cp_N2 = mxGetPr(CP_N2); 
    rho_sol = mxGetScalar(RHO_SOL);
    E = mxGetScalar(EE);
    ug0 = mxGetScalar(UG0);
    Pprofile = mxGetPr(PPROFILE);
    MW_inert = mxGetScalar(MW_INERT);
    step_grid_drying = mxGetScalar(STEP_GRID_DRYING);
    vl_crit = mxGetPr(VL_CRIT);
    vl_eq = mxGetPr(VL_EQ);
    coeff_antoine = mxGetPr(COEFF_ANTOINE);
    MW_components = mxGetPr(MW_COMPONENTS);
    h_M = mxGetPr(H_M);
    a_V = mxGetScalar(A_V);
    RhoLComp = mxGetPr(RHOLCOMP);
    LatHeat = mxGetPr(LATHEAT);
    CpLComp = mxGetPr(CPLCOMP);
    cp_s = mxGetScalar(CP_S);
    h_T_j = mxGetScalar(H_T_J);

    /* Create a matrix for the return argument */
    DXDT = mxCreateDoubleMatrix((mwSize) number_nodes*(1+2*number_volatile_components), (mwSize) 1, mxREAL);
    dxdt = mxGetPr(DXDT);


    drying_pde_model(t, x, number_nodes,
                     number_volatile_components, epsL_non_vol, Tinlet_drying, cp_N2,
                     rho_sol, E, ug0, Pprofile, MW_inert, step_grid_drying,
                     vl_crit, vl_eq, coeff_antoine, MW_components, h_M, a_V,
                     RhoLComp, LatHeat, CpLComp, cp_s, dxdt, h_T_j);

    return;
}

/* LocalWords:  yp maxlhs
 */
