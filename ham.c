#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tb_util.c"
#include "eigen.c"

int n4;

/* For periodic boundary condition ********************************************/
double RegionH[3];
double SignR(double v, double x) {if (x>0) return v; else return -v;}

void Hamiltonian(){
    double **h; /* Hamiltonian matrix */
    double *d; /* Eigenvalues */
    double *e; /* Work array for matrix diagonalization */
    n4 = 4*nAtom;
    h = dmatrix(1,n4,1,n4);
    d = dvector(1,n4);
    e = dvector(1,n4);
    double r0=2.360352;
    double Es=-5.25,Ep=1.20;
    double h_lambda0[]={-2.038,1.745,2.75,-1.075};
    double h_lambda[4];
    double n_lambda[]={9.5,8.5,7.5,7.5};
    double r_lambda[]={3.4,3.55,3.7,3.7};
    double unit_r_relative[3];
    double length_r_relative;
    
    int i,j,k;
    
    for(i=1;i<=n4;i++){
        for(j=1;j<=n4;j++){
            h[i][j]=0;
        }
    }
    for(i=0;i<nAtom;i++){
        h[4*i+1][4*i+1]=Es;
        h[4*i+2][4*i+2]=Ep;
        h[4*i+3][4*i+3]=Ep;
        h[4*i+4][4*i+4]=Ep;
    }
    for(i=0;i<nAtom;i++){
        for(j=i+1;j<nAtom;j++){
            length_r_relative=0;
            for(k=0;k<3;k++){
                RegionH[k] = 0.5*LCNS*InitUcell[k];
                unit_r_relative[k]=r[i][k]-r[j][k];
                //Boundary condition ---->
                unit_r_relative[k] = unit_r_relative[k]
                -SignR(RegionH[k],unit_r_relative[k]-RegionH[k])
                -SignR(RegionH[k],unit_r_relative[k]+RegionH[k]);
                // <----
                length_r_relative+=pow(unit_r_relative[k],2);
            }
            length_r_relative=sqrt(length_r_relative);
            for(k=0;k<3;k++){
                unit_r_relative[k]/=length_r_relative;
            }
            //At this point unit_r_relative is the unit vector from atom i to j.
            
            for(k=0;k<4;k++){
                h_lambda[k]=h_lambda0[k]*pow(r0/length_r_relative,2)*
                exp(2*(-pow(length_r_relative/r_lambda[k],n_lambda[k])+ pow(r0/r_lambda[k],n_lambda[k])));
            }
            
            h[4*i+1][4*j+1]=h_lambda[0];
            h[4*i+1][4*j+2]=unit_r_relative[0]*h_lambda[1];
            h[4*i+1][4*j+3]=unit_r_relative[1]*h_lambda[1];
            h[4*i+1][4*j+4]=unit_r_relative[2]*h_lambda[1];
            
            h[4*i+2][4*j+1]=-unit_r_relative[0]*h_lambda[1];
            h[4*i+2][4*j+2]=unit_r_relative[0]*unit_r_relative[0]*h_lambda[2]
            + (1-unit_r_relative[0]*unit_r_relative[0])*h_lambda[3];
            h[4*i+2][4*j+3]=unit_r_relative[0]*unit_r_relative[1]*(h_lambda[2]-h_lambda[3]);
            h[4*i+2][4*j+4]=unit_r_relative[0]*unit_r_relative[2]*(h_lambda[2]-h_lambda[3]);
            
            h[4*i+3][4*j+1]=-unit_r_relative[1]*h_lambda[1];
            h[4*i+3][4*j+2]= unit_r_relative[0]*unit_r_relative[1]*(h_lambda[2]-h_lambda[3]);
            h[4*i+3][4*j+3]=unit_r_relative[1]*unit_r_relative[1]*h_lambda[2]
            + (1-unit_r_relative[1]*unit_r_relative[1])*h_lambda[3];
            h[4*i+3][4*j+4]=unit_r_relative[1]*unit_r_relative[2]*(h_lambda[2]-h_lambda[3]);
            
            h[4*i+4][4*j+1]=-unit_r_relative[2]*h_lambda[1];
            h[4*i+4][4*j+2]=unit_r_relative[0]*unit_r_relative[2]*(h_lambda[2]-h_lambda[3]);
            h[4*i+4][4*j+3]=unit_r_relative[1]*unit_r_relative[2]*(h_lambda[2]-h_lambda[3]);
            h[4*i+4][4*j+4]=unit_r_relative[2]*unit_r_relative[2]*h_lambda[2]
            + (1-unit_r_relative[2]*unit_r_relative[2])*h_lambda[3];
            
        }
    }
    
    for(i=1;i<n4;i++){
        for(j=i+1;j<n4;j++){
            h[j][i]=h[i][j];
        }
    }
    tred2(h,n4,d,e);
    
    tqli(d,e,n4,h);
    
    for(i=1;i<=n4;i++){
        printf("%le ",d[i]);
    }
    
}

int main()
{
    InitConf();
    Hamiltonian();
    return 0;
}
