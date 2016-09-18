/*******************************************************************************
 (Spectral method)
 
 Quantum dynamics (QD) simulation of an electron in one dimension.
 
 USAGE
 
 %cc -o sqd sqd.c four1.c -lm
 %sqd < sqd.in (see qd1.h for the input-file format)
 *******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "sqd.h"

int main(int argc, char **argv) {
    int step; /* Simulation loop iteration index */
    
    init_param();  /* Read input parameters */
    init_prop();   /* Initialize the kinetic & potential propagators */
    init_wavefn(); /* Initialize the electron wave function */
    
    for (step=1; step<=NSTEP; step++) {
        single_step(); /* Time propagation for one step, DT */
        if (step%NECAL==0) {
            calc_energy();
            printf("%le %le %le %le\n",DT*step,ekin,epot,etot);
        }
    }
    
    calc_den();
    
    return 0;
}

/*----------------------------------------------------------------------------*/
void init_param() {
    /*------------------------------------------------------------------------------
     Initializes parameters by reading them from standard input.
     ------------------------------------------------------------------------------*/
    /* Read control parameters */
    scanf("%le",&LX);
    scanf("%le",&DT);
    scanf("%d",&NSTEP);
    scanf("%d",&NECAL);
    scanf("%le%le%le",&X0,&S0,&E0);
    scanf("%le%le",&BH,&BW);
    scanf("%le",&EH);
    
    /* Calculate the mesh size */
    dx = LX/NX;
}

/*----------------------------------------------------------------------------*/
void init_prop() {
    /*------------------------------------------------------------------------------
     Initializes the kinetic & potential propagators.
     ------------------------------------------------------------------------------*/
    int i;
    double x;
    int m;
    
    for (m=0; m<NX; m++) {
        if (m<NX/2)
            km[m]=2*M_PI*m/LX;
        else
            km[m] = 2*M_PI*(m-NX)/LX;
        ut[2*m+0] = cos(-0.5*DT*km[m]*km[m]);
        ut[2*m+1] = sin(-0.5*DT*km[m]*km[m]);
    }
    
    
    /* Set up potential propagator */
    for (i=0; i<=NX-1; i++) {
        x = dx*i;
        /* Construct the edge potential */
        if (i==0 || i==NX-1)
            v[i] = EH;
        /* Construct the barrier potential */
        else if (0.5*(LX-BW)<x && x<0.5*(LX+BW))
            v[i] = BH;
        else
            v[i] = 0.0;
        /* Half-step potential propagator */
        uv[2*i+0] = cos(-0.5*DT*v[i]);
        uv[2*i+1] = sin(-0.5*DT*v[i]);
    }
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
    /*------------------------------------------------------------------------------
     Initializes the wave function as a traveling Gaussian wave packet.
     ------------------------------------------------------------------------------*/
    int sx,s;
    double x,gauss,psisq,norm_fac;
    
    /* Calculate the the wave function value mesh point-by-point */
    for (sx=0; sx<=NX-1; sx++) {
        x = dx*sx-X0;
        gauss = exp(-0.25*x*x/(S0*S0));
        psi[2*sx+0] = gauss*cos(sqrt(2.0*E0)*x);
        psi[2*sx+1] = gauss*sin(sqrt(2.0*E0)*x);
    }
    
    /* Normalize the wave function */
    psisq=0.0;
    for (sx=0; sx<NX; sx++)
        for (s=0; s<2; s++)
            psisq += psi[2*sx+s]*psi[2*sx+s];
    psisq *= dx;
    norm_fac = 1.0/sqrt(psisq);
    for (sx=0; sx<NX; sx++)
        for (s=0; s<2; s++)
            psi[2*sx+s] *= norm_fac;
}

/*----------------------------------------------------------------------------*/
void single_step() {
    /*------------------------------------------------------------------------------
     Propagates the electron wave function for a unit time step, DT.
     ------------------------------------------------------------------------------*/
    int j;
    
    prop(uv);
    four1(psi-1, (unsigned long) NX, -1);
    for (j=0; j<2*NX; j++)
        psi[j]/= NX;
    
    prop(ut);
    four1(psi-1, (unsigned long) NX, 1);
    
    prop(uv);
}

/*----------------------------------------------------------------------------*/
void prop(double u[]) {
    /*------------------------------------------------------------------------------
     Potential propagator for a half time step, DT/2.
     ------------------------------------------------------------------------------*/
    int sx;
    double wr,wi;
    
    for (sx=0; sx<NX; sx++) {
        wr = u[2*sx+0]*psi[2*sx+0]-u[2*sx+1]*psi[2*sx+1];
        wi = u[2*sx+0]*psi[2*sx+1]+u[2*sx+1]*psi[2*sx+0];
        psi[2*sx+0]=wr;
        psi[2*sx+1]=wi;
    }
}

/*----------------------------------------------------------------------------*/
void calc_energy() {
    /*------------------------------------------------------------------------------
     Calculates the kinetic, potential & total energies, EKIN, EPOT & ETOT.
     ------------------------------------------------------------------------------*/
    int sx;
    int m;
    int j;
    
    
    four1(psi-1, (unsigned long) NX, -1);
    for (j=0; j<2*NX; j++)
        psi[j]/= NX;
    
    /* Kinetic energy = <PSI|(-1/2)Laplacian|PSI> = <PSI|WRK> */
    ekin = 0.0;
    for (m=0; m<NX; m++)
        ekin += 0.5*km[m]*km[m]*(psi[2*m+0]*psi[2*m+0]+psi[2*m+1]*psi[2*m+1]);
    ekin *= dx*NX;
    
    four1(psi-1, (unsigned long) NX, 1);
    
    /* Potential energy */
    epot = 0.0;
    for (sx=0; sx<NX; sx++)
        epot += v[sx]*(psi[2*sx+0]*psi[2*sx+0]+psi[2*sx+1]*psi[2*sx+1]);
    epot *= dx;
    
    /* Total energy */
    etot = ekin+epot;
}


/* Calculates the probability density at each position at the final time */

void calc_den() {
    int k;
    
    for (k=0; k<NX; k++) {
        real = psi[2*k];
        imag = psi[2*k+1];
        den = psi[2*k]*psi[2*k]+psi[2*k+1]*psi[2*k+1];
        printf("%le %le %le %le\n",dx*k,den,real,imag);
    }
    
}

