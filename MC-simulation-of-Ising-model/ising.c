/* Monte Carlo simulation of the Ising model */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define L 20 //Lattice size
int s[L][L]; //Spins s[i][j]
int s_new, s_neighbor, k, l;
int Sta_step = 2000000; //no. of iterations
double exp_dV[2][5];
double JdivT; //J/kBT
double HdivT = 0.0; //H/kBT
double dV;

void table_set() {
    for (k=0; k<2; k++) {
        s_new = 2*k-1;
        for (l=0; l<5; l++) {
            s_neighbor = 2*l-4;
            exp_dV[k][l] = exp(2.0*s_new*(JdivT*s_neighbor+HdivT));
        }
    }
} //pre-compute the acceptance probability

int main() {
    double x, pi, sigM, sumM = 0.0, sumM2 = 0.0, avgM, hist[2*L*L+1], exp_val, runM;
    int i, j, im, ip, jm, jp, step;
    
    printf("Input JdivT\n");
    scanf("%le",&JdivT);
    
    table_set();
    
    for (i=0; i<L; i++) {
        for (j=0; j<L; j++) {
            s[i][j] = 1; //cold start
        }
    }
    
    runM = 1.0*L*L;
    
    for (i=0; i<2*L*L+1; i++) {
        hist[i] = 0;
    } //initialize histogram
    
    srand((unsigned)time((NULL)));
    for (step=1; step<=Sta_step; step++) {
        //randomly select a grid point (i,j)
        i = rand()%L;
        j = rand()%L;
        //compute the change in potential eneergy dV with a single spin
        s_new = -s[i][j]; //flip
        im = (i+L-1)%L; //up
        ip = (i+1)%L; //down
        jm = (j+L-1)%L; //left
        jp = (j+1)%L; //right
        s_neighbor = s[im][j]+s[ip][j]+s[i][jm]+s[i][jp];
        k = (1+s_new)/2; //hash functions
        l = (4+s_neighbor)/2;
        exp_val = exp_dV[k][l];
        
        if (exp_val > 1.0) {
            s[i][j] = s_new; //update spin
            runM += 2.0*s_new; //update magnetization
        } //unconditionally accept
        
        else if (rand()/(double)RAND_MAX <= exp_val) {
            s[i][j] = s_new;
            runM += 2.0*s_new;
        } //conditionally accept
        
        sumM += runM;
        sumM2 += runM*runM;
        ++hist[(int)runM+L*L];
    }
    
    avgM = fabs(sumM/Sta_step); //absolute value of the mean magnetization
    sigM = sqrt(sumM2/Sta_step-avgM*avgM);
    
    printf("Magnetization = %e %e\n",avgM,sigM);
    for (i=0; i<2*L*L+1; i++) {
        printf("%d %f\n", i-L*L, hist[i]);
    }
    return 0;
}

