/***** Brownian Dynamic Simulation of chain in 2D *****/
/* Parallel Implementation with openMP */

/* Author: Ziyi Zhu
 * Contact: wazzytrent@gmail.com
 * Last_updated: March. 9th, 2017
 */

// EDITED: steve kuei (kuei.steve@rice.edu)
// Last Updated: 080417

/* References
1. Ermak and McCammon JCP 1978
2. Rotne and Prager JCP 1969
3. Li and Biswal Soft Matter 2010
4. Du, Li and Biswal Soft Matter 2013
5. Du, Toffo and Biswal PRE 2014
6. Allison Macromolecule 1986
7. Gauger and Stark PRE 2006
8. Crocker PRL 1996
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "vector.h"
#include "utility.h"
#include "BDMP_sk_openMP.h"

#include "in_namelist.c"

/* BROWNIAN DYNAMICS SIMULATION Main Function */
/* Ermak-McCammon scheme */

int main(int argc, char **argv) {
    /*** Step 1: Read in Namelist Parameters from Input File ***/
    GetNameList(argc, argv);
    PrintNameList(stdout);
    if (ibm != 0 && ibm != 1 && ibm != 2) {
        printf("Invalid brownian motion calculator!");
        exit (1);
    }
    if (mag != 1 && mag != 2) {
        printf("Invalid magnetic force calculator!");
        exit(1);
    }


    /*** Step 2: Specify Parameters ***/
    // Constant parameters
    double kb = 1.38e-23;   //Boltzman constant

    // Magnetic Force Parameters
    double mu0 = pi * 4.0e-7;   //permeability in vacuum
    double Ht = Bt/mu0;         //field strength in A/m

    // Bending constant
    double g_new = Lp*kb*T / (2*a+L0) / (2*a+L0); // g/L0 in Gauger's presentation

    // Other Parameters
    int kf = (int)(1.0/dt/Nf+0.5);   //frame parameter
    int Nup = 2;                 //moment updating factor

    /*** Step 3: Declaration of variables and arrays ***/
    FILE* fp;     //file pointer for output
    double t;     //simulation timer
    int istep;    //counter for time step
    int i, j, k;  //loop variables

    double R[N][N];
    VecR2 u[N][N];
    VecR2 r[N], DF[N], dR[N];
    VecR2 Fmagc[N][N], Fstret[N], Fbend[N], Ftot[N];
    VecR2 Dx[N][N], Dy[N][N], Dz[N][N], eff_f[N];
    VecR2 Mdp, H0, Hind, rv;
    VecR2 tmpstret, tmpsrfrpl;

    // Hydrodynamic Tensor Parameters
    double rpt_coeff1, rpt_coeff2, rpt_coeff3;
    // Random Displacement Parameters
    double Dif, rand_coeff1, rand_coeff2;
    //VecR2 fj;
    // Magnetic Parameter
    double Mdp2, Mdpt, rnorm;
    // Surface repulsive force parameter
    double Rsfc; //Surface to surface distance in nm


    /*** Step 4: Initialization ***/
    t = 0;
    istep = 0;

    Dif = kb*T/(6*pi*eta*a);
    rand_coeff1 = sqrt(2*kb*T/dt)*sqrt(6*pi*eta*a);  //sqrt(2*(kb*T)^2/(D*dt))
    rand_coeff2 = Dif*dt/(kb*T);

    Mdpt = 4.0/3.0*pi*Cub(a)*Chi;

    // initial positions: straight line
    double rxtmp = 0;
    for (i=0; i<N; i++) {
        VSet2(r[i], i*(2*a+L0), 0);
        rxtmp += r[i].x;
    }
    rxtmp = rxtmp/N;
    // move the mass center back to origin
    for (i=0; i<N; i++) {
        r[i].x -= rxtmp;
    }

    // open the output file to write
    fp = fopen(argv[2],"w");
    if (fp == NULL) {
        printf("can't open file");
        exit(1);
    }
    setlinebuf(fp);


    /*** Step 5: BD Simulation Loop ***/
    srand((unsigned) time(NULL));   //system time seed

    while (t < tmax) {
        if (istep % kf == 0) {
            for (i=0; i<N; i++) {
                fprintf(fp, "%15.4e\t", r[i].x);
            }
            for (i=0; i<N; i++) {
                fprintf(fp, "%15.4e\t", r[i].y);
            }
            fprintf(fp,"\n");
        }

        t = t + dt;
        istep = istep + 1;

        //TODO sk: modify mag force. currently constant in X direction
        // VSet2(H0, Ht, 0);
        VSet2(H0, (Ht*cos(f * 2.0 * pi * (t - ti))), (Ht*sin(f * 2.0 * pi * (t - ti))));
        VSCopy2(Mdp, Mdpt, H0);
        Mdp2 = VDot2(Mdp,Mdp);
        // END MAG FORCE


        #pragma omp parallel num_threads(Td)
        {
        /*******************************************/

        #pragma omp for private(i) nowait
        for (i=0; i<N; i++) {
            VSet2(Ftot[i], 0, 0);
            VSet2(DF[i], 0, 0);
            VSet2(eff_f[i], 0, 0);
            if (ibm == 0) {
	           VSet2(dR[i], 0, 0);
            } else {
		VSet2(dR[i], normrand(), normrand());
            }
        }

        /*******************************************/

        #pragma omp for private(i,j) collapse(2)
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                // Relative Distance
                R[i][j] = sqrt(Sqr(r[j].x-r[i].x)+Sqr(r[j].y-r[i].y));
                // Unit Vector
                if (i==j) {
                    VSet2(u[i][j], 0, 0);
                } else {
                    u[i][j].x = (r[j].x-r[i].x)/R[i][j];
                    u[i][j].y = (r[j].y-r[i].y)/R[i][j];
                }
            }
        }

        /*******************************************/

        #pragma omp for private(i,j) //collapse(2)
        for (i=0; i<N-1; i++) {
            for (j=i+1; j<N; j++) {
                Fmagc[i][j].x = -3*mu0/(4*pi*pow(R[i][j],4)) * ((Mdp2-5.0*Sqr(VDot2(Mdp,u[i][j])))*u[i][j].x + VDot2(Mdp,u[i][j])*Mdp.x + VDot2(Mdp,u[i][j])*Mdp.x);
                Fmagc[i][j].y = -3*mu0/(4*pi*pow(R[i][j],4)) * ((Mdp2-5.0*Sqr(VDot2(Mdp,u[i][j])))*u[i][j].y + VDot2(Mdp,u[i][j])*Mdp.y + VDot2(Mdp,u[i][j])*Mdp.y);
                Fmagc[j][i].x = -Fmagc[i][j].x;
                Fmagc[j][i].y = -Fmagc[i][j].y;
            }
        }

        /*******************************************/

        #pragma omp for private(i, j, Rsfc, tmpstret, tmpsrfrpl)
        for (i=0; i<N; i++) {
            // Stretching and Surface Repulsive Forces
            if (i == 0) {
                VSCopy2(Fstret[0], (k0*(R[0][1] - 2 * a - L0)), u[0][1]);
                // surface repuslive forces
                //Rsfc = (R[0][1]-2*a)*1.0e9; //Surface to surface distance in nm, Rsfc[i-1][i]
                //VSCopy2(Fsrfrpl[0], ((-C1*log(Rsfc)/Sqr(Rsfc)+C1/Sqr(Rsfc)-C2/Sqr(Rsfc)+2*C3*Rsfc)*1.0e-12), u[1][0]);
            } else if (i == N-1) {
                VSCopy2(Fstret[N-1], (-k0*(R[N-2][N-1] - 2 * a - L0)), u[N-2][N-1]);
                // surface repuslive forces
                //Rsfc = (R[N-2][N-1]-2*a)*1.0e9;
                //VSCopy2(Fsrfrpl[N-1], ((-C1*log(Rsfc)/Sqr(Rsfc)+C1/Sqr(Rsfc)-C2/Sqr(Rsfc)+2*C3*Rsfc)*1.0e-12), u[N-2][N-1]);
            } else {
                VSCopy2(tmpstret, (-k0*(R[i-1][i] - 2 * a - L0)), u[i - 1][i]);
                VSAdd2(Fstret[i], tmpstret, (k0*(R[i][i+1] - 2 * a - L0)), u[i][i+1]);
                // surface repuslive forces
                //Rsfc = (R[i-1][i]-2*a)*1.0e9;
                //VSCopy2(tmpsrfrpl, ((-C1*log(Rsfc)/Sqr(Rsfc)+C1/Sqr(Rsfc)-C2/Sqr(Rsfc)+2*C3*Rsfc)*1.0e-12), u[i-1][i]);
                //Rsfc = (R[i][i+1]-2*a)*1.0e9;
                //VSAdd2(Fsrfrpl[i], tmpsrfrpl, ((-C1*log(Rsfc)/Sqr(Rsfc)+C1/Sqr(Rsfc)-C2/Sqr(Rsfc)+2*C3*Rsfc)*1.0e-12), u[i+1][i]);
            }

            // Bending
            if (i == 0) {
                Fbend[0].x = (VDot2(u[0][1], u[1][2]))*u[0][1].x - u[1][2].x;
                Fbend[0].y = (VDot2(u[0][1], u[1][2]))*u[0][1].y - u[1][2].y;
            } else if (i == 1) {
                Fbend[1].x = -(1 + VDot2(u[0][1],u[1][2]))*u[0][1].x + (1 + VDot2(u[0][1],u[1][2]) + VDot2(u[1][2],u[2][3]))*u[1][2].x - u[2][3].x;
                Fbend[1].y = -(1 + VDot2(u[0][1],u[1][2]))*u[0][1].y + (1 + VDot2(u[0][1],u[1][2]) + VDot2(u[1][2],u[2][3]))*u[1][2].y - u[2][3].y;
            } else if (i == N-2) {
                Fbend[N-2].x = u[N-4][N-3].x - (VDot2(u[N-4][N-3],u[N-3][N-2]) + 1 + VDot2(u[N-3][N-2],u[N-2][N-1])) * u[N-3][N-2].x + (VDot2(u[N-3][N-2],u[N-2][N-1]) + 1) * u[N-2][N-1].x;
                Fbend[N-2].y = u[N-4][N-3].y - (VDot2(u[N-4][N-3],u[N-3][N-2]) + 1 + VDot2(u[N-3][N-2],u[N-2][N-1])) * u[N-3][N-2].y + (VDot2(u[N-3][N-2],u[N-2][N-1]) + 1) * u[N-2][N-1].y;
            } else if (i == N-1) {
                Fbend[N-1].x = u[N-3][N-2].x - VDot2(u[N-3][N-2],u[N-2][N-1]) * u[N-2][N-1].x;
                Fbend[N-1].y = u[N-3][N-2].y - VDot2(u[N-3][N-2],u[N-2][N-1]) * u[N-2][N-1].y;
            } else {
                Fbend[i].x = u[i-2][i-1].x - (VDot2(u[i-2][i-1],u[i-1][i]) + 1 + VDot2(u[i-1][i],u[i][i+1])) * u[i-1][i].x + (VDot2(u[i-1][i],u[i][i+1]) + 1 + VDot2(u[i][i+1],u[i+1][i+2])) * u[i][i+1].x - u[i+1][i+2].x;
                Fbend[i].y = u[i-2][i-1].y - (VDot2(u[i-2][i-1],u[i-1][i]) + 1 + VDot2(u[i-1][i],u[i][i+1])) * u[i-1][i].y + (VDot2(u[i-1][i],u[i][i+1]) + 1 + VDot2(u[i][i+1],u[i+1][i+2])) * u[i][i+1].y - u[i+1][i+2].y;
            }
            VSCopy2(Fbend[i], g_new, Fbend[i]);

            // Total Force
            VAdd2(Ftot[i], Ftot[i], Fstret[i]);
            //VAdd2(Ftot[i], Ftot[i], Fsrfrpl[i]);
            VAdd2(Ftot[i], Ftot[i], Fbend[i]);
            for (j=0; j<N; j++) {
                VAdd2(Ftot[i], Ftot[i], Fmagc[i][j]);
            }
        }

        /*******************************************/

        #pragma omp for private(i, j, rpt_coeff1, rpt_coeff2, rpt_coeff3)
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                if (i==j) {
                    // Compute Hydrodynamics Tensor
                    VSet2(Dx[i][j], Dif, 0);
                    VSet2(Dy[i][j], 0, Dif);
                    // Compute deterministic displacement
                    DF[i].x += Dif*Ftot[j].x;
                    DF[i].y += Dif*Ftot[j].y;
                } else {
                    // Compute Hydrodynamics Tensor
                    rpt_coeff1 = kb*T/(8.0*pi*eta*R[i][j]);
                    rpt_coeff2 = 1+2.0/3.0*Sqr(a)/Sqr(R[i][j]);
                    rpt_coeff3 = 1-2.0*Sqr(a)/Sqr(R[i][j]);

                    Dx[i][j].x = rpt_coeff1*(rpt_coeff2+rpt_coeff3*Sqr(u[i][j].x));
                    Dx[i][j].y = rpt_coeff1*(rpt_coeff3*u[i][j].x*u[i][j].y);
                    Dy[i][j].x = rpt_coeff1*(rpt_coeff3*u[i][j].y*u[i][j].x);
                    Dy[i][j].y = rpt_coeff1*(rpt_coeff2+rpt_coeff3*Sqr(u[i][j].y));

                    // Compute deterministic displacement from hydro-coupling
                    DF[i].x += VDot2(Dx[i][j], Ftot[j]);
                    DF[i].y += VDot2(Dy[i][j], Ftot[j]);
                }

                if (ibm == 2) {
                    // Compute effective random force
                    VecR2 fj;
                    fj.x = rand_coeff1 * dR[j].x;
                    fj.y = rand_coeff1 * dR[j].y;
                    if (i==j) {
                        eff_f[i].x += fj.x;
                        eff_f[i].y += fj.y;
                    } else {
                        eff_f[i].x += (Dx[i][j].x/Dif)*fj.x + (Dx[i][j].y/Dif)*fj.y;
                        eff_f[i].y += (Dy[i][j].x/Dif)*fj.x + (Dy[i][j].y/Dif)*fj.y;
                    }
                }
            }

            // Compute Random Displacement
            if (ibm == 2) {
                // with Hydrodynamics
                VSCopy2(dR[i], rand_coeff2, eff_f[i]);
            } else if (ibm == 1) {
                // without Hydrodynamics
                VSCopy2(dR[i], sqrt(2.0*Dif*dt), dR[i]);
            } else {}

            // Update position
            r[i].x += DF[i].x*dt/kb/T + dR[i].x;
            r[i].y += DF[i].y*dt/kb/T + dR[i].y;
        }

        /*******************************************/
        }

        // The end of one time step (while loop)
    }

    /* The end of BD simulation */
    fclose(fp);
    return 0;
}
