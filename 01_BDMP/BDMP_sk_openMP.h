#include "in_namelist.h"

// Parameters for General Simulation
double a, eta, T, Salt, dt, tmax, Nf;
int ibm, mag, N;
/* a: Radius of beads in m
 * eta: Solvent viscosity in Pa*s
 * T: Temperature of the system in K
 * Salt: Concentration of salt
 * dt: Time of each time step in sec
 * tmax: Maximum simulation time in sec
 * Nf: Number of frame captured per sec
 * ibm: 0,1,2 - non-Brownian, Brownian, Hydrodynamic-coupled Brownian
 * N: Number of beads in the system
 */



// Magnetic Parameters
double Bt, Chi;
/* Bt: Magnetic field strength in Tesla
 * Chi: Particle susceptibility
 */



// Chain Parameters
double Lp, L0, k0, C1, C2, C3;
/* Lp: Persistence length in m
 * L0: Relaxed length of Hookesian spring between beads in m
 * k0: Hookesian ppring constant
 * C1, C2, C3: fitting variables from relationship between surface to surface
 * distance D(nm) and force Fsrfrpl(pN) by Fsrfrpl = -C1*ln(D)/D^2+C1/D^2-C2/D^2+2*C3*D.
 */


// Buckling Parameters
double ti, trex;
/* ti: Initial Time
 * trex: Time of Relaxation of field (Rinaldi JAP 2014)
 */


// OpenMP Parameters
int Td;
/* Number of threads used for Simulation */



NameList nameList[] = {
    NameR(a),         // Particle Raidus in m
    NameR(eta),       // Solvent viscosity in Pa*s
    NameR(T),         // Temperature in K
    NameR(Nf),        // Number of frames recorded per second
    NameR(dt),        // Time step in sec, usually 1-2 us
    NameR(tmax),      // Total running time in sec
    NameR(Lp),        // Persistence Length in m
    NameR(L0),        // Relaxed length between particles in m
    NameR(k0),        // Hooke's spring constant of stiffness N/m
    NameR(C1),        // Fitting variables from relationship between...
    NameR(C2),        // surface to surface distance D(nm) and force Fsrfrpl(pN)...
    NameR(C3),        // by Fsrfrpl = -C1*ln(D)/D^2+C1/D^2-C2/D^2+2*C3*D.
    NameR(Bt),        // Magnetic field strength in Tesla
    NameR(Chi),       // Particle susceptibility
    NameR(ti),        // Initial time in sec
    NameR(trex),      // Time of relaxation of field (Rinaldi JAP 2014)
    NameI(ibm),       // Include brownian motion or not
    NameI(mag),       // Magnetic force calculation: 1 for DM and 2 for MDM
    NameI(N),         // Number of particles
    NameI(Td)         // Number of threads for openMP
};
