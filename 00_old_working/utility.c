/* Utility Module */

#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "utility.h"

/****************************************************************/
/* Random Number Generator Module **
 ***********************************/
double nrand()
/*uniform distribution,(0,1]*/
{
    return (rand()+1.0)/(RAND_MAX+1.0);
}

double normrand()
/*normal distribution, centered on 0, std dev 1, standard norm distribution*/
// Box-Muller Transform
{
    return sqrt(-2*log(nrand())) * cos(2*pi*nrand());
}


/****************************************************************/
/** Position Initialization Module **
 ************************************/




/****************************************************************/
/** Output Module **
 *******************/



/****************************************************************/
/** Force & Hydrodynamics Calc Module **
 ***************************************/



/****************************************************************/
/** Collision Handler Module **
 ******************************/

