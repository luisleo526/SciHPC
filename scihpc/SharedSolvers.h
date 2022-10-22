//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_SHAREDSOLVERS_H
#define SCIHPC_SHAREDSOLVERS_H

#include "ccd_solver.h"
#include "uccd_solver.h"
#include "weno_solver.h"
#include "second_order_solver.h"
#include "multigrid.h"

struct SharedSolvers {
    ccd_solver *ccd;
    uccd_solver *uccd;
    weno_solver *weno;
    second_order_solver *secSol;
    multigrid *mg;
};

SharedSolvers *SharedSolvers_alloc(scalar_data *f, structured_grid *geo);

#endif //SCIHPC_SHAREDSOLVERS_H
