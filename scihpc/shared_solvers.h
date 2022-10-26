//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_SHARED_SOLVERS_H
#define SCIHPC_SHARED_SOLVERS_H

#include "ccd_solver.h"
#include "uccd_solver.h"
#include "weno_solver.h"
#include "second_order_solver.h"
#include "multigrid.h"

struct shared_solvers {
    ccd_solver *ccd;
    uccd_solver *uccd;
    weno_solver *weno;
    second_order_solver *secSol;
    multigrid *mg;
};

shared_solvers *SharedSolvers_alloc(scalar_data *f, structured_grid *geo);

void shared_solvers_mg_init_Neumann(shared_solvers *solvers);

void shared_solvers_mg_init_Dirichlet(shared_solvers *solvers);

#endif //SCIHPC_SHARED_SOLVERS_H
