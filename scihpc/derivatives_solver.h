//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_DERIVATIVES_SOLVER_H
#define SCIHPC_DERIVATIVES_SOLVER_H

#include "ccd_solver.h"
#include "uccd_solver.h"
#include "weno_solver.h"

struct derivatives_solver{
    ccd_solver *ccd;
    uccd_solver *uccd;
    weno_solver *weno;
};

derivatives_solver* derivatives_solver_alloc(scalar_data* f, structured_grid* geo);

#endif //SCIHPC_DERIVATIVES_SOLVER_H
