//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_DERIVATIVES_SOLVER_H
#define SCIHPC_DERIVATIVES_SOLVER_H

#include "ccd_solver.h"
#include "uccd_solver.h"

struct solvers_ptr{
    ccd_solver *ccd;
    uccd_solver *uccd;
};

#endif //SCIHPC_DERIVATIVES_SOLVER_H
