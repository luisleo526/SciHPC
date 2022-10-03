//
// Created by leo on 10/4/22.
//

#include "derivatives_solver.h"

derivatives_solver *derivatives_solver_alloc(scalar_data *f, structured_grid *geo) {
    auto solvers = new derivatives_solver;
    solvers->ccd = new ccd_solver(f, geo);
    solvers->uccd = new uccd_solver(f, geo);
    solvers->weno = new weno_solver(f, geo);
    return solvers;
}
