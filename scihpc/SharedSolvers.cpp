//
// Created by leo on 10/4/22.
//

#include "SharedSolvers.h"

SharedSolvers *SharedSolvers_alloc(scalar_data *f, structured_grid *geo) {
    auto solvers = new SharedSolvers;
    solvers->ccd = new ccd_solver(f, geo);
    solvers->uccd = new uccd_solver(f, geo);
    solvers->weno = new weno_solver(f, geo);
    solvers->secSol = new second_order_solver(geo->dx, geo->dy, geo->dz);
    solvers->mg = new multigrid(f, geo);
    return solvers;
}
