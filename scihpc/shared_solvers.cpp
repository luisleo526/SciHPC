//
// Created by leo on 10/4/22.
//

#include "shared_solvers.h"

shared_solvers *SharedSolvers_alloc(scalar_data *f, structured_grid *geo) {
    auto solvers = new shared_solvers;
    solvers->ccd = new ccd_solver(f, geo);
    solvers->uccd = new uccd_solver(f, geo);
    solvers->weno = new weno_solver(f, geo);
    solvers->secSol = new second_order_solver(geo->dx, geo->dy, geo->dz);
    solvers->mg = new multigrid(f, geo);
    return solvers;
}

void shared_solvers_mg_init_Neumann(shared_solvers *solvers) {
    for (int i = 0; i < solvers->mg->level_num; ++i) {
        solvers->mg->at[i]->init_NeumannBC();
    }
}

void shared_solvers_mg_init_Dirichlet(shared_solvers *solvers) {
    for (int i = 0; i < solvers->mg->level_num; ++i) {
        solvers->mg->at[i]->init_DirichletBC();
        solvers->mg->at[i]->no_compatibility = true;
    }
}

