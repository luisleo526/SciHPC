//
// Created by leo on 10/3/22.
//

#ifndef SCIHPC_WRAPPER_H
#define SCIHPC_WRAPPER_H

#include "global.h"
#include "derivatives_solver.h"
#include "scalar_data.h"
#include "vector_data.h"
#include "dummy_data.h"
#include "velocities_bc.h"

class wrapper {
public:
    bool is_scalar;
    scalar_data* scalar;
    vector_data* vector;
    explicit wrapper(scalar_data* _scalar);
    explicit wrapper(vector_data* _vector);
    derivatives_solver* solvers;
    problem_parameters* params;
    dummy_data* dummy;
    velocities_bc* velbc;

    void link_solvers(derivatives_solver* _solvers);
    void link_params(problem_parameters* _params);
    void link_dummy(dummy_data* _dummy);
    void link_bc(velocities_bc* _bc);

};


#endif //SCIHPC_WRAPPER_H
