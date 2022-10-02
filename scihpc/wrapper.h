//
// Created by leo on 10/3/22.
//

#ifndef SCIHPC_WRAPPER_H
#define SCIHPC_WRAPPER_H

#include "global.h"
#include "derivatives_solver.h"
#include "scalar_data.h"
#include "vector_data.h"

class wrapper {
public:
    scalar_data* scalar;
    vector_data* vector;
    bool is_scalar;
    explicit wrapper(scalar_data* _scalar);
    explicit wrapper(vector_data* _vector);
    solvers_ptr* solvers;
    problem_parameters* params;

    void link_solvers(solvers_ptr* _solvers);
    void link_params(problem_parameters* _params);
};


#endif //SCIHPC_WRAPPER_H
