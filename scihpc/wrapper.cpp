//
// Created by leo on 10/3/22.
//

#include "wrapper.h"

wrapper::wrapper(scalar_data *_scalar) {
    scalar = _scalar;
    vector = nullptr;
    is_scalar = true;
}

wrapper::wrapper(vector_data *_vector) {
    scalar = nullptr;
    vector = _vector;
    is_scalar = false;
}

void wrapper::link_solvers(derivatives_solver *_solvers) {
    solvers = _solvers;
}

void wrapper::link_params(problem_parameters *_params) {
    params = _params;
}

void wrapper::link_dummy(dummy_data *_dummy) {
    dummy = _dummy;
}

void wrapper::link_bc(velocities_bc *_bc) {
    velbc = _bc;
}
