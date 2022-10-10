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
#include "bc_factory.h"
#include <string>

class wrapper {
public:
    bool is_scalar;
    scalar_data *scalar;
    vector_data *vector;

    wrapper(bool as_scalar, structured_grid* _geo, bc_info xlbc_type, bc_info xrbc_type);
    wrapper(bool as_scalar, structured_grid* _geo, bc_info xlbc_type, bc_info xrbc_type,
            bc_info ylbc_type, bc_info yrbc_type);
    wrapper(bool as_scalar, structured_grid* _geo, bc_info xlbc_type, bc_info xrbc_type,
            bc_info ylbc_type, bc_info yrbc_type, bc_info zlbc_type, bc_info zrbc_type);
    bc_factory *bcFactory;

    //shared data
    derivatives_solver *solvers;
    problem_parameters *params;
    dummy_data *dummy;
    structured_grid* geo;

    void link_solvers(derivatives_solver *_solvers);
    void link_params(problem_parameters *_params);
    void link_dummy(dummy_data *_dummy);

    void apply_vel_bc();
    void apply_nvel_bc();
    void apply_scalar_bc();
    void apply_vel_x_bc();
    void apply_vel_y_bc();
    void apply_vel_z_bc();


};


#endif //SCIHPC_WRAPPER_H
