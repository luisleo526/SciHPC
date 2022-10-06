//
// Created by 溫晧良 on 2022/10/5.
//

#ifndef SCIHPC_PROJECTION_METHOD_H
#define SCIHPC_PROJECTION_METHOD_H

#include "wrapper_func.h"
#include "no_slip.h"
#include "flux.h"
#include "source.h"
#include "boundary_condition.h"

class projection_method {
public:
    explicit projection_method(scalar_data* f);
    DataType ***u_src, ***v_src, ***w_src;
    DataType ***u_src_old, ***v_src_old, ***w_src_old;
    DataType ***CC, ***CR, ***CL, ***CU, ***CD, ***CF, ***CB, ***RHS;
    void add_stress_x(wrapper* nvel, wrapper* lsf);
    void add_stress_y(wrapper* nvel, wrapper* lsf);
    void add_stress_z(wrapper* nvel, wrapper* lsf);
    void find_source(wrapper *vel, wrapper *nvel, wrapper *lsf, structured_grid *geo);
    void find_intermediate_velocity(wrapper *vel);
    void solve_ppe(wrapper *pressure, wrapper *lsf, wrapper *vel, structured_grid *geo);
    void find_final_velocity(wrapper *vel, wrapper* pressure, wrapper *lsf, structured_grid *geo);
    void solve(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf, structured_grid *geo);
};


#endif //SCIHPC_PROJECTION_METHOD_H
