//
// Created by 溫晧良 on 2022/10/5.
//

#ifndef SCIHPC_PROJECTION_METHOD_H
#define SCIHPC_PROJECTION_METHOD_H

#include "wrapper_func.h"
#include "flux.h"
#include "source.h"
#include "structured_grid.h"

class projection_method {
public:
    explicit projection_method(scalar_data *f);

    DataType ***u_src, ***v_src, ***w_src;
    DataType ***u_src_old, ***v_src_old, ***w_src_old;
    DataType ***CC, ***CR, ***CL, ***CU, ***CD, ***CF, ***CB, ***RHS;

    void add_stress_x(wrapper *vel, wrapper *lsf, wrapper *nvel) const;

    void add_stress_y(wrapper *vel, wrapper *lsf, wrapper *nvel) const;

    void add_stress_z(wrapper *vel, wrapper *lsf, wrapper *nvel) const;

    void add_body_force(wrapper *lsf) const;

    void add_surface_force(wrapper *lsf) const;

    void find_source(wrapper *vel, wrapper *nvel, wrapper *lsf) const;

    void find_source_sec(wrapper *vel, wrapper *nvel, wrapper *lsf) const;

    void find_intermediate_velocity(wrapper *vel) const;

    void solve_ppe(wrapper *pressure, wrapper *lsf, wrapper *vel) const;

    static void find_final_velocity(wrapper *vel, wrapper *pressure, wrapper *lsf);

    void ab_solve(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf) const;

    void ab_solve_sec(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf) const;
};


#endif //SCIHPC_PROJECTION_METHOD_H
