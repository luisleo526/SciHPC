//
// Created by 溫晧良 on 2022/10/5.
//

#ifndef SCIHPC_PROJECTION_METHOD_H
#define SCIHPC_PROJECTION_METHOD_H

#include "wrapper_func.h"
#include "flux.h"
#include "source.h"
#include "structured_grid.h"
#include <Eigen/Sparse>

class projection_method {
private:
    int nx, ny, nz;
    int pos(int i, int j, int k);
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

    void projection(wrapper *pressure, wrapper *lsf, wrapper *vel) const;

    static void fast_pressure_correction(wrapper *pressure, wrapper *lsf, wrapper *vel);

    void ab_solve(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf);

    void ab_solve_sec(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf);
};


#endif //SCIHPC_PROJECTION_METHOD_H
