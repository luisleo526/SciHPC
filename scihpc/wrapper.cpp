//
// Created by leo on 10/3/22.
//

#include "wrapper.h"
#include <cassert>
#include <utility>

void wrapper::link_solvers(derivatives_solver *_solvers) {
    solvers = _solvers;
}

void wrapper::link_params(problem_parameters *_params) {
    params = _params;
}

void wrapper::link_dummy(dummy_data *_dummy) {
    dummy = _dummy;
}

void wrapper::apply_vel_bc() {
    assert(!is_scalar);
    bcFactory->apply_node_staggered_x(&vector->x);
    bcFactory->apply_node_staggered_y(&vector->y);
    bcFactory->apply_node_staggered_z(&vector->z);
}

void wrapper::apply_nvel_bc() {
    assert(!is_scalar);
    bcFactory->apply_node_centered(&vector->x);
    bcFactory->apply_node_centered(&vector->y);
    bcFactory->apply_node_centered(&vector->z);
}

void wrapper::apply_scalar_bc() {
    assert(is_scalar);
    bcFactory->apply_node_centered(scalar);
}

void wrapper::apply_vel_x_bc() {
    assert(!is_scalar);
    bcFactory->apply_node_staggered_x(&vector->x);
    bcFactory->apply_node_staggered_x(&vector->y);
    bcFactory->apply_node_staggered_x(&vector->z);
}

void wrapper::apply_vel_y_bc() {
    assert(!is_scalar);
    bcFactory->apply_node_staggered_y(&vector->x);
    bcFactory->apply_node_staggered_y(&vector->y);
    bcFactory->apply_node_staggered_y(&vector->z);
}

void wrapper::apply_vel_z_bc() {
    assert(!is_scalar);
    bcFactory->apply_node_staggered_z(&vector->x);
    bcFactory->apply_node_staggered_z(&vector->y);
    bcFactory->apply_node_staggered_z(&vector->z);
}

wrapper::wrapper(bool as_scalar, structured_grid *_geo, bc_info xlbc_type, bc_info xrbc_type) {
    is_scalar = as_scalar;
    geo = _geo;
    if (is_scalar) {
        scalar = new scalar_data(geo->x_axis.n);
        vector = nullptr;
    } else {
        vector = new vector_data(geo->x_axis.n);
        scalar = nullptr;
    }
    bcFactory = new bc_factory(geo, xlbc_type, xrbc_type);
}

wrapper::wrapper(bool as_scalar, structured_grid *_geo, bc_info xlbc_type, bc_info xrbc_type, bc_info ylbc_type,
                 bc_info yrbc_type) {
    is_scalar = as_scalar;
    geo = _geo;
    if (is_scalar) {
        scalar = new scalar_data(geo->x_axis.n, geo->y_axis.n);
        vector = nullptr;
    } else {
        vector = new vector_data(geo->x_axis.n, geo->y_axis.n);
        scalar = nullptr;
    }
    bcFactory = new bc_factory(geo,
                               xlbc_type, xrbc_type,
                               ylbc_type, yrbc_type);

}

wrapper::wrapper(bool as_scalar, structured_grid *_geo, bc_info xlbc_type, bc_info xrbc_type, bc_info ylbc_type,
                 bc_info yrbc_type, bc_info zlbc_type, bc_info zrbc_type) {
    is_scalar = as_scalar;
    geo = _geo;
    if (is_scalar) {
        scalar = new scalar_data(geo->x_axis.n, geo->y_axis.n, geo->z_axis.n);
        vector = nullptr;
    } else {
        vector = new vector_data(geo->x_axis.n, geo->y_axis.n, geo->z_axis.n);
        scalar = nullptr;
    }
    bcFactory = new bc_factory(geo,
                               xlbc_type, xrbc_type,
                               ylbc_type, yrbc_type,
                               zlbc_type, zrbc_type);
}
