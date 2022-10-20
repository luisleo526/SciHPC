//
// Created by leo on 10/3/22.
//

#include "wrapper.h"
#include <cassert>
#include <utility>

void wrapper::link_solvers(SharedSolvers *_solvers) {
    solvers = _solvers;
}

void wrapper::link_params(problem_parameters *_params) {
    params = _params;
}

void wrapper::link_dummy(dummy_data *_dummy) {
    dummy = _dummy;
}

void wrapper::apply_vel_bc() const {
    assert(!is_scalar);
    bcFactoryU->apply_node_staggered_x(&vector->x);
    bcFactoryV->apply_node_staggered_y(&vector->y);
    bcFactoryW->apply_node_staggered_z(&vector->z);
}

void wrapper::apply_nvel_bc() const {
    assert(!is_scalar);
    bcFactoryU->apply_node_centered(&vector->x);
    bcFactoryV->apply_node_centered(&vector->y);
    bcFactoryW->apply_node_centered(&vector->z);
}

void wrapper::apply_scalar_bc() const {
    assert(is_scalar);
    bcFactory->apply_node_centered(scalar);
}

void wrapper::apply_vel_x_bc() const {
    assert(!is_scalar);
    bcFactoryU->apply_node_staggered_x(&vector->x);
    bcFactoryV->apply_node_staggered_x(&vector->y);
    bcFactoryW->apply_node_staggered_x(&vector->z);
}

void wrapper::apply_vel_y_bc() const {
    assert(!is_scalar);
    bcFactoryU->apply_node_staggered_y(&vector->x);
    bcFactoryV->apply_node_staggered_y(&vector->y);
    bcFactoryW->apply_node_staggered_y(&vector->z);
}

void wrapper::apply_vel_z_bc() const {
    assert(!is_scalar);
    bcFactoryU->apply_node_staggered_z(&vector->x);
    bcFactoryV->apply_node_staggered_z(&vector->y);
    bcFactoryW->apply_node_staggered_z(&vector->z);
}

wrapper::wrapper(bool as_scalar, structured_grid *_geo,
                 const bc_info &xlbc_type, const bc_info &xrbc_type) {
    is_scalar = as_scalar;
    geo = _geo;
    if (is_scalar) {
        scalar = new scalar_data(geo->x_axis.n);
        vector = nullptr;
        bcFactory = new bc_factory(geo, xlbc_type, xrbc_type);
    } else {
        vector = new vector_data(geo->x_axis.n);
        scalar = nullptr;
        bcFactoryU = new bc_factory(geo, xlbc_type, xrbc_type);
        bcFactoryV = new bc_factory(geo, xlbc_type, xrbc_type);
        bcFactoryW = new bc_factory(geo, xlbc_type, xrbc_type);
    }

    if (xlbc_type.type == SLIP) {
        bcFactoryU->xlbc.type = NO_SLIP;
    }

    if (xrbc_type.type == SLIP) {
        bcFactoryU->xrbc.type = NO_SLIP;
    }
}

wrapper::wrapper(bool as_scalar, structured_grid *_geo,
                 const bc_info &xlbc_type, const bc_info &xrbc_type,
                 const bc_info &ylbc_type, const bc_info &yrbc_type) {
    is_scalar = as_scalar;
    geo = _geo;
    if (is_scalar) {
        scalar = new scalar_data(geo->x_axis.n, geo->y_axis.n);
        vector = nullptr;
        bcFactory = new bc_factory(geo,
                                   xlbc_type, xrbc_type,
                                   ylbc_type, yrbc_type);
    } else {
        vector = new vector_data(geo->x_axis.n, geo->y_axis.n);
        scalar = nullptr;
        bcFactoryU = new bc_factory(geo,
                                    xlbc_type, xrbc_type,
                                    ylbc_type, yrbc_type);
        bcFactoryV = new bc_factory(geo,
                                    xlbc_type, xrbc_type,
                                    ylbc_type, yrbc_type);
        bcFactoryW = new bc_factory(geo,
                                    xlbc_type, xrbc_type,
                                    ylbc_type, yrbc_type);
    }

    if (xlbc_type.type == SLIP) {
        bcFactoryU->xlbc.type = NO_SLIP;
    }

    if (xrbc_type.type == SLIP) {
        bcFactoryU->xrbc.type = NO_SLIP;
    }

    if (ylbc_type.type == SLIP) {
        bcFactoryV->ylbc.type = NO_SLIP;
    }

    if (yrbc_type.type == SLIP) {
        bcFactoryV->yrbc.type = NO_SLIP;
    }

}

wrapper::wrapper(bool as_scalar, structured_grid *_geo,
                 const bc_info &xlbc_type, const bc_info &xrbc_type,
                 const bc_info &ylbc_type, const bc_info &yrbc_type,
                 const bc_info &zlbc_type, const bc_info &zrbc_type) {
    is_scalar = as_scalar;
    geo = _geo;
    if (is_scalar) {
        scalar = new scalar_data(geo->x_axis.n, geo->y_axis.n, geo->z_axis.n);
        vector = nullptr;
        bcFactory = new bc_factory(geo,
                                   xlbc_type, xrbc_type,
                                   ylbc_type, yrbc_type,
                                   zlbc_type, zrbc_type);
    } else {
        vector = new vector_data(geo->x_axis.n, geo->y_axis.n, geo->z_axis.n);
        scalar = nullptr;
        bcFactoryU = new bc_factory(geo,
                                    xlbc_type, xrbc_type,
                                    ylbc_type, yrbc_type,
                                    zlbc_type, zrbc_type);
        bcFactoryV = new bc_factory(geo,
                                    xlbc_type, xrbc_type,
                                    ylbc_type, yrbc_type,
                                    zlbc_type, zrbc_type);
        bcFactoryW = new bc_factory(geo,
                                    xlbc_type, xrbc_type,
                                    ylbc_type, yrbc_type,
                                    zlbc_type, zrbc_type);
    }

    if (xlbc_type.type == SLIP) {
        bcFactoryU->xlbc.type = NO_SLIP;
    }

    if (xrbc_type.type == SLIP) {
        bcFactoryU->xrbc.type = NO_SLIP;
    }

    if (ylbc_type.type == SLIP) {
        bcFactoryV->ylbc.type = NO_SLIP;
    }

    if (yrbc_type.type == SLIP) {
        bcFactoryV->yrbc.type = NO_SLIP;
    }

    if (zlbc_type.type == SLIP) {
        bcFactoryW->zlbc.type = NO_SLIP;
    }

    if (zrbc_type.type == SLIP) {
        bcFactoryW->zrbc.type = NO_SLIP;
    }
}
