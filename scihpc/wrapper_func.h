//
// Created by leo on 10/4/22.
//

#ifndef SCIHPC_WRAPPER_FUNC_H
#define SCIHPC_WRAPPER_FUNC_H

#include "wrapper.h"
#include "lsm.h"
#include "structured_grid.h"

void find_heavyside(wrapper *lsf);

void find_delta(wrapper *lsf);

void find_sign(wrapper *lsf);

void find_gradient(wrapper *lsf);

void find_curvature(wrapper *lsf);

void find_density(wrapper *lsf);

void find_viscosity(wrapper *lsf);

void store_tmp(wrapper *f);

DataType l2norm(wrapper *f);

DataType divergence(wrapper *vel);

void all_to_face_x(wrapper *ref, wrapper *tgt);

void all_to_face_y(wrapper *ref, wrapper *tgt);

void all_to_face_z(wrapper *ref, wrapper *tgt);

void node_from_face(wrapper *ref, wrapper *tgt);

void find_dt(wrapper* vel, wrapper* lsf);

#endif //SCIHPC_WRAPPER_FUNC_H
