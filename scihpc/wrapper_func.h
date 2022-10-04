//
// Created by leo on 10/4/22.
//

#ifndef SCIHPC_WRAPPER_FUNC_H
#define SCIHPC_WRAPPER_FUNC_H

#include "wrapper.h"
#include "lsm.h"

void find_heavyside(wrapper* lsf);
void find_delta(wrapper* lsf);
void find_sign(wrapper* lsf);
void find_gradient(wrapper* lsf);
void store(wrapper* f);
DataType l2norm(wrapper* f);
void integrate_a(wrapper* f);
void integrate_b(wrapper* f);
#endif //SCIHPC_WRAPPER_FUNC_H
