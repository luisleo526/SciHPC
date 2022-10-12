//
// Created by leo on 10/4/22.
//

#ifndef SCIHPC_GODUNOV_GRADIENT_H
#define SCIHPC_GODUNOV_GRADIENT_H

#include "global.h"
#include "wrapper.h"
#include "structured_grid.h"

DataType weno5_for_godunov(DataType a, DataType b, DataType c, DataType d);

void godunov_gradient(wrapper *f);

void stabilized_upon_gradient(wrapper *f);

#endif //SCIHPC_GODUNOV_GRADIENT_H
