//
// Created by leo on 10/4/22.
//

#ifndef SCIHPC_GODUNOV_GRADIENT_H
#define SCIHPC_GODUNOV_GRADIENT_H

#include "global.h"
#include "wrapper.h"
#include "structured_grid.h"
#include <cmath>

DataType godunov_limiter_p(DataType fp, DataType fm);
DataType godunov_limiter_m(DataType fp, DataType fm);

void godunov_gradient(wrapper* f, structured_grid* geo);
void stabilized_upon_gradient(wrapper* f, structured_grid* geo);

#endif //SCIHPC_GODUNOV_GRADIENT_H
