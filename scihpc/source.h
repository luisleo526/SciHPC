//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_SOURCE_H
#define SCIHPC_SOURCE_H

#include "global.h"
#include "structured_grid.h"
#include "derivatives_solver.h"
#include "wrapper.h"
#include "lsm.h"

void convection(wrapper *f, wrapper *vel, structured_grid *geo, DataType ***s,
                void (*flux)(scalar_data *, vector_data *));

void Hamilton_Jacobi(wrapper *f, wrapper *vel, structured_grid *geo, DataType ***s,
                     void (*flux)(scalar_data *, vector_data *));

void mpls(wrapper *phi, wrapper *vel, structured_grid *geo, DataType ***s,
          void (*flux)(scalar_data *, vector_data *));

void lsf_init(wrapper *phi, wrapper *vel, structured_grid *geo, DataType ***s,
              void (*flux)(scalar_data *, vector_data *));

#endif //SCIHPC_SOURCE_H
