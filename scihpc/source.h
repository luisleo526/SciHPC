//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_SOURCE_H
#define SCIHPC_SOURCE_H


#include "wrapper_func.h"
#include "godunov_gradient.h"

void convection(wrapper *f, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *));
void convection_sec(wrapper *f, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *));
void Hamilton_Jacobi(wrapper *f, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *));
void mpls(wrapper *phi, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *));
void lsf_redistance_no_lambda(wrapper *phi, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *));
void lsf_redistance_lambda(wrapper *phi, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *));

#endif //SCIHPC_SOURCE_H
