//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_BOUNDARY_CONDITION_H
#define SCIHPC_BOUNDARY_CONDITION_H

#include "scalar_data.h"
#include "structured_grid.h"
#include "vector_data.h"

void periodic_x(scalar_data *f);
void periodic_y(scalar_data *f);
void periodic_z(scalar_data *f);
void periodic(scalar_data *f);
void periodic(vector_data *f);

void zero_order_extrapolation_x(scalar_data *f);
void zero_order_extrapolation_y(scalar_data *f);
void zero_order_extrapolation_z(scalar_data *f);
void zero_order_extrapolation(scalar_data *f);
void zero_order_extrapolation(vector_data *f);


void first_order_extrapolation_x(scalar_data *f);
void first_order_extrapolation_y(scalar_data *f);
void first_order_extrapolation_z(scalar_data *f);
void first_order_extrapolation(scalar_data *f);
void first_order_extrapolation(vector_data *f);

#endif //SCIHPC_BOUNDARY_CONDITION_H
