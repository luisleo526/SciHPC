//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_BOUNDARY_CONDITION_H
#define SCIHPC_BOUNDARY_CONDITION_H

#include "scalar_data.h"
#include "structured_grid.h"

void periodic_x(scalar_data *f, structured_grid *geo);

void periodic_y(scalar_data *f, structured_grid *geo);

void periodic_z(scalar_data *f, structured_grid *geo);

void periodic(scalar_data *f, structured_grid *geo);


#endif //SCIHPC_BOUNDARY_CONDITION_H
