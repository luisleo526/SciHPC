//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_SOURCE_H
#define SCIHPC_SOURCE_H
#include "global.h"
#include "scalar_data.h"
#include "structured_grid.h"
#include "vector_data.h"
#include "ccd.h"

void convection(scalar_data *f, vector_data *vel, structured_grid *geo, DataType ***s,
                void(*flux)(scalar_data *,vector_data *pData));

#endif //SCIHPC_SOURCE_H
