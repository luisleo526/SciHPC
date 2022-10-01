//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_LSM_H
#define SCIHPC_LSM_H

#include "global.h"
#include "scalar_data.h"
#include <cmath>

DataType ***Heaviside(scalar_data* f,params* param);
DataType ***Delta(scalar_data* f,params* param);
DataType ***Sign(scalar_data* f,params* param);

DataType lsf_mass(scalar_data* f,params* param);
DataType lsf_volume(scalar_data* f,params* param);

#endif //SCIHPC_LSM_H
