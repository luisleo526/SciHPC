//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_LSM_H
#define SCIHPC_LSM_H

#include "global.h"
#include <cmath>
#include "wrapper.h"

DataType Heaviside(DataType x, DataType ls_width);
DataType Delta(DataType x, DataType ls_width);
DataType Sign(const DataType& x, const DataType& ls_width);

DataType lsf_mass(wrapper *f);
DataType lsf_volume(wrapper *f);

#endif //SCIHPC_LSM_H
