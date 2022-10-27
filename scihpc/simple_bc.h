//
// Created by 溫晧良 on 2022/10/27.
//

#ifndef SCIHPC_SIMPLE_BC_H
#define SCIHPC_SIMPLE_BC_H

#include "global.h"

void zero_order_extrapolation(DataType*** f, int nx, int ny, int nz, int ghc);
void first_order_extrapolation(DataType*** f, int nx, int ny, int nz, int ghc);


#endif //SCIHPC_SIMPLE_BC_H
