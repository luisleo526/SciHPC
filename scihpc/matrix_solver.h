//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_MATRIX_SOLVER_H
#define SCIHPC_MATRIX_SOLVER_H

#include "global.h"

void tdma(DataType *a, DataType *b, DataType *c, DataType *x, int n);

void twin_dec(DataType **a, DataType **b, DataType **aa, DataType **bb, int n);

void twin_bks(DataType **a, DataType **b, DataType **aa, DataType **bb, DataType *s, DataType *ss, int n);

#endif //SCIHPC_MATRIX_SOLVER_H
