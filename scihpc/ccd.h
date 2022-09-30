//
// Created by 溫晧良 on 2022/9/30.
//

#ifndef SCIHPC_CCD_H
#define SCIHPC_CCD_H

#include "global.h"
#include "matrix_solver.h"
#include "scalar_data.h"

void ccd_coefficient_boundary_condition(DataType ***coeff, int n, DataType dx);

DataType*** ccd_coefficient_matrix(int n, DataType dx);

DataType** ccd_src_matrix(const DataType *f, long d, int n, DataType dx);

void ccd_find_fx(scalar_data* f, DataType dx);

void ccd_find_fy(scalar_data* f, DataType dy);

void ccd_find_fz(scalar_data* f, DataType dz);

#endif //SCIHPC_CCD_H
