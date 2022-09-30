//
// Created by 溫晧良 on 2022/9/30.
//

#ifndef SCIHPC_CCD_H
#define SCIHPC_CCD_H

#include "global.h"
#include "matrix_solver.h"
#include "scalar_data.h"
#include "structured_grid.h"
#include "vector_data.h"

void ccd_coefficient_boundary_condition(DataType ***coeff, int n, DataType dx);

DataType ***ccd_coefficient_matrix(int n, DataType h);
DataType ****uccd_coefficient_matrix(int n, DataType h);

void ccd_find_fx(scalar_data *f, structured_grid *geo);
void ccd_find_fy(scalar_data *f, structured_grid *geo);
void ccd_find_fz(scalar_data *f, structured_grid *geo);
void ccd_find_derivatives(scalar_data *f, structured_grid *geo);

void uccd_find_fx(scalar_data *f, structured_grid *geo, vector_data *nvel);
void uccd_find_fy(scalar_data *f, structured_grid *geo, vector_data *nvel);
void uccd_find_fz(scalar_data *f, structured_grid *geo, vector_data *nvel);

#endif //SCIHPC_CCD_H
