//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_UCCD_SOLVER_H
#define SCIHPC_UCCD_SOLVER_H

#include "global.h"
#include "matrix_solver.h"
#include "uccd_base.h"
#include "scalar_data.h"
#include "vector_data.h"
#include "structured_grid.h"

class uccd_solver {
public:
    uccd_solver(int _nx, int _ny, int _nz, DataType _dx, DataType _dy, DataType _dz);
    uccd_solver(scalar_data* f, structured_grid* geo);
    uccd_base *x, *y, *z;

    void find_fx(scalar_data *f, vector_data *vel) const;
    void find_fy(scalar_data *f, vector_data *vel) const;
    void find_fz(scalar_data *f, vector_data *vel) const;
    void find_derivatives(scalar_data *f, vector_data *vel) const;
};


#endif //SCIHPC_UCCD_SOLVER_H
