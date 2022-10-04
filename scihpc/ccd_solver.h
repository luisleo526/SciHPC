//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_CCD_SOLVER_H
#define SCIHPC_CCD_SOLVER_H

#include "global.h"
#include "matrix_solver.h"
#include "structured_grid.h"
#include "scalar_data.h"
#include "vector_data.h"
#include "ccd_base.h"

class ccd_solver {
public:
    ccd_solver(scalar_data* f, structured_grid* geo);

    ccd_base *x, *y, *z;

    void find_fx(scalar_data *f) const;
    void find_fy(scalar_data *f) const;
    void find_fz(scalar_data *f) const;

    void mixed_xy(scalar_data *f) const;
    void mixed_yz(scalar_data *f) const;
    void mixed_zx(scalar_data *f) const;

    void find_derivatives(scalar_data *f) const;
    void find_derivatives(vector_data *f) const;

    void find_derivatives_all(scalar_data *f) const;
    void find_derivatives_all(vector_data *f) const;

};


#endif //SCIHPC_CCD_SOLVER_H
