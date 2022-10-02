//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_UCCD_BASE_H
#define SCIHPC_UCCD_BASE_H

#include "ccd_base.h"
#include "global.h"
#include "matrix_solver.h"

class uccd_base {
private:
    void init_coefficients() const;
public:
    ccd_base *upwind, *downwind;
    DataType h;
    int n;
    uccd_base(int _n, DataType _h);
};


#endif //SCIHPC_UCCD_BASE_H
