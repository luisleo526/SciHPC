//
// Created by leo on 10/2/22.
//

#ifndef SCIHPC_CCD_BASE_H
#define SCIHPC_CCD_BASE_H

#include "global.h"
#include "matrix_solver.h"

class ccd_base {
private:
    void init_coefficients() const;

public:
    ccd_base(int _n, DataType _h);

    DataType **a, **b, **aa, **bb, *s, *ss;
    DataType h;
    int n;
};

#endif //SCIHPC_CCD_BASE_H
