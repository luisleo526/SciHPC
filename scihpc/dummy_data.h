//
// Created by leo on 10/3/22.
//

#ifndef SCIHPC_DUMMY_DATA_H
#define SCIHPC_DUMMY_DATA_H

#include "global.h"
#include "scalar_data.h"

struct dummy_data {
    // for lsf
    DataType ***heaviside;
    DataType ***sign;
    DataType ***delta;
    DataType ***grad;
    DataType ***curvature;
    // for l2 norm
    DataType ***tmp;
    DataType ***u_tmp, ***v_tmp, ***w_tmp;
    // for domain integration
    DataType ***a, ***b;
    DataType ***a_int, ***b_int;
};

dummy_data* dummy_data_alloc(scalar_data* f);

#endif //SCIHPC_DUMMY_DATA_H
