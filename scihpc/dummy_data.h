//
// Created by leo on 10/3/22.
//

#ifndef SCIHPC_DUMMY_DATA_H
#define SCIHPC_DUMMY_DATA_H

#include "global.h"

struct dummy_data {
    DataType ***heaviside;
    DataType ***sign;
    DataType ***delta;
    DataType ***tmp;
    DataType ***u_tmp, ***v_tmp, ***w_tmp;
    DataType ***a, ***b, ***c;
};

dummy_data* dummy_data_alloc(int nx, int ny, int nz);

#endif //SCIHPC_DUMMY_DATA_H
