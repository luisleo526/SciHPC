//
// Created by 溫晧良 on 2022/10/1.
//

#include "lsm.h"
#include "wrapper_func.h"

DataType Heaviside(DataType x, DataType ls_width) {

    if (x < -ls_width) {
        return 0;
    } else if (x > ls_width) {
        return 1;
    } else {
        return 0.5 * (1 + x / ls_width + sin(pi * x / ls_width) / pi);
    }

}

DataType Delta(DataType x, DataType ls_width) {

    if (x < -ls_width) {
        return 0;
    } else if (x > ls_width) {
        return 0;
    } else {
        return 0.5 * (1 + cos(pi * x / ls_width)) / ls_width;
    }

}


DataType Sign(const DataType &x, const DataType &ls_width) {

    if (x < -ls_width) {
        return -1;
    } else if (x > ls_width) {
        return 1;
    } else {
        return x / sqrt(x * x + ls_width * ls_width);
    }

}








