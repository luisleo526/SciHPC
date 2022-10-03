//
// Created by 溫晧良 on 2022/10/1.
//

#include "lsm.h"

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


DataType Sign(const DataType& x, const DataType& ls_width) {

    if (x < -ls_width) {
        return -1;
    } else if (x > ls_width) {
        return 1;
    } else {
        return x / sqrt(x * x + ls_width * ls_width);
    }

}

DataType lsf_mass(wrapper *f) {
    DataType mass = 0.0;

#pragma omp parallel for default(none) shared(f) reduction(+:mass) collapse(3)
    for (int i = 0; i < f->scalar->nx; ++i) {
        for (int j = 0; j < f->scalar->ny; ++j) {
            for (int k = 0; k < f->scalar->nz; ++k) {
                auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                auto h = Heaviside(f->scalar->data[index.i][index.j][index.k], f->params->ls_width);
                auto rho = h + (1.0-h)*f->params->density_ratio;
                mass += rho * h;
            }
        }
    }
    return mass;
}

DataType lsf_volume(wrapper *f) {
    DataType volume = 0;

#pragma omp parallel for default(none) shared(f) reduction(+:volume) collapse(3)
    for (int i = 0; i < f->scalar->nx; ++i) {
        for (int j = 0; j < f->scalar->ny; ++j) {
            for (int k = 0; k < f->scalar->nz; ++k) {
                auto index = f->scalar->index_mapping(i + 1, j + 1, k + 1);
                auto h = Heaviside(f->scalar->data[index.i][index.j][index.k], f->params->ls_width);
                volume += h;
            }
        }
    }
    return volume;
}








