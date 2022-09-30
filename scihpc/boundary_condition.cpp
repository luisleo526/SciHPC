//
// Created by 溫晧良 on 2022/9/28.
//

#include "boundary_condition.h"

void periodic_x(scalar_data *f) {
    for (int j = 1; j < f->ny + 1; ++j) {
        for (int k = 1; k < f->nz + 1; ++k) {
            indices index_l = f->index_mapping(1, j, k);
            indices index_r = f->index_mapping(f->nx, j, k);
            for (int i = 0; i < f->ghc; ++i) {
                f->data[index_l.i - i - 1][index_l.j][index_l.k] = f->data[index_r.i - i][index_r.j][index_r.k];
                f->data[index_r.i + i + 1][index_r.j][index_r.k] = f->data[index_l.i + i][index_l.j][index_l.k];
            }
        }
    }
}

void periodic_y(scalar_data *f) {
    for (int i = 1; i < f->nx + 1; ++i) {
        for (int k = 1; k < f->nz + 1; ++k) {
            indices index_l = f->index_mapping(i, 1, k);
            indices index_r = f->index_mapping(i, f->ny, k);
            for (int j = 0; j < f->ghc; ++j) {
                f->data[index_l.i][index_l.j - j - 1][index_l.k] = f->data[index_r.i][index_r.j - j][index_r.k];
                f->data[index_r.i][index_r.j + j + 1][index_r.k] = f->data[index_l.i][index_l.j + j][index_l.k];
            }
        }
    }
}

void periodic_z(scalar_data *f) {
    for (int i = 1; i < f->nx + 1; ++i) {
        for (int j = 1; j < f->ny + 1; ++j) {
            indices index_l = f->index_mapping(i, j, 1);
            indices index_r = f->index_mapping(i, j, f->nz);
            for (int k = 0; k < f->ghc; ++k) {
                f->data[index_l.i][index_l.j][index_l.k - k - 1] = f->data[index_r.i][index_r.j][index_r.k - k];
                f->data[index_r.i][index_r.j][index_r.k + k + 1] = f->data[index_l.i][index_l.j][index_l.k + k];
            }
        }
    }
}

void periodic(scalar_data *f) {
    periodic_x(f);
    if (f->ndim > 1) {
        periodic_y(f);
    }
    if (f->ndim > 2) {
        periodic_z(f);
    }
}

void first_order_extrapolation_x(scalar_data *f) {
    for (int j = 1; j < f->ny + 1; ++j) {
        for (int k = 1; k < f->nz + 1; ++k) {
            indices index_l = f->index_mapping(1, j, k);
            indices index_r = f->index_mapping(f->nx, j, k);
            for (int i = 0; i < f->ghc; ++i) {
                f->data[index_l.i - i - 1][index_l.j][index_l.k] =
                        2.0 * f->data[index_l.i - i + 1][index_l.j][index_l.k] -
                        f->data[index_l.i - i][index_l.j][index_l.k];
                f->data[index_r.i + i + 1][index_r.j][index_r.k] =
                        2.0 * f->data[index_r.i + i - 1][index_r.j][index_r.k] -
                        f->data[index_r.i + i][index_r.j][index_r.k];
            }
        }
    }
}

void first_order_extrapolation_y(scalar_data *f) {
    for (int i = 1; i < f->nx + 1; ++i) {
        for (int k = 1; k < f->nz + 1; ++k) {
            indices index_l = f->index_mapping(i, 1, k);
            indices index_r = f->index_mapping(i, f->ny, k);
            for (int j = 0; j < f->ghc; ++j) {
                f->data[index_l.i][index_l.j - j - 1][index_l.k] =
                        2.0 * f->data[index_l.i][index_l.j - j + 1][index_l.k] -
                        f->data[index_l.i][index_l.j - j][index_l.k];
                f->data[index_r.i][index_r.j + j + 1][index_r.k] =
                        2.0 * f->data[index_r.i][index_r.j + j - 1][index_r.k] -
                        f->data[index_r.i][index_r.j + j][index_r.k];
            }
        }
    }
}

void first_order_extrapolation_z(scalar_data *f) {
    for (int i = 1; i < f->nx + 1; ++i) {
        for (int j = 1; j < f->ny + 1; ++j) {
            indices index_l = f->index_mapping(i, j, 1);
            indices index_r = f->index_mapping(i, j, f->nz);
            for (int k = 0; k < f->ghc; ++k) {
                f->data[index_l.i][index_l.j][index_l.k - k - 1] =
                        2.0 * f->data[index_l.i][index_l.j][index_l.k - k + 1] -
                        f->data[index_l.i][index_l.j][index_l.k - k];
                f->data[index_r.i][index_r.j][index_r.k + k + 1] =
                        2.0 * f->data[index_r.i][index_r.j][index_r.k + k - 1] -
                        f->data[index_r.i][index_r.j][index_r.k + k];
            }
        }
    }
}

void first_order_extrapolation(scalar_data *f) {
    first_order_extrapolation_x(f);
    if (f->ndim > 1) {
        first_order_extrapolation_y(f);
    }
    if (f->ndim > 2) {
        first_order_extrapolation_z(f);
    }
}
