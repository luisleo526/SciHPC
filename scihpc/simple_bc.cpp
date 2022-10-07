//
// Created by 溫晧良 on 2022/9/28.
//

#include "simple_bc.h"

void periodic_x(scalar_data *f) {
#pragma omp parallel for default(none) shared(f)
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
#pragma omp parallel for default(none) shared(f)
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
#pragma omp parallel for default(none) shared(f)
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

void periodic(vector_data *f) {
    periodic_x(&f->x);
    periodic_x(&f->y);
    periodic_x(&f->z);
    if (f->x.ndim > 1) {
        periodic_y(&f->x);
        periodic_y(&f->y);
        periodic_y(&f->z);
    }
    if (f->x.ndim > 2) {
        periodic_z(&f->x);
        periodic_z(&f->x);
        periodic_z(&f->x);
    }
}

void first_order_extrapolation_x(scalar_data *f) {
#pragma omp parallel for default(none) shared(f)
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
#pragma omp parallel for default(none) shared(f)
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
#pragma omp parallel for default(none) shared(f)
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

void first_order_extrapolation(vector_data *f) {
    first_order_extrapolation_x(&f->x);
    first_order_extrapolation_x(&f->y);
    first_order_extrapolation_x(&f->z);
    if (f->x.ndim > 1) {
        first_order_extrapolation_y(&f->x);
        first_order_extrapolation_y(&f->y);
        first_order_extrapolation_y(&f->z);
    }
    if (f->x.ndim > 2) {
        first_order_extrapolation_z(&f->x);
        first_order_extrapolation_z(&f->x);
        first_order_extrapolation_z(&f->x);
    }
}

void zero_order_extrapolation_x(scalar_data *f) {
#pragma omp parallel for  default(none) shared(f)
    for (int j = 0; j < f->ny; ++j) {
        for (int k = 0; k < f->nz; ++k) {
            indices index_l = f->index_mapping(1, j+1, k+1);
            indices index_r = f->index_mapping(f->nx, j+1, k+1);
            for (int i = 0; i < f->ghc; ++i) {
                f->data[index_l.i - i - 1][index_l.j][index_l.k] = f->data[index_l.i][index_l.j][index_l.k];
                f->data[index_r.i + i + 1][index_r.j][index_r.k] = f->data[index_r.i][index_r.j][index_r.k];
            }
        }
    }
}

void zero_order_extrapolation_y(scalar_data *f) {
#pragma omp parallel for default(none) shared(f)
    for (int i = 0; i < f->nx; ++i) {
        for (int k = 0; k < f->nz; ++k) {
            indices index_l = f->index_mapping(i+1, 1, k+1);
            indices index_r = f->index_mapping(i+1, f->ny, k+1);
            for (int j = 0; j < f->ghc; ++j) {
                f->data[index_l.i][index_l.j - j - 1][index_l.k] = f->data[index_l.i][index_l.j][index_l.k];
                f->data[index_r.i][index_r.j + j + 1][index_r.k] = f->data[index_r.i][index_r.j][index_r.k];
            }
        }
    }
}

void zero_order_extrapolation_z(scalar_data *f) {
#pragma omp parallel for default(none) shared(f)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            indices index_l = f->index_mapping(i+1, j+1, 1);
            indices index_r = f->index_mapping(i+1, j+1, f->nz);
            for (int k = 0; k < f->ghc; ++k) {
                f->data[index_l.i][index_l.j][index_l.k - k - 1] = f->data[index_l.i][index_l.j][index_l.k];
                f->data[index_r.i][index_r.j][index_r.k + k + 1] = f->data[index_r.i][index_r.j][index_r.k];
            }
        }
    }
}

void zero_order_extrapolation(scalar_data *f) {
    zero_order_extrapolation_x(f);
    if (f->ndim > 1) {
        zero_order_extrapolation_y(f);
    }
    if (f->ndim > 2) {
        zero_order_extrapolation_z(f);
    }
}

void zero_order_extrapolation(vector_data *f) {
    zero_order_extrapolation_x(&f->x);
    zero_order_extrapolation_x(&f->y);
    zero_order_extrapolation_x(&f->z);
    if (f->x.ndim > 1) {
        zero_order_extrapolation_y(&f->x);
        zero_order_extrapolation_y(&f->y);
        zero_order_extrapolation_y(&f->z);
    }
    if (f->x.ndim > 2) {
        zero_order_extrapolation_z(&f->x);
        zero_order_extrapolation_z(&f->x);
        zero_order_extrapolation_z(&f->x);
    }
}

