//
// Created by 溫晧良 on 2022/10/9.
//

#include "dirichlet.h"

void Dirichlet_node_xr(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(f->nx, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i + i][index.j][index.k] = U;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i - 1][index.j][index.k] / 3.0 + 2.0 * U / 3.0;
        }
    }
}

void Dirichlet_node_xl(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i - i][index.j][index.k] = U;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i + 1][index.j][index.k] / 3.0 + 2.0 * U / 3.0;
        }
    }
}

void Dirichlet_node_yr(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, f->ny, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j + j][index.k] = U;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j - 1][index.k] / 3.0 + 2.0 * U / 3.0;
        }
    }
}

void Dirichlet_node_yl(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j - j][index.k] = U;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j + 1][index.k] / 3.0 + 2.0 * U / 3.0;
        }
    }
}

void Dirichlet_node_zr(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, f->nz);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k + k] = U;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j][index.k - 1] / 3.0 + 2.0 * U / 3.0;
        }
    }
}

void Dirichlet_node_zl(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k - k] = U;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j][index.k + 1] / 3.0 + 2.0 * U / 3.0;
        }
    }
}

void Dirichlet_face_xr(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(f->nx, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i + i][index.j][index.k] = U;
            }
            f->data[index.i][index.j][index.k] = U;
        }
    }
}

void Dirichlet_face_xl(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i - i][index.j][index.k] = U;
            }
        }
    }
}

void Dirichlet_face_yr(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, f->ny, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j + j][index.k] = U;
            }
            f->data[index.i][index.j][index.k] = U;
        }
    }
}

void Dirichlet_face_yl(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j - j][index.k] = U;
            }
        }
    }
}

void Dirichlet_face_zr(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, f->nz);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k + k] = U;
            }
            f->data[index.i][index.j][index.k] = U;
        }
    }
}

void Dirichlet_face_zl(scalar_data *f, DataType U) {
#pragma omp parallel for default(none) shared(f, U) collapse(2)
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k - k] = U;
            }
        }
    }
}