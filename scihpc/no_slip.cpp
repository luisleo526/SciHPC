//
// Created by 溫晧良 on 2022/10/5.
//

#include "no_slip.h"

void no_slip_node_xl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int j = 0; j < vel->ny; j++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= vel->ghc; ++i) {
                vel->data[index.i - i][index.j][index.k] = -vel->data[index.i - 1 + i][index.j][index.k];
            }
        }
    }
}

void no_slip_node_xr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int j = 0; j < vel->ny; j++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(vel->nx, j + 1, k + 1);
            for (int i = 1; i <= vel->ghc; ++i) {
                vel->data[index.i + i][index.j][index.k] = -vel->data[index.i + 1 - i][index.j][index.k];
            }
        }
    }
}

void no_slip_node_yl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= vel->ghc; ++j) {
                vel->data[index.i][index.j - j][index.k] = -vel->data[index.i][index.j - 1 + j][index.k];
            }
        }
    }
}

void no_slip_node_yr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(i + 1, vel->ny, k + 1);
            for (int j = 1; j <= vel->ghc; ++j) {
                vel->data[index.i][index.j + j][index.k] = -vel->data[index.i][index.j + 1 - j][index.k];
            }
        }
    }
}

void no_slip_node_zl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int j = 0; j < vel->ny; j++) {
            auto index = vel->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= vel->ghc; ++k) {
                vel->data[index.i][index.j][index.k - k] = -vel->data[index.i][index.j][index.k - 1 + k];
            }
        }
    }
}

void no_slip_node_zr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int j = 0; j < vel->ny; j++) {
            auto index = vel->index_mapping(i + 1, j + 1, vel->nz);
            for (int k = 1; k <= vel->ghc; ++k) {
                vel->data[index.i][index.j][index.k + k] = -vel->data[index.i][index.j][index.k + 1 - k];
            }
        }
    }
}

void no_slip_face_xl(scalar_data *u) {
#pragma omp parallel for default(none) shared(u) collapse(2)
    for (int j = 0; j < u->ny; j++) {
        for (int k = 0; k < u->nz; k++) {
            auto index = u->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= u->ghc; ++i) {
                u->data[index.i - i][index.j][index.k] = -u->data[index.i - 2 + i][index.j][index.k];
            }
            u->data[index.i - 1][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_face_xr(scalar_data *u) {
#pragma omp parallel for default(none) shared(u) collapse(2)
    for (int j = 0; j < u->ny; j++) {
        for (int k = 0; k < u->nz; k++) {
            auto index = u->index_mapping(u->nx, j + 1, k + 1);
            for (int i = 1; i <= u->ghc; ++i) {
                u->data[index.i + i][index.j][index.k] = -u->data[index.i - i][index.j][index.k];
            }
            u->data[index.i][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_face_yl(scalar_data *v) {
#pragma omp parallel for default(none) shared(v) collapse(2)
    for (int i = 0; i < v->nx; i++) {
        for (int k = 0; k < v->nz; k++) {
            auto index = v->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= v->ghc; ++j) {
                v->data[index.i][index.j - j][index.k] = -v->data[index.i][index.j - 2 + j][index.k];
            }
            v->data[index.i][index.j - 1][index.k] = 0.0;
        }
    }
}

void no_slip_face_yr(scalar_data *v) {
#pragma omp parallel for default(none) shared(v) collapse(2)
    for (int i = 0; i < v->nx; i++) {
        for (int k = 0; k < v->nz; k++) {
            auto index = v->index_mapping(i + 1, v->ny, k + 1);
            for (int j = 1; j <= v->ghc; ++j) {
                v->data[index.i][index.j + j][index.k] = -v->data[index.i][index.j - j][index.k];
            }
            v->data[index.i][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_face_zl(scalar_data *w) {
#pragma omp parallel for default(none) shared(w) collapse(2)
    for (int i = 0; i < w->nx; i++) {
        for (int j = 0; j < w->ny; j++) {
            auto index = w->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= w->ghc; ++k) {
                w->data[index.i][index.j][index.k - k] = -w->data[index.i][index.j][index.k - 2 + k];
            }
            w->data[index.i][index.j][index.k - 1] = 0.0;
        }
    }
}

void no_slip_face_zr(scalar_data *w) {
#pragma omp parallel for default(none) shared(w) collapse(2)
    for (int i = 0; i < w->nx; i++) {
        for (int j = 0; j < w->ny; j++) {
            auto index = w->index_mapping(i + 1, j + 1, w->nz);
            for (int k = 1; k <= w->ghc; ++k) {
                w->data[index.i][index.j][index.k + k] = -w->data[index.i][index.j][index.k - k];
            }
            w->data[index.i][index.j][index.k] = 0.0;
        }
    }
}







