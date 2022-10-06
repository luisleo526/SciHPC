//
// Created by 溫晧良 on 2022/10/5.
//

#include "no_slip.h"


void no_slip_node_xr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int j = 0; j < vel->ny; j++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(vel->nx, j + 1, k + 1);
            for (int i = 0; i < vel->ghc; ++i) {
                vel->data[index.i + 1 + i][index.j][index.k] = -vel->data[index.i - i][index.j][index.k];
            }
        }
    }
}

void no_slip_face_xr(scalar_data *u) {
#pragma omp parallel for default(none) shared(u) collapse(2)
    for (int j = 0; j < u->ny; j++) {
        for (int k = 0; k < u->nz; k++) {
            auto index = u->index_mapping(u->nx, j + 1, k + 1);
            for (int i = 0; i < u->ghc; ++i) {
                u->data[index.i + 1 + i][index.j][index.k] = -u->data[index.i - 1 - i][index.j][index.k];
            }
            u->data[index.i][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_node_xl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int j = 0; j < vel->ny; j++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(1, j + 1, k + 1);
            for (int i = 0; i < vel->ghc; ++i) {
                vel->data[index.i - 1 - i][index.j][index.k] = -vel->data[index.i + i][index.j][index.k];
            }
        }
    }
}

void no_slip_face_xl(scalar_data *u) {
#pragma omp parallel for default(none) shared(u) collapse(2)
    for (int j = 0; j < u->ny; j++) {
        for (int k = 0; k < u->nz; k++) {
            auto index = u->index_mapping(1, j + 1, k + 1);
            for (int i = 0; i < u->ghc - 1; ++i) {
                u->data[index.i - 2 - i][index.j][index.k] = -u->data[index.i + i][index.j][index.k];
            }
            u->data[index.i - 1][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_node_yr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(i + 1, vel->ny, k + 1);
            for (int j = 0; j < vel->ghc; ++j) {
                vel->data[index.i][index.j + 1 + j][index.k] = -vel->data[index.i][index.j - j][index.k];
            }
        }
    }
}

void no_slip_face_yr(scalar_data *v) {
#pragma omp parallel for default(none) shared(v) collapse(2)
    for (int i = 0; i < v->nx; i++) {
        for (int k = 0; k < v->nz; k++) {
            auto index = v->index_mapping(i + 1, v->ny, k + 1);
            for (int j = 0; j < v->ghc; ++j) {
                v->data[index.i][index.j + 1 + j][index.k] = -v->data[index.i][index.j - 1 - j][index.k];
            }
            v->data[index.i][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_node_yl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(i + 1, 1, k + 1);
            for (int j = 0; j < vel->ghc; ++j) {
                vel->data[index.i][index.j - 1 - j][index.k] = -vel->data[index.i][index.j + j][index.k];
            }
        }
    }
}

void no_slip_face_yl(scalar_data *v) {
#pragma omp parallel for default(none) shared(v) collapse(2)
    for (int i = 0; i < v->nx; i++) {
        for (int k = 0; k < v->nz; k++) {
            auto index = v->index_mapping(i + 1, 1, k + 1);
            for (int j = 0; j < v->ghc - 1; ++j) {
                v->data[index.i][index.j - 2 - j][index.k] = -v->data[index.i][index.j + j][index.k];
            }
            v->data[index.i][index.j - 1][index.k] = 0.0;
        }
    }
}

void no_slip_node_zr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int j = 0; j < vel->ny; j++) {
            auto index = vel->index_mapping(i + 1, j + 1, vel->nz);
            for (int k = 0; k < vel->ghc; ++k) {
                vel->data[index.i][index.j][index.k + 1 + k] = -vel->data[index.i][index.j][index.k - k];
            }
        }
    }
}

void no_slip_face_zr(scalar_data *w) {
#pragma omp parallel for default(none) shared(w) collapse(2)
    for (int i = 0; i < w->nx; i++) {
        for (int j = 0; j < w->ny; j++) {
            auto index = w->index_mapping(i + 1, j + 1, w->nz);
            for (int k = 0; k < w->ghc; ++k) {
                w->data[index.i][index.j][index.k + 1 + k] = -w->data[index.i][index.j][index.k - 1 - k];
            }
            w->data[index.i][index.j][index.k] = 0.0;
        }
    }
}

void no_slip_node_zl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int j = 0; j < vel->ny; j++) {
            auto index = vel->index_mapping(i + 1, j + 1, 1);
            for (int k = 0; k < vel->ghc; ++k) {
                vel->data[index.i][index.j][index.k - 1 - k] = -vel->data[index.i][index.j][index.k + k];
            }
        }
    }
}

void no_slip_face_zl(scalar_data *w) {
#pragma omp parallel for default(none) shared(w) collapse(2)
    for (int i = 0; i < w->nx; i++) {
        for (int j = 0; j < w->ny; j++) {
            auto index = w->index_mapping(i + 1, j + 1, 1);
            for (int k = 0; k < w->ghc - 1; ++k) {
                w->data[index.i][index.j][index.k - 2 - k] = -w->data[index.i][index.j][index.k + k];
            }
            w->data[index.i][index.j][index.k - 1] = 0.0;
        }
    }
}

void no_slip_node_xr(vector_data *vel) {
    no_slip_node_xr(&vel->x);
    if (vel->x.ndim > 1) {
        no_slip_node_xr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_xr(&vel->z);
    }
}

void no_slip_node_xl(vector_data *vel) {
    no_slip_node_xl(&vel->x);
    if (vel->x.ndim > 1) {
        no_slip_node_xl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_xl(&vel->z);
    }
}

void no_slip_node_yr(vector_data *vel) {
    if (vel->x.ndim > 1) {
        no_slip_node_yr(&vel->x);
        no_slip_node_yr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_yr(&vel->z);
    }
}

void no_slip_node_yl(vector_data *vel) {
    if (vel->x.ndim > 1) {
        no_slip_node_yl(&vel->x);
        no_slip_node_yl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_yl(&vel->z);
    }
}

void no_slip_node_zr(vector_data *vel) {
    if (vel->x.ndim > 2) {
        no_slip_node_zr(&vel->x);
        no_slip_node_zr(&vel->y);
        no_slip_node_zr(&vel->z);
    }
}

void no_slip_node_zl(vector_data *vel) {
    if (vel->x.ndim > 2) {
        no_slip_node_zl(&vel->x);
        no_slip_node_zl(&vel->y);
        no_slip_node_zl(&vel->z);
    }
}

void no_slip_face_xr(vector_data *vel) {
    no_slip_face_xr(&vel->x);
    if (vel->x.ndim > 1) {
        no_slip_node_xr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_xr(&vel->z);
    }
}

void no_slip_face_xl(vector_data *vel) {
    no_slip_face_xl(&vel->x);
    if (vel->x.ndim > 1) {
        no_slip_node_xl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_xl(&vel->z);
    }
}

void no_slip_face_yr(vector_data *vel) {
    if (vel->x.ndim > 1) {
        no_slip_node_yr(&vel->x);
        no_slip_face_yr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_yr(&vel->z);
    }
}

void no_slip_face_yl(vector_data *vel) {
    if (vel->x.ndim > 1) {
        no_slip_node_yl(&vel->x);
        no_slip_face_yl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        no_slip_node_yl(&vel->z);
    }
}

void no_slip_face_zr(vector_data *vel) {
    if (vel->x.ndim > 2) {
        no_slip_node_zr(&vel->x);
        no_slip_node_zr(&vel->y);
        no_slip_face_zr(&vel->z);
    }
}

void no_slip_face_zl(vector_data *vel) {
    if (vel->x.ndim > 2) {
        no_slip_node_zl(&vel->x);
        no_slip_node_zl(&vel->y);
        no_slip_face_zl(&vel->z);
    }
}

void no_slip_node(vector_data *nvel) {
    no_slip_node_xr(nvel);
    no_slip_node_xl(nvel);
    no_slip_node_yr(nvel);
    no_slip_node_yl(nvel);
    no_slip_node_zr(nvel);
    no_slip_node_zl(nvel);
}

void no_slip_face(vector_data *vel) {
    no_slip_face_xr(vel);
    no_slip_face_xl(vel);
    no_slip_face_yr(vel);
    no_slip_face_yl(vel);
    no_slip_face_zr(vel);
    no_slip_face_zl(vel);
}

void no_slip_face_x(vector_data *vel) {
    if (vel->x.ndim == 1) {
        no_slip_face_xr(&vel->x);
        no_slip_face_xl(&vel->x);
    } else if (vel->x.ndim == 2) {
        no_slip_face_xr(&vel->x);
        no_slip_face_xr(&vel->y);
        no_slip_face_xl(&vel->x);
        no_slip_face_xl(&vel->y);
        no_slip_node_yr(&vel->x);
        no_slip_node_yr(&vel->y);
        no_slip_node_yl(&vel->x);
        no_slip_node_yl(&vel->y);
    } else if (vel->x.ndim == 3) {
        no_slip_face_xr(&vel->x);
        no_slip_face_xr(&vel->y);
        no_slip_face_xr(&vel->z);
        no_slip_face_xl(&vel->x);
        no_slip_face_xl(&vel->y);
        no_slip_face_xl(&vel->z);
        no_slip_node_yr(&vel->x);
        no_slip_node_yr(&vel->y);
        no_slip_node_yr(&vel->z);
        no_slip_node_yl(&vel->x);
        no_slip_node_yl(&vel->y);
        no_slip_node_yl(&vel->z);
        no_slip_node_zr(&vel->x);
        no_slip_node_zr(&vel->y);
        no_slip_node_zr(&vel->z);
        no_slip_node_zl(&vel->x);
        no_slip_node_zl(&vel->y);
        no_slip_node_zl(&vel->z);
    }
}

void no_slip_face_y(vector_data *vel) {
    if (vel->x.ndim == 2) {
        no_slip_node_xr(&vel->x);
        no_slip_node_xr(&vel->y);
        no_slip_node_xl(&vel->x);
        no_slip_node_xl(&vel->y);
        no_slip_face_yr(&vel->x);
        no_slip_face_yr(&vel->y);
        no_slip_face_yl(&vel->x);
        no_slip_face_yl(&vel->y);
    } else if (vel->x.ndim == 3) {
        no_slip_node_xr(&vel->x);
        no_slip_node_xr(&vel->y);
        no_slip_node_xr(&vel->z);
        no_slip_node_xl(&vel->x);
        no_slip_node_xl(&vel->y);
        no_slip_node_xl(&vel->z);
        no_slip_face_yr(&vel->x);
        no_slip_face_yr(&vel->y);
        no_slip_face_yr(&vel->z);
        no_slip_face_yl(&vel->x);
        no_slip_face_yl(&vel->y);
        no_slip_face_yl(&vel->z);
        no_slip_node_zr(&vel->x);
        no_slip_node_zr(&vel->y);
        no_slip_node_zr(&vel->z);
        no_slip_node_zl(&vel->x);
        no_slip_node_zl(&vel->y);
        no_slip_node_zl(&vel->z);
    }
}

void no_slip_face_z(vector_data *vel) {
    if (vel->x.ndim == 3) {
        no_slip_node_xr(&vel->x);
        no_slip_node_xr(&vel->y);
        no_slip_node_xr(&vel->z);
        no_slip_node_xl(&vel->x);
        no_slip_node_xl(&vel->y);
        no_slip_node_xl(&vel->z);
        no_slip_node_yr(&vel->x);
        no_slip_node_yr(&vel->y);
        no_slip_node_yr(&vel->z);
        no_slip_node_yl(&vel->x);
        no_slip_node_yl(&vel->y);
        no_slip_node_yl(&vel->z);
        no_slip_face_zr(&vel->x);
        no_slip_face_zr(&vel->y);
        no_slip_face_zr(&vel->z);
        no_slip_face_zl(&vel->x);
        no_slip_face_zl(&vel->y);
        no_slip_face_zl(&vel->z);
    }
}








