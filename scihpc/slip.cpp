//
// Created by leo on 10/7/22.
//

#include "slip.h"

void slip_node_xr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int j = 0; j < vel->ny; j++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(vel->nx, j + 1, k + 1);
            for (int i = 1; i <= vel->ghc; ++i) {
                vel->data[index.i + i][index.j][index.k] = vel->data[index.i + 1 - i][index.j][index.k];
            }
        }
    }
}

void slip_node_xl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int j = 0; j < vel->ny; j++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= vel->ghc; ++i) {
                vel->data[index.i - i][index.j][index.k] = vel->data[index.i - 1 + i][index.j][index.k];
            }
        }
    }
}

void slip_node_yr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(i + 1, vel->ny, k + 1);
            for (int j = 1; j <= vel->ghc; ++j) {
                vel->data[index.i][index.j + j][index.k] = vel->data[index.i][index.j + 1 - j][index.k];
            }
        }
    }
}

void slip_node_yl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int k = 0; k < vel->nz; k++) {
            auto index = vel->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= vel->ghc; ++j) {
                vel->data[index.i][index.j - j][index.k] = vel->data[index.i][index.j - 1 + j][index.k];
            }
        }
    }
}

void slip_node_zr(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int j = 0; j < vel->ny; j++) {
            auto index = vel->index_mapping(i + 1, j + 1, vel->nz);
            for (int k = 1; k <= vel->ghc; ++k) {
                vel->data[index.i][index.j][index.k + k] = vel->data[index.i][index.j][index.k + 1 - k];
            }
        }
    }
}

void slip_node_zl(scalar_data *vel) {
#pragma omp parallel for default(none) shared(vel) collapse(2)
    for (int i = 0; i < vel->nx; i++) {
        for (int j = 0; j < vel->ny; j++) {
            auto index = vel->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= vel->ghc; ++k) {
                vel->data[index.i][index.j][index.k - k] = vel->data[index.i][index.j][index.k - 1 + k];
            }
        }
    }
}

void slip_face_xr(scalar_data *u) {
#pragma omp parallel for default(none) shared(u) collapse(2)
    for (int j = 0; j < u->ny; j++) {
        for (int k = 0; k < u->nz; k++) {
            auto index = u->index_mapping(u->nx, j + 1, k + 1);
            for (int i = 1; i <= u->ghc; ++i) {
                u->data[index.i + i][index.j][index.k] = u->data[index.i - i][index.j][index.k];
            }
        }
    }
}

void slip_face_xl(scalar_data *u) {
#pragma omp parallel for default(none) shared(u) collapse(2)
    for (int j = 0; j < u->ny; j++) {
        for (int k = 0; k < u->nz; k++) {
            auto index = u->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= u->ghc; ++i) {
                u->data[index.i - i][index.j][index.k] = u->data[index.i - 2 + i][index.j][index.k];
            }
        }
    }
}

void slip_face_yr(scalar_data *v) {
#pragma omp parallel for default(none) shared(v) collapse(2)
    for (int i = 0; i < v->nx; i++) {
        for (int k = 0; k < v->nz; k++) {
            auto index = v->index_mapping(i + 1, v->ny, k + 1);
            for (int j = 1; j <= v->ghc; ++j) {
                v->data[index.i][index.j + j][index.k] = v->data[index.i][index.j - j][index.k];
            }
        }
    }
}

void slip_face_yl(scalar_data *v) {
#pragma omp parallel for default(none) shared(v) collapse(2)
    for (int i = 0; i < v->nx; i++) {
        for (int k = 0; k < v->nz; k++) {
            auto index = v->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= v->ghc; ++j) {
                v->data[index.i][index.j - j][index.k] = v->data[index.i][index.j - 2 + j][index.k];
            }
        }
    }
}

void slip_face_zr(scalar_data *w) {
#pragma omp parallel for default(none) shared(w) collapse(2)
    for (int i = 0; i < w->nx; i++) {
        for (int j = 0; j < w->ny; j++) {
            auto index = w->index_mapping(i + 1, j + 1, w->nz);
            for (int k = 1; k <= w->ghc; ++k) {
                w->data[index.i][index.j][index.k + k] = w->data[index.i][index.j][index.k - k];
            }
        }
    }
}

void slip_face_zl(scalar_data *w) {
#pragma omp parallel for default(none) shared(w) collapse(2)
    for (int i = 0; i < w->nx; i++) {
        for (int j = 0; j < w->ny; j++) {
            auto index = w->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= w->ghc; ++k) {
                w->data[index.i][index.j][index.k - k] = w->data[index.i][index.j][index.k - 2 + k];
            }
        }
    }
}

void slip_node_xr(vector_data *vel) {
    slip_node_xr(&vel->x);
    if (vel->x.ndim > 1) {
        slip_node_xr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_xr(&vel->z);
    }
}

void slip_node_xl(vector_data *vel) {
    slip_node_xl(&vel->x);
    if (vel->x.ndim > 1) {
        slip_node_xl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_xl(&vel->z);
    }
}

void slip_node_yr(vector_data *vel) {
    if (vel->x.ndim > 1) {
        slip_node_yr(&vel->x);
        slip_node_yr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_yr(&vel->z);
    }
}

void slip_node_yl(vector_data *vel) {
    if (vel->x.ndim > 1) {
        slip_node_yl(&vel->x);
        slip_node_yl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_yl(&vel->z);
    }
}

void slip_node_zr(vector_data *vel) {
    if (vel->x.ndim > 2) {
        slip_node_zr(&vel->x);
        slip_node_zr(&vel->y);
        slip_node_zr(&vel->z);
    }
}

void slip_node_zl(vector_data *vel) {
    if (vel->x.ndim > 2) {
        slip_node_zl(&vel->x);
        slip_node_zl(&vel->y);
        slip_node_zl(&vel->z);
    }
}

void slip_face_xr(vector_data *vel) {
    slip_face_xr(&vel->x);
    if (vel->x.ndim > 1) {
        slip_node_xr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_xr(&vel->z);
    }
}

void slip_face_xl(vector_data *vel) {
    slip_face_xl(&vel->x);
    if (vel->x.ndim > 1) {
        slip_node_xl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_xl(&vel->z);
    }
}

void slip_face_yr(vector_data *vel) {
    if (vel->x.ndim > 1) {
        slip_node_yr(&vel->x);
        slip_face_yr(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_yr(&vel->z);
    }
}

void slip_face_yl(vector_data *vel) {
    if (vel->x.ndim > 1) {
        slip_node_yl(&vel->x);
        slip_face_yl(&vel->y);
    }
    if (vel->x.ndim > 2) {
        slip_node_yl(&vel->z);
    }
}

void slip_face_zr(vector_data *vel) {
    if (vel->x.ndim > 2) {
        slip_node_zr(&vel->x);
        slip_node_zr(&vel->y);
        slip_face_zr(&vel->z);
    }
}

void slip_face_zl(vector_data *vel) {
    if (vel->x.ndim > 2) {
        slip_node_zl(&vel->x);
        slip_node_zl(&vel->y);
        slip_face_zl(&vel->z);
    }
}

void slip_node(vector_data *nvel) {
    slip_node_xr(nvel);
    slip_node_xl(nvel);
    slip_node_yr(nvel);
    slip_node_yl(nvel);
    slip_node_zr(nvel);
    slip_node_zl(nvel);
}

void slip_face(vector_data *vel) {
    slip_face_xr(vel);
    slip_face_xl(vel);
    slip_face_yr(vel);
    slip_face_yl(vel);
    slip_face_zr(vel);
    slip_face_zl(vel);
}

void slip_face_x(vector_data *vel) {
    if (vel->x.ndim == 1) {
        slip_face_xr(&vel->x);
        slip_face_xl(&vel->x);
    } else if (vel->x.ndim == 2) {
        slip_face_xr(&vel->x);
        slip_face_xr(&vel->y);
        slip_face_xl(&vel->x);
        slip_face_xl(&vel->y);
        slip_node_yr(&vel->x);
        slip_node_yr(&vel->y);
        slip_node_yl(&vel->x);
        slip_node_yl(&vel->y);
    } else if (vel->x.ndim == 3) {
        slip_face_xr(&vel->x);
        slip_face_xr(&vel->y);
        slip_face_xr(&vel->z);
        slip_face_xl(&vel->x);
        slip_face_xl(&vel->y);
        slip_face_xl(&vel->z);
        slip_node_yr(&vel->x);
        slip_node_yr(&vel->y);
        slip_node_yr(&vel->z);
        slip_node_yl(&vel->x);
        slip_node_yl(&vel->y);
        slip_node_yl(&vel->z);
        slip_node_zr(&vel->x);
        slip_node_zr(&vel->y);
        slip_node_zr(&vel->z);
        slip_node_zl(&vel->x);
        slip_node_zl(&vel->y);
        slip_node_zl(&vel->z);
    }
}


















