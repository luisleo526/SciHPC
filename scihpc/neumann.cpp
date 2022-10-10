//
// Created by 溫晧良 on 2022/10/10.
//

#include "neumann.h"

void Neumann_node_xr(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(f->nx, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i + i][index.j][index.k] =
                        f->data[index.i + i - (2 * i - 1)][index.j][index.k] + (2 * i - 1) * df * h;
            }
        }
    }
}

void Neumann_node_xl(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i - i][index.j][index.k] =
                        f->data[index.i - i + (2 * i - 1)][index.j][index.k] - (2 * i - 1) * df * h;
            }
        }
    }
}

void Neumann_node_yr(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, f->ny, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j + j][index.k] =
                        f->data[index.i][index.j + j - (2 * j - 1)][index.k] + (2 * j - 1) * df * h;
            }
        }
    }
}

void Neumann_node_yl(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j - j][index.k] =
                        f->data[index.i][index.j - j + (2 * j - 1)][index.k] - (2 * j - 1) * df * h;
            }
        }
    }
}

void Neumann_node_zr(scalar_data *f, DataType df, DataType h) {
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, f->nz);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k + k] =
                        f->data[index.i][index.j][index.k + k - (2 * k - 1)] + (2 * k - 1) * df * h;
            }
        }
    }
}

void Neumann_node_zl(scalar_data *f, DataType df, DataType h) {
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k - k] =
                        f->data[index.i][index.j][index.k - k + (2 * k - 1)] - (2 * k - 1) * df * h;
            }
        }
    }
}

void Neumann_face_xr(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(f->nx, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i + i][index.j][index.k] =
                        f->data[index.i - i][index.j][index.k] + 2 * i * df * h;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i - 1][index.j][index.k] + df * h;

        }
    }
}

void Neumann_face_xl(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(0, j + 1, k + 1);
            for (int i = 1; i <= f->ghc - 1; ++i) {
                f->data[index.i - i][index.j][index.k] = f->data[index.i + i][index.j][index.k] - 2 * i * df * h;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i + 1][index.j][index.k] - df * h;
        }
    }
}

void Neumann_face_yr(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, f->ny, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j + j][index.k] =
                        f->data[index.i][index.j - j][index.k] + 2 * j * df * h;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j - 1][index.k] + df * h;
        }
    }
}

void Neumann_face_yl(scalar_data *f, DataType df, DataType h) {
    for (int k = 0; k < f->nz; k++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, 0, k + 1);
            for (int j = 1; j <= f->ghc - 1; ++j) {
                f->data[index.i][index.j - j][index.k] = f->data[index.i][index.j + j][index.k] - 2 * j * df * h;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j + 1][index.k] - df * h;
        }
    }
}

void Neumann_face_zr(scalar_data *f, DataType df, DataType h) {
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, f->nz);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k + k] =
                        f->data[index.i][index.j][index.k - k] + 2 * k * df * h;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j][index.k - 1] + df * h;
        }
    }
}

void Neumann_face_zl(scalar_data *f, DataType df, DataType h) {
    for (int j = 0; j < f->ny; j++) {
        for (int i = 0; i < f->nx; i++) {
            auto index = f->index_mapping(i + 1, j + 1, 0);
            for (int k = 1; k <= f->ghc - 1; ++k) {
                f->data[index.i][index.j][index.k - k] = f->data[index.i][index.j][index.k + k] - 2 * k * df * h;
            }
            f->data[index.i][index.j][index.k] = f->data[index.i][index.j][index.k + 1] - df * h;
        }
    }
}




