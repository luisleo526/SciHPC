//
// Created by 溫晧良 on 2022/10/10.
//

#include "periodic.h"

void periodic_node_xr(scalar_data *f) {
    for (int j = 0; j < f->ny; j++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(f->nx, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i + i][index.j][index.k] = f->data[f->ghc + (i - 1)][index.j][index.k];
            }
        }
    }
}

void periodic_node_xl(scalar_data *f) {
    for (int j = 0; j < f->ny; j++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i - i][index.j][index.k] = f->data[f->Nx - f->ghc - 1 - (i - 1)][index.j][index.k];
            }
        }
    }
}

void periodic_node_yr(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(i + 1, f->ny, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j + j][index.k] = f->data[index.i][f->ghc + (j - 1)][index.k];
            }
        }
    }
}

void periodic_node_yl(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j - j][index.k] = f->data[index.i][f->Ny - f->ghc - 1 - (j - 1)][index.k];
            }
        }
    }
}

void periodic_node_zr(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(i + 1, j + 1, f->nz);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k + k] = f->data[index.i][index.j][f->ghc + (k - 1)];
            }
        }
    }
}

void periodic_node_zl(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k - k] = f->data[index.i][index.j][f->Nz - f->ghc - 1 - (k - 1)];
            }
        }
    }
}

//  ...                                                            ...
//  ...    ╔═══════╦═══════╦═══════╗ | ╔═══════╦═══════╦═══════╗   ...
//  ...    ║       ║       ║       ║ | ║       ║       ║       ║   ...
//  ...    ║       ║       ║       ║ | ║       ║       ║       ║   ...
//  ...    ╠═══════╬═══════╬═══════╣ | ╠═══════╬═══════╬═══════╣   ...
//  ...    ║       ║       ║       ║ | ║       ║       ║       ║   ...
//  ...    ║ nx-2  ║ nx-1  ║  nx   ║ | ║   1   ║   2   ║   3   ║   ...
//  ...    ╠═══════╬═══════╬═══════╣ | ╠═══════╬═══════╬═══════╣   ...
//  ...    ║       ║       ║       ║ | ║       ║       ║       ║   ...
//  ...    ║       ║       ║       ║ | ║       ║       ║       ║   ...
//  ...    ╚═══════╩═══════╩═══════╝ | ╚═══════╩═══════╩═══════╝   ...

void periodic_face_xr(scalar_data *f) {
    for (int j = 0; j < f->ny; j++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(f->nx, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i + i][index.j][index.k] = f->data[f->ghc + (i - 1)][index.j][index.k];
            }
        }
    }
}

void periodic_face_xl(scalar_data *f) {
    for (int j = 0; j < f->ny; j++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(1, j + 1, k + 1);
            for (int i = 1; i <= f->ghc; ++i) {
                f->data[index.i - i][index.j][index.k] = f->data[f->Nx - f->ghc - 1 - (i - 1)][index.j][index.k];
            }
        }
    }
}

void periodic_face_yr(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(i + 1, f->ny, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j + j][index.k] = f->data[index.i][f->ghc + (j - 1)][index.k];
            }
        }
    }
}

void periodic_face_yl(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int k = 0; k < f->nz; k++) {
            auto index = f->index_mapping(i + 1, 1, k + 1);
            for (int j = 1; j <= f->ghc; ++j) {
                f->data[index.i][index.j - j][index.k] = f->data[index.i][f->Ny - f->ghc - 1 - (j - 1)][index.k];
            }
        }
    }
}

void periodic_face_zr(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(i + 1, j + 1, f->nz);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k + k] = f->data[index.i][index.j][f->ghc + (k - 1)];
            }
        }
    }
}

void periodic_face_zl(scalar_data *f) {
    for (int i = 0; i < f->nx; i++) {
        for (int j = 0; j < f->ny; j++) {
            auto index = f->index_mapping(i + 1, j + 1, 1);
            for (int k = 1; k <= f->ghc; ++k) {
                f->data[index.i][index.j][index.k - k] = f->data[index.i][index.j][f->Nz - f->ghc - 1 - (k - 1)];
            }
        }
    }
}