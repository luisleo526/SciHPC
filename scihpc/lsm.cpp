//
// Created by 溫晧良 on 2022/10/1.
//

#include "lsm.h"

DataType ***Heaviside(scalar_data *f) {
    auto heaviside = init_array(f->Nx, f->Ny, f->Nz);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                if (f->data[i][j][k] < -f->params->ls_width) {
                    heaviside[i][j][k] = 0;
                } else if (f->data[i][j][k] > f->params->ls_width) {
                    heaviside[i][j][k] = 1;
                } else {
                    heaviside[i][j][k] = 0.5 * (1 + f->data[i][j][k] / f->params->ls_width +
                                        sin(pi * f->data[i][j][k] / f->params->ls_width) / pi);
                }
            }
        }
    }
    return heaviside;
}

DataType ***Delta(scalar_data *f) {
    auto delta = init_array(f->Nx, f->Ny, f->Nz);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                if (f->data[i][j][k] < -f->params->ls_width) {
                    delta[i][j][k] = 0;
                } else if (f->data[i][j][k] > f->params->ls_width) {
                    delta[i][j][k] = 0;
                } else {
                    delta[i][j][k] = 0.5 * (1 + cos(pi * f->data[i][j][k] / f->params->ls_width)) / f->params->ls_width;
                }
            }
        }
    }
    return delta;
}

DataType ***Sign(scalar_data *f) {
    auto sign = init_array(f->Nx, f->Ny, f->Nz);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                if (f->data[i][j][k] < -f->params->ls_width) {
                    sign[i][j][k] = -1;
                } else if (f->data[i][j][k] > f->params->ls_width) {
                    sign[i][j][k] = 1;
                } else {
                    sign[i][j][k] = f->data[i][j][k] /
                            sqrt(f->data[i][j][k] * f->data[i][j][k] + f->params->ls_width * f->params->ls_width);
                }
            }
        }
    }
    return sign;
}

DataType lsf_mass(scalar_data *f) {
    DataType mass = 0;
    auto h = Heaviside(f);
#pragma omp parallel for reduction(+:mass)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i+1, j+1, k+1);
                auto rho = h[index.i][index.j][index.k] + (1.0-h[index.i][index.j][index.k]) * f->params->density_ratio;
                mass += rho * h[index.i][index.j][index.k];
            }
        }
    }
    delete3d(h, f->Nx, f->Ny);
    return mass;
}

DataType lsf_volume(scalar_data *f) {
    DataType volume = 0;
    auto h = Heaviside(f);
#pragma omp parallel for reduction(+:volume)
    for (int i = 0; i < f->nx; ++i) {
        for (int j = 0; j < f->ny; ++j) {
            for (int k = 0; k < f->nz; ++k) {
                auto index = f->index_mapping(i+1, j+1, k+1);
                volume += h[index.i][index.j][index.k];
            }
        }
    }
    delete3d(h, f->Nx, f->Ny);
    return volume;
}




