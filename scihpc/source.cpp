//
// Created by 溫晧良 on 2022/10/1.
//

#include "source.h"

void convection(scalar_data *f, vector_data *vel, structured_grid *geo, DataType ***s,
                void(*flux)(scalar_data *, vector_data *)) {
    flux(f, vel);
    uccd_find_derivatives(f, geo, vel);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                s[i][j][k] = -vel->x.data[i][j][k] * f->fx[i][j][k];
            }
        }
    }
    if (f->ndim > 1) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= vel->y.data[i][j][k] * f->fy[i][j][k];
                }
            }
        }
    }
    if (f->ndim > 2) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= vel->z.data[i][j][k] * f->fz[i][j][k];
                }
            }
        }
    }
}

void Hamilton_Jacobi(scalar_data *f, vector_data *vel, structured_grid *geo, DataType ***s,
                     void (*flux)(scalar_data *, vector_data *)) {
    flux(f, vel);
    uccd_find_derivatives(f, geo, vel);
#pragma omp parallel for
    for (int i = 0; i < f->Nx; ++i) {
        for (int j = 0; j < f->Ny; ++j) {
            for (int k = 0; k < f->Nz; ++k) {
                s[i][j][k] = -f->fx[i][j][k];
            }
        }
    }
    if (f->ndim > 1) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= f->fy[i][j][k];
                }
            }
        }
    }
    if (f->ndim > 2) {
#pragma omp parallel for
        for (int i = 0; i < f->Nx; ++i) {
            for (int j = 0; j < f->Ny; ++j) {
                for (int k = 0; k < f->Nz; ++k) {
                    s[i][j][k] -= f->fz[i][j][k];
                }
            }
        }
    }

}

void mpls(scalar_data *phi, vector_data *vel, structured_grid *geo, DataType ***s,
          void (*flux)(scalar_data *, vector_data *)) {

    flux(phi, vel);
    ccd_find_derivatives(phi, geo);
    auto mass = lsf_mass(phi);
    auto heavy = Heaviside(phi);
    auto delta = Delta(phi);

    DataType eta = 0.0;

#pragma omp parallel for reduction(+:eta)
    for (int i = 0; i < phi->nx; ++i) {
        for (int j = 0; j < phi->ny; ++j) {
            for (int k = 0; k < phi->nz; ++k) {

                auto index = phi->index_mapping(i + 1, j + 1, k + 1);

                auto gradient = pow(phi->fx[index.i][index.j][index.k], 2)
                                     + pow(phi->fy[index.i][index.j][index.k], 2);
                if (phi->ndim > 2) {
                    gradient += pow(phi->fz[index.i][index.j][index.k], 2);
                }
                gradient = sqrt(gradient);

                auto g = delta[index.i][index.j][index.k] *
                         (2.0 * (1.0 - phi->params->density_ratio) * heavy[index.i][index.j][index.k] + phi->params->density_ratio);

                eta += g * delta[index.i][index.j][index.k] * gradient;

            }
        }
    }

    eta = (phi->params->lsf_mass0 - mass) / (eta * phi->params->dt);

#pragma omp parallel for
    for (int i = 0; i < phi->Nx; ++i) {
        for (int j = 0; j < phi->Ny; ++j) {
            for (int k = 0; k < phi->Nz; ++k) {

                auto gradient = pow(phi->fx[i][j][k], 2)
                                      + pow(phi->fy[i][j][k], 2);
                if (phi->ndim > 2) {
                    gradient += pow(phi->fz[i][j][k], 2);
                }
                gradient = sqrt(gradient);

                s[i][j][k] = eta * delta[i][j][k] * gradient;
            }
        }
    }

    delete3d(heavy, phi->nx, phi->ny);
    delete3d(delta, phi->nx, phi->ny);

}
