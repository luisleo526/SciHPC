//
// Created by 溫晧良 on 2022/10/1.
//

#include "source.h"

void convection(wrapper *f, wrapper *vel, structured_grid *geo, DataType ***s,
                void (*flux)(scalar_data *, vector_data *)) {
    flux(f->scalar, vel->vector);
    f->solvers->uccd->find_derivatives(f->scalar, vel->vector);
#pragma omp parallel for default(none) shared(f, vel, geo, s)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                s[i][j][k] = -vel->vector->x.data[i][j][k] * f->scalar->fx[i][j][k];
            }
        }
    }
    if (f->scalar->ndim > 1) {
#pragma omp parallel for default(none) shared(f, vel, geo, s)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= vel->vector->y.data[i][j][k] * f->scalar->fy[i][j][k];
                }
            }
        }
    }
    if (f->scalar->ndim > 2) {
#pragma omp parallel for default(none) shared(f, vel, geo, s)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= vel->vector->z.data[i][j][k] * f->scalar->fz[i][j][k];
                }
            }
        }
    }
}

void Hamilton_Jacobi(wrapper *f, wrapper *vel, structured_grid *geo, DataType ***s,
                     void (*flux)(scalar_data *, vector_data *)) {
    flux(f->scalar, vel->vector);
    f->solvers->uccd->find_derivatives(f->scalar, vel->vector);
#pragma omp parallel for default(none) shared(f, vel, geo, s)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                s[i][j][k] = -f->scalar->fx[i][j][k];
            }
        }
    }
    if (f->scalar->ndim > 1) {
#pragma omp parallel for default(none) shared(f, vel, geo, s)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= f->scalar->fy[i][j][k];
                }
            }
        }
    }
    if (f->scalar->ndim > 2) {
#pragma omp parallel for default(none) shared(f, vel, geo, s)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= f->scalar->fz[i][j][k];
                }
            }
        }
    }

}

void mpls(wrapper *phi, wrapper *vel, structured_grid *geo, DataType ***s,
          void (*flux)(scalar_data *, vector_data *)) {

    flux(phi->scalar, vel->vector);
    phi->solvers->ccd->find_derivatives(phi->scalar);
    auto mass = lsf_mass(phi);

    DataType eta = 0.0;

#pragma omp parallel for default(none) shared(phi) reduction(+:eta)
    for (int i = 0; i < phi->scalar->nx; ++i) {
        for (int j = 0; j < phi->scalar->ny; ++j) {
            for (int k = 0; k < phi->scalar->nz; ++k) {

                auto index = phi->scalar->index_mapping(i + 1, j + 1, k + 1);

                auto gradient = pow(phi->scalar->fx[index.i][index.j][index.k], 2)
                                + pow(phi->scalar->fy[index.i][index.j][index.k], 2);
                if (phi->scalar->ndim > 2) {
                    gradient += pow(phi->scalar->fz[index.i][index.j][index.k], 2);
                }
                gradient = sqrt(gradient);

                auto heavy = Heaviside(phi->scalar->data[index.i][index.j][index.k], phi->params->ls_width);
                auto delta = Delta(phi->scalar->data[index.i][index.j][index.k], phi->params->ls_width);
                auto g = delta * (2.0 * (1.0 - phi->params->density_ratio) * heavy + phi->params->density_ratio);

                eta += g * delta * gradient;

            }
        }
    }

    eta = (phi->params->lsf_mass0 - mass) / (eta * phi->params->dt);

#pragma omp parallel for default(none) shared(s, phi, eta)
    for (int i = 0; i < phi->scalar->nx; ++i) {
        for (int j = 0; j < phi->scalar->ny; ++j) {
            for (int k = 0; k < phi->scalar->nz; ++k) {

                auto gradient = pow(phi->scalar->fx[i][j][k], 2)
                                + pow(phi->scalar->fy[i][j][k], 2);
                if (phi->scalar->ndim > 2) {
                    gradient += pow(phi->scalar->fz[i][j][k], 2);
                }
                gradient = sqrt(gradient);

                auto delta = Delta(phi->scalar->data[i][j][k], phi->params->ls_width);

                s[i][j][k] = eta * delta * gradient;
            }
        }
    }

}
