//
// Created by 溫晧良 on 2022/10/1.
//

#include "source.h"
#include "wrapper_func.h"

void convection(wrapper *f, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *)) {
    flux(f->scalar, vel->vector);
    f->solvers->uccd->find_derivatives(f->scalar, vel->vector);
//    f->solvers->weno->weno5_find_derivatives(f->scalar, vel->vector);
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                s[i][j][k] = -vel->vector->x.data[i][j][k] * f->scalar->fx[i][j][k];
            }
        }
    }
    if (f->scalar->ndim > 1) {
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= vel->vector->y.data[i][j][k] * f->scalar->fy[i][j][k];
                }
            }
        }
    }
    if (f->scalar->ndim > 2) {
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= vel->vector->z.data[i][j][k] * f->scalar->fz[i][j][k];
                }
            }
        }
    }
}

void convection_sec(wrapper *f, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *)) {
    flux(f->scalar, vel->vector);
    f->solvers->secSol->find_derivatives(f->scalar, vel->vector);
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                s[i][j][k] = -vel->vector->x.data[i][j][k] * f->scalar->fx[i][j][k];
            }
        }
    }
    if (f->scalar->ndim > 1) {
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= vel->vector->y.data[i][j][k] * f->scalar->fy[i][j][k];
                }
            }
        }
    }
    if (f->scalar->ndim > 2) {
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= vel->vector->z.data[i][j][k] * f->scalar->fz[i][j][k];
                }
            }
        }
    }
}

void Hamilton_Jacobi(wrapper *f, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *)) {
    flux(f->scalar, vel->vector);
    f->solvers->uccd->find_derivatives(f->scalar, vel->vector);
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
    for (int i = 0; i < f->scalar->Nx; ++i) {
        for (int j = 0; j < f->scalar->Ny; ++j) {
            for (int k = 0; k < f->scalar->Nz; ++k) {
                s[i][j][k] = -f->scalar->fx[i][j][k];
            }
        }
    }
    if (f->scalar->ndim > 1) {
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= f->scalar->fy[i][j][k];
                }
            }
        }
    }
    if (f->scalar->ndim > 2) {
#pragma omp parallel for default(none) shared(f, vel, s) collapse(3)
        for (int i = 0; i < f->scalar->Nx; ++i) {
            for (int j = 0; j < f->scalar->Ny; ++j) {
                for (int k = 0; k < f->scalar->Nz; ++k) {
                    s[i][j][k] -= f->scalar->fz[i][j][k];
                }
            }
        }
    }

}

void mpls(wrapper *phi, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *)) {

    flux(phi->scalar, vel->vector);
    identity_with_extrapolation(phi->scalar);
    phi->solvers->ccd->find_derivatives(phi->scalar);
    auto mass = lsf_mass(phi);

    find_heavyside(phi);
    find_delta(phi);
    find_gradient(phi);

    DataType eta = 0.0;

#pragma omp parallel for default(none) shared(phi) reduction(+:eta) collapse(3)
    for (int i = 0; i < phi->scalar->nx; ++i) {
        for (int j = 0; j < phi->scalar->ny; ++j) {
            for (int k = 0; k < phi->scalar->nz; ++k) {

                auto index = phi->scalar->index_mapping(i + 1, j + 1, k + 1);

                auto gradient = phi->dummy->grad[index.i][index.j][index.k];
                auto heavy = phi->dummy->heaviside[index.i][index.j][index.k];
                auto delta = phi->dummy->delta[index.i][index.j][index.k];
                auto g = delta * (2.0 * (1.0 - phi->params->density_ratio) * heavy + phi->params->density_ratio);
                if (!phi->params->positive_ref) {
                    g = delta * (2.0 * (phi->params->density_ratio - 1.0) * heavy + 1.0);
                }
                eta += g * delta * gradient;
            }
        }
    }

    eta = (phi->params->lsf_mass0 - mass) / (eta * phi->geo->dv * phi->params->dt);

#pragma omp parallel for default(none) shared(s, phi, eta) collapse(3)
    for (int i = 0; i < phi->scalar->Nx; ++i) {
        for (int j = 0; j < phi->scalar->Ny; ++j) {
            for (int k = 0; k < phi->scalar->Nz; ++k) {
                auto gradient = phi->dummy->grad[i][j][k];
                auto delta = phi->dummy->delta[i][j][k];
                s[i][j][k] = eta * delta * gradient;
            }
        }
    }

}

void
lsf_redistance_no_lambda(wrapper *phi, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *)) {

    godunov_gradient(phi);
#pragma omp parallel for default(none) shared(phi, s) collapse(3)
    for (int i = 0; i < phi->scalar->Nx; ++i) {
        for (int j = 0; j < phi->scalar->Ny; ++j) {
            for (int k = 0; k < phi->scalar->Nz; ++k) {
                s[i][j][k] = -phi->dummy->sign[i][j][k] * (phi->dummy->grad[i][j][k] - 1.0) * phi->params->rdt /
                             phi->params->dt;
            }
        }
    }
}

void lsf_redistance_lambda(wrapper *phi, wrapper *vel, DataType ***s, void (*flux)(scalar_data *, vector_data *)) {

    lsf_redistance_no_lambda(phi, vel, s, flux);

#pragma omp parallel for default(none) shared(phi) collapse(3)
    for (int i = 0; i < phi->scalar->Nx; ++i) {
        for (int j = 0; j < phi->scalar->Ny; ++j) {
            for (int k = 0; k < phi->scalar->Nz; ++k) {
                phi->dummy->a[i][j][k] =
                        phi->dummy->sign[i][j][k] * phi->dummy->delta[i][j][k] * (phi->dummy->grad[i][j][k] - 1.0);
                phi->dummy->b[i][j][k] =
                        phi->dummy->grad[i][j][k] * phi->dummy->delta[i][j][k] * phi->dummy->delta[i][j][k];
            }
        }
    }

    zero_order_extrapolation(phi->dummy->a,
                             phi->scalar->Nx, phi->scalar->Ny, phi->scalar->Nz,
                             phi->scalar->ghc);
    zero_order_extrapolation(phi->dummy->b,
                             phi->scalar->Nx, phi->scalar->Ny, phi->scalar->Nz,
                             phi->scalar->ghc);

#pragma omp parallel for default(none) shared(phi, s) collapse(3)
    for (int I = 0; I < phi->scalar->nx; ++I) {
        for (int J = 0; J < phi->scalar->ny; ++J) {
            for (int K = 0; K < phi->scalar->nz; ++K) {
                auto index = phi->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto a = 0.0;
                auto b = 0.0;
                if (phi->scalar->ndim == 2) {
                    a = 16.0 * phi->dummy->a[i][j][k];
                    b = 16.0 * phi->dummy->b[i][j][k];
                    for (int ii = -1; ii < 2; ++ii) {
                        for (int jj = -1; jj < 2; ++jj) {
                            if (ii != 0 && jj != 0) {
                                a += phi->dummy->a[i + ii][j + jj][k];
                                b += phi->dummy->b[i + ii][j + jj][k];
                            }
                        }
                    }
                } else if (phi->scalar->ndim == 3) {
                    a = 51.0 * phi->dummy->a[i][j][k];
                    b = 51.0 * phi->dummy->b[i][j][k];
                    for (int ii = -1; ii < 2; ++ii) {
                        for (int jj = -1; jj < 2; ++jj) {
                            for (int kk = -1; kk < 2; ++kk) {
                                if (ii != 0 && jj != 0 && kk != 0) {
                                    a += phi->dummy->a[i + ii][j + jj][k + kk];
                                    b += phi->dummy->b[i + ii][j + jj][k + kk];
                                }
                            }
                        }
                    }
                }

                if (fabs(b) > 1e-8) {
                    s[i][j][k] += a / b *
                                  phi->dummy->delta[i][j][k] * phi->dummy->grad[i][j][k] * phi->params->rdt /
                                  phi->params->dt;
                }

            }
        }
    }

}
