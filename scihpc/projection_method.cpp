//
// Created by 溫晧良 on 2022/10/5.
//

#include "projection_method.h"
#include <iostream>
#include <cassert>
#include "wrapper_func.h"

projection_method::projection_method(scalar_data *f) {
    u_src = init_array(f->Nx, f->Ny, f->Nz);
    v_src = init_array(f->Nx, f->Ny, f->Nz);
    w_src = init_array(f->Nx, f->Ny, f->Nz);

    u_src_old = init_array(f->Nx, f->Ny, f->Nz);
    v_src_old = init_array(f->Nx, f->Ny, f->Nz);
    w_src_old = init_array(f->Nx, f->Ny, f->Nz);

    // For PPE iterative solver
    CC = init_array(f->Nx, f->Ny, f->Nz);
    CR = init_array(f->Nx, f->Ny, f->Nz);
    CL = init_array(f->Nx, f->Ny, f->Nz);
    CU = init_array(f->Nx, f->Ny, f->Nz);
    CD = init_array(f->Nx, f->Ny, f->Nz);
    CF = init_array(f->Nx, f->Ny, f->Nz);
    CB = init_array(f->Nx, f->Ny, f->Nz);
    RHS = init_array(f->Nx, f->Ny, f->Nz);

}

void projection_method::add_stress_x(wrapper *vel, wrapper *lsf, wrapper *nvel) const {

    auto ref = 1.0;
    if (!lsf->params->positive_ref) {
        ref = -1.0;
    }

#pragma omp parallel for default(none) shared(vel, lsf, nvel, ref) collapse(3)
    for (int I = 0; I < lsf->scalar->nx; ++I) {
        for (int J = 0; J < lsf->scalar->ny; ++J) {
            for (int K = 0; K < lsf->scalar->nz; ++K) {

                auto index = lsf->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i + 1][j][k]);
                auto mu = 0.5 * (lsf->dummy->viscosity[i][j][k] + lsf->dummy->viscosity[i + 1][j][k]);
                auto rho = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i + 1][j][k]);

                auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i + 1][j][k]);
                auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i + 1][j][k]);

                auto uxx = vel->vector->x.fxx[i][j][k];
                auto uyy = vel->vector->x.fyy[i][j][k];
                auto ux = vel->vector->x.fx[i][j][k];
                auto uy = vel->vector->x.fy[i][j][k];

                auto vx = nvel->vector->y.fx[i][j][k];

                auto stress_part1 = uxx + uyy;
                auto stress_part2 = 2.0 * phix * ux + phiy * (uy + vx);

                if (lsf->scalar->ndim > 2) {
                    auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i + 1][j][k]);
                    auto uzz = vel->vector->x.fzz[i][j][k];
                    auto uz = vel->vector->x.fz[i][j][k];
                    auto wx = nvel->vector->z.fx[i][j][k];
                    stress_part1 += uzz;
                    stress_part2 += phiz * (uz + wx);
                }

                stress_part1 *= mu;
                stress_part2 *= ref * delta * (1.0 - lsf->params->viscosity_ratio);

                u_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;
            }
        }
    }
}

void projection_method::add_stress_y(wrapper *vel, wrapper *lsf, wrapper *nvel) const {

    auto ref = 1.0;
    if (!lsf->params->positive_ref) {
        ref = -1.0;
    }

#pragma omp parallel for default(none) shared(vel, lsf, nvel, ref) collapse(3)
    for (int I = 0; I < lsf->scalar->nx; ++I) {
        for (int J = 0; J < lsf->scalar->ny; ++J) {
            for (int K = 0; K < lsf->scalar->nz; ++K) {

                auto index = lsf->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j + 1][k]);
                auto mu = 0.5 * (lsf->dummy->viscosity[i][j][k] + lsf->dummy->viscosity[i][j + 1][k]);
                auto rho = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j + 1][k]);

                auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i][j + 1][k]);
                auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i][j + 1][k]);

                auto vxx = vel->vector->y.fxx[i][j][k];
                auto vyy = vel->vector->y.fyy[i][j][k];
                auto vx = vel->vector->y.fx[i][j][k];
                auto vy = vel->vector->y.fy[i][j][k];

                auto uy = nvel->vector->x.fy[i][j][k];

                auto stress_part1 = vxx + vyy;
                auto stress_part2 = 2.0 * phiy * vy + phix * (vx + uy);

                if (lsf->scalar->ndim > 2) {
                    auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i][j + 1][k]);
                    auto vzz = vel->vector->y.fzz[i][j][k];
                    auto vz = vel->vector->y.fz[i][j][k];
                    auto wy = nvel->vector->z.fy[i][j][k];
                    stress_part1 += vzz;
                    stress_part2 += phiz * (vz + wy);
                }

                stress_part1 *= mu;
                stress_part2 *= ref * delta * (1.0 - lsf->params->viscosity_ratio);

                v_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;

            }
        }
    }

}

void projection_method::add_stress_z(wrapper *vel, wrapper *lsf, wrapper *nvel) const {

    auto ref = 1.0;
    if (!lsf->params->positive_ref) {
        ref = -1.0;
    }

#pragma omp parallel for default(none) shared(vel, lsf, nvel, ref) collapse(3)
    for (int I = 0; I < lsf->scalar->nx; ++I) {
        for (int J = 0; J < lsf->scalar->ny; ++J) {
            for (int K = 0; K < lsf->scalar->nz; ++K) {

                auto index = lsf->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j][k + 1]);
                auto mu = 0.5 * (lsf->dummy->viscosity[i][j][k] + lsf->dummy->viscosity[i][j][k + 1]);
                auto rho = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j][k + 1]);

                auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i][j][k + 1]);
                auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i][j][k + 1]);
                auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i][j][k + 1]);

                auto wxx = vel->vector->z.fxx[i][j][k];
                auto wyy = vel->vector->z.fyy[i][j][k];
                auto wx = vel->vector->z.fx[i][j][k];
                auto wy = vel->vector->z.fy[i][j][k];
                auto wz = vel->vector->z.fz[i][j][k];

                auto uz = nvel->vector->x.fz[i][j][k];
                auto vz = nvel->vector->y.fz[i][j][k];

                auto stress_part1 = wxx + wyy;
                auto stress_part2 = 2.0 * phiz * wz + phix * (wx + uz) + phiy * (wy + vz);

                stress_part1 *= mu;
                stress_part2 *= ref * delta * (1.0 - lsf->params->viscosity_ratio);

                w_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;

            }
        }
    }
}

void projection_method::add_body_force(wrapper *lsf) const {

    if (lsf->params->Froude_number > 0.0) {
        if (lsf->scalar->ndim == 2) {
            for (int i = 0; i < lsf->scalar->Nx; ++i) {
                for (int j = 0; j < lsf->scalar->Ny; ++j) {
                    for (int k = 0; k < lsf->scalar->Nz; ++k) {
                        v_src[i][j][k] += -1.0 / lsf->params->Froude_number;
                    }
                }
            }
        } else if (lsf->scalar->ndim == 3) {
            for (int i = 0; i < lsf->scalar->Nx; ++i) {
                for (int j = 0; j < lsf->scalar->Ny; ++j) {
                    for (int k = 0; k < lsf->scalar->Nz; ++k) {
                        w_src[i][j][k] += -1.0 / lsf->params->Froude_number;
                    }
                }
            }
        }
    }
}

void projection_method::add_surface_force(wrapper *lsf) const {

    if (lsf->params->Weber_number > 0.0) {

        find_curvature(lsf);
        for (int i = 0; i < lsf->scalar->Nx - 1; ++i) {
            for (int j = 0; j < lsf->scalar->Ny - 1; ++j) {
                for (int k = 0; k < lsf->scalar->Nz; ++k) {
                    auto curv = 0.5 * (lsf->dummy->curvature[i][j][k] + lsf->dummy->curvature[i + 1][j][k]);
                    auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i + 1][j][k]);
                    auto phix = 0.5 * (lsf->scalar->fx[i][j][k] + lsf->scalar->fx[i + 1][j][k]);
                    auto rho = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i + 1][j][k]);

                    u_src[i][j][k] -= curv * delta * phix / (lsf->params->Weber_number * rho);

                    curv = 0.5 * (lsf->dummy->curvature[i][j][k] + lsf->dummy->curvature[i][j + 1][k]);
                    delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j + 1][k]);
                    auto phiy = 0.5 * (lsf->scalar->fy[i][j][k] + lsf->scalar->fy[i][j + 1][k]);
                    rho = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j + 1][k]);

                    v_src[i][j][k] -= curv * delta * phiy / (lsf->params->Weber_number * rho);

                    if (lsf->scalar->ndim > 2 && k < lsf->scalar->Nz - 1) {
                        curv = 0.5 * (lsf->dummy->curvature[i][j][k] + lsf->dummy->curvature[i][j][k + 1]);
                        delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j][k + 1]);
                        auto phiz = 0.5 * (lsf->scalar->fz[i][j][k] + lsf->scalar->fz[i][j][k + 1]);
                        rho = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j][k + 1]);

                        w_src[i][j][k] -= curv * delta * phiz / (lsf->params->Weber_number * rho);
                    }

                }
            }
        }
    }
}

void projection_method::find_source(wrapper *vel, wrapper *nvel, wrapper *lsf) const {

#pragma omp parallel for default(none) shared(vel, nvel, lsf) collapse(3)
    for (int i = 0; i < vel->vector->x.Nx; ++i) {
        for (int j = 0; j < vel->vector->x.Ny; ++j) {
            for (int k = 0; k < vel->vector->x.Nz; ++k) {
                u_src_old[i][j][k] = u_src[i][j][k];
                v_src_old[i][j][k] = v_src[i][j][k];
                w_src_old[i][j][k] = w_src[i][j][k];
            }
        }
    }

    // prepare level set function
    find_delta(lsf);
    find_heavyside(lsf);
    find_density(lsf);
    find_viscosity(lsf);
    identity_with_extrapolation(lsf->scalar);
    lsf->solvers->ccd->find_derivatives_all(lsf->scalar);

    // Compute velocities at face x for convection part
    all_to_face_x(vel, nvel);
    nvel->apply_vel_x_bc();
    // Compute convection source
    vel->scalar = &vel->vector->x;
    convection(vel, nvel, u_src, identity_with_extrapolation_face_x);
    // Compute derivatives of face x velocities for stress tensor at x direction
    identity_with_extrapolation_face_x(nvel->vector);
    nvel->solvers->uccd->find_fx(&nvel->vector->y, nvel->vector);
    if (nvel->vector->x.ndim > 2) {
        nvel->solvers->uccd->find_fx(&nvel->vector->z, nvel->vector);
    }
    // Compute stress tensor at x direction
    add_stress_x(vel, lsf, nvel);

    // Compute velocities at face y for convection part
    all_to_face_y(vel, nvel);
    nvel->apply_vel_y_bc();
    // Compute convection source
    vel->scalar = &vel->vector->y;
    convection(vel, nvel, v_src, identity_with_extrapolation_face_y);
    // Compute derivatives of face y velocities for stress tensor at y direction
    identity_with_extrapolation_face_y(nvel->vector);
    nvel->solvers->uccd->find_fy(&nvel->vector->x, nvel->vector);
    if (nvel->vector->x.ndim > 2) {
        nvel->solvers->uccd->find_fy(&nvel->vector->z, nvel->vector);
    }
    // Compute stress tensor at y direction
    add_stress_y(vel, lsf, nvel);

    if (lsf->scalar->ndim > 2) {
        // Compute velocities at face z for convection part
        all_to_face_z(vel, nvel);
        nvel->apply_vel_z_bc();
        // Compute convection source
        vel->scalar = &vel->vector->z;
        convection(vel, nvel, w_src, identity_with_extrapolation_face_z);
        // Compute derivatives of face z velocities for stress tensor at z direction
        identity_with_extrapolation_face_z(nvel->vector);
        nvel->solvers->uccd->find_fz(&nvel->vector->x, nvel->vector);
        nvel->solvers->uccd->find_fz(&nvel->vector->y, nvel->vector);
        // Compute stress tensor at z direction
        add_stress_z(vel, lsf, nvel);
    }

    add_body_force(lsf);
    add_surface_force(lsf);
}

void projection_method::find_source_sec(wrapper *vel, wrapper *nvel, wrapper *lsf) const {

#pragma omp parallel for default(none) shared(vel, nvel, lsf) collapse(3)
    for (int i = 0; i < vel->vector->x.Nx; ++i) {
        for (int j = 0; j < vel->vector->x.Ny; ++j) {
            for (int k = 0; k < vel->vector->x.Nz; ++k) {
                u_src_old[i][j][k] = u_src[i][j][k];
                v_src_old[i][j][k] = v_src[i][j][k];
                w_src_old[i][j][k] = w_src[i][j][k];
            }
        }
    }

    // prepare level set function
    find_delta(lsf);
    find_heavyside(lsf);
    find_density(lsf);
    find_viscosity(lsf);
    identity_with_extrapolation(lsf->scalar);
    lsf->solvers->secSol->find_derivatives_all(lsf->scalar);

    // Compute velocities at face x for convection part
    all_to_face_x(vel, nvel);
    nvel->apply_vel_x_bc();
    // Compute convection source
    vel->scalar = &vel->vector->x;
    convection_sec(vel, nvel, u_src, identity_with_extrapolation_face_x);
    // Compute derivatives of face x velocities for stress tensor at x direction
    identity_with_extrapolation_face_x(nvel->vector);
    nvel->solvers->secSol->find_fx(&nvel->vector->y, nvel->vector);
    if (nvel->vector->x.ndim > 2) {
        nvel->solvers->secSol->find_fx(&nvel->vector->z, nvel->vector);
    }
    // Compute stress tensor at x direction
    add_stress_x(vel, lsf, nvel);

    // Compute velocities at face y for convection part
    all_to_face_y(vel, nvel);
    nvel->apply_vel_y_bc();
    // Compute convection source
    vel->scalar = &vel->vector->y;
    convection_sec(vel, nvel, v_src, identity_with_extrapolation_face_y);
    // Compute derivatives of face y velocities for stress tensor at y direction
    identity_with_extrapolation_face_y(nvel->vector);
    nvel->solvers->secSol->find_fy(&nvel->vector->x, nvel->vector);
    if (nvel->vector->x.ndim > 2) {
        nvel->solvers->secSol->find_fy(&nvel->vector->z, nvel->vector);
    }
    // Compute stress tensor at y direction
    add_stress_y(vel, lsf, nvel);

    if (lsf->scalar->ndim > 2) {
        // Compute velocities at face z for convection part
        all_to_face_z(vel, nvel);
        nvel->apply_vel_z_bc();
        // Compute convection source
        vel->scalar = &vel->vector->z;
        convection_sec(vel, nvel, w_src, identity_with_extrapolation_face_z);
        // Compute derivatives of face z velocities for stress tensor at z direction
        identity_with_extrapolation_face_z(nvel->vector);
        nvel->solvers->secSol->find_fz(&nvel->vector->x, nvel->vector);
        nvel->solvers->secSol->find_fz(&nvel->vector->y, nvel->vector);
        // Compute stress tensor at z direction
        add_stress_z(vel, lsf, nvel);
    }

    add_body_force(lsf);
    add_surface_force(lsf);
}

void projection_method::find_intermediate_velocity(wrapper *vel) const {

#pragma omp parallel for default(none) shared(vel) collapse(3)
    for (int i = 0; i < vel->vector->x.Nx; ++i) {
        for (int j = 0; j < vel->vector->x.Ny; ++j) {
            for (int k = 0; k < vel->vector->x.Nz; ++k) {
                vel->vector->x.data[i][j][k] +=
                        vel->params->dt * (1.5 * u_src[i][j][k] - 0.5 * u_src_old[i][j][k]);
                vel->vector->y.data[i][j][k] +=
                        vel->params->dt * (1.5 * v_src[i][j][k] - 0.5 * v_src_old[i][j][k]);
                vel->vector->z.data[i][j][k] +=
                        vel->params->dt * (1.5 * w_src[i][j][k] - 0.5 * w_src_old[i][j][k]);
            }
        }
    }

    vel->apply_vel_bc();

}

void projection_method::init_ppe_proj_coefficients(wrapper *pressure) {

#pragma omp parallel for default(none) shared(pressure) collapse(3)
    for (int I = 0; I < pressure->scalar->nx; ++I) {
        for (int J = 0; J < pressure->scalar->ny; ++J) {
            for (int K = 0; K < pressure->scalar->nz; ++K) {

                auto index = pressure->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto R = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i + 1][j][k]);
                auto L = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i - 1][j][k]);
                auto U = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j + 1][k]);
                auto D = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j - 1][k]);

                CR[i][j][k] = 1.0 / R / pressure->geo->dx / pressure->geo->dx;
                CL[i][j][k] = 1.0 / L / pressure->geo->dx / pressure->geo->dx;
                CU[i][j][k] = 1.0 / U / pressure->geo->dy / pressure->geo->dy;
                CD[i][j][k] = 1.0 / D / pressure->geo->dy / pressure->geo->dy;

                if (I == 0) {
                    CL[i][j][k] = 0.0;
                }

                if (I == pressure->scalar->nx - 1) {
                    CR[i][j][k] = 0.0;
                }

                if (J == 0) {
                    CD[i][j][k] = 0.0;
                }

                if (J == pressure->scalar->ny - 1) {
                    CU[i][j][k] = 0.0;
                }

                if (pressure->scalar->ndim > 2) {
                    auto F = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j][k + 1]);
                    auto B = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j][k - 1]);

                    CF[i][j][k] = 1.0 / F / pressure->geo->dz / pressure->geo->dz;
                    CB[i][j][k] = 1.0 / B / pressure->geo->dz / pressure->geo->dz;
                    if (K == 0) {
                        CB[i][j][k] = 0.0;
                    }
                    if (K == pressure->scalar->nz - 1) {
                        CF[i][j][k] = 0.0;
                    }
                }

                CC[i][j][k] = -(CR[i][j][k] + CL[i][j][k] + CU[i][j][k] + CD[i][j][k]);

                if (pressure->scalar->ndim > 2) {
                    CC[i][j][k] -= (CF[i][j][k] + CB[i][j][k]);
                }

            }
        }
    }
}

void projection_method::init_ppe_fpc_coefficients(wrapper *pressure) {

    auto rho12 = pressure->params->density_ratio;

#pragma omp parallel for default(none) shared(pressure, rho12) collapse(3)
    for (int I = 0; I < pressure->scalar->nx; ++I) {
        for (int J = 0; J < pressure->scalar->ny; ++J) {
            for (int K = 0; K < pressure->scalar->nz; ++K) {

                auto index = pressure->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto R = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i + 1][j][k]);
                auto L = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i - 1][j][k]);
                auto U = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j + 1][k]);
                auto D = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j - 1][k]);

                CR[i][j][k] = (1.0 - rho12 / R) / pressure->geo->dx / pressure->geo->dx;
                CL[i][j][k] = (1.0 - rho12 / L) / pressure->geo->dx / pressure->geo->dx;
                CU[i][j][k] = (1.0 - rho12 / U) / pressure->geo->dy / pressure->geo->dy;
                CD[i][j][k] = (1.0 - rho12 / D) / pressure->geo->dy / pressure->geo->dy;

                if (I == 0) {
                    CL[i][j][k] = 0.0;
                }

                if (I == pressure->scalar->nx - 1) {
                    CR[i][j][k] = 0.0;
                }

                if (J == 0) {
                    CD[i][j][k] = 0.0;
                }

                if (J == pressure->scalar->ny - 1) {
                    CU[i][j][k] = 0.0;
                }

                if (pressure->scalar->ndim > 2) {
                    auto F = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j][k + 1]);
                    auto B = 0.5 * (pressure->dummy->density[i][j][k] + pressure->dummy->density[i][j][k - 1]);

                    CF[i][j][k] = (1.0 - rho12 / F) / pressure->geo->dz / pressure->geo->dz;
                    CB[i][j][k] = (1.0 - rho12 / B) / pressure->geo->dz / pressure->geo->dz;
                    if (K == 0) {
                        CB[i][j][k] = 0.0;
                    }
                    if (K == pressure->scalar->nz - 1) {
                        CF[i][j][k] = 0.0;
                    }
                }

                CC[i][j][k] = -(CR[i][j][k] + CL[i][j][k] + CU[i][j][k] + CD[i][j][k]);

                if (pressure->scalar->ndim > 2) {
                    CC[i][j][k] -= (CF[i][j][k] + CB[i][j][k]);
                }

                RHS[i][j][k] = CC[i][j][k] * pressure->scalar->flux[i][j][k] +
                               CR[i][j][k] * pressure->scalar->flux[i + 1][j][k] +
                               CL[i][j][k] * pressure->scalar->flux[i - 1][j][k] +
                               CU[i][j][k] * pressure->scalar->flux[i][j + 1][k] +
                               CD[i][j][k] * pressure->scalar->flux[i][j - 1][k];

                if (pressure->scalar->ndim > 2) {
                    RHS[i][j][k] += CF[i][j][k] * pressure->scalar->flux[i][j][k + 1] +
                                    CB[i][j][k] * pressure->scalar->flux[i][j][k - 1];
                }

            }
        }
    }
}


void projection_method::projection(wrapper *pressure, wrapper *lsf, wrapper *vel) const {

#pragma omp parallel for default(none) shared(pressure, lsf, vel) collapse(3)
    for (int I = 0; I < lsf->scalar->nx; ++I) {
        for (int J = 0; J < lsf->scalar->ny; ++J) {
            for (int K = 0; K < lsf->scalar->nz; ++K) {

                auto index = lsf->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                RHS[i][j][k] = (vel->vector->x.data[i][j][k] - vel->vector->x.data[i - 1][j][k]) /
                               pressure->geo->dx +
                               (vel->vector->y.data[i][j][k] - vel->vector->y.data[i][j - 1][k]) /
                               pressure->geo->dy;
                if (lsf->scalar->ndim > 2) {
                    RHS[i][j][k] += (vel->vector->z.data[i][j][k] - vel->vector->z.data[i][j][k - 1]) /
                                    pressure->geo->dz;
                }
                RHS[i][j][k] /= lsf->params->dt;

            }
        }
    }

    // solve pressure poisson equation

    int iter = 0;
    DataType error, sump;
    do {
        store_tmp(pressure);
        sump = 0.0;
        error = 0.0;
#pragma omp parallel for default(none) shared(pressure, lsf) reduction(+:sump, error) collapse(3)
        for (int i = 0; i < pressure->scalar->nx; ++i) {
            for (int j = 0; j < pressure->scalar->ny; ++j) {
                for (int k = 0; k < pressure->scalar->nz; ++k) {
                    auto index = pressure->scalar->index_mapping(i + 1, j + 1, k + 1);
                    auto I = index.i;
                    auto J = index.j;
                    auto K = index.k;

                    auto pnew = (RHS[I][J][K] -
                                 (CL[I][J][K] * pressure->scalar->data[I - 1][J][K] +
                                  CR[I][J][K] * pressure->scalar->data[I + 1][J][K] +
                                  CD[I][J][K] * pressure->scalar->data[I][J - 1][K] +
                                  CU[I][J][K] * pressure->scalar->data[I][J + 1][K]));
                    if (lsf->scalar->ndim > 2) {
                        pnew -= (CB[I][J][K] * pressure->scalar->data[I][J][K - 1] +
                                 CF[I][J][K] * pressure->scalar->data[I][J][K + 1]);
                    }

                    pnew /= CC[I][J][K];

                    pressure->scalar->data[I][J][K] = pressure->params->ppe_omega * pnew +
                                                      (1.0 - pressure->params->ppe_omega) *
                                                      pressure->scalar->data[I][J][K];
                    sump += pressure->scalar->data[I][J][K];
                    error += fabs(pnew * CC[I][J][K] - CC[I][J][K] * pressure->scalar->data[I][J][K]);
                }
            }
        }
        sump = sump / (pressure->scalar->nx * pressure->scalar->ny);
        if (lsf->scalar->ndim > 2) {
            sump = sump / pressure->scalar->nz;
        }

#pragma omp parallel for default(none) shared(pressure, sump) collapse(3)
        for (int i = 0; i < pressure->scalar->nx; ++i) {
            for (int j = 0; j < pressure->scalar->ny; ++j) {
                for (int k = 0; k < pressure->scalar->nz; ++k) {
                    auto index = pressure->scalar->index_mapping(i + 1, j + 1, k + 1);
                    pressure->scalar->data[index.i][index.j][index.k] -= sump;
                }
            }
        }
        pressure->apply_scalar_bc();
        if (++iter % 10000 == 0) {
            std::cout << "iter: " << iter << " error: " << error << " linfity-p: " << linfity(pressure) << std::endl;
        }
    } while (iter < pressure->params->ppe_max_iter and error > pressure->params->ppe_tol and
             linfity(pressure) > pressure->params->ppe_tol2);

    // update velocity

#pragma omp parallel for default(none) shared(vel, pressure, lsf) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx - 1; ++i) {
        for (int j = 0; j < lsf->scalar->Ny - 1; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {

                auto density = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i + 1][j][k]);
                vel->vector->x.data[i][j][k] -= pressure->params->dt / pressure->geo->dx / density *
                                                (pressure->scalar->data[i + 1][j][k] -
                                                 pressure->scalar->data[i][j][k]);

                density = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j + 1][k]);
                vel->vector->y.data[i][j][k] -= pressure->params->dt / pressure->geo->dy / density *
                                                (pressure->scalar->data[i][j + 1][k] -
                                                 pressure->scalar->data[i][j][k]);

                if (lsf->scalar->ndim > 2 && k < lsf->scalar->Nz - 1) {
                    density = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j][k + 1]);
                    vel->vector->z.data[i][j][k] -= pressure->params->dt / pressure->geo->dz / density *
                                                    (pressure->scalar->data[i][j][k + 1] -
                                                     pressure->scalar->data[i][j][k]);
                }

            }
        }
    }
    vel->apply_vel_bc();

}

void projection_method::fast_pressure_correction(wrapper *pressure, wrapper *lsf, wrapper *vel) {

    auto dx = pressure->geo->dx;
    auto dy = pressure->geo->dy;
    auto dz = pressure->geo->dz;
    auto dt = pressure->params->dt;
    auto rho12 = pressure->params->density_ratio;

#pragma omp parallel for default(none) shared(pressure, lsf, vel, dx, dy, dz, dt, rho12) collapse(3)
    for (int I = 0; I < lsf->scalar->nx; ++I) {
        for (int J = 0; J < lsf->scalar->ny; ++J) {
            for (int K = 0; K < lsf->scalar->nz; ++K) {

                auto index = lsf->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;

                auto pos = pressure->solvers->mg->at[0]->of(I, J, K);

                pressure->solvers->mg->at[0]->rhs[pos] =
                        (vel->vector->x.data[i][j][k] - vel->vector->x.data[i - 1][j][k]) / dx +
                        (vel->vector->y.data[i][j][k] - vel->vector->y.data[i][j - 1][k]) / dy;
                if (lsf->scalar->ndim > 2) {
                    pressure->solvers->mg->at[0]->rhs[pos] +=
                            (vel->vector->z.data[i][j][k] - vel->vector->z.data[i][j][k - 1]) / dz;
                }
                pressure->solvers->mg->at[0]->rhs[pos] *= rho12 / lsf->params->dt;

                pressure->solvers->mg->at[0]->rhs[pos] += RHS[i][j][k];

            }
        }
    }

    auto error = pressure->solvers->mg->at[0]->residual();
    int iter = 0;

    while (error > pressure->params->ppe_tol && iter < 100) {
        iter++;
        store_tmp(pressure);
        pressure->solvers->mg->full_cycle();
        error = pressure->solvers->mg->at[0]->residual();
        if (iter % 10 == 0) {
            std::cout << "PPE residual from MG: " << error << ","
                      << pressure->solvers->mg->at[pressure->solvers->mg->level_num - 1]->residual() << std::endl;
        }
    }

#pragma omp parallel for default(none) shared(pressure) collapse(3)
    for (int I = 0; I < pressure->scalar->nx; ++I) {
        for (int J = 0; J < pressure->scalar->ny; ++J) {
            for (int K = 0; K < pressure->scalar->nz; ++K) {
                auto index = pressure->scalar->index_mapping(I + 1, J + 1, K + 1);
                auto i = index.i;
                auto j = index.j;
                auto k = index.k;
                auto pos = pressure->solvers->mg->at[0]->of(I, J, K);
                pressure->scalar->data[i][j][k] = pressure->solvers->mg->at[0]->sol[pos];
            }
        }
    }
    pressure->apply_scalar_bc();

#pragma omp parallel for default(none) shared(vel, pressure, lsf, dx, dy, dz, dt, rho12) collapse(3)
    for (int i = 0; i < lsf->scalar->Nx - 1; ++i) {
        for (int j = 0; j < lsf->scalar->Ny - 1; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {

                auto px = (pressure->scalar->data[i + 1][j][k] - pressure->scalar->data[i][j][k]) / dx;
                auto px2 = (pressure->scalar->flux[i + 1][j][k] - pressure->scalar->flux[i][j][k]) / dx;
                auto density = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i + 1][j][k]);

                vel->vector->x.data[i][j][k] -= dt * ((px - px2) / rho12 + px2 / density);

                auto py = (pressure->scalar->data[i][j + 1][k] - pressure->scalar->data[i][j][k]) / dy;
                auto py2 = (pressure->scalar->flux[i][j + 1][k] - pressure->scalar->flux[i][j][k]) / dy;
                density = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j + 1][k]);

                vel->vector->y.data[i][j][k] -= dt * ((py - py2) / rho12 + py2 / density);

                if (lsf->scalar->ndim > 2 && k < lsf->scalar->Nz - 1) {
                    auto pz = (pressure->scalar->data[i][j][k + 1] - pressure->scalar->data[i][j][k]) / dz;
                    auto pz2 = (pressure->scalar->flux[i][j][k + 1] - pressure->scalar->flux[i][j][k]) / dz;
                    density = 0.5 * (lsf->dummy->density[i][j][k] + lsf->dummy->density[i][j][k + 1]);

                    vel->vector->z.data[i][j][k] -= dt * ((pz - pz2) / rho12 + pz2 / density);
                }

            }
        }
    }
    vel->apply_vel_bc();

}

void projection_method::ab_solve(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf) {
    pressure->scalar->store();
    find_source(vel, nvel, lsf);
    find_intermediate_velocity(vel);

    for (int cnt = 0; cnt < pressure->params->ppe_initer; ++cnt) {
        if (pressure->params->iter < 5) {
            if (cnt == 0) {
                init_ppe_proj_coefficients(pressure);
            }
            projection(pressure, lsf, vel);
        } else {
            if (cnt == 0) {
                pressure_guess_solution(pressure);
                init_ppe_fpc_coefficients(pressure);
            }
            fast_pressure_correction(pressure, lsf, vel);
        }
    }

    node_from_face(vel, nvel);
    nvel->apply_nvel_bc();
}

void projection_method::ab_solve_sec(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf) {
    pressure->scalar->store();
    find_source_sec(vel, nvel, lsf);
    find_intermediate_velocity(vel);

    for (int cnt = 0; cnt < pressure->params->ppe_initer; ++cnt) {
        if (pressure->params->iter < 5) {
            if (cnt == 0) {
                init_ppe_proj_coefficients(pressure);
            }
            projection(pressure, lsf, vel);
        } else {
            if (cnt == 0) {
                pressure_guess_solution(pressure);
                init_ppe_fpc_coefficients(pressure);
            }
            fast_pressure_correction(pressure, lsf, vel);
        }
    }

    node_from_face(vel, nvel);
    nvel->apply_nvel_bc();
}

void projection_method::pressure_guess_solution(wrapper *pressure) {

#pragma omp parallel for default(none) shared(pressure) collapse(3)
    for (int i = 0; i < pressure->scalar->Nx; ++i) {
        for (int j = 0; j < pressure->scalar->Ny; ++j) {
            for (int k = 0; k < pressure->scalar->Nz; ++k) {
                auto dpdt =
                        (pressure->scalar->old[i][j][k] - pressure->scalar->old2[i][j][k]) / pressure->params->dt_old2;
                pressure->scalar->flux[i][j][k] = pressure->scalar->old[i][j][k] + pressure->params->dt_old * dpdt;
            }
        }
    }
}
