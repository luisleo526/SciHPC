//
// Created by 溫晧良 on 2022/10/5.
//

#include "projection_method.h"
#include <iostream>

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


#pragma omp parallel for default(none) shared(vel, lsf, geo, nvel, u_src) collapse(3)
    for (int I = 1; I <= lsf->scalar->nx; ++I) {
        for (int J = 1; J <= lsf->scalar->ny; ++J) {
            for (int k = 0; k < lsf->scalar->nz; ++k) {

                auto index = lsf->scalar->index_mapping(I, J, 1);
                auto i = index.i;
                auto j = index.j;

                auto h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i + 1][j][k]);
                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i + 1][j][k]);
                auto mu = h + (1.0 - h) * lsf->params->viscosity_ratio;
                auto rho = h + (1.0 - h) * lsf->params->density_ratio;

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
                stress_part2 *= delta * (1.0 - lsf->params->viscosity_ratio);

                u_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;
            }
        }
    }
}

void projection_method::add_stress_y(wrapper *vel, wrapper *lsf, wrapper *nvel) const {


#pragma omp parallel for default(none) shared(vel, lsf, geo, nvel, v_src) collapse(3)
    for (int I = 1; I <= lsf->scalar->nx; ++I) {
        for (int J = 1; J <= lsf->scalar->ny; ++J) {
            for (int k = 0; k < lsf->scalar->nz; ++k) {

                auto index = lsf->scalar->index_mapping(I, J, 1);
                auto i = index.i;
                auto j = index.j;

                auto h = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j + 1][k]);
                auto delta = 0.5 * (lsf->dummy->delta[i][j][k] + lsf->dummy->delta[i][j + 1][k]);
                auto mu = h + (1.0 - h) * lsf->params->viscosity_ratio;
                auto rho = h + (1.0 - h) * lsf->params->density_ratio;

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
                stress_part2 *= delta * (1.0 - lsf->params->viscosity_ratio);

                v_src[i][j][k] += (stress_part1 + stress_part2) / lsf->params->Reynolds_number / rho;
//                                  - 1.0 / lsf->params->Froude_number;

            }
        }
    }

}

void projection_method::add_stress_z(wrapper *vel, wrapper *lsf, wrapper *nvel) const {

}

void
projection_method::find_source(wrapper *vel, wrapper *nvel, wrapper *lsf) const {

#pragma omp parallel for default(none) collapse(3) shared(u_src, v_src, w_src, u_src_old, v_src_old, w_src_old, vel)
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
    find_curvature(lsf);

    all_to_face_x(vel, nvel);
    nvel->apply_vel_x_bc();
    vel->scalar = &vel->vector->x;
    convection(vel, nvel, u_src, identity_flux);
    identity_flux(nvel->vector);
    nvel->solvers->uccd->find_derivatives(&nvel->vector->y, nvel->vector);
    add_stress_x(vel, lsf, nvel);

    all_to_face_y(vel, nvel);
    nvel->apply_vel_y_bc();
    vel->scalar = &vel->vector->y;
    convection(vel, nvel, v_src, identity_flux);
    identity_flux(nvel->vector);
    nvel->solvers->uccd->find_derivatives(&nvel->vector->y, nvel->vector);
    add_stress_y(vel, lsf, nvel);

    if (lsf->scalar->ndim > 2) {
        all_to_face_z(vel, nvel);
        nvel->apply_vel_z_bc();
        vel->scalar = &vel->vector->z;
        convection(vel, nvel, w_src, identity_flux);
        identity_flux(nvel->vector);
        nvel->solvers->uccd->find_derivatives(&nvel->vector->z, nvel->vector);
        add_stress_z(vel, lsf, nvel);
    }

}

void projection_method::find_intermediate_velocity(wrapper *vel) const {

#pragma omp parallel for default(none) shared(vel, u_src, v_src, w_src) collapse(3)
    for (int i = 0; i < vel->vector->x.Nx; ++i) {
        for (int j = 0; j < vel->vector->x.Ny; ++j) {
            for (int k = 0; k < vel->vector->x.Nz; ++k) {
                vel->vector->x.data[i][j][k] += vel->params->dt * (1.5 * u_src[i][j][k] - 0.5 * u_src_old[i][j][k]);
                vel->vector->y.data[i][j][k] += vel->params->dt * (1.5 * v_src[i][j][k] - 0.5 * v_src_old[i][j][k]);
                vel->vector->z.data[i][j][k] += vel->params->dt * (1.5 * w_src[i][j][k] - 0.5 * w_src_old[i][j][k]);
            }
        }
    }

    vel->apply_vel_bc();

}

void projection_method::solve_ppe(wrapper *pressure, wrapper *lsf, wrapper *vel) const {

#pragma omp parallel for default(none) shared(lsf, geo, vel, CR, CL, CU, CD, CC, RHS) collapse(3)
    for (int i = 1; i < lsf->scalar->Nx - 1; ++i) {
        for (int j = 1; j < lsf->scalar->Ny - 1; ++j) {
            for (int k = 0; k < lsf->scalar->Nz; ++k) {
                auto R = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i + 1][j][k]);
                auto L = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i - 1][j][k]);
                auto U = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j + 1][k]);
                auto D = 0.5 * (lsf->dummy->heaviside[i][j][k] + lsf->dummy->heaviside[i][j - 1][k]);

                R = R + (1.0 - R) * lsf->params->density_ratio;
                L = L + (1.0 - L) * lsf->params->density_ratio;
                U = U + (1.0 - U) * lsf->params->density_ratio;
                D = D + (1.0 - D) * lsf->params->density_ratio;

                CR[i][j][k] = 1.0 / R / pressure->geo->dx / pressure->geo->dx;
                CL[i][j][k] = 1.0 / L / pressure->geo->dx / pressure->geo->dx;
                CU[i][j][k] = 1.0 / U / pressure->geo->dy / pressure->geo->dy;
                CD[i][j][k] = 1.0 / D / pressure->geo->dy / pressure->geo->dy;

                CC[i][j][k] = -(CR[i][j][k] + CL[i][j][k] + CU[i][j][k] + CD[i][j][k]);

                RHS[i][j][k] = (vel->vector->x.data[i][j][k] - vel->vector->x.data[i - 1][j][k]) / pressure->geo->dx +
                               (vel->vector->y.data[i][j][k] - vel->vector->y.data[i][j - 1][k]) / pressure->geo->dy;
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
#pragma omp parallel for default(none) shared(lsf, geo, pressure) collapse(2) reduction(+:sump) reduction(max:error)
        for (int i = 0; i < pressure->scalar->nx; ++i) {
            for (int j = 0; j < pressure->scalar->ny; ++j) {
                auto index = pressure->scalar->index_mapping(i + 1, j + 1, 1);
                auto I = index.i;
                auto J = index.j;
                auto K = index.k;

                auto pnew = (RHS[I][J][K] -
                             (CL[I][J][K] * pressure->scalar->data[I - 1][J][K] +
                              CR[I][J][K] * pressure->scalar->data[I + 1][J][K] +
                              CD[I][J][K] * pressure->scalar->data[I][J - 1][K] +
                              CU[I][J][K] * pressure->scalar->data[I][J + 1][K])) / CC[I][J][K];

                pressure->scalar->data[I][J][K] = 1.5 * pnew - 0.5 * pressure->scalar->data[I][J][K];
                sump += pressure->scalar->data[I][J][K];
                error = fmax(error, fabs(pnew * CC[I][J][K] - CC[I][J][K] * pressure->scalar->data[I][J][K]));
            }
        }
        sump = sump / (pressure->scalar->nx * pressure->scalar->ny);

#pragma omp parallel for default(none) shared(pressure, sump) collapse(2)
        for (int i = 0; i < pressure->scalar->nx; ++i) {
            for (int j = 0; j < pressure->scalar->ny; ++j) {
                auto index = pressure->scalar->index_mapping(i + 1, j + 1, 1);
                pressure->scalar->data[index.i][index.j][index.k] -= sump;
            }
        }
        pressure->apply_scalar_bc();
        if (++iter % 5000 == 0) {
            std::cout << "iter: " << iter << " error: " << error << " sump: " << sump << std::endl;
        }
    } while (iter < pressure->params->ppe_max_iter and error > pressure->params->ppe_tol);

}

void projection_method::find_final_velocity(wrapper *vel, wrapper *pressure, wrapper *lsf) {

#pragma omp parallel for default(none) shared(vel, pressure, geo, lsf) collapse(2)
    for (int i = 0; i < vel->vector->x.Nx - 1; ++i) {
        for (int j = 0; j < vel->vector->x.Ny - 1; ++j) {
            auto h = 0.5 * (lsf->dummy->heaviside[i][j][0] + lsf->dummy->heaviside[i + 1][j][0]);
            auto density = (1 - h) * vel->params->density_ratio + h;
            vel->vector->x.data[i][j][0] -= pressure->params->dt / pressure->geo->dx / density *
                                            (pressure->scalar->data[i + 1][j][0] - pressure->scalar->data[i][j][0]);

            h = 0.5 * (lsf->dummy->heaviside[i][j][0] + lsf->dummy->heaviside[i][j + 1][0]);
            density = (1 - h) * vel->params->density_ratio + h;

            vel->vector->y.data[i][j][0] -= pressure->params->dt / pressure->geo->dy / density *
                                            (pressure->scalar->data[i][j + 1][0] - pressure->scalar->data[i][j][0]);
        }
    }
    vel->apply_vel_bc();
}

void projection_method::solve(wrapper *vel, wrapper *nvel, wrapper *pressure, wrapper *lsf) const {
    find_source(vel, nvel, lsf);
    find_intermediate_velocity(vel);
    solve_ppe(pressure, lsf, vel);
    find_final_velocity(vel, pressure, lsf);
    node_from_face(vel, nvel);
    nvel->apply_nvel_bc();
}
