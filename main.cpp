#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <ctime>

#include "scihpc/wrapper.h"
#include "scihpc/structured_grid.h"
#include "scihpc/flux.h"
#include "scihpc/wrapper_func.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/projection_method.h"
#include "scihpc/vtkWriter.h"
#include "scihpc/bc_factory.h"

int main() {

    auto geo = structured_grid(axis{0.0, 5.0, 100},
                               axis{0.0, 2.0, 40});

    auto phi = wrapper(true, &geo,
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN});
    auto vel = wrapper(false, &geo,
                       bc_info{NO_SLIP}, bc_info{NO_SLIP},
                       bc_info{NO_SLIP}, bc_info{NO_SLIP});
    auto nvel = wrapper(false, &geo,
                        bc_info{NO_SLIP}, bc_info{NO_SLIP},
                        bc_info{NO_SLIP}, bc_info{NO_SLIP});
    auto pressure = wrapper(true, &geo,
                            bc_info{NEUMANN}, bc_info{NEUMANN},
                            bc_info{NEUMANN}, bc_info{NEUMANN});

    auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);
    auto flow_solver = projection_method(phi.scalar);
    auto vtk = vtkWriter(&geo, "dambreak");

    auto param = new problem_parameters{};
    auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);

    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            auto index = phi.scalar->index_mapping(i + 1, j + 1, 1);
            if (geo.xc[i] < 1.0 and geo.yc[j] < 1.0) {
                phi.scalar->data[index.i][index.j][index.k] = 1.0;
            } else {
                phi.scalar->data[index.i][index.j][index.k] = -1.0;
            }
        }
    }
    phi.apply_scalar_bc();

    for (int i = 0; i < phi.scalar->Nx; ++i) {
        for (int j = 0; j < phi.scalar->Ny; ++j) {
            vel.vector->x.data[i][j][0] = 0.0;
            vel.vector->y.data[i][j][0] = 0.0;
            vel.vector->z.data[i][j][0] = 0.0;
            nvel.vector->x.data[i][j][0] = 0.0;
            nvel.vector->y.data[i][j][0] = 0.0;
            nvel.vector->z.data[i][j][0] = 0.0;
            pressure.scalar->data[i][j][0] = 0.0;
        }
    }

    phi.link_params(param);
    phi.link_solvers(deri_solvers);
    phi.link_dummy(dummy);

    pressure.link_params(param);
    pressure.link_solvers(deri_solvers);
    pressure.link_dummy(dummy);

    vel.link_params(param);
    vel.link_solvers(deri_solvers);
    vel.link_dummy(dummy);

    nvel.link_params(param);
    nvel.link_solvers(deri_solvers);
    nvel.link_dummy(dummy);

    param->ls_width = 1.5 * geo.h;
    param->rdt = 0.1 * geo.h;
    param->dt = 0.01 * geo.h;
    param->viscosity_ratio = 0.01;
    param->density_ratio = 0.001;
    param->Froude_number = 1.0;
    param->Reynolds_number = 45000.0;

    int step, instep;

    // Init phi
    find_sign(&phi);
    stabilized_upon_gradient(&phi, &geo);
    step = 0;
    do {
        store_tmp(&phi);
        solver.tvd_rk3(&phi, &vel, identity_flux, lsf_redistance_lambda);
    } while (++step * param->rdt < 5.0 and l2norm(&phi) > 1e-6);

    flow_solver.find_source(&vel, &nvel, &phi);
    param->lsf_mass0 = lsf_mass(&phi);

    vtk.create(0);
    vtk.add_scalar(phi.scalar, "phi");
    vtk.add_vector(nvel.vector, "nvel");
    vtk.close();

    step = 0;
    int pltid = 1;
    do {
        solver.tvd_rk3(&phi, &nvel, &identity_flux, &convection);

        do {
            solver.euler(&phi, &nvel, &identity_flux, &mpls);
        } while (fabs(1.0 - lsf_mass(&phi) / param->lsf_mass0) > 1e-10);

        flow_solver.solve(&vel, &nvel, &pressure, &phi);

        if (++step % 10 == 0) {
            std::cout << "----------------------------------------" << std::endl;
            std::cout << " time: " << step * param->dt << std::endl;
            std::cout << " mass loss ratio: " << fabs(1.0 - lsf_mass(&phi) / param->lsf_mass0) << std::endl;
            std::cout << " div: " << divergence(&vel, &geo) << std::endl;
            std::cout << " l2norm: " << l2norm(&pressure) << std::endl;

            instep = 0;
            while (instep * param->rdt < 2.0 * param->ls_width and step % 20 == 0) {
                if (instep == 0) {
                    find_sign(&phi);
                }
                instep++;
                solver.tvd_rk3(&phi, &nvel, identity_flux, lsf_redistance_lambda);
            };
        }

        if (step % 250 == 0) {
            vtk.create(pltid++);
            vtk.add_scalar(phi.scalar, "phi");
            vtk.add_vector(nvel.vector, "nvel");
            vtk.close();
        }
    } while (step * param->dt < 5.0);

    return 0;
}
