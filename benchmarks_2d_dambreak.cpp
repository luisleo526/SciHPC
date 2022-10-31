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

    auto geo = structured_grid(axis{0.0, 5.0, 128},
                               axis{0.0, 5.0, 128});

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

    auto param = set_air_water(0.0571);
    auto deri_solvers = SharedSolvers_alloc(phi.scalar, &geo);
    shared_solvers_mg_init_Neumann(deri_solvers);
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
    param->rdt = 0.5 * geo.h;
    param->dt = 0.01 * geo.h;
    param->max_CFL = 0.05;
    param->Weber_number = -1.0;
    param->ppe_tol = 1e-8;
    param->ppe_initer = 1;
    param->ppe_tol2 = 1e-10;

    std::cout << "Reynolds number: " << param->Reynolds_number << std::endl;
    std::cout << "Weber number: " << param->Weber_number << std::endl;
    std::cout << "Fround number: " << param->Froude_number << std::endl;

    int step, instep;

    // Init phi
    find_sign(&phi);
    stabilized_upon_gradient(&phi);
    step = 0;
    do {
        store_tmp(&phi);
        solver.tvd_rk3(&phi, &vel, identity_flux, lsf_redistance_lambda);
    } while (++step * param->rdt < 5.0 and l2norm(&phi) > 1e-6);

    flow_solver.find_source(&vel, &nvel, &phi);
    param->lsf_mass0 = lsf_mass(&phi);

    vtk.create(0);
    vtk.add_scalar(phi.scalar, "phi");
    vtk.add_scalar(pressure.scalar, "pressure");
    vtk.add_vector(nvel.vector, "nvel");
    vtk.close();

    step = 0;
    int pltid = 1;
    int reinit_id = 1;
    do {

//        find_dt(&vel);
        param->t += param->dt;

        solver.tvd_rk3(&phi, &nvel, &identity_flux, &convection);

        do {
            solver.euler(&phi, &nvel, &identity_flux, &mpls);
        } while (fabs(1.0 - lsf_mass(&phi) / param->lsf_mass0) > 1e-10);

        flow_solver.ab_solve(&vel, &nvel, &pressure, &phi);

        if (++step % 10 == 0) {
            std::cout << "----------------------------------------" << std::endl;
            std::cout << " time: " << param->t << std::endl;
            std::cout << " stable CFL: " << param->stable_CFL << std::endl;
            std::cout << " mass loss ratio(%): " << fabs(1.0 - lsf_mass(&phi) / param->lsf_mass0) * 100 << std::endl;
            std::cout << " div: " << divergence(&vel) << std::endl;
            std::cout << " l2norm: " << l2norm(&pressure) << std::endl;
        }

        instep = 0;
        while (instep * param->rdt < 2.0 * param->ls_width and param->t > reinit_id * 20 * param->dt) {
            if (instep == 0) {
                find_sign(&phi);
                reinit_id++;
            }
            instep++;
            solver.tvd_rk3(&phi, &nvel, identity_flux, lsf_redistance_lambda);
        };

        if (param->t > pltid * 0.1) {
            vtk.create(pltid++);
            vtk.add_scalar(phi.scalar, "phi");
            vtk.add_scalar(pressure.scalar, "pressure");
            vtk.add_vector(nvel.vector, "nvel");
            vtk.close();
        }
    } while (param->t < 10.0);

    return 0;
}
