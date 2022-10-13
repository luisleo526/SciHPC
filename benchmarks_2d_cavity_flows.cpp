//
// Created by 溫晧良 on 2022/10/11.
//
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


int main() {

    auto geo = structured_grid(axis{0.0, 1.0, 64},
                               axis{0.0, 1.0, 64});

    auto phi = wrapper(true, &geo,
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN});
    auto vel = wrapper(false, &geo,
                       bc_info{NO_SLIP}, bc_info{NO_SLIP},
                       bc_info{NO_SLIP}, bc_info{DIRICHLET});
    vel.bcFactoryU->yrbc.value = 1.0;
    vel.bcFactoryV->yrbc.value = 0.0;
    auto nvel = wrapper(false, &geo,
                        bc_info{NO_SLIP}, bc_info{NO_SLIP},
                        bc_info{NO_SLIP}, bc_info{DIRICHLET});
    nvel.bcFactoryU->yrbc.value = 1.0;
    nvel.bcFactoryV->yrbc.value = 0.0;
    auto pressure = wrapper(true, &geo,
                            bc_info{NEUMANN}, bc_info{NEUMANN},
                            bc_info{NEUMANN}, bc_info{NEUMANN});

    auto flow_solver = projection_method(phi.scalar);
    auto vtk = vtkWriter(&geo, "lid_driven_cavity");

    auto param = new problem_parameters{};
    auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);

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
    param->dt = 0.1 * geo.h;
    param->viscosity_ratio = 1.0;
    param->density_ratio = 1.0;
    param->Reynolds_number = 100.0;
    param->Froude_number = -1.0;
    param->ppe_omega = 1.5;
    param->ppe_tol2 = 1e-6;

    flow_solver.find_source(&vel, &nvel, &phi);

    vtk.create(0);
    vtk.add_scalar(pressure.scalar, "pressure");
    vtk.add_vector(nvel.vector, "nvel");
    vtk.close();

    DataType error = 0.0;
    int step = 0, pltid = 1;
    do {
        step++;
        store_tmp(&vel);
        flow_solver.ab_solve(&vel, &nvel, &pressure, &phi);
        error = l2norm(&vel);
        std::cout << "------------------------------------" << std::endl;
        std::cout << "time : " << step * param->dt << std::endl;
        std::cout << "error: " << error << std::endl;
        std::cout << "div: " << divergence(&vel) << std::endl;
        if (step % 50 == 0) {
            vtk.create(pltid++);
            vtk.add_scalar(pressure.scalar, "pressure");
            vtk.add_vector(nvel.vector, "nvel");
            vtk.close();
        }
    } while (error > 1.0e-10);

    return 0;
}