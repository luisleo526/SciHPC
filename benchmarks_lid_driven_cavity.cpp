#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <ctime>

#include "wrapper.h"
#include "structured_grid.h"
#include "scihpc/flux.h"
#include "wrapper_func.h"
#include "runge_kutta.h"
#include "scihpc/projection_method.h"
#include "scihpc/vtkWriter.h"


int main() {

    auto geo = structured_grid(axis{0.0, 5.0, 128},
                               axis{0.0, 2.0, 128});

    auto phi = wrapper(true, &geo,
                       bc_info{NEUMANN}, bc_info{NEUMANN},
                       bc_info{NEUMANN}, bc_info{NEUMANN});
    auto vel = wrapper(false, &geo,
                       bc_info{NO_SLIP}, bc_info{NO_SLIP},
                       bc_info{NO_SLIP}, bc_info{NO_SLIP});
    vel.bcFactoryV->yrbc.type = DIRICHLET;
    vel.bcFactoryV->yrbc.value = 1.0;
    auto nvel = wrapper(false, &geo,
                        bc_info{NO_SLIP}, bc_info{NO_SLIP},
                        bc_info{NO_SLIP}, bc_info{NO_SLIP});
    nvel.bcFactoryV->yrbc.type = DIRICHLET;
    nvel.bcFactoryV->yrbc.value = 1.0;
    auto pressure = wrapper(true, &geo,
                            bc_info{NEUMANN}, bc_info{NEUMANN},
                            bc_info{NEUMANN}, bc_info{NEUMANN});

    auto flow_solver = projection_method(phi.scalar);
    auto vtk = vtkWriter(&geo, "lid_driven_cavity");

    auto param = new problem_parameters{};
    auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);

    for (int i = 0; i < phi.scalar->Nx; ++i) {
        for (int j = 0; j < phi.scalar->Ny; ++j) {
            vel.vector->x.data[i][j][0] = 0.0;
            vel.vector->y.data[i][j][0] = 0.0;
            vel.vector->z.data[i][j][0] = 0.0;
            nvel.vector->x.data[i][j][0] = 0.0;
            nvel.vector->y.data[i][j][0] = 0.0;
            nvel.vector->z.data[i][j][0] = 0.0;
            pressure.scalar->data[i][j][0] = 0.0;
            phi.scalar->data[i][j][0] = 1.0;
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
    param->dt = 0.01 * geo.h;
    param->viscosity_ratio = 1.0;
    param->density_ratio = 1.0;
    param->Reynolds_number = 100.0;
    param->ppe_max_iter = 5000;

    flow_solver.find_source(&vel, &nvel, &phi);

    vtk.create(0);
    vtk.add_vector(nvel.vector, "nvel");
    vtk.close();

    DataType error = 0.0;
    do {
        store_tmp(&vel);
        flow_solver.solve(&vel, &nvel, &pressure, &phi);
        error = l2norm(&vel);
        std::cout << error << std::endl;
    } while (error > 1.0e-10);

    return 0;
}