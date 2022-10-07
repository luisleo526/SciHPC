#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <ctime>

#include "wrapper.h"
#include "structured_grid.h"
#include "scihpc/flux.h"
#include "simple_bc.h"
#include "wrapper_func.h"
#include "runge_kutta.h"
#include "scihpc/projection_method.h"
#include "scihpc/vtkWriter.h"


int main() {

    auto phi = wrapper(new scalar_data(128, 128));
    auto pressure = wrapper(new scalar_data(phi.scalar->nx, phi.scalar->ny));
    auto vel = wrapper(new vector_data(phi.scalar->nx, phi.scalar->ny));
    auto nvel = wrapper(new vector_data(phi.scalar->nx, phi.scalar->ny));
    auto geo = structured_grid(axis{0.0, 1.0, phi.scalar->nx},
                               axis{0.0, 1.0, phi.scalar->ny});

    auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Ny);
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
    param->Reynolds_number = 1000.0;

    flow_solver.find_source(&vel, &nvel, &phi, &geo);

    vtk.create(0);
    vtk.add_vector(nvel.vector, "nvel");
    vtk.close();

    DataType error = 0.0;
    do {
        store_tmp(&vel);
        flow_solver.solve(&vel, &nvel, &pressure, &phi, &geo);
        error = l2norm(&vel);
    } while (error > 1.0e-5);

    return 0;
}