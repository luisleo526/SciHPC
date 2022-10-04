#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include "scalar_data.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include "scihpc/vector_data.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/source.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/flux.h"
#include "scihpc/vtkWriter.h"
#include "scihpc/derivatives_solver.h"
#include "scihpc/wrapper.h"
#include "scihpc/lsm.h"

int main() {

    auto phi = wrapper(new scalar_data(64, 64));
    auto vel = wrapper(new vector_data(phi.scalar->nx, phi.scalar->ny));
    auto geo = structured_grid(axis{0.0, 1.0, phi.scalar->nx},
                               axis{0.0, 1.0, phi.scalar->ny});

    auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);

    for (int i = 0; i < phi.scalar->nx; ++i) {
        for (int j = 0; j < phi.scalar->ny; ++j) {
            auto index = phi.scalar->index_mapping(i + 1, j + 1, 1);
            phi.scalar->data[index.i][index.j][index.k] = -sqrt(pow(geo.xc[i] - 0.5, 2) +
                                                                pow(geo.yc[j] - 0.75, 2)) + 0.15;
            vel.vector->x.data[index.i][index.j][index.k] =
                    sin(pi * geo.xc[i]) * sin(pi * geo.xc[i]) * sin(2.0 * pi * geo.yc[j]);
            vel.vector->y.data[index.i][index.j][index.k] =
                    -sin(pi * geo.yc[j]) * sin(pi * geo.yc[j]) * sin(2.0 * pi * geo.xc[i]);
        }
    }

    zero_order_extrapolation(phi.scalar);
    zero_order_extrapolation(vel.vector);

    auto vtk = vtkWriter(&geo, "vortex_deformation");
    vtk.create(0);
    vtk.add_scalar_data(phi.scalar, "phi");
    vtk.add_vector_data(vel.vector, "vel");
    vtk.close();

    auto param = new problem_parameters{};
    param->density_ratio = 1.0;
    param->viscosity_ratio = 1.0;
    param->ls_width = 1.5 * geo.h;
    param->dt = 0.1 * geo.h;
    auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
    auto dummy = dummy_data_alloc(phi.scalar);
    phi.link_params(param);
    phi.link_solvers(deri_solvers);
    phi.link_dummy(dummy);


    param->lsf_mass0 = lsf_mass(&phi);

    DataType period = 16.0;
    int cnt = 0;
    int plt_id = 1;

    auto start = std::chrono::high_resolution_clock::now();
    while (cnt * param->dt < period) {

        cnt++;

        solver.tvd_rk3(&phi, &vel, &geo, &identity_flux, &zero_order_extrapolation, &convection);
        do {
            solver.euler(&phi, &vel, &geo, &identity_flux, &zero_order_extrapolation, &mpls);
        } while ( fabs(1.0 - lsf_mass(&phi) / param->lsf_mass0) > 1e-10);

        for (int i = 0; i < phi.scalar->nx; ++i) {
            for (int j = 0; j < phi.scalar->ny; ++j) {
                auto index = phi.scalar->index_mapping(i + 1, j + 1, 1);
                vel.vector->x.data[index.i][index.j][index.k] =
                        sin(pi * geo.xc[i]) * sin(pi * geo.xc[i]) * sin(2.0 * pi * geo.yc[j]) *
                        cos(pi * cnt * param->dt / period);
                vel.vector->y.data[index.i][index.j][index.k] =
                        -sin(pi * geo.yc[j]) * sin(pi * geo.yc[j]) * sin(2.0 * pi * geo.xc[i]) *
                        cos(pi * cnt * param->dt / period);
            }
        }
        zero_order_extrapolation(vel.vector);

        if (cnt * param->dt >= plt_id * period / 8.0) {
            vtk.create(plt_id);
            vtk.add_scalar_data(phi.scalar, "phi");
            vtk.add_vector_data(vel.vector, "vel");
            vtk.close();
            plt_id++;
        }

        DataType mass = lsf_mass(&phi);

        if (cnt % 5 == 0) {
            std::cout << "Time: " << cnt * param->dt << " Mass Error: " << (1.0 - mass / param->lsf_mass0) * 100.0
                      << "%"
                      << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "CPU time" << duration.count() / 1e6 << "s" << std::endl;

}