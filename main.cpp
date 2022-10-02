#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include "scihpc/scalar_data.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include "scihpc/vector_data.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/source.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/flux.h"
#include "scihpc/vtkWriter.h"
#include "scihpc/lsm.h"

int main() {

    auto phi = scalar_data(128, 64);
    auto vel = vector_data(phi.nx, phi.ny);
    auto geo = structured_grid(axis{0.0, 1.0, phi.nx},
                               axis{0.0, 1.0, phi.ny});

    auto solver = runge_kutta(phi.Nx, phi.Ny, phi.Nz);
    auto param = new problem_parameters{};
    param->ls_width = 1.5 * geo.h;
    phi.link(param);

    for (int i = 0; i < phi.nx; ++i) {
        for (int j = 0; j < phi.ny; ++j) {
            auto index = phi.index_mapping(i + 1, j + 1, 1);
            phi.data[index.i][index.j][index.k] = -sqrt(pow(geo.xc[i] - 0.5, 2) +
                                                        pow(geo.yc[j] - 0.75, 2)) + 0.15;
            vel.x.data[index.i][index.j][index.k] =
                    sin(pi * geo.xc[i]) * sin(pi * geo.xc[i]) * sin(2.0 * pi * geo.yc[j]);
            vel.y.data[index.i][index.j][index.k] =
                    -sin(pi * geo.yc[j]) * sin(pi * geo.yc[j]) * sin(2.0 * pi * geo.xc[i]);
        }
    }

    zero_order_extrapolation(&phi);
    zero_order_extrapolation(&vel);

    auto vtk = vtkWriter(&geo, "vortex_deformation");
    vtk.create(0);
    vtk.add_scalar_data(&phi, "phi");
    vtk.add_vector_data(&vel, "vel");
    vtk.close();

    auto mass0 = lsf_mass(&phi);
    DataType period = 4.0;
    auto dt = 0.1 * geo.h;
    int cnt = 0;
    int plt_id = 1;


    auto start = std::chrono::high_resolution_clock::now();
    while (cnt * dt < period) {

        cnt++;

        solver.tvd_rk3(dt, &phi, &vel, &geo, &identity_flux, &zero_order_extrapolation, &convection);

        for (int i = 0; i < phi.nx; ++i) {
            for (int j = 0; j < phi.ny; ++j) {
                auto index = phi.index_mapping(i + 1, j + 1, 1);
                vel.x.data[index.i][index.j][index.k] =
                        sin(pi * geo.xc[i]) * sin(pi * geo.xc[i]) * sin(2.0 * pi * geo.yc[j]) *
                        cos(pi * cnt * dt / period);
                vel.y.data[index.i][index.j][index.k] =
                        -sin(pi * geo.yc[j]) * sin(pi * geo.yc[j]) * sin(2.0 * pi * geo.xc[i]) *
                        cos(pi * cnt * dt / period);
            }
        }
        zero_order_extrapolation(&vel);

        if (cnt * dt >= plt_id * period / 4.0) {
            vtk.create(plt_id);
            vtk.add_scalar_data(&phi, "phi");
            vtk.add_vector_data(&vel, "vel");
            vtk.close();
            plt_id++;
        }

        auto mass = lsf_mass(&phi);

        if (cnt % 5 == 0) {
            std::cout << "Time: " << cnt * dt << " Mass Error: " << (mass0 - mass) / mass0 * 100.0 << "%" << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);


}