#include <iostream>
#include <cmath>
#include <fstream>

#include "scihpc/scalar_data.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include "scihpc/vector_data.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/source.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/flux.h"

int main() {

    for (int cnt = 0; cnt < 5; ++cnt) {

        auto phi = scalar_data(16 * static_cast<int>(pow(2, cnt)));
        auto nvel = vector_data(phi.nx);
        auto geo = structured_grid(axis{-1.0, 1.0, phi.nx});
        auto solver = runge_kutta(phi.Nx, phi.Ny, phi.Nz);
        for (int i = 0; i < phi.nx; ++i) {
            auto index = phi.index_mapping(i + 1, 1, 1);
            phi.data[index.i][index.j][index.k] = sin(pi * geo.xc[i]);
            nvel.x.data[index.i][index.j][index.k] = 1.0;
        }

        periodic(&phi);

        int tcnt = 0;
        int plt_id = 0;
        DataType cfl = 0.01;
        auto dt = cfl * geo.dx;

        while (tcnt * dt < 2.0) {
//            if (tcnt * dt >= plt_id * 0.1) {
//                std::ofstream fout;
//                fout.open("./phi_" + std::to_string(plt_id) + ".dat");
//                for (int i = 0; i < phi.Nx; ++i) {
//                    auto index = phi.index_mapping(1, 1, 1);
//                    fout << phi.data[i][index.j][index.k] << std::endl;
//                }
//                fout.close();
//                plt_id++;
//            }
            solver.tvd_rk3(dt, &phi, &nvel, &geo, &identity_flux, &periodic, &convection);
            tcnt++;
        }

        DataType error = 0.0;
        for (int i = 0; i < phi.nx; ++i) {
            auto index = phi.index_mapping(i + 1, 1, 1);
            error += pow(phi.data[index.i][index.j][index.k] - sin(pi * (geo.xc[i] - tcnt * dt)), 2);
        }
        std::cout << sqrt(error / phi.nx) << std::endl;

    }

}