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

int main() {

    std::cout.precision(4);
    DataType prev_error = 0.0;
    DataType prev_h = 0.0;

    for (int cnt = 0; cnt < 5; ++cnt) {

        auto phi = scalar_data(16 * static_cast<int>(pow(2, cnt)));
        auto vel = vector_data(phi.nx);
        auto geo = structured_grid(axis{-1.0, 1.0, phi.nx});
        auto solver = runge_kutta(phi.Nx, phi.Ny, phi.Nz);
        for (int i = 0; i < phi.nx; ++i) {
            auto index = phi.index_mapping(i + 1, 1, 1);
            phi.data[index.i][index.j][index.k] = sin(pi * geo.xc[i]);
            vel.x.data[index.i][index.j][index.k] = 1.0;
        }

        periodic(&phi);

        int tcnt = 0;
        DataType cfl = 0.01;
        auto dt = cfl * geo.dx;

        auto begin = std::chrono::high_resolution_clock::now();
        while (tcnt * dt < 2.0) {
            solver.tvd_rk3(dt, &phi, &vel, &geo, &identity_flux, &periodic, &convection);
            tcnt++;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        // calculate l2 norm
        DataType error = 0.0;
        for (int i = 0; i < phi.nx; ++i) {
            auto index = phi.index_mapping(i + 1, 1, 1);
            error += pow(phi.data[index.i][index.j][index.k] - sin(pi * (geo.xc[i] - tcnt * dt)), 2);
        }
        error = sqrt(error / phi.nx);

        // output
        std::cout << std::setw(5) << phi.nx << " | " << std::scientific << error;
        if (cnt > 0) {
            std::cout << " | " << std::fixed << std::setw(7)
                      << (log10(error) - log10(prev_error)) / (log10(geo.dx) - log10(prev_h));
        } else {
            std::cout << " | " << std::setw(7) << "-";
        }
        std::cout << " | " << std::fixed << std::setw(7) << elapsed.count() << std::endl;

        prev_error = error;
        prev_h = geo.dx;

    }

}