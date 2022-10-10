#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

#include "scihpc/scalar_data.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include "scihpc/vector_data.h"
#include "scihpc/runge_kutta.h"
#include "scihpc/source.h"
#include "scihpc/flux.h"
#include "scihpc/derivatives_solver.h"
#include "scihpc/wrapper.h"
#include "scihpc/bc_factory.h"

int main() {

    std::cout.precision(4);
    DataType prev_error = 0.0;
    DataType prev_h = 0.0;

    for (int cnt = 0; cnt < 8; ++cnt) {

        auto geo = structured_grid(axis{-1.0, 1.0, 16 * static_cast<int>(pow(2, cnt))});

        auto phi = wrapper(true, &geo, bc_info{PERIODIC}, bc_info{PERIODIC});
        auto vel = wrapper(false, &geo, bc_info{PERIODIC}, bc_info{PERIODIC});
        auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);
        for (int i = 0; i < phi.scalar->nx; ++i) {
            auto index = phi.scalar->index_mapping(i + 1, 1, 1);
            phi.scalar->data[index.i][index.j][index.k] = sin(pi * geo.xc[i]);
            vel.vector->x.data[index.i][index.j][index.k] = 1.0;
        }
        phi.apply_scalar_bc();
        vel.apply_nvel_bc();

        auto params = new problem_parameters;
        params->dt = 0.01 * geo.h;

        auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);

        phi.link_solvers(deri_solvers);
        phi.link_params(params);

        int tcnt = 0;
        auto begin = std::chrono::high_resolution_clock::now();
        while (tcnt * params->dt < 2.0) {
            solver.tvd_rk3(&phi, &vel, &identity_flux, &convection);
            tcnt++;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        // calculate l2 norm
        DataType error = 0.0;
        for (int i = 0; i < phi.scalar->nx; ++i) {
            auto index = phi.scalar->index_mapping(i + 1, 1, 1);
            error += pow(phi.scalar->data[index.i][index.j][index.k] - sin(pi * (geo.xc[i] - tcnt * params->dt)), 2);
        }
        error = sqrt(error / phi.scalar->nx);

        // output
        std::cout << std::setw(5) << phi.scalar->nx << " | " << std::scientific << error;
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