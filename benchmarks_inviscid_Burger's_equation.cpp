#include <iostream>
#include <iomanip>
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
#include "scihpc/derivatives_solver.h"
#include "scihpc/wrapper.h"

int main() {

    std::cout.precision(4);
    DataType prev_error = 0.0;
    DataType prev_h = 0.0;

    const int cnt_max = 5;

    for (int cnt = 0; cnt <= cnt_max; ++cnt) {

        DataType shock_formation = 2.0;

        auto phi = wrapper(new scalar_data(16 * static_cast<int>(pow(2, cnt))));
        auto vel = wrapper(new vector_data(phi.scalar->nx));
        auto geo = structured_grid(axis{0.0, 1.0, phi.scalar->nx});
        auto solver = runge_kutta(phi.scalar->Nx, phi.scalar->Ny, phi.scalar->Nz);
        for (int i = 0; i < phi.scalar->nx; ++i) {
            auto index = phi.scalar->index_mapping(i + 1, 1, 1);
            phi.scalar->data[index.i][index.j][index.k] = sin(2.0 * pi * geo.xc[i]) / (2.0 * pi * shock_formation);
        }

        periodic(phi.scalar);

        auto params = new problem_parameters;
        params->dt = 0.01 * geo.h;

        auto deri_solvers = new solvers_ptr;
        deri_solvers->ccd = new ccd_solver(phi.scalar, &geo);
        deri_solvers->uccd = new uccd_solver(phi.scalar, &geo);

        phi.link_solvers(deri_solvers);
        phi.link_params(params);

        int tcnt = 0;
        auto begin = std::chrono::high_resolution_clock::now();
        while (tcnt * params->dt < shock_formation * 0.75) {
            for (int i = 0; i < phi.scalar->Nx; ++i) {
                for (int j = 0; j < phi.scalar->Ny; ++j) {
                    for (int k = 0; k < phi.scalar->Nz; ++k) {
                        vel.vector->x.data[i][j][k] = phi.scalar->data[i][j][k];
                    }
                }
            }
            solver.tvd_rk3(&phi, &vel, &geo, &burgers_flux, &periodic, &Hamilton_Jacobi);
            tcnt++;
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        DataType phi_exact[phi.scalar->nx];
        for (int i = 0; i < phi.scalar->nx; ++i) {
            phi_exact[i] = sin(2.0 * pi * geo.xc[i]) / (2.0 * pi * shock_formation);
        }

        DataType error = 0.0;
        do {
            error = 0.0;
            for (int i = 0; i < phi.scalar->nx; ++i) {
                auto new_phi = sin(2.0 * pi * (geo.xc[i] - phi_exact[i] * tcnt * params->dt)) / (2.0 * pi * shock_formation);
                error += pow(new_phi - phi_exact[i], 2);
                phi_exact[i] = new_phi;
            }
            error = sqrt(error / phi.scalar->nx);
        } while (error > 1.0e-14);

        // calculate l2 norm
        error = 0.0;
        for (int i = 0; i < phi.scalar->nx; ++i) {
            auto index = phi.scalar->index_mapping(i + 1, 1, 1);
            error += pow(phi.scalar->data[index.i][index.j][index.k] - phi_exact[i], 2);
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