#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <ctime>

#include "scihpc/wrapper.h"
#include "scihpc/structured_grid.h"
#include "scihpc/flux.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/wrapper_func.h"

int main() {

    int level_max = 4;
    DataType prev_error = 0.0;
    DataType prev_h = 0.0;

    for (int level = 0; level < level_max; ++level) {

        auto n = 32 * static_cast<int>(pow(2, level));
        auto phi = wrapper(new scalar_data(n, n, n));
        auto vel = wrapper(new vector_data(phi.scalar->nx, phi.scalar->ny));
        auto geo = structured_grid(axis{-1.0, 1.0, phi.scalar->nx},
                                   axis{-1.0, 1.0, phi.scalar->ny},
                                   axis{-1.0, 1.0, phi.scalar->nz});
        auto param = new problem_parameters{};
        auto deri_solvers = derivatives_solver_alloc(phi.scalar, &geo);
        auto dummy = dummy_data_alloc(phi.scalar);
        phi.link_params(param);
        phi.link_solvers(deri_solvers);
        phi.link_dummy(dummy);

        param->ls_width = 1.5 * geo.h;

        for (int i = 0; i < phi.scalar->nx; ++i) {
            for (int j = 0; j < phi.scalar->ny; ++j) {
                for (int k = 0; k < phi.scalar->nz; ++k) {
                    auto index = phi.scalar->index_mapping(i + 1, j + 1, k + 1);
                    phi.scalar->data[index.i][index.j][index.k] =
                            -sqrt(geo.xc[i] * geo.xc[i] + geo.yc[j] * geo.yc[j] + geo.zc[k] * geo.zc[k]) + 0.5;
                }
            }
        }
        zero_order_extrapolation(phi.scalar);

        identity_flux(phi.scalar, vel.vector);
        find_curvature(&phi);
        find_delta(&phi);

        auto error = 0.0;
        int cnt = 0;
        for (int i = 0; i < phi.scalar->nx; ++i) {
            for (int j = 0; j < phi.scalar->ny; ++j) {
                for (int k = 0; k < phi.scalar->nz; ++k) {
                    auto index = phi.scalar->index_mapping(i + 1, j + 1, k + 1);
                    if (phi.dummy->delta[index.i][index.j][index.k] > 0.0) {
                        error += pow(phi.dummy->curvature[index.i][index.j][index.k] - 2.0, 2);
                        cnt++;
                    }
                }
            }
        }

        error = sqrt(error / cnt);

        // output
        std::cout << std::setw(5) << phi.scalar->nx << " | " << std::scientific << error;
        if (level > 0) {
            std::cout << " | " << std::fixed << std::setw(7)
                      << (log10(error) - log10(prev_error)) / (log10(geo.dx) - log10(prev_h));
        } else {
            std::cout << " | " << std::setw(7) << "-";
        }
        std::cout << std::endl;

        prev_error = error;
        prev_h = geo.dx;

    }

    return 0;
}
