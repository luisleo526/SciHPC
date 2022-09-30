#include <iostream>
#include "scihpc/scalar_data.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include <cmath>
#include "scihpc/ccd.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/vector_data.h"

const DataType pi = acos(static_cast<DataType>(-1.0));

int main() {

    for (int cnt = 0; cnt < 10; ++cnt) {

        auto phi = scalar_data(16 * static_cast<int>(pow(2, cnt)));
        auto nvel = vector_data(phi.nx);
        auto geo = structured_grid(axis{0, pi, phi.nx});

        //init phi
        for (int i = 1; i <= phi.nx; ++i) {
            auto index = phi.index_mapping(i, 1, 1);
            phi.data[index.i][index.j][index.k] = sin(geo.xc[i - 1]) + 0.1*sin(2 * geo.xc[i - 1]);
            nvel.x.data[index.i][index.j][index.k] = phi.data[index.i][index.j][index.k];
        }

        periodic(&phi);
        periodic(&(nvel.x));

//        uccd_find_fx(&phi, &geo, &nvel);
        ccd_find_fx(&phi, &geo);

        DataType err_fx = 0.0, err_fxx = 0.0;
        for (int i = 1; i <= phi.nx; ++i) {
            auto index = phi.index_mapping(i, 1, 1);
            auto err1 = phi.fx[index.i][index.j][index.k] - cos(geo.xc[i - 1]) - 0.2 * cos(2 * geo.xc[i - 1]);
            auto err2 = phi.fxx[index.i][index.j][index.k] + sin(geo.xc[i - 1]) + 0.4 * sin(2 * geo.xc[i - 1]);
            err_fx += err1 * err1;
            err_fxx += err2 * err2;
        }

        std::cout << phi.nx << "," << sqrtl(err_fx / phi.nx) << "," << sqrtl(err_fxx / phi.nx) << std::endl;

    }

}