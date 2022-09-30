#include <iostream>
#include "scihpc/scalar_data.h"
#include "scihpc/boundary_condition.h"
#include "scihpc/global.h"
#include "scihpc/structured_grid.h"
#include <cmath>
#include <fstream>

const DataType pi = acos(static_cast<DataType>(-1.0));

int main() {

    auto phi = scalar_data(100);
    auto geo = structured_grid(axis{-1.0, 1.0, phi.nx});

    //init phi
    for (int i = 0; i < phi.nx; ++i) {
        indices index = phi.index_mapping(i + 1, 1, 1);
        phi.data[index.i][index.j][index.k] = sinl(pi * geo.xc[i]);
    }

    periodic(&phi, &geo);

    std::ofstream file;
    file.open("./test-periodic.txt");

    for (int i = 0; i < phi.Nx; ++i) {
//        indices index = phi.index_mapping(i, 1, 1);
//        file << (geo.x[i] + geo.x[i - 1]) / 2.0 << "," << phi.data[index.i][index.j][index.k] << std::endl;
        file << phi.data[i][0][0] << std::endl;
    }


    file.close();

}