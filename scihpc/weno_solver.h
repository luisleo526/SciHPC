//
// Created by leo on 10/3/22.
//

#ifndef SCIHPC_WENO_SOLVER_H
#define SCIHPC_WENO_SOLVER_H

#include "global.h"
#include "structured_grid.h"
#include "scalar_data.h"
#include "vector_data.h"

struct weno_data {
    DataType x0, x1, x2;
};

class weno_solver {
private:
    static weno_data indicators_p(DataType a, DataType b, DataType c, DataType d, DataType e);

    static weno_data indicators_m(DataType a, DataType b, DataType c, DataType d, DataType e);

    static weno_data interpolation_p(DataType a, DataType b, DataType c, DataType d, DataType e);

    static weno_data interpolation_m(DataType a, DataType b, DataType c, DataType d, DataType e);

    static weno_data wenojs_ceofficients(DataType b0, DataType b1, DataType b2);

public:
    weno_solver(scalar_data *f, structured_grid *geo);

    DataType dx, dy, dz;
    DataType ***fp, ***fm, ***fh;

    void weno5_flux_x(scalar_data *f);

    void weno5_flux_y(scalar_data *f);

    void weno5_flux_z(scalar_data *f);

    void weno5_find_fx(scalar_data *f, vector_data *vel);

    void weno5_find_fy(scalar_data *f, vector_data *vel);

    void weno5_find_fz(scalar_data *f, vector_data *vel);

    void weno5_find_derivatives(scalar_data *f, vector_data *vel);
};


#endif //SCIHPC_WENO_SOLVER_H
