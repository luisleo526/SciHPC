//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_RUNGE_KUTTA_H
#define SCIHPC_RUNGE_KUTTA_H

#include "global.h"
#include "scalar_data.h"
#include "vector_data.h"
#include "structured_grid.h"

class runge_kutta {
public:
    runge_kutta(int nx, int ny, int nz);

    DataType ***s1, ***s2, ***s3;

    void tvd_rk3(DataType dt, scalar_data *f, vector_data *vel, structured_grid *geo,
                 void(*flux)(scalar_data *, vector_data *),
                 void (*bc)(scalar_data *),
                 void (*rhs)(scalar_data *, vector_data *, structured_grid *, DataType ***,
                             void (*)(scalar_data *, vector_data *)));

};


#endif //SCIHPC_RUNGE_KUTTA_H
