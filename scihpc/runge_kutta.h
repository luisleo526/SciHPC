//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_RUNGE_KUTTA_H
#define SCIHPC_RUNGE_KUTTA_H

#include "global.h"
#include "scalar_data.h"
#include "wrapper.h"

class runge_kutta {
public:
    runge_kutta(int nx, int ny, int nz);

    DataType ***s1, ***s2, ***s3;

    void tvd_rk3(wrapper *f, wrapper *vel, void (*flux)(scalar_data *, vector_data *),
                 void (*rhs)(wrapper *, wrapper *, DataType ***, void (*)(scalar_data *, vector_data *))) const;

    void euler(wrapper *f, wrapper *vel, void (*flux)(scalar_data *, vector_data *),
               void (*rhs)(wrapper *, wrapper *, DataType ***, void (*)(scalar_data *, vector_data *))) const;
};


#endif //SCIHPC_RUNGE_KUTTA_H
