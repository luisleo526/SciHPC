//
// Created by 溫晧良 on 2022/10/11.
//

#ifndef SCIHPC_SECOND_ORDER_SOLVER_H
#define SCIHPC_SECOND_ORDER_SOLVER_H

#include "scalar_data.h"
#include "vector_data.h"

class second_order_solver {
public:
    DataType dx, dy, dz;

    void find_fx(scalar_data *f) const;

    void find_fy(scalar_data *f) const;

    void find_fz(scalar_data *f) const;

    void find_fx(scalar_data *f, vector_data *vel) const;

    void find_fy(scalar_data *f, vector_data *vel) const;

    void find_fz(scalar_data *f, vector_data *vel) const;

    void find_derivatives(scalar_data *f) const;

    void find_derivatives(scalar_data *f, vector_data *vel) const;

    second_order_solver(DataType _dx, DataType _dy, DataType _dz);
};


#endif //SCIHPC_SECOND_ORDER_SOLVER_H
