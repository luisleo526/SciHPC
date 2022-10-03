//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_STRUCTURED_GRID_H
#define SCIHPC_STRUCTURED_GRID_H

#include "global.h"

class structured_grid {
public:
    explicit structured_grid(axis x_data);

    structured_grid(axis x_data, axis y_data);

    structured_grid(axis x_data, axis y_data, axis z_data);

    DataType *x, *xc;
    DataType *y, *yc;
    DataType *z, *zc;
    DataType dx, dy, dz, h, dv;
    struct {
        DataType xstart, ystart, zstart;
        int nx, ny, nz;
    } vtk_info;
};


#endif //SCIHPC_STRUCTURED_GRID_H
