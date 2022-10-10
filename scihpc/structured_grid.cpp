//
// Created by 溫晧良 on 2022/9/28.
//

#include "structured_grid.h"

structured_grid::structured_grid(axis x_data) {
    x = new DataType[x_data.n + 1];
    y = new DataType[1];
    z = new DataType[1];

    xc = new DataType[x_data.n];
    yc = new DataType[1];
    zc = new DataType[1];

    dx = (x_data.end - x_data.start) / static_cast<DataType>(x_data.n);
    dy = static_cast<DataType>(0.0);
    dz = static_cast<DataType>(0.0);

    h = dx;
    dv = dx;

    for (int i = 0; i < x_data.n + 1; ++i) {
        x[i] = x_data.start + static_cast<DataType>(i * dx);
        if (i < x_data.n) {
            xc[i] = x[i] + static_cast<DataType>(0.5 * dx);
        }
    }

    vtk_info.nx = x_data.n + 1;
    vtk_info.ny = 1;
    vtk_info.nz = 1;

    vtk_info.xstart = x_data.start;
    vtk_info.ystart = 0.0;
    vtk_info.zstart = 0.0;

    ndim = 1;
    x_axis = x_data;

}

structured_grid::structured_grid(axis x_data, axis y_data) {
    x = new DataType[x_data.n + 1];
    y = new DataType[y_data.n + 1];
    z = new DataType[1];

    xc = new DataType[x_data.n];
    yc = new DataType[y_data.n];
    zc = new DataType[1];

    dx = (x_data.end - x_data.start) / static_cast<DataType>(x_data.n);
    dy = (y_data.end - y_data.start) / static_cast<DataType>(y_data.n);
    dz = static_cast<DataType>(0.0);

    h = fmin(dx, dy);
    dv = dx * dy;

    for (int i = 0; i < x_data.n + 1; ++i) {
        x[i] = x_data.start + static_cast<DataType>(i * dx);
        if (i < x_data.n) {
            xc[i] = x[i] + static_cast<DataType>(0.5 * dx);
        }
    }
    for (int j = 0; j < y_data.n + 1; ++j) {
        y[j] = y_data.start + static_cast<DataType>(j * dy);
        if (j < y_data.n) {
            yc[j] = y[j] + static_cast<DataType>(0.5 * dy);
        }
    }

    vtk_info.nx = x_data.n + 1;
    vtk_info.ny = y_data.n + 1;
    vtk_info.nz = 1;

    vtk_info.xstart = x_data.start;
    vtk_info.ystart = y_data.start;
    vtk_info.zstart = 0.0;

    ndim = 2;
    x_axis = x_data;
    y_axis = y_data;

}

structured_grid::structured_grid(axis x_data, axis y_data, axis z_data) {
    x = new DataType[x_data.n + 1];
    y = new DataType[y_data.n + 1];
    z = new DataType[z_data.n + 1];

    xc = new DataType[x_data.n];
    yc = new DataType[y_data.n];
    zc = new DataType[z_data.n];

    dx = (x_data.end - x_data.start) / static_cast<DataType>(x_data.n);
    dy = (y_data.end - y_data.start) / static_cast<DataType>(y_data.n);
    dz = (z_data.end - z_data.start) / static_cast<DataType>(z_data.n);

    h = fmin(fmin(dx, dy), dz);
    dv = dx * dy * dz;

    for (int i = 0; i < x_data.n + 1; ++i) {
        x[i] = x_data.start + static_cast<DataType>(i * dx);
        if (i < x_data.n) {
            xc[i] = x[i] + static_cast<DataType>(0.5 * dx);
        }
    }
    for (int j = 0; j < y_data.n + 1; ++j) {
        y[j] = y_data.start + static_cast<DataType>(j * dy);
        if (j < y_data.n) {
            yc[j] = y[j] + static_cast<DataType>(0.5 * dy);
        }
    }
    for (int k = 0; k < z_data.n + 1; ++k) {
        z[k] = z_data.start + static_cast<DataType>(k * dz);
        if (k < z_data.n) {
            zc[k] = z[k] + static_cast<DataType>(0.5 * dz);
        }
    }

    vtk_info.nx = x_data.n + 1;
    vtk_info.ny = y_data.n + 1;
    vtk_info.nz = z_data.n + 1;

    vtk_info.xstart = x_data.start;
    vtk_info.ystart = y_data.start;
    vtk_info.zstart = z_data.start;

    ndim = 3;
    x_axis = x_data;
    y_axis = y_data;
    z_axis = z_data;

}
