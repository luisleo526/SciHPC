//
// Created by 溫晧良 on 2022/9/28.
//

#include "structured_grid.h"

structured_grid::structured_grid(axis x_data) {
    x = new DataType[x_data.n + 1];
    y = new DataType[1];
    z = new DataType[1];
    dx = (x_data.end - x_data.start) / static_cast<DataType>(x_data.n);
    dy = static_cast<DataType>(0.0);
    dz = static_cast<DataType>(0.0);
    for (int i = 0; i < x_data.n + 1; ++i) {
        x[i] = static_cast<DataType>(i * dx);
    }
}

structured_grid::structured_grid(axis x_data, axis y_data) {
    x = new DataType[x_data.n + 1];
    y = new DataType[y_data.n + 1];
    z = new DataType[1];
    dx = (x_data.end - x_data.start) / static_cast<DataType>(x_data.n);
    dy = (y_data.end - y_data.start) / static_cast<DataType>(y_data.n);
    dz = static_cast<DataType>(0.0);
    for (int i = 0; i < x_data.n + 1; ++i) {
        x[i] = static_cast<DataType>(i * dx);
    }
    for (int j = 0; j < y_data.n + 1; ++j) {
        y[j] = static_cast<DataType>(j * dy);
    }
}

structured_grid::structured_grid(axis x_data, axis y_data, axis z_data) {
    x = new DataType[x_data.n + 1];
    y = new DataType[y_data.n + 1];
    z = new DataType[z_data.n + 1];
    dx = (x_data.end - x_data.start) / static_cast<DataType>(x_data.n);
    dy = (y_data.end - y_data.start) / static_cast<DataType>(y_data.n);
    dz = (z_data.end - z_data.start) / static_cast<DataType>(z_data.n);
    for (int i = 0; i < x_data.n + 1; ++i) {
        x[i] = static_cast<DataType>(i * dx);
    }
    for (int j = 0; j < y_data.n + 1; ++j) {
        y[j] = static_cast<DataType>(j * dy);
    }
    for (int k = 0; k < z_data.n + 1; ++k) {
        z[k] = static_cast<DataType>(k * dz);
    }
}
