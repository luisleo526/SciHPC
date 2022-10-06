//
// Created by 溫晧良 on 2022/10/1.
//

#include "vtkWriter.h"

void SwapEnd(DataType &var) {
    char *varArray = reinterpret_cast<char *>(&var);
    for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
        std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

vtkWriter::vtkWriter(structured_grid *_geo, std::string _filename) {
    grid = _geo;
    filename = std::move(_filename);
}

void vtkWriter::create(unsigned int id) {
    file = std::ofstream(filename + std::to_string(id) + std::string(".vtk"), std::ios::binary);
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "vtk output" << std::endl;
    file << "BINARY" << std::endl;
    file << "DATASET STRUCTURED_POINTS" << std::endl;
    file << "DIMENSIONS " << grid->vtk_info.nx << " " << grid->vtk_info.ny << " " << grid->vtk_info.nz << std::endl;
    file << "ORIGIN " << grid->vtk_info.xstart << " " << grid->vtk_info.ystart << " " << grid->vtk_info.zstart
         << std::endl;
    file << "SPACING " << grid->dx << " " << grid->dy << " " << grid->dz << std::endl;
}

void vtkWriter::add_scalar(scalar_data *data, const std::string &name) {
    file << "CELL_DATA " << data->nx * data->ny * data->nz << std::endl;
    file << "SCALARS " << name << " double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < data->nz; ++k) {
        for (int j = 0; j < data->ny; ++j) {
            for (int i = 0; i < data->nx; ++i) {
                auto index = data->index_mapping(i + 1, j + 1, k + 1);
                auto value = data->data[index.i][index.j][index.k];
                SwapEnd(value);
                file.write(reinterpret_cast<char *>(&value), sizeof(DataType));
            }
        }
    }
}

void vtkWriter::add_vector(vector_data *data, const std::string &name) {

    file << "VECTORS " << name << " double" << std::endl;
    for (int k = 0; k < data->x.nz; ++k) {
        for (int j = 0; j < data->x.ny; ++j) {
            for (int i = 0; i < data->x.nx; ++i) {
                auto index = data->x.index_mapping(i + 1, j + 1, k + 1);
                auto value = data->x.data[index.i][index.j][index.k];
                SwapEnd(value);
                file.write(reinterpret_cast<char *>(&value), sizeof(DataType));
                value = data->y.data[index.i][index.j][index.k];
                SwapEnd(value);
                file.write(reinterpret_cast<char *>(&value), sizeof(DataType));
                value = data->z.data[index.i][index.j][index.k];
                SwapEnd(value);
                file.write(reinterpret_cast<char *>(&value), sizeof(DataType));
            }
        }
    }
}

void vtkWriter::close() {
    file.close();
}
