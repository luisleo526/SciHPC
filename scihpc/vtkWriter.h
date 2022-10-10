//
// Created by 溫晧良 on 2022/10/1.
//

#ifndef SCIHPC_VTKWRITER_H
#define SCIHPC_VTKWRITER_H

#include <string>
#include <fstream>
#include "global.h"
#include "structured_grid.h"
#include "vector_data.h"
#include "scalar_data.h"

class vtkWriter {

public:
    vtkWriter(structured_grid *_geo, std::string _filename);

    structured_grid *grid;
    std::string filename;
    std::ofstream file;

    void create(unsigned int id);

    void add_scalar(scalar_data *data, const std::string &name);

    void add_scalar(DataType ***_data, scalar_data *data, const std::string &name);

    void add_vector(vector_data *data, const std::string &name);

    void add_vector(DataType ***u, DataType ***v, DataType ***w, vector_data *data, const std::string &name);

    void close();
};


#endif //SCIHPC_VTKWRITER_H
