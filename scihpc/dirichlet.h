//
// Created by 溫晧良 on 2022/10/9.
//

#ifndef SCIHPC_DIRICHLET_H
#define SCIHPC_DIRICHLET_H

#include "vector_data.h"
#include "global.h"

void Dirichlet_node_xr(scalar_data *f, DataType U);
void Dirichlet_node_xl(scalar_data *f, DataType U);
void Dirichlet_node_yr(scalar_data *f, DataType U);
void Dirichlet_node_yl(scalar_data *f, DataType U);
void Dirichlet_node_zr(scalar_data *f, DataType U);
void Dirichlet_node_zl(scalar_data *f, DataType U);

void Dirichlet_face_xr(scalar_data *f, DataType U);
void Dirichlet_face_xl(scalar_data *f, DataType U);
void Dirichlet_face_yr(scalar_data *f, DataType U);
void Dirichlet_face_yl(scalar_data *f, DataType U);
void Dirichlet_face_zr(scalar_data *f, DataType U);
void Dirichlet_face_zl(scalar_data *f, DataType U);

#endif //SCIHPC_DIRICHLET_H
