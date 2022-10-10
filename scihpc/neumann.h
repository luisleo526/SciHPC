//
// Created by 溫晧良 on 2022/10/10.
//

#ifndef SCIHPC_NEUMANN_H
#define SCIHPC_NEUMANN_H

#include "scalar_data.h"

void Neumann_node_xr(scalar_data *f, DataType df, DataType h);
void Neumann_node_xl(scalar_data *f, DataType df, DataType h);
void Neumann_node_yr(scalar_data *f, DataType df, DataType h);
void Neumann_node_yl(scalar_data *f, DataType df, DataType h);
void Neumann_node_zr(scalar_data *f, DataType df, DataType h);
void Neumann_node_zl(scalar_data *f, DataType df, DataType h);

void Neumann_face_xr(scalar_data *f, DataType df, DataType h);
void Neumann_face_xl(scalar_data *f, DataType df, DataType h);
void Neumann_face_yr(scalar_data *f, DataType df, DataType h);
void Neumann_face_yl(scalar_data *f, DataType df, DataType h);
void Neumann_face_zr(scalar_data *f, DataType df, DataType h);
void Neumann_face_zl(scalar_data *f, DataType df, DataType h);

#endif //SCIHPC_NEUMANN_H
