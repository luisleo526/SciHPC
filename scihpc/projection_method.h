//
// Created by 溫晧良 on 2022/10/5.
//

#ifndef SCIHPC_PROJECTION_METHOD_H
#define SCIHPC_PROJECTION_METHOD_H

#include "wrapper.h"

class projection_method {
public:
    projection_method(scalar_data* f);
    DataType ***u_src, ***v_src, ***w_src;
    DataType ***u_src_old, ***v_src_old, ***w_src_old;
};


#endif //SCIHPC_PROJECTION_METHOD_H
