//
// Created by 溫晧良 on 2022/9/28.
//

#ifndef SCIHPC_GLOBAL_H
#define SCIHPC_GLOBAL_H
typedef long double DataType;

struct indices {
    long i, j, k;
};

struct axis {
    float start, end;
    int n;
};
#endif //SCIHPC_GLOBAL_H
