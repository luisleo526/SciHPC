//
// Created by 溫晧良 on 2022/9/28.
//

#include "matrix_solver.h"

void tdma(DataType *a, DataType *b, DataType *c, DataType *x, int n) {
    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    for (int i = 1; i < n; ++i) {
        c[i] = c[i] / (b[i] - a[i] * c[i - 1]);
        x[i] = (x[i] - a[i] * x[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }
    for (int i = n - 2; i > -1; --i) {
        x[i] = x[i] - c[i] * x[i + 1];
    }
}
