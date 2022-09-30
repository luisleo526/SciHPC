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


void twin_dec(DataType **a, DataType **b, DataType **aa, DataType **bb, int n) {

    a[0][0] = static_cast<DataType>(0.0);
    aa[0][0] = static_cast<DataType>(0.0);
    b[0][0] = static_cast<DataType>(1.0);
    bb[0][0] = static_cast<DataType>(1.0);

    a[2][n - 1] = static_cast<DataType>(0.0);
    aa[2][n - 1] = static_cast<DataType>(0.0);
    b[2][n - 1] = static_cast<DataType>(1.0);
    bb[2][n - 1] = static_cast<DataType>(1.0);

    DataType den, sc1, sc2;

    for (int i = 1; i < n; ++i) {
        den = a[1][i - 1] * bb[1][i - 1] - aa[1][i - 1] * b[1][i - 1];

        sc1 = -(aa[1][i - 1] * b[0][i] - a[0][i] * bb[1][i - 1]) / den;
        sc2 = (b[0][i] * a[1][i - 1] - b[1][i - 1] * a[0][i]) / den;

        a[0][i] = sc1;
        a[1][i] = a[1][i] - (sc1 * a[2][i - 1] + sc2 * aa[2][i - 1]);
        b[0][i] = sc2;
        b[1][i] = b[1][i] - (sc1 * b[2][i - 1] + sc2 * bb[2][i - 1]);

        sc1 = -(aa[1][i - 1] * bb[0][i] - aa[0][i] * bb[1][i - 1]) / den;
        sc2 = (bb[0][i] * a[1][i - 1] - b[1][i - 1] * aa[0][i]) / den;

        aa[0][i] = sc1;
        aa[1][i] = aa[1][i] - (sc1 * a[2][i - 1] + sc2 * aa[2][i - 1]);
        bb[0][i] = sc2;
        bb[1][i] = bb[1][i] - (sc1 * b[2][i - 1] + sc2 * bb[2][i - 1]);

    }
}

void twin_bks(DataType **a, DataType **b, DataType **aa, DataType **bb, DataType *s, DataType *ss, int n) {

    for (int i = 1; i < n; ++i) {
        s[i] = s[i] - (a[0][i] * s[i - 1] + b[0][i] * ss[i - 1]);
        ss[i] = ss[i] - (aa[0][i] * s[i - 1] + bb[0][i] * ss[i - 1]);
    }

    DataType den, sols, solss;

    den = a[1][n - 1] * bb[1][n - 1] - aa[1][n - 1] * b[1][n - 1];
    sols = -(b[1][n - 1] * ss[n - 1] - bb[1][n - 1] * s[n - 1]) / den;
    solss = (a[1][n - 1] * ss[n - 1] - aa[1][n - 1] * s[n - 1]) / den;
    s[n - 1] = sols;
    ss[n - 1] = solss;

    for (int i = n - 2; i > -1; --i) {

        s[i] = s[i] - (a[2][i] * s[i + 1] + b[2][i] * ss[i + 1]);
        ss[i] = ss[i] - (aa[2][i] * s[i + 1] + bb[2][i] * ss[i + 1]);

        den = a[1][i] * bb[1][i] - aa[1][i] * b[1][i];
        sols = -(b[1][i] * ss[i] - bb[1][i] * s[i]) / den;
        solss = (a[1][i] * ss[i] - aa[1][i] * s[i]) / den;

        s[i] = sols;
        ss[i] = solss;
    }
}
