//
// Created by leo on 10/2/22.
//

#include "ccd_base.h"

ccd_base::ccd_base(int _n, DataType _h) {

    h = _h;
    n = _n;

    a = new DataType *[3];
    b = new DataType *[3];
    aa = new DataType *[3];
    bb = new DataType *[3];
    s = new DataType[n];
    ss = new DataType[n];

    for (int i = 0; i < 3; ++i) {
        a[i] = new DataType[n];
        b[i] = new DataType[n];
        aa[i] = new DataType[n];
        bb[i] = new DataType[n];
        for (int j = 0; j < n; ++j) {
            a[i][j] = 0;
            b[i][j] = 0;
            aa[i][j] = 0;
            bb[i][j] = 0;
        }
    }

    for (int i = 0; i < n; ++i) {
        s[i] = 0;
        ss[i] = 0;
    }

    if (h > 0.0) {
        init_coefficients();
    }

}

void ccd_base::init_coefficients() const {

    for (int i = 0; i < n; ++i) {
        a[0][i] = 7.0 / 16.0;
        a[1][i] = 1.0;
        a[2][i] = 7.0 / 16.0;

        b[0][i] = h / 16.0;
        b[1][i] = 0.0;
        b[2][i] = -h / 16.0;

        aa[0][i] = -9.0 / 8.0 / h;
        aa[1][i] = 0.0;
        aa[2][i] = 9.0 / 8.0 / h;

        bb[0][i] = -1.0 / 8.0;
        bb[1][i] = 1.0;
        bb[2][i] = -1.0 / 8.0;
    }

    a[0][0] = 0.0;
    b[0][0] = 0.0;
    aa[0][0] = 0.0;
    bb[0][0] = 0.0;

    a[1][0] = 1.0;
    b[1][0] = 0.0;
    aa[1][0] = 0.0;
    bb[1][0] = 1.0;

    a[2][0] = 2.0;
    b[2][0] = -h;
    aa[2][0] = -2.5 / h;
    bb[2][0] = 8.5;

    a[0][n - 1] = 2.0;
    b[0][n - 1] = h;
    aa[0][n - 1] = 2.5 / h;
    bb[0][n - 1] = 8.5;

    a[1][n - 1] = 1.0;
    b[1][n - 1] = 0.0;
    aa[1][n - 1] = 0.0;
    bb[1][n - 1] = 1.0;

    a[2][n - 1] = 0.0;
    b[2][n - 1] = 0.0;
    aa[2][n - 1] = 0.0;
    bb[2][n - 1] = 0.0;

    twin_dec(a, b, aa, bb, n);

}
