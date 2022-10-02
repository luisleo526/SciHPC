//
// Created by leo on 10/2/22.
//

#include "uccd_base.h"

uccd_base::uccd_base(int _n, DataType _h) {
    h = _h;
    n = _n;
    upwind = new ccd_base(_n, _h);
    downwind = new ccd_base(_n, _h);
    if (h > 0.0) {
        init_coefficients();
    }
}

void uccd_base::init_coefficients() const {

    DataType a1, b1, b2, b3;

    a1 = 0.875;
    b1 = 0.1251282341599089;
    b2 = -0.2487176584009104;
    b3 = 0.0001282341599089;

    for (int i = 0; i < n; ++i) {

        upwind->a[0][i] = a1;
        upwind->a[1][i] = 1.0;
        upwind->a[2][i] = 0.0;
        downwind->a[0][i] = 0.0;
        downwind->a[1][i] = 1.0;
        downwind->a[2][i] = a1;

        upwind->b[0][i] = b1 * h;
        upwind->b[1][i] = b2 * h;
        upwind->b[2][i] = b3 * h;
        downwind->b[0][i] = -b3 * h;
        downwind->b[1][i] = -b2 * h;
        downwind->b[2][i] = -b1 * h;

        upwind->aa[0][i] = -9.0 / 8.0 / h;
        upwind->aa[1][i] = 0.0;
        upwind->aa[2][i] = 9.0 / 8.0 / h;
        downwind->aa[0][i] = -9.0 / 8.0 / h;
        downwind->aa[1][i] = 0.0;
        downwind->aa[2][i] = 9.0 / 8.0 / h;

        upwind->bb[0][i] = -1.0 / 8.0;
        upwind->bb[1][i] = 1.0;
        upwind->bb[2][i] = -1.0 / 8.0;
        downwind->bb[0][i] = -1.0 / 8.0;
        downwind->bb[1][i] = 1.0;
        downwind->bb[2][i] = -1.0 / 8.0;

    }

    upwind->a[0][0] = 0.0;
    upwind->a[1][0] = 1.0;
    upwind->a[2][0] = 2.0;
    upwind->b[0][0] = 0.0;
    upwind->b[1][0] = 0.0;
    upwind->b[2][0] = -h;
    upwind->aa[0][0] = 0.0;
    upwind->aa[1][0] = 0.0;
    upwind->aa[2][0] = -2.5 / h;
    upwind->bb[0][0] = 0.0;
    upwind->bb[1][0] = 1.0;
    upwind->bb[2][0] = 8.5;
    upwind->a[0][n - 1] = 2.0;
    upwind->a[1][n - 1] = 1.0;
    upwind->a[2][n - 1] = 0.0;
    upwind->b[0][n - 1] = h;
    upwind->b[1][n - 1] = 0.0;
    upwind->b[2][n - 1] = 0.0;
    upwind->aa[0][n - 1] = 2.5 / h;
    upwind->aa[1][n - 1] = 0.0;
    upwind->aa[2][n - 1] = 0.0;
    upwind->bb[0][n - 1] = 8.5;
    upwind->bb[1][n - 1] = 1.0;
    upwind->bb[2][n - 1] = 0.0;

    downwind->a[0][0] = 0.0;
    downwind->a[1][0] = 1.0;
    downwind->a[2][0] = 2.0;
    downwind->b[0][0] = 0.0;
    downwind->b[1][0] = 0.0;
    downwind->b[2][0] = -h;
    downwind->aa[0][0] = 0.0;
    downwind->aa[1][0] = 0.0;
    downwind->aa[2][0] = -2.5 / h;
    downwind->bb[0][0] = 0.0;
    downwind->bb[1][0] = 1.0;
    downwind->bb[2][0] = 8.5;
    downwind->a[0][n - 1] = 2.0;
    downwind->a[1][n - 1] = 1.0;
    downwind->a[2][n - 1] = 0.0;
    downwind->b[0][n - 1] = h;
    downwind->b[1][n - 1] = 0.0;
    downwind->b[2][n - 1] = 0.0;
    downwind->aa[0][n - 1] = 2.5 / h;
    downwind->aa[1][n - 1] = 0.0;
    downwind->aa[2][n - 1] = 0.0;
    downwind->bb[0][n - 1] = 8.5;
    downwind->bb[1][n - 1] = 1.0;
    downwind->bb[2][n - 1] = 0.0;

    twin_dec(upwind->a, upwind->b, upwind->aa, upwind->bb, upwind->n);
    twin_dec(downwind->a, downwind->b, downwind->aa, downwind->bb, downwind->n);

}
