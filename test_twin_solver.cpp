//
// Created by 溫晧良 on 2022/9/30.
//

#include "scihpc/global.h"
#include "scihpc/matrix_solver.h"
#include <random>
#include <iostream>
#include <ctime>

using namespace std;

const int n = 10;

int main() {

    default_random_engine generator(time(NULL));
    uniform_real_distribution<DataType> distribution(-1.0, 1.0);

    DataType **a, **b, **aa, **bb, *s, *ss, *fx, *fxx;
    a = new DataType *[3];
    b = new DataType *[3];
    aa = new DataType *[3];
    bb = new DataType *[3];
    s = new DataType[n];
    ss = new DataType[n];
    fx = new DataType[n];
    fxx = new DataType[n];

    for (int i = 0; i < 3; ++i) {
        a[i] = new DataType[n];
        b[i] = new DataType[n];
        aa[i] = new DataType[n];
        bb[i] = new DataType[n];
    }

    for (int i = 0; i < n; ++i) {
        fx[i] = distribution(generator);
        fxx[i] = distribution(generator);
        for (int j = 0; j < 3; ++j) {
            a[j][i] = distribution(generator);
            b[j][i] = distribution(generator);
            aa[j][i] = distribution(generator);
            bb[j][i] = distribution(generator);
        }
    }

    for (int i = 0; i < n; ++i) {
        s[i] = a[1][i] * fx[i] + b[1][i] * fxx[i];
        ss[i] = aa[1][i] * fx[i] + bb[1][i] * fxx[i];
        if (i > 0) {
            s[i] += a[0][i] * fx[i - 1] + b[0][i] * fxx[i - 1];
            ss[i] += aa[0][i] * fx[i - 1] + bb[0][i] * fxx[i - 1];
        }
        if (i < n - 1) {
            s[i] += a[2][i] * fx[i + 1] + b[2][i] * fxx[i + 1];
            ss[i] += aa[2][i] * fx[i + 1] + bb[2][i] * fxx[i + 1];
        }
    }

    twin_dec(a, b, aa, bb, n);
    twin_bks(a, b, aa, bb, s, ss, n);

    DataType err=0.0;
    for (int i = 0; i < n; ++i) {
        err += abs(s[i] - fx[i]) + abs(ss[i] - fxx[i]);
    }

    cout << err / n << endl;

}
