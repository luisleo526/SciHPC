//
// Created by 溫晧良 on 2022/9/28.
//

#include <iostream>
#include "data_array.h"
#include "matrix_solver.h"
#include <random>
#include <ctime>
#include <cmath>

const int n = 100;

using namespace std;

int main() {

    DataType a[n], b[n], c[n], d[n], x[n];

    default_random_engine generator(time(NULL));
    uniform_real_distribution<DataType> unif(-1.0, 1.0);
    for (int i = 0; i < n; ++i) {
        a[i] = unif(generator);
        b[i] = -2.0;
        c[i] = unif(generator);
        x[i] = unif(generator);
    }

    for (int i = 0; i < n; ++i) {
        d[i] = b[i] * x[i];
        if (i > 0) {
            d[i] += a[i] * x[i - 1];
        }

        if (i < n - 1) {
            d[i] += c[i] * x[i + 1];
        }
    }

    tdma(a, b, c, d, n);

    DataType error = 0.0;
    for (int i = 0; i < n; ++i) {
        error += (d[i] - x[i]) * (d[i] - x[i]);
    }
    cout << sqrtl(error / static_cast<DataType>(n)) <<endl;

    return 0;

}