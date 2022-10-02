#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <omp.h>

#include "scalar_data.h"
#include "scihpc/global.h"


DataType ***test();

int main() {

    do {
        auto arr = test();
        delete3d(arr, 100, 100);
        std::cout << arr << std::endl;
    } while (true);

    return 0;
}
