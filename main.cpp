#include <iostream>
#include "data_array.h"

int main() {
    auto f = data_array(10, 10, 10);
    for (int i = 0; i < f.Nx; ++i) {
        for (int j = 0; j < f.Ny; ++j) {
            for (int k = 0; k < f.Nz; ++k) {
                f.data[i][j][k] = static_cast<DataType>(i + 10 * j + 100 * k);
            }
        }
    }
    auto ptr = &f.data[5][4][5];
    std::cout << *ptr << "," << *(ptr+f.dis.i)<< "," << *(ptr+f.dis.j)<< "," << *(ptr+f.dis.k) << std::endl;
    std::cout << sizeof(ptr);
    return 0;
}
