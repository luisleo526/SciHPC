#include <iostream>
#include "Data.h"
int main() {
    auto f = Data(10);
    std::cout << f.index_mapping(f.Nx()+3, 1, 1).i << f.ndim ;
    return 0;
}
