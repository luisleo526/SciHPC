#include <iostream>
#include <Eigen/Core>
#include "scihpc/global.h"

using namespace Eigen;
int main()
{

    Matrix<double, Dynamic, Dynamic> a;

    a = Matrix<double, Dynamic, Dynamic>(3, 3);

    std::cout << a << std::endl;
}
