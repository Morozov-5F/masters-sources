#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include <algorithm>

#include "point.hpp"

#include "foundation.hpp"

int main(int argc, const char * argv[])
{
    std::cout << "---Foundation 1---" << std::endl;
    auto f1 = Foundation(1, 45, 90, 0, 0);
    std::cout << "---Foundation 2---" << std::endl;
    auto f2 = Foundation(1, 60, 60, 0, 0);
    std::cout << "---Foundation 3---" << std::endl;
    auto f3 = Foundation(1, 120, 30, 0, 0);
    std::cout << "---Foundation 4---" << std::endl;
    auto f4 = Foundation(1, 90, 45, 1, 0.1);
    std::cout << "---Foundation 5---" << std::endl;
    auto f5 = Foundation(3, 30, 120, 0.2, 0.1);

    return 0;
}
