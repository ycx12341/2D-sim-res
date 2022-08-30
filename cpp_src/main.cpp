#include <iostream>

#include "algo/sim2d.h"

int main() {
    set_seed(SEED);
    auto scc = *Sim_2D_Factory::SCC(DEFAULT_N_DIMS);

//    Matrix<int> *a = new MatrixDy<int>(3, 4, 1);
//    a->print("%d ");
//    (*a)(1, 1) = 2;
//    a->print("%d ");

    for (int i = 4; i < 5; ++i) {
//    for (int i = 0; i < 10000; ++i) {
        std::cout << i << std::endl;
        scc.calculate_sse(i);
    }

    return 0;
}
