#include <iostream>

#include "algo/calculate_sse.h"
#include "collection/matrix.h"

int main() {
    set_seed(SEED);
    Sim_2D sim2D = Sim_2D::sim_scc(DEFAULT_N_DIMS);

    for (int i = 4; i < 5; ++i) {
        sim2D.calculate_sse(i);
    }

    return 0;
}
