#include "calculate_sse.h"

void Sim_2D::calculate_sse(const int idx) {
    this->IDX = idx;
    generate_pattern();
}