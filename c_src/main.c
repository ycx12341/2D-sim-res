#include <stdio.h>
#include <math.h>

#include "util/math/math_ext.h"
#include "util/matrix/matrix.h"
#include "algorithm/scc.h"
#include "util/collection/collection.h"

int main() {
    set_seed(SEED);

    sse_pars_t pars;
    init_pars(&pars);
    for (int i = 4; i < 5; ++i) {
        calculate_scc(&pars, i);
    }

//    double a[5] = {0.1, 0.1, 0.2, 0.3, 0.4};
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));
//    printf("%d\n", sample_prob1(5, a));

//    arraylist_t *l = new_arraylist();
//    for (int    i  = 0; i < 10; ++i) {
//        node_t n;
//        n._int = i;
//        arraylist_append(l, n);
//    }
//    for (int    i  = 0; i < l->size; ++i) {
//        printf("%d, ", arraylist_get(l, i)._int);
//    }
//
//    int d[3] = {2, 3, 2};
//
//    arraylist_remove_many(l, 3, d);
//    printf("\n");
//    for (int i = 0; i < l->size; ++i) {
//        printf("%d, ", arraylist_get(l, i)._int);
//    }
//    printf("\n");
//
//    double   a[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
//    int      an    = double_array_delete_many(10, a, 3, d);
//    for (int i     = 0; i < an; ++i) {
//        printf("%d, ", (int) a[i]);
//    }

    return 0;
}
