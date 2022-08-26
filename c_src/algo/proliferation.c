#include "algo.h"
#include "../seed/seed.h"
#include "../matrix/matrix.h"

#include <assert.h>
#include <stdio.h>

/**
 * If the cell has more than two neighbouring positions which are
 * not occupied, it will proliferate. The original cell will vanish
 * and split into two daughter cells, which will be randomly
 * distributed into two unoccupied neighbouring locations.
 * @param nbr_num number of neighbouring positions
 */
void cell_proliferate(
        const int nbr_num,
        const int nbr_temp[nbr_num],
        const int *nbr_coord[2],
        double ind_position[Y_LEN][X_LEN],
        const int x_coord,
        const int y_coord
) {
    const int zeros = int_array_count(nbr_num, nbr_temp, 0);
    if (zeros >= 2) {
        int_arr_2_t sample_idx = unif_index2(zeros);

        const int *sample_a = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[0] + 1)];
        const int *sample_b = nbr_coord[int_array_find(nbr_num, nbr_temp, 0, sample_idx.arr[1] + 1)];

        ind_position[x_coord][y_coord]         = 0;     // Original cell vanishes.
        ind_position[sample_a[0]][sample_a[1]] = 1;     // Daughter cells being allocated.
        ind_position[sample_b[0]][sample_b[1]] = 1;
    }
}

void proliferation(
        const int PROF_CELLS_NUM,
        double prof_cells[PROF_CELLS_NUM],
        double ind_position[Y_LEN][X_LEN],
        arraylist_t *coord
) {
    for (int i = 0; i < PROF_CELLS_NUM; ++i) {
        node_t    cell_pos = arraylist_get(coord, (int) prof_cells[i]);
        const int x_coord  = cell_pos._intPair[0];
        const int y_coord  = cell_pos._intPair[1];

        /* Possible locations for daughter cells (8 surrounding points.) */
        const int right_position[2]      = {x_coord, y_coord + 1};
        const int right_down_position[2] = {x_coord + 1, y_coord + 1};
        const int down_position[2]       = {x_coord + 1, y_coord};
        const int up_right_position[2]   = {x_coord - 1, y_coord + 1};
        const int up_position[2]         = {x_coord - 1, y_coord};
        const int up_left_position[2]    = {x_coord - 1, y_coord - 1};
        const int left_position[2]       = {x_coord, y_coord - 1};
        const int left_down_position[2]  = {x_coord + 1, y_coord - 1};

        const int right      = (int) ind_position[right_position[0]][right_position[1]];
        const int right_down = (int) ind_position[right_down_position[0]][right_down_position[1]];
        const int down       = (int) ind_position[down_position[0]][down_position[1]];
        const int up_right   = (int) ind_position[up_right_position[0]][up_right_position[1]];
        const int up         = (int) ind_position[up_position[0]][up_position[1]];
        const int up_left    = (int) ind_position[up_left_position[0]][up_left_position[1]];
        const int left       = (int) ind_position[left_position[0]][left_position[1]];
        const int left_down  = (int) ind_position[left_down_position[0]][left_down_position[1]];

        if (x_coord == 0) {
            if (y_coord == 0) {                    /* Special case: top left corner */
                /* Possible directions to move, check if these points are occupied. */
                const int neighbouring_temp[3]   = {right, right_down, down};
                const int *neighbouring_coord[3] = {right_position, right_down_position, down_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else if (y_coord == X_LEN - 1) {     /* Special case: top right corner */
                const int neighbouring_temp[3]   = {left, left_down, down};
                const int *neighbouring_coord[3] = {left_position, left_down_position, down_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else {
                const int neighbouring_temp[5]   = {left, right, down, left_down, right_down};
                const int *neighbouring_coord[5] = {left_position, right_position, down_position,
                                                    left_down_position, right_down_position};
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            }
        } else if (x_coord == Y_LEN - 1) {         /* Lower boundary */
            if (y_coord == 0) {                    /* Special case: bottom left corner */
                const int neighbouring_temp[3]   = {up, up_right, right};
                const int *neighbouring_coord[3] = {up_position, up_right_position, right_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else if (y_coord == X_LEN - 1) {     /* Special case: bottom right corner */
                const int neighbouring_temp[3]   = {left, up, up_left};
                const int *neighbouring_coord[3] = {left_position, up_position, up_left_position};
                cell_proliferate(3, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            } else {                                        // WARNING: BE CAREFUL with the order (keep for consistency)
                const int neighbouring_temp[5]   = {left, up_left, up, up_right, right};
                const int *neighbouring_coord[5] = {left_position, up_left_position, up_position,
                                                    up_right_position, right_position};
                cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
            }
        } else if (y_coord == 0 && x_coord != Y_LEN - 1) {              /* Left boundary */
            const int neighbouring_temp[5]   = {up, up_right, right, right_down, down};
            const int *neighbouring_coord[5] = {
                    up_position, up_right_position, right_position, right_down_position, down_position
            };
            cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
        } else if (y_coord == X_LEN - 1 && x_coord != Y_LEN - 1) {      /* Right boundary */
            const int neighbouring_temp[5]   = {up, up_left, left, left_down, down};
            const int *neighbouring_coord[5] = {
                    up_position, up_left_position, left_position, left_down_position,
                    down_position
            };
            cell_proliferate(5, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
        } else {                                                        /* Other cells that are not at the boundary */
            const int neighbouring_temp[8]   = {
                    left, right, up, down, up_left, up_right, left_down, right_down
            };
            const int *neighbouring_coord[8] = {
                    left_position, right_position, up_position, down_position,
                    up_left_position, up_right_position, left_down_position, right_down_position
            };
            cell_proliferate(8, neighbouring_temp, neighbouring_coord, ind_position, x_coord, y_coord);
        }
    }

    arraylist_clear(coord);
    node_t   node;
    for (int i = 0; i < Y_LEN; ++i) {
        for (int j = 0; j < X_LEN; ++j) {
            if (ind_position[i][j] == 1) {
                node._intPair[0] = i;
                node._intPair[1] = j;
                arraylist_append(coord, node);
            }
        }
    }
    assert(coord->size > 0);

    MATRIX_INIT(ind_position, Y_LEN, X_LEN, 0)
    for (int i = 0; i < coord->size; ++i) {
        node = arraylist_get(coord, i);
        ind_position[node._intPair[0]][node._intPair[1]] = 1;
    }
}