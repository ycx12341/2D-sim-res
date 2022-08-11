#ifndef C_SRC_COLLECTION_H
#define C_SRC_COLLECTION_H

/**
 * Generate regular sequences.
 * @param buf   A buffer to store the generated sequence.
 * @param from  Lower bound of the sequence.
 * @param to    Upper bound of the sequence.
 * @param num   Length of sequence.
 * @return The number of sequence.
 */
int seq(double *buf, double from, double to, int num);

typedef struct {
    double _double;
    int    _int;
} node_t;

typedef struct {
    node_t x;
    node_t y;
} pair_t;

#define ARRAYLIST_INITIAL_CAPACITY 8

typedef struct {
    void         **buf;
    unsigned int buf_size;
    unsigned int len;
} arraylist_t;

arraylist_t *new_arraylist();

#endif //C_SRC_COLLECTION_H