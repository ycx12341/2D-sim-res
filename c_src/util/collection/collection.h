#ifndef C_SRC_COLLECTION_H
#define C_SRC_COLLECTION_H

#include <stdlib.h>
#include <stdbool.h>

/**
 * Generate regular sequences.
 * @param buf   A buffer to store the generated sequence.
 * @param from  Lower bound of the sequence.
 * @param to    Upper bound of the sequence.
 * @param num   Length of sequence.
 * @return The number of sequence.
 */
int seq(double *buf, double from, double to, int num);

typedef struct pair pair_t;
typedef struct node node_t;

struct node {
    double _double;
    int    _int;
    int    _intPair[2];
    void   *_ptr;
};

struct pair {
    node_t x;
    node_t y;
};

#define ARRAYLIST_INITIAL_CAPACITY 4

typedef struct {
    unsigned int size;          // Count of items currently in list
    unsigned int capacity;      // Allocated memory size, in items
    node_t       *body;         // Pointer to allocated memory for items (of size capacity * sizeof(void*))
    int          need_free;     // Set to 1 if the element is dynamically allocated

} arraylist_t;

arraylist_t *new_arraylist(bool need_free);

void arraylist_free(arraylist_t *l);

void arraylist_append(arraylist_t *l, node_t item);

node_t arraylist_get(arraylist_t *l, unsigned int index);

#endif //C_SRC_COLLECTION_H