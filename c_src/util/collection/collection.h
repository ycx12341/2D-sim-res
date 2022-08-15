#ifndef C_SRC_COLLECTION_H
#define C_SRC_COLLECTION_H

#include <stdlib.h>
#include <stdbool.h>

/**
 * Generate regular sequences.
 * @param buf        A buffer to store the generated sequence.
 * @param from       Lower bound of the sequence.
 * @param to         Upper bound of the sequence.
 * @param length_out Desired length of the sequence.
 * @return The number of sequence.
 */
int seq_length_out(double *buf, double from, double to, int length_out);

/**
 * Generate regular sequences.
 * @param buf        A buffer to store the generated sequence.
 * @param from       Lower bound of the sequence.
 * @param to         Upper bound of the sequence.
 * @param by         Increment of the sequence.
 * @return The number of sequence.
 */
int seq_by(double *buf, double from, double to, double by);

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
    int    size;                // Count of items currently in list
    int    capacity;            // Allocated memory size, in items
    node_t *body;               // Pointer to allocated memory for items (of size capacity * sizeof(void*))
    int    need_free;           // Set to 1 if the element is dynamically allocated

} arraylist_t;

arraylist_t *new_arraylist(bool need_free);

void arraylist_free(arraylist_t *l);

void arraylist_append(arraylist_t *l, node_t item);

node_t arraylist_get(arraylist_t *l, int index);

node_t arraylist_remove(arraylist_t *l, int index);

/**
 * Find the minimum value in the array and the count of them.
 * @return Pair(min_value, count_of_min_values)
 */
pair_t array_min(int len, const double arr[len]);

/**
 * Find the n_th element in the array with `element == value`
 * @return the index of found element
 */
int array_find(int len, const double arr[len], double val, int n);

#endif //C_SRC_COLLECTION_H