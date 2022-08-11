#include "collection.h"

#include <stdlib.h>
#include <assert.h>

arraylist_t *new_arraylist() {
    arraylist_t *list = malloc(sizeof(arraylist_t));
    list->len      = 0;
    list->buf_size = ARRAYLIST_INITIAL_CAPACITY;
    list->buf      = malloc(sizeof(void *) * ARRAYLIST_INITIAL_CAPACITY);
    assert(list->buf);
    return list;
}