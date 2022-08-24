#include "collection.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

arraylist_t *new_arraylist() {
    arraylist_t *new = malloc(sizeof(arraylist_t));
    assert(new);
    new->size     = 0;
    new->body     = malloc(sizeof(node_t) * ARRAYLIST_INITIAL_CAPACITY);
    new->capacity = ARRAYLIST_INITIAL_CAPACITY;
    assert(new->body);
    return new;
}

void arraylist_free(arraylist_t *l) {
    if (l->body) {
        free(l->body);
        l->body = NULL;
    }
    free(l);
    l = NULL;
}

void arraylist_clear(arraylist_t *l) {
    l->size = 0;
}

void arraylist_allocate(arraylist_t *l, int size) {
    assert(size > 0);
    if (size > l->capacity) {
        int new_capacity = l->capacity;
        while (new_capacity < size) {
            new_capacity *= 2;
        }
        l->body = realloc(l->body, sizeof(node_t) * new_capacity);
        assert(l->body);
        l->capacity = new_capacity;
    }
}

void arraylist_append(arraylist_t *l, node_t item) {
    arraylist_allocate(l, l->size + 1);
    l->body[l->size++] = item;
}

node_t arraylist_get(arraylist_t *l, int index) {
    assert(index < l->size);
    return l->body[index];
}

#define arraylist_memshift(s, offset, length) memmove((s) + (offset), s, (length) * sizeof(node_t));

node_t arraylist_remove(arraylist_t *l, const int index) {
    node_t value = l->body[index];
    arraylist_memshift(l->body + index + 1, -1, l->size - index)
    l->size--;
    return value;
}

void arraylist_remove_many(arraylist_t *l, const int num, const int indices[num]) {
    node_t *new = malloc(sizeof(node_t) * l->capacity);
    for (int i = 0, len = l->size; i < len; ++i) {
        if (!int_array_contain(num, indices, i)) {
            new[i] = arraylist_get(l, i);
        }
    }
    free(l->body);
    l->body = new;
}

void arraylist_set(arraylist_t *l, const unsigned int index, const node_t value) {
    assert(index < l->size);
    l->body[index] = value;
}