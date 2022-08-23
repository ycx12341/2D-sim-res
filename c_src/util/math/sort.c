#include "math_ext.h"

void swap(double *a, double *b) {
    double t = *a;
    *a = *b;
    *b = t;
}

int partition(double arr[], const int low, const int high, const bool ascending) {
    double pivot = arr[high];
    int    i     = (low - 1);

    for (int j = low; j < high; j++) {
        if (ascending && arr[j] > pivot || !ascending && arr[j] < pivot) { continue; }
        i++;
        swap(&arr[i], &arr[j]);
    }
    swap(&arr[i + 1], &arr[high]);
    return i + 1;
}

void quickSort_recur(double array[], const int low, const int high, const bool ascending) {
    if (low < high) {
        int pi = partition(array, low, high, ascending);
        quickSort_recur(array, low, pi - 1, ascending);
        quickSort_recur(array, pi + 1, high, ascending);
    }
}

void quickSort(const int len, double array[], bool ascending) {
    quickSort_recur(array, 0, len - 1, ascending);
}