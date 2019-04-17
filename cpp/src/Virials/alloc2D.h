//
// Created by arpit on 4/17/2019.
//

#ifndef ZENO_ALLOC2D_H
#define ZENO_ALLOC2D_H

#include <stdlib.h>

inline static void** malloc2D(int rows, int cols, size_t s) {
    void* raw = malloc(s*rows*cols);
    void** array = (void**)malloc(rows*sizeof(void*));
    for (int i=0; i<rows; i++) {
        array[i] = ((char*)raw)+s*cols*i;
    }
    return array;
}

inline static void** realloc2D(void** array, int rows, int cols, size_t s) {
    if (!array) return malloc2D(rows, cols, s);
    void* newRaw = realloc(*array, s*rows*cols);
    void** newArray = (void**)realloc(array, rows*sizeof(void*));
    for (int i=0; i<rows; i++) {
        newArray[i] = ((char*)newRaw)+s*cols*i;
    }
    return newArray;
}

inline static void free2D(void** array) {
    if (!array) return;
    free(*array);
    free(array);
}
#endif //ZENO_ALLOC2D_H
