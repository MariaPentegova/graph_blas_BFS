#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "GraphBLAS.h"
#include "LAGraph.h"

int generate_random_graph(GrB_Matrix *A, int num_nodes, float density) {
    if (GrB_Matrix_new(A, GrB_BOOL, num_nodes, num_nodes) != GrB_SUCCESS) {
        return -1;
    }

    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes; j++) {
            if (i != j && ((float)rand() / RAND_MAX) < density) {
                if (GrB_Matrix_setElement(*A, true, i, j) != GrB_SUCCESS) {
                    return -1;
                }
            }
        }
    }
    return 0;
}

int initialize_vectors(GrB_Vector *level, GrB_Vector *parent, GrB_Index n) {
    if (GrB_Vector_new(level, GrB_INT32, n) != GrB_SUCCESS) return -1;
    if (GrB_Vector_new(parent, GrB_INT32, n) != GrB_SUCCESS) return -1;

    for (GrB_Index i = 0; i < n; i++) {
        GrB_Vector_setElement(*level, -1, i);
        GrB_Vector_setElement(*parent, -1, i);
    }
    return 0;
}

int main() {
    GrB_init(GrB_NONBLOCKING);

    srand((unsigned)time(NULL));

    int num_nodes = 10;
    float density = 0.3f;

    GrB_Matrix A = NULL;
    GrB_Vector level = NULL;
    GrB_Vector parent = NULL;

    generate_random_graph(&A, num_nodes, density);

    initialize_vectors(&level, &parent, num_nodes);

    LAGr_BreadthFirstSearch(&level, &parent, A, 0, "Test BFS");

    printf("Вершина : Уровень\n");
    for (GrB_Index i = 0; i < num_nodes; i++) {
        int level_value;
        if (GrB_Vector_getElement(&level_value, level, i) == GrB_SUCCESS) {
            printf("%ld : %d\n", i, level_value);
        } else {
            printf("%ld : unreachable\n", i);
        }
    }

    return 0;
}
