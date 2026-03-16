#include "GraphBLAS.h"
#include <stdio.h>

#define CHECK(status) do { if (status != GrB_SUCCESS) { printf("Error at line %d\n", __LINE__); return -1; } } while (0)

int initialize_vectors(GrB_Vector *levels, GrB_Vector *parents, GrB_Index n) {
    if (GrB_Vector_new(levels, GrB_INT32, n) != GrB_SUCCESS) return -1;
    if (GrB_Vector_new(parents, GrB_INT32, n) != GrB_SUCCESS) {
        GrB_Vector_free(levels);
        return -1;
    }
    for (GrB_Index i = 0; i < n; i++) {
        GrB_Vector_setElement(*levels, -1, i);
        GrB_Vector_setElement(*parents, -1, i);
    }
    return 0;
}


int MultiSource_BFS(GrB_Vector *levels, GrB_Vector *parents, GrB_Matrix A, GrB_Vector sources) {
    GrB_Index n;
    CHECK(GrB_Matrix_nrows(&n, A));

    if (initialize_vectors(levels, parents, n) != 0) return -1;

    for (GrB_Index i = 0; i < n; i++) {
        bool has_source;
        if (GrB_Vector_extractElement(&has_source, sources, i) == GrB_SUCCESS && has_source) {
            CHECK(GrB_Vector_setElement(*levels, 0, i));
            CHECK(GrB_Vector_setElement(*parents, i, i));
        }
    }

    GrB_Vector frontier = NULL, new_frontier = NULL;
    CHECK(GrB_Vector_new(&frontier, GrB_BOOL, n));
    CHECK(GrB_Vector_new(&new_frontier, GrB_BOOL, n));

    CHECK(GrB_assign(frontier, NULL, NULL, true, sources, NULL));

    int current_level = 0;

    while (true) {
        bool frontier_exists = false;
        CHECK(GrB_reduce(&frontier_exists, NULL, GrB_LOR_MONOID_BOOL, frontier, NULL));
        if (!frontier_exists) break;

        CHECK(GrB_vxm(&new_frontier, NULL, NULL, GrB_LOR_SEMIRING_BOOL, A, frontier, NULL));

        for (GrB_Index i = 0; i < n; i++) {
            bool reached;
            if (GrB_Vector_extractElement(&reached, new_frontier, i) == GrB_SUCCESS && reached) {
                int32_t curr_level;
                if (GrB_Vector_extractElement(&curr_level, *levels, i) != GrB_SUCCESS || curr_level == -1) {
                    CHECK(GrB_Vector_setElement(*levels, current_level + 1, i));
                }

                int32_t parent_value;
                if (GrB_Vector_extractElement(&parent_value, *parents, i) != GrB_SUCCESS || parent_value == -1) {
                    int index;
                    if (GrB_Vector_extractElement(&index, frontier, i) == GrB_SUCCESS) {
                        CHECK(GrB_Vector_setElement(*parents, index, i));
                    }
                }
            }
        }

        GrB_Vector tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        CHECK(GrB_Vector_clear(new_frontier));
        current_level++;
    }

    if (frontier) GrB_Vector_free(&frontier);
    if (new_frontier) GrB_Vector_free(&new_frontier);

    return 0;
}


int Parent_BFS(GrB_Vector *levels, GrB_Vector *parents, GrB_Matrix A, GrB_Index start_node, const char *label) {
    GrB_Index n;
    CHECK(GrB_Matrix_nrows(&n, A));

    if (initialize_vectors(levels, parents, n) != 0) return -1;

    CHECK(GrB_Vector_setElement(*levels, 0, start_node));
    CHECK(GrB_Vector_setElement(*parents, start_node, start_node));

    GrB_Vector frontier = NULL, new_frontier = NULL;
    CHECK(GrB_Vector_new(&frontier, GrB_BOOL, n));
    CHECK(GrB_Vector_new(&new_frontier, GrB_BOOL, n));

    CHECK(GrB_Vector_setElement(frontier, true, start_node));

    int current_level = 0;

    while (true) {

        bool frontier_exists = false;
        CHECK(GrB_reduce(&frontier_exists, NULL, GrB_LOR_MONOID_BOOL, frontier, NULL));
        if (!frontier_exists) break;

        CHECK(GrB_vxm(&new_frontier, NULL, NULL, GrB_LOR_SEMIRING_BOOL, A, frontier, NULL));

        for (GrB_Index i = 0; i < n; i++) {
            bool reached;
            if (GrB_Vector_extractElement(&reached, new_frontier, i) == GrB_SUCCESS && reached) {
                int32_t curr_level;
                if (GrB_Vector_extractElement(&curr_level, *levels, i) != GrB_SUCCESS || curr_level == -1) {
                    CHECK(GrB_Vector_setElement(*levels, current_level + 1, i));
                }

                int parent_index;
                if (GrB_Vector_extractElement(&parent_index, frontier, i) == GrB_SUCCESS) {
                    CHECK(GrB_Vector_setElement(*parents, parent_index, i));
                }
            }
        }

        GrB_Vector tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        CHECK(GrB_Vector_clear(new_frontier));

        current_level++;
    }

    if (frontier) GrB_Vector_free(&frontier);
    if (new_frontier) GrB_Vector_free(&new_frontier);

    return 0;
}
