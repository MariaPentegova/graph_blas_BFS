// от нескольких источников 

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <GraphBLAS.h>
#include <LAGraph.h>

int bfs_multi_source(GrB_Matrix A, GrB_Vector dist, GrB_Index *starts, GrB_Index num_starts) {
    GrB_Info info;
    GrB_Index n;

    GrB_Matrix_nrows(&n, A);

    GrB_Vector visited = NULL, frontier = NULL, new_frontier = NULL;
    GrB_Vector_new(&visited, GrB_BOOL, n);
    GrB_Vector_new(&frontier, GrB_BOOL, n);
    GrB_Vector_new(&new_frontier, GrB_BOOL, n);

    for (GrB_Index i = 0; i < n; i++) {
        GrB_Vector_setElement_INT32(dist, -1, i);
    }

    for (GrB_Index i = 0; i < num_starts; i++) {
        GrB_Vector_setElement_BOOL(frontier, true, starts[i]);
        GrB_Vector_setElement_BOOL(visited, true, starts[i]);
        GrB_Vector_setElement_INT32(dist, 0, starts[i]);
    }

    int level = 0;

    while (true) {
        GrB_mxv(new_frontier, visited, NULL, LAGraph_LorLand_BOOL, A, frontier, NULL);

        GrB_Vector not_visited = NULL;
        GrB_Vector_new(&not_visited, GrB_BOOL, n);
        GrB_Vector_assign_BOOL(not_visited, NULL, NULL, true, GrB_ALL, n, NULL);
        GrB_eWiseAdd_BinaryOp(not_visited, NULL, NULL, GrB_LAND, not_visited, visited, NULL);
        GrB_apply(not_visited, NULL, NULL, GrB_LNOT, not_visited, NULL);

        GrB_eWiseMult_BinaryOp(new_frontier, NULL, NULL, GrB_LAND, new_frontier, not_visited, NULL);

        GrB_Vector_free(&not_visited);

        GrB_Index nvals;
        GrB_Vector_nvals(&nvals, new_frontier);
        if (nvals == 0) {
            break;
        }

        GrB_Index *indices = malloc(nvals * sizeof(GrB_Index));
        bool *values = malloc(nvals * sizeof(bool));
        GrB_Vector_extractTuples_BOOL(indices, values, &nvals, new_frontier);

        for (GrB_Index i = 0; i < nvals; i++) {
            GrB_Vector_setElement_INT32(dist, level + 1, indices[i]);
            GrB_Vector_setElement_BOOL(visited, true, indices[i]);
        }

        free(indices);
        free(values);

        GrB_Vector_clear(frontier);
        GrB_Vector_dup(&frontier, new_frontier);

        level++;
    }

    GrB_Vector_free(&visited);
    GrB_Vector_free(&frontier);
    GrB_Vector_free(&new_frontier);

    return 0;
}


// с построением дерева обхода

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <GraphBLAS.h>
#include <LAGraph.h>  

int bfs_tree(GrB_Matrix A, GrB_Vector dist, GrB_Vector parent, GrB_Index start) {
    GrB_Info info;
    GrB_Index n;

    GrB_Matrix_nrows(&n, A);

    GrB_Vector visited = NULL, frontier = NULL, new_frontier = NULL;
    GrB_Vector_new(&visited, GrB_BOOL, n);
    GrB_Vector_new(&frontier, GrB_BOOL, n);
    GrB_Vector_new(&new_frontier, GrB_BOOL, n);

    GrB_assign(dist, NULL, NULL, -1, GrB_ALL, n, NULL);
    GrB_assign(parent, NULL, NULL, -1, GrB_ALL, n, NULL);

    GrB_Vector_setElement(frontier, true, start);
    GrB_Vector_setElement(visited, true, start);
    GrB_Vector_setElement(dist, 0, start);

    int level = 0;

    GrB_Index *rows = malloc(n * sizeof(GrB_Index));
    for (GrB_Index i = 0; i < n; i++) {
        rows[i] = i;
    }

    while (true) {
        GrB_mxv(new_frontier, visited, NULL, LAGraph_LorLand_BOOL, A, frontier, NULL);

        GrB_Vector not_visited = NULL;
        GrB_Vector_new(&not_visited, GrB_BOOL, n);
        GrB_assign(not_visited, NULL, NULL, true, GrB_ALL, n, NULL);
        GrB_eWiseAdd_BinaryOp(not_visited, NULL, NULL, GrB_LAND, not_visited, visited, NULL);
        GrB_apply(not_visited, NULL, NULL, GrB_LNOT, not_visited, NULL);

        GrB_eWiseMult_BinaryOp(new_frontier, NULL, NULL, GrB_LAND, new_frontier, not_visited, NULL);
        GrB_Vector_free(&not_visited);

        GrB_Index nvals;
        GrB_Vector_nvals(&nvals, new_frontier);
        if (nvals == 0) {
            break;  
        }

        GrB_Index *indices = malloc(nvals * sizeof(GrB_Index));
        bool *values = malloc(nvals * sizeof(bool));
        GrB_Vector_extractTuples_BOOL(indices, values, &nvals, new_frontier);

        for (GrB_Index i = 0; i < nvals; i++) {
            GrB_Index v = indices[i];

            GrB_Vector_setElement(dist, level + 1, v);
            GrB_Vector_setElement(visited, true, v);

            GrB_Vector col_v;
            GrB_Vector_new(&col_v, GrB_BOOL, n);
            info = GrB_extract(col_v, NULL, NULL, A, rows, n, &v, 1, NULL);
            if (info != GrB_SUCCESS) {
                fprintf(stderr, "Ошибка при извлечении столбца %lu\n", v);
                free(indices);
                free(values);
                free(rows);
                GrB_Vector_free(&col_v);
                return info;
            }

            GrB_Vector parent_candidates;
            GrB_Vector_new(&parent_candidates, GrB_BOOL, n);
            GrB_eWiseMult_BinaryOp(parent_candidates, NULL, NULL, GrB_LAND, col_v, frontier, NULL);

            GrB_Index pc_nvals;
            GrB_Vector_nvals(&pc_nvals, parent_candidates);
            if (pc_nvals > 0) {
                GrB_Index *pc_indices = malloc(pc_nvals * sizeof(GrB_Index));
                bool *pc_values = malloc(pc_nvals * sizeof(bool));
                GrB_Vector_extractTuples_BOOL(pc_indices, pc_values, &pc_nvals, parent_candidates);

                GrB_Vector_setElement(parent, pc_indices[0], v);
                
                free(pc_indices);
                free(pc_values);
            } else {
                GrB_Vector_setElement(parent, -1, v);
            }

            GrB_Vector_free(&col_v);
            GrB_Vector_free(&parent_candidates);
        }

        free(indices);
        free(values);

        GrB_Vector_clear(frontier);
        GrB_Vector_dup(&frontier, new_frontier);

        level++;
    }
