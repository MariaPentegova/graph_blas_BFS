#include <stdio.h>
#include "GraphBLAS.h"
#include "LAGraph.h"

int main() {
    GrB_Info info = GrB_init(GrB_NONBLOCKING);
    if (info != GrB_SUCCESS) {
        printf("Ошибка инициализации GraphBLAS\n");
        return 1;
    }

    printf("Успех! GraphBLAS и LAGraph подключены.\n");

    // Завершаем работу
    GrB_finalize();
    return 0;
}