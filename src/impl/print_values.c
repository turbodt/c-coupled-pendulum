#include "../printers.h"
#include <stdio.h>


void print_values(
    size_t step,
    double t,
    double y[],
    SystemParameters const *params
) {
    (void) step;
    printf("%.5e", t);
    printf(" %.5e %.5e", y[0], y[1]);
    printf("\n");
};
