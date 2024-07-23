#include "../printers.h"
#include <stdio.h>
#include <math.h>


void print_headers(SystemParameters const * params) {
    (void) params;
    printf("step time");
    printf(" theta_1 theta_2 dot_theta_1 dot_theta_2");
    printf(" x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6");
    printf(" T U_k U_h E");
    printf("\n");
};


void print_values(
    size_t step,
    double t,
    double y[],
    SystemParameters const *params
) {
    (void) step;
    ProblemParamaters const * pp = &params->pp;

    double const theta1 = y[0];
    double const theta2 = y[1];
    double const theta_dot_1 = y[2];
    double const theta_dot_2 = y[3];
    double const theta_diff = theta2 - theta1;

    double const s1 = sin(theta1);
    double const s2 = sin(theta2);

    double const c_diff = cos(theta_diff);
    double const diff_s = s2 - s1;

    double const B =
        pp->d * pp->d
        + 2 * pp->l3 * pp->l3 * (1 - c_diff)
        + 2 * pp->d * pp->l3 * diff_s;
    double const D = sqrt(B);

    double const x1 = 0, y1 = 0, x2 = pp->d, y2 = 0;
    double const x3 = x1 + pp->l3 * sin(theta1), y3 = y1 - pp->l3 * cos(theta1);
    double const x4 = x2 + pp->l3 * sin(theta2), y4 = y2 - pp->l3 * cos(theta2);
    double const x5 = x1 + pp->l1 * sin(theta1), y5 = y1 - pp->l1 * cos(theta1);
    double const x6 = x2 + pp->l2 * sin(theta2), y6 = y2 - pp->l2 * cos(theta2);

    double const T = 0.5 * (
        pp->m1 * pp->l1 * pp->l1 * theta_dot_1 * theta_dot_1
        + pp->m2 * pp->l2 * pp->l2 * theta_dot_2 * theta_dot_2
    );

    double const dD = pp->d - D;
    double const Uk = 0.5 * pp->k * dD * dD;
    double const Uh_0 = -pp->g * (pp->m1 * pp->l1 + pp->m2 *pp->l2);
    double const Uh = pp->g * (pp->m1 * y5 + pp->m2 * y6) - Uh_0;
    double const E = T + Uk + Uh;

    printf("%zu %.5e", step, t);
    printf(" %.5e %.5e %.5e %.5e", theta1, theta2, theta_dot_1, theta_dot_2);
    printf(" %.5e %.5e %.5e %.5e", x1, y1, x2, y2);
    printf(" %.5e %.5e %.5e %.5e", x3, y3, x4, y4);
    printf(" %.5e %.5e %.5e %.5e", x5, y5, x6, y6);
    printf(" %.5e %.5e %.5e %.5e", T, Uk, Uh, E);
    printf("\n");
};
