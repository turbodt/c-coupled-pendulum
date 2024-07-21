#include <stdio.h>
#include <memory.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "./system.h"
#include "./printers.h"


#define H_START 1e-6
#define EPS_ABS 1e-10
#define EPS_REL 0.0


int main(void) {

    SystemParameters sys_params = {
        .step_count = 100,
        .t0 = 0.0,
        .y0 = {1.0, 0.0},
        .t1 = 10.0,
        .pp = {
            .d = 1.0,
            .l1 = 1.0,
            .l2 = 1.0,
            .l3 = 0.5,
            .m1 = 1.0,
            .m2 = 1.0,
            .k = 1.0,
            .g = 9.8,
        }
    };
    SystemEquations const * eq = get_system_eq();

    gsl_odeiv2_system system = {
        eq->function,
        eq->jacobian,
        SYSTEM_DIM,
        &sys_params
    };

    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(
        &system,
        gsl_odeiv2_step_rk8pd,
        H_START,
        EPS_ABS,
        EPS_REL
    );

    double t = sys_params.t0;
    double y[SYSTEM_DIM];

    memcpy(y, sys_params.y0, sizeof(sys_params.y0));

    print_headers(&sys_params);
    print_values(0, t, y, &sys_params);

    for (size_t i = 0; i < sys_params.step_count; i++) {
        double ti = (i + 1) * sys_params.t1 / ((double) sys_params.step_count);
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf("Error, return value=%d\n", status);
            break;
        }

        print_values(i + 1, t, y, &sys_params);
    }

    gsl_odeiv2_driver_free(driver);

    return 0;
}
