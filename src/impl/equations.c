#include "../system.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>


static int function (double, const double[], double [], void *);
static int jacobian (double, const double[], double *, double[], void *);
static SystemEquations const system_equations = {
    .function=function,
    .jacobian=jacobian,
};


SystemEquations const * get_system_eq(void) {
    return &system_equations;
};


int function (double t, const double y[], double dydt[], void * p) {
    (void)(t);
    SystemParameters const * params = p;
    dydt[0] = y[1];
    dydt[1] = -y[0];
    return GSL_SUCCESS;
};


int jacobian (
    double t,
    const double y[],
    double * dfdy,
    double dfdt[],
    void * p
) {
    (void)(t);
    SystemParameters * params = (SystemParameters *) p;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, -1.0);
    gsl_matrix_set(m, 1, 1, 0.0);

    return GSL_SUCCESS;
}
