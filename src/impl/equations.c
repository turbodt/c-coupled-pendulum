#include "../system.h"
#include <math.h>
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
    ProblemParamaters const * pp = &params->pp;

    double const g1 = pp->g / pp->l1;
    double const g2 = pp->g / pp->l2;
    double const kml1 = pp->k / pp->m1 / pp->l1 / pp->l1;
    double const kml2 = pp->k / pp->m2 / pp->l2 / pp->l2;
    double const theta1 = y[0];
    double const theta2 = y[1];
    double const theta_dot_1 = y[2];
    double const theta_dot_2 = y[3];
    double const theta_diff = theta2 - theta1;

    double const c1 = cos(theta1);
    double const c2 = cos(theta2);
    double const s1 = sin(theta1);
    double const s2 = sin(theta2);

    double const c_diff = cos(theta_diff);
    double const s_diff = sin(theta_diff);
    double const diff_s = s2 - s1;

    double const B =
        pp->d * pp->d
        + 2 * pp->l3 * pp->l3 * (1 - c_diff)
        + 2 * pp->d * pp->l3 * diff_s;
    double const D = sqrt(B);

    double const dBd1 = -2 * pp->l3 * (pp->d * c1 + pp->l3 * s_diff);
    double const dBd2 = 2 * pp->l3 * (pp->d * c2 + pp->l3 * s_diff);

    double const dDd1 = dBd1 / D / 2;
    double const dDd2 = dBd2 / D / 2;

    double const Dd = D - pp->d;

    dydt[0] = theta_dot_1;
    dydt[1] = theta_dot_2;
    dydt[2] = - kml1 * Dd * dDd1 - g1 * s1 - pp->mu * theta_dot_1;
    dydt[3] = - kml2 * Dd * dDd2 - g2 * s2 - pp->mu * theta_dot_2;

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
    ProblemParamaters const * pp = &params->pp;

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    dfdt[4] = 0.0;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
    gsl_matrix * m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, 0.0);

    gsl_matrix_set(m, 0, 2, 1.0);
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 1, 2, 0.0);
    gsl_matrix_set(m, 1, 3, 1.0);

    double const g1 = pp->g / pp->l1;
    double const g2 = pp->g / pp->l2;
    double const kml1 = pp->k / pp->m1 / pp->l1 / pp->l1;
    double const kml2 = pp->k / pp->m2 / pp->l2 / pp->l2;
    double const theta1 = y[0];
    double const theta2 = y[1];
    double const theta_diff = theta2 - theta1;

    double const c1 = cos(theta1);
    double const c2 = cos(theta2);
    double const s1 = sin(theta1);
    double const s2 = sin(theta2);

    double const c_diff = cos(theta_diff);
    double const s_diff = sin(theta_diff);
    double const diff_s = s2 - s1;

    double const B =
        pp->d * pp->d
        + 2 * pp->l3 * pp->l3 * (1 - c_diff)
        + 2 * pp->d * pp->l3 * diff_s;
    double const D = sqrt(B);
    double const D3 = B*D;

    double const dBd1 = -2 * pp->l3 * (pp->d * c1 + pp->l3 * s_diff);
    double const dBd2 = 2 * pp->l3 * (pp->d * c2 + pp->l3 * s_diff);
    double const d2Bd1d1 = 2 * pp->l3 * (pp->d * s1 + pp->l3 * c_diff);
    double const d2Bd2d2 = -2 * pp->l3 * (pp->d * s2 - pp->l3 * c_diff);
    double const d2Bd1d2 = -2 * pp->l3 * pp->l3 * c_diff;

    double const dDd1 = dBd1 / D / 2;
    double const dDd2 = dBd2 / D / 2;
    double const d2Dd1d1 = d2Bd1d1 / D / 2 - dBd1 * dBd1 / D3 / 4;
    double const d2Dd1d2 = d2Bd1d2 / D / 2 - dBd1 * dBd2 / D3 / 4;
    double const d2Dd2d2 = d2Bd2d2 / D / 2 - dBd2 * dBd2 / D3 / 4;

    double const Dd = D - pp->d;

    double const a[4] = {
        -kml1 * (dDd1 * dDd1 + Dd * d2Dd1d1) - g1 * c1,
        -kml1 * (dDd1 * dDd2 + Dd * d2Dd1d2),
        -kml2 * (dDd1 * dDd2 + Dd * d2Dd1d2),
        -kml2 * (dDd2 * dDd2 + Dd * d2Dd2d2) - g2 * c2,
    };

    gsl_matrix_set(m, 2, 0, a[0]);
    gsl_matrix_set(m, 2, 1, a[1]);
    gsl_matrix_set(m, 3, 0, a[2]);
    gsl_matrix_set(m, 3, 1, a[3]);

    gsl_matrix_set(m, 2, 2, -pp->mu);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -pp->mu);

    return GSL_SUCCESS;
}
