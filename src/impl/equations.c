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

static double calcule_dTkdtheta1(double const [], SystemParameters const *);
static double calcule_dTkdtheta2(double const [], SystemParameters const *);
static double calcule_d2Tkdtheta11(double const [], SystemParameters const *);
static double calcule_d2Tkdtheta12(double const [], SystemParameters const *);
static double calcule_d2Tkdtheta22(double const [], SystemParameters const *);


SystemEquations const * get_system_eq(void) {
    return &system_equations;
};


int function (double t, const double y[], double dydt[], void * p) {
    (void)(t);
    SystemParameters const * params = p;
    ProblemParamaters const * pp = &params->pp;

    double const theta1 = y[0];
    double const theta2 = y[1];
    double const l1_2 = pp->l1*pp->l1;
    double const l2_2 = pp->l2*pp->l2;
    double const dTkdtheta1 = calcule_dTkdtheta1(y, params);
    double const dTkdtheta2 = calcule_dTkdtheta2(y, params);

    dydt[0] = y[2];
    dydt[1] = y[3];
    dydt[2] = - pp->g / pp->l1 * sin(theta1) - dTkdtheta1 / pp->m1 / l1_2;
    dydt[2] = - pp->g / pp->l2 * sin(theta2) - dTkdtheta2 / pp->m2 / l2_2;

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

    gsl_matrix_set(m, 2, 2, 0.0);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, 0.0);

    double const theta1 = y[0];
    double const theta2 = y[1];
    double const l1_2 = pp->l1*pp->l1;
    double const l2_2 = pp->l2*pp->l2;
    double const d2Tkdtheta11 = calcule_d2Tkdtheta11(y, params);
    double const d2Tkdtheta12 = calcule_d2Tkdtheta12(y, params);
    double const d2Tkdtheta22 = calcule_d2Tkdtheta22(y, params);

    double const a[4] = {
        - pp->g / pp->l1 * cos(theta1) - d2Tkdtheta11 / pp->m1 / l1_2,
        - d2Tkdtheta12 / pp->m1 / l1_2,
        - d2Tkdtheta12 / pp->m2 / l2_2,
        - pp->g / pp->l2 * cos(theta1) - d2Tkdtheta22 / pp->m2 / l2_2,
    };

    gsl_matrix_set(m, 2, 0, a[0]);
    gsl_matrix_set(m, 2, 1, a[1]);
    gsl_matrix_set(m, 3, 0, a[2]);
    gsl_matrix_set(m, 3, 1, a[3]);

    return GSL_SUCCESS;
}


double calcule_dTkdtheta1(double const y[], SystemParameters const *params) {
    ProblemParamaters const * pp = &params->pp;

    double const theta1 = y[0];
    double const theta2 = y[1];
    double const theta_diff = theta2 - theta1;

    double const c1 = cos(theta1);
    double const c2 = cos(theta2);
    double const s1 = sin(theta1);
    double const s2 = sin(theta2);

    double const c_diff = cos(theta_diff);
    double const s_diff = sin(theta_diff);
    double const diff_c = c2 - c1;
    double const diff_s = s2 - s1;

    double const d_ = pp->d / pp->l3;

    double const rad = d_ * d_ + 2 * (1 - c_diff) + 2 * d_ * diff_s;
    double const root = sqrt(rad);

    double const A = d_ * c1 + s_diff;
    return - pp->k * pp->l3 * pp->l3 * (A + d_ * A / root);
};


double calcule_dTkdtheta2(double const y[], SystemParameters const *params) {
    ProblemParamaters const * pp = &params->pp;

    double const theta1 = y[0];
    double const theta2 = y[1];
    double const theta_diff = theta2 - theta1;

    double const c1 = cos(theta1);
    double const c2 = cos(theta2);
    double const s1 = sin(theta1);
    double const s2 = sin(theta2);

    double const c_diff = cos(theta_diff);
    double const s_diff = sin(theta_diff);
    double const diff_c = c2 - c1;
    double const diff_s = s2 - s1;

    double const d_ = pp->d / pp->l3;

    double const rad = d_ * d_ + 2 * (1 - c_diff) + 2 * d_ * diff_s;
    double const root = sqrt(rad);

    double const A = d_ * c2 + s_diff;
    return pp->k * pp->l3 * pp->l3 * (A + d_ * A / root);
};
