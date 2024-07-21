#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


#define H_START 1e-6
#define EPS_ABS 1e-10
#define EPS_REL 0.0


typedef struct SystemParameters {
    size_t n;
    double t0;
    double t1;
} SystemParameters;


static void output_values(size_t, double, double[], SystemParameters const *);


int function (double t, const double y[], double dydt[], void * p) {
    (void)(t);
    SystemParameters * params = (SystemParameters *) p;
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


int main(void) {

    SystemParameters params = {
        .n = 100,
        .t0 = 0.0,
        .t1 = 10.0,
    };

    gsl_odeiv2_system system = {
        function,
        jacobian,
        2,
        &params
    };

    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(
        &system,
        gsl_odeiv2_step_rk8pd,
        H_START,
        EPS_ABS,
        EPS_REL
    );

    double t = params.t0;
    double y[2] = {1, 0};
    output_values(0, t, y, &params);

    for (size_t i = 0; i < params.n; i++) {
        double ti = (i + 1) * params.t1 / ((double) params.n);
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf ("error, return value=%d\n", status);
            break;
        }

        output_values(i, t, y, &params);
    }

    gsl_odeiv2_driver_free(driver);

    return 0;
}


void output_values(
    size_t step,
    double t,
    double y[],
    SystemParameters const *params
) {
    printf("%.5e", t);
    printf(" %.5e %.5e", y[0], y[1]);
    printf("\n");
};
