#ifndef CP_SYSTEM_H
#define CP_SYSTEM_H


#include <stdlib.h>


#define SYSTEM_DIM 4


typedef struct SystemEquations {
    int (* const function) (double, const double[], double[], void *);
    int (* const jacobian) (double, const double[], double *, double[], void *);
} SystemEquations;


typedef struct ProblemParamaters {
    double d;
    double l1;
    double l2;
    double l3;
    double m1;
    double m2;
    double k;
    double g;
} ProblemParamaters;


typedef struct SystemParameters {
    size_t step_count;
    double t0;
    double y0[SYSTEM_DIM];
    double t1;
    ProblemParamaters pp;
} SystemParameters;


SystemEquations const * get_system_eq(void);


#endif
