
typedef struct FitData{
    double *xvals;
    double *yvals;
    size_t sample_size;
    double *params;
    size_t params_size;
}

void gen_fit(const double *x,const double *y);


void fit_conjugate_gradient(const double *x,const double *y);

void fit_steepest_descent(const double *x,const double *y);
void fit_gauss_newton(const double *x,const double *y);
void fit_levenberg_marquardt(const double *x,const double *y);

void fit_linear_system(const double *x,const double *y);
