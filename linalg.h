
typedef struct Matr{
    size_t m_width;
    size_t m_height;
    double *data;
} Matrix;

void print_matr(const double* matr);
void print_colvec(const double* colvec);
void print_rowvec(const double* rowvec);
double veclen(const double* vec);
void multiply(const double* L,const double* R, double* Out);
void gen_reflection(int kpiv, const double* V, double* H);
void gen_reflection_by_matr(const double* M, double* H);
void invert_triangle(const double *R, double *Inv);
int do_linalg();
void invert(const double *M,const double *Minv);
void mulcolvec(const double *M, const double *V, double *Res);
