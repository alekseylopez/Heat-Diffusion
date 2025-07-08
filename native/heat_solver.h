#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <jni.h>

#ifdef __cplusplus
extern "C"
{
#endif

// Structure to represent a 2D grid
typedef struct
{
    int nx, ny;
    double dx, dy, alpha;
    double* data;
} Grid2D;

// Structure for sparse matrix (CSR format)
typedef struct
{
    int rows, cols;
    int* row_ptr;
    int* col_idx;
    double* vals;
    int nnz;  // number of non-zeros
} SparseMatrix;

// grid operations
Grid2D* grid_create(int nx, int ny, double dx, double dy, double alpha);
void grid_destroy(Grid2D* grid);
void grid_set_value(Grid2D* grid, int i, int j, double value);
double grid_get_value(Grid2D* grid, int i, int j);
void grid_apply_boundary_conditions(Grid2D* grid, double time);

// sparse matrix operations
SparseMatrix* sparse_matrix_create(int rows, int cols);
void sparse_matrix_destroy(SparseMatrix* matrix);
void sparse_matrix_add_entry(SparseMatrix* matrix, int row, int col, double value);
void sparse_matrix_finalize(SparseMatrix* matrix);
void sparse_matrix_multiply(const SparseMatrix* A, const double* x, double* y);

// conjugate gradient solver
int conjugate_gradient_solve(const SparseMatrix* A, const double* b, double* x, double tolerance, int max_iterations);

// time stepping methods
void explicit_step(Grid2D* grid, double dt);
void implicit_step(Grid2D* grid, double dt, double theta);

// JNI interface functions
JNIEXPORT jlong JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_createGrid
  (JNIEnv *, jobject, jint, jint, jdouble, jdouble, jdouble);

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_destroyGrid
  (JNIEnv *, jobject, jlong);

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_setGridValue
  (JNIEnv *, jobject, jlong, jint, jint, jdouble);

JNIEXPORT jdouble JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_getGridValue
  (JNIEnv *, jobject, jlong, jint, jint);

JNIEXPORT jdoubleArray JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_getGridData
  (JNIEnv *, jobject, jlong);

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_explicitStep
  (JNIEnv *, jobject, jlong, jdouble);

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_implicitStep
  (JNIEnv *, jobject, jlong, jdouble, jdouble);

#ifdef __cplusplus
}
#endif

#endif // HEAT_SOLVER_H