#include "heat_solver.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

// grid operations
Grid2D* grid_create(int nx, int ny, double dx, double dy, double alpha)
{
    Grid2D* grid = (Grid2D*) malloc(sizeof(Grid2D));

    grid->nx = nx;
    grid->ny = ny;
    grid->dx = dx;
    grid->dy = dy;
    grid->alpha = alpha;
    grid->data = (double*) calloc(nx * ny, sizeof(double));

    return grid;
}

void grid_destroy(Grid2D* grid)
{
    if(grid)
    {
        free(grid->data);
        free(grid);
    }
}

void grid_set_value(Grid2D* grid, int i, int j, double value)
{
    if(i >= 0 && i < grid->nx && j >= 0 && j < grid->ny)
        grid->data[i + j * grid->nx] = value;
}

double grid_get_value(Grid2D* grid, int i, int j)
{
    if(i >= 0 && i < grid->nx && j >= 0 && j < grid->ny)
        return grid->data[i + j * grid->nx];

    return 0.0;
}

void grid_apply_boundary_conditions(Grid2D* grid, double time)
{
    // apply zero boundary conditions

    for(int i = 0; i < grid->nx; i++)
    {
        grid_set_value(grid, i, 0, 0.0); // bottom
        grid_set_value(grid, i, grid->ny - 1, 0.0); // top
    }

    for (int j = 0; j < grid->ny; j++)
    {
        grid_set_value(grid, 0, j, 0.0); // left
        grid_set_value(grid, grid->nx - 1, j, 0.0); // right
    }
}

// sparse matrix operations
SparseMatrix* sparse_matrix_create(int rows, int cols)
{
    SparseMatrix* matrix = (SparseMatrix*) malloc(sizeof(SparseMatrix));

    matrix->rows = rows;
    matrix->cols = cols;
    matrix->row_ptr = (int*) calloc(rows + 1, sizeof(int));
    matrix->col_idx = NULL;
    matrix->vals = NULL;
    matrix->nnz = 0;

    return matrix;
}

void sparse_matrix_destroy(SparseMatrix* matrix)
{
    if(matrix)
    {
        free(matrix->row_ptr);
        free(matrix->col_idx);
        free(matrix->vals);
        free(matrix);
    }
}

void sparse_matrix_multiply(const SparseMatrix* A, const double* x, double* y)
{
    for(int i = 0; i < A->rows; i++)
    {
        y[i] = 0.0;

        for(int idx = A->row_ptr[i]; idx < A->row_ptr[i + 1]; idx++)
            y[i] += A->vals[idx] * x[A->col_idx[idx]];
    }
}

// conjugate gradient solver
int conjugate_gradient_solve(const SparseMatrix* A, const double* b, double* x, double tolerance, int max_iterations)
{
    int n = A->rows;
    double* r = (double*)malloc(n * sizeof(double));
    double* p = (double*)malloc(n * sizeof(double));
    double* Ap = (double*)malloc(n * sizeof(double));
    
    // r = b - A*x
    sparse_matrix_multiply(A, x, r);
    for(int i = 0; i < n; i++)
    {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }
    
    double rsold = 0.0;
    for(int i = 0; i < n; i++)
        rsold += r[i] * r[i];
    
    int iter;
    for(iter = 0; iter < max_iterations; iter++)
    {
        if(sqrt(rsold) < tolerance)
            break;
        
        sparse_matrix_multiply(A, p, Ap);
        
        double pAp = 0.0;
        for(int i = 0; i < n; i++)
            pAp += p[i] * Ap[i];
        
        double alpha = rsold / pAp;
        
        for(int i = 0; i < n; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        double rsnew = 0.0;
        for(int i = 0; i < n; i++)
            rsnew += r[i] * r[i];
        
        double beta = rsnew / rsold;
        for(int i = 0; i < n; i++)
            p[i] = r[i] + beta * p[i];
        
        rsold = rsnew;
    }
    
    free(r);
    free(p);
    free(Ap);
    
    return iter;
}

// explicit time stepping
void explicit_step(Grid2D* grid, double dt)
{
    int nx = grid->nx;
    int ny = grid->ny;
    double dx2 = grid->dx * grid->dx;
    double dy2 = grid->dy * grid->dy;
    double alpha = grid->alpha;
    
    double* new_data = (double*) malloc(nx * ny * sizeof(double));
    memcpy(new_data, grid->data, nx * ny * sizeof(double));
    
    // update interior points
    for(int i = 1; i < nx - 1; i++)
    {
        for(int j = 1; j < ny - 1; j++)
        {
            double u_ij = grid_get_value(grid, i, j);
            double u_ip1 = grid_get_value(grid, i + 1, j);
            double u_im1 = grid_get_value(grid, i - 1, j);
            double u_jp1 = grid_get_value(grid, i, j + 1);
            double u_jm1 = grid_get_value(grid, i, j - 1);
            
            double d2u_dx2 = (u_ip1 - 2.0 * u_ij + u_im1) / dx2;
            double d2u_dy2 = (u_jp1 - 2.0 * u_ij + u_jm1) / dy2;
            
            new_data[i + j * nx] = u_ij + dt * alpha * (d2u_dx2 + d2u_dy2);
        }
    }
    
    free(grid->data);
    grid->data = new_data;
    
    grid_apply_boundary_conditions(grid, 0.0);
}

// implicit time stepping
void implicit_step(Grid2D* grid, double dt, double theta)
{
    int nx = grid->nx;
    int ny = grid->ny;
    int n = nx * ny;
    
    double dx2 = grid->dx * grid->dx;
    double dy2 = grid->dy * grid->dy;
    double alpha = grid->alpha;
    
    // TODO: for now, just do a simple backward Euler step - full implicit later
    double* b = (double*) malloc(n * sizeof(double));
    double* x = (double*) malloc(n * sizeof(double));
    
    // copy current solution as initial guess
    memcpy(x, grid->data, n * sizeof(double));
    memcpy(b, grid->data, n * sizeof(double));
    
    // simple iterative solver for demonstration
    for(int iter = 0; iter < 10; iter++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            for(int j = 1; j < ny - 1; j++)
            {
                int idx = i + j * nx;
                double coeff = 1.0 + dt * alpha * theta * (2.0 / dx2 + 2.0 / dy2);
                
                double rhs = b[idx];
                if(i > 0)
                    rhs += dt * alpha * theta * x[idx - 1] / dx2;
                if(i < nx - 1)
                    rhs += dt * alpha * theta * x[idx + 1] / dx2;
                if(j > 0)
                    rhs += dt * alpha * theta * x[idx - nx] / dy2;
                if(j < ny - 1)
                    rhs += dt * alpha * theta * x[idx + nx] / dy2;
                
                x[idx] = rhs / coeff;
            }
        }
    }
    
    memcpy(grid->data, x, n * sizeof(double));
    grid_apply_boundary_conditions(grid, 0.0);
    
    free(b);
    free(x);
}

// JNI interface implementations
JNIEXPORT jlong JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_createGrid
  (JNIEnv *env, jobject obj, jint nx, jint ny, jdouble dx, jdouble dy, jdouble alpha)
{
    Grid2D* grid = grid_create(nx, ny, dx, dy, alpha);
    return (jlong)grid;
}

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_destroyGrid
  (JNIEnv *env, jobject obj, jlong gridPtr)
{
    Grid2D* grid = (Grid2D*) gridPtr;
    grid_destroy(grid);
}

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_setGridValue
  (JNIEnv *env, jobject obj, jlong gridPtr, jint i, jint j, jdouble value)
{
    Grid2D* grid = (Grid2D*) gridPtr;
    grid_set_value(grid, i, j, value);
}

JNIEXPORT jdouble JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_getGridValue
  (JNIEnv *env, jobject obj, jlong gridPtr, jint i, jint j)
{
    Grid2D* grid = (Grid2D*) gridPtr;
    return grid_get_value(grid, i, j);
}

JNIEXPORT jdoubleArray JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_getGridData
  (JNIEnv *env, jobject obj, jlong gridPtr)
{
    Grid2D* grid = (Grid2D*) gridPtr;
    int size = grid->nx * grid->ny;
    
    jdoubleArray result = (*env)->NewDoubleArray(env, size);
    if(result == NULL)
        return NULL;
    
    (*env)->SetDoubleArrayRegion(env, result, 0, size, grid->data);
    return result;
}

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_explicitStep
  (JNIEnv *env, jobject obj, jlong gridPtr, jdouble dt)
{
    Grid2D* grid = (Grid2D*) gridPtr;
    explicit_step(grid, dt);
}

JNIEXPORT void JNICALL Java_com_alekseylopez_heatdiffusion_nativebridge_HeatSolver_implicitStep
  (JNIEnv *env, jobject obj, jlong gridPtr, jdouble dt, jdouble theta)
{
    Grid2D* grid = (Grid2D*) gridPtr;
    implicit_step(grid, dt, theta);
}