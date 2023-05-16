#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#define code_OK 0
#define code_ERROR 1  // Error, incorrect matrix
#define code_ERR_CULC 2  // Calculation error (mismatched matrix sizes; matrix for which calculations cannot be performed, etc.)
#define SUCCESS 1
#define FAILURE 0

#define eps 1e-6
#define eps_eq 1e-7

typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
} matrix_t;


// --------------- Main functions ---------------
int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

// --------------- Helper ---------------
int check_exist(matrix_t *A);
int check_match_size(matrix_t *A, matrix_t *B);

int matr_del_ij(matrix_t *A, int k, int l, matrix_t *result);
double minr(matrix_t *A, int k, int l);

void diag_init(matrix_t *A, double a);
void index_init(matrix_t *A);
void zero_init(matrix_t *A);
void up_init(matrix_t *A, double a);
void low_init(matrix_t *A, double a);

// void print_matr(matrix_t *A);

#endif // SRC_S21_MATRIX_H_
