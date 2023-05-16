#include "s21_matrix.h"
// ----------------------------- Checks -----------------------------
// 1 - OK, 0 - memory was not allocated or the pointer is empty
int check_exist(matrix_t *A) {
  int error = 1;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = 0;
  } else {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i] == NULL) {
        error = 0;
        break;
      }
    }
  }
  return error;
}

int check_match_size(matrix_t *A, matrix_t *B) {
  int error = 1;
  if (A->rows != B->rows || A->columns != B->columns) error = 0;
  return error;
}

// ----------------------------- Calculations -----------------------------
// Deletes row and column for element A[k][l]
// 1 - OK, 0 - memory was not allocated or the pointer is empty
int matr_del_ij(matrix_t *A, int k, int l, matrix_t *result) {
  int res = 1;
  if (!check_exist(A) || (k < 0 || l < 0 || k >= A->rows || l >= A->columns) ||
      A->rows != A->columns || A->rows == 1 || A->columns == 1) {
    res = 0;
  } else if (s21_create_matrix(A->rows - 1, A->columns - 1, result) !=
             code_OK) {
    res = 0;
  } else {
    int p, q;
    for (int i = 0; i < result->rows; i++)
      for (int j = 0; j < result->columns; j++) {
        p = i >= k ? i + 1 : i;
        q = j >= l ? j + 1 : j;
        result->matrix[i][j] = A->matrix[p][q];
      }
  }
  return res;
}

// Minor of element A[k][l]
double minr(matrix_t *A, int k, int l) {
  double res = 0.0;
  if (check_exist(A) && (k >= 0 && l >= 0 && k < A->rows && l < A->columns) &&
      A->rows == A->columns) {
    matrix_t M;
    if (matr_del_ij(A, k, l, &M)) {
      if (s21_determinant(&M, &res) != code_OK) {
        res = 0.0;
      }
    }
    s21_remove_matrix(&M);
  }
  return res;
}

// ----------------------------- Init for tests -----------------------------
// Initializes matrix by 0
void zero_init(matrix_t *A) {
  if (check_exist(A)) {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        A->matrix[i][j] = 0;
      }
  }
}

// Initializes the matrix by element (a) on the diagonal
void diag_init(matrix_t *A, double a) {
  if (check_exist(A)) {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        if (i == j)
          A->matrix[i][j] = a;
        else
          A->matrix[i][j] = 0;
      }
  }
}

// Initializes the matrix with numbers with element numbers (counting from 1)
void index_init(matrix_t *A) {
  if (check_exist(A))
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        A->matrix[i][j] = (i + 1) * 10 + (j + 1);
      }
}

// Initializes an upper triangular matrix with a number (a)
void up_init(matrix_t *A, double a) {
  if (check_exist(A)) {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        if (i <= j) A->matrix[i][j] = a;
        // A->matrix[i][j] = (i + 1) * 10 + (j + 1);
        else
          A->matrix[i][j] = 0;
      }
  }
}

// Initializes a lower triangular matrix with number (a)
void low_init(matrix_t *A, double a) {
  if (check_exist(A)) {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        if (i >= j)
          A->matrix[i][j] = a;
        else
          A->matrix[i][j] = 0;
      }
  }
}

// ----------------------------- Print -----------------------------
// void print_matr(matrix_t *A) {
//     if (check_exist(A)) {
//         for (int i = 0; i < A->rows; i++) {
//             for (int j = 0; j < A->columns; j++) {
//                 printf("%15lf", A->matrix[i][j]);
//             }
//             printf("\n");
//         }
//     }
// }
