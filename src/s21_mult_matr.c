#include "s21_matrix.h"

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  if (check_exist(A) && check_exist(B)) {
    if (A->columns != B->rows) {
      error = code_ERR_CULC;
    } else if (s21_create_matrix(A->rows, B->columns, result) == code_OK) {
      for (int i = 0; i < result->rows; i++)
        for (int j = 0; j < result->columns; j++) {
          result->matrix[i][j] = 0;
          for (int k = 0; k < A->columns; k++)
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      error = code_OK;
    }
  }
  return error;
}
