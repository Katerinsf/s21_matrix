#include "s21_matrix.h"

int s21_transpose(matrix_t *A, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  if (check_exist(A) &&
      s21_create_matrix(A->columns, A->rows, result) == code_OK) {
    for (int i = 0; i < result->rows; i++)
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    error = code_OK;
  }
  return error;
}
