#include "s21_matrix.h"

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  if (check_exist(A) &&
      s21_create_matrix(A->rows, A->columns, result) == code_OK) {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    error = code_OK;
  }
  return error;
}
