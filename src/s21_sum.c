#include "s21_matrix.h"

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  if (check_exist(A) && check_exist(B)) {
    if (!check_match_size(A, B)) {
      error = code_ERR_CULC;
    } else if (s21_create_matrix(A->rows, A->columns, result) == code_OK) {
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      error = code_OK;
    }
  }
  return error;
}
