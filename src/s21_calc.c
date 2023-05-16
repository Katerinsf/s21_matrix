#include "s21_matrix.h"

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error (mismatched matrix sizes; matrix for which
  // calculations cannot be performed, etc.)
  int error = code_ERROR;
  if (check_exist(A)) {
    if (A->rows != A->columns || (A->rows == 1 && A->rows == A->columns)) {
      error = code_ERR_CULC;
    } else if (s21_create_matrix(A->rows, A->columns, result) == code_OK) {
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = pow(-1, i + j) * minr(A, i, j);
        }
      error = code_OK;
    }
  }
  return error;
}
