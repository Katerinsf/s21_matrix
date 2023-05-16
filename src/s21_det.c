#include "s21_matrix.h"

int s21_determinant(matrix_t *A, double *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  if (check_exist(A)) {
    if (A->rows != A->columns) {
      error = code_ERR_CULC;
    } else if (A->rows == 1) {
      *result = A->matrix[0][0];
      error = code_OK;
    } else {
      *result = 0.0;
      for (int i = 0; i < A->rows; i++)
        *result += pow(-1, i) * A->matrix[i][0] * minr(A, i, 0);
      error = code_OK;
    }
  }
  return error;
}
