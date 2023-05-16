#include "s21_matrix.h"

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  matrix_t M, M_T;
  double det;
  if (check_exist(A)) {
    if (A->rows != A->columns || (s21_determinant(A, &det) == code_OK &&
                                  fabs(det) < eps)) {  // Проверяем det = 0
      error = code_ERR_CULC;
    } else if (A->rows == 1) {
      if (s21_create_matrix(A->rows, A->columns, result) == code_OK) {
        error = code_OK;
        result->matrix[0][0] = 1 / A->matrix[0][0];
      }
    } else if (s21_calc_complements(A, &M) == code_OK) {
      if (s21_transpose(&M, &M_T) == code_OK) {
        if (s21_mult_number(&M_T, 1 / det, result) == code_OK) {
          error = code_OK;
        }
        s21_remove_matrix(&M_T);
      }
      s21_remove_matrix(&M);
    }
  }
  return error;
}
