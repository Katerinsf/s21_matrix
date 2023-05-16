#include <float.h>

#include "s21_matrix.h"

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (!check_exist(A) || !check_exist(B) || !check_match_size(A, B))
    return FAILURE;
  else
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > eps_eq) {
          return FAILURE;
        }
  return SUCCESS;
}
