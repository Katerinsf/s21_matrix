#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  // 0 - OK
  // 1 - Error, incorrect matrix
  // 2 - Calculation error
  int error = code_ERROR;
  if (rows > 0 && columns > 0 && result != NULL) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    int i;
    if (result->matrix != NULL)
      for (i = 0; i < rows; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
        if (result->matrix[i] != NULL)
          error = code_OK;
        else {
          for (int j = 0; j < i; j++) free(result->matrix[j]);
          free(result->matrix);
          break;
        }
      }
  }
  return error;
}
