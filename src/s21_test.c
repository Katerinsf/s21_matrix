#include <check.h>

#include "s21_matrix.h"
// #include "check.h"

// ------------------- create ------------------
START_TEST(create_ok) {
  matrix_t A;
  int res;
  res = s21_create_matrix(2, 3, &A);
  ck_assert_int_eq(res, code_OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(create_nm_err) {
  matrix_t A;
  int res;
  res = s21_create_matrix(2, -3, &A);
  ck_assert_int_eq(res, code_ERROR);
}
END_TEST

// ------------------- eq ------------------
START_TEST(eq_diag_success) {
  matrix_t A, B;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    diag_init(&A, (i - 0.123) * 7.9 / (i + 11.88));
    diag_init(&B, (i - 0.123) * 7.9 / (i + 11.88));
    res = s21_eq_matrix(&A, &B);

    ck_assert_int_eq(res, SUCCESS);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(A.matrix[j][k], B.matrix[j][k], eps_eq);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
  }
}
END_TEST

START_TEST(eq_index_success) {
  matrix_t A, B;
  int res;
  for (int i = 1; i < 7; i++)
    for (int j = 1; j < 8; j++) {
      s21_create_matrix(i, j, &A);
      s21_create_matrix(i, j, &B);
      index_init(&A);
      index_init(&B);
      res = s21_eq_matrix(&A, &B);

      ck_assert_int_eq(res, SUCCESS);
      for (int k = 0; k < i; k++)
        for (int l = 0; l < j; l++) {
          ck_assert_double_eq_tol(A.matrix[k][l], B.matrix[k][l], eps_eq);
        }
      s21_remove_matrix(&A);
      s21_remove_matrix(&B);
    }
}
END_TEST

START_TEST(eq_nm_failure) {
  matrix_t A, B;
  int res;
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(3, 4, &B);
  zero_init(&A);
  zero_init(&B);
  res = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res, FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(eq_elem_failure) {
  matrix_t A, B;
  int res;
  s21_create_matrix(5, 3, &A);
  s21_create_matrix(5, 3, &B);
  index_init(&A);
  zero_init(&B);
  B.matrix[4][2] = 53.01;
  res = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(res, FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

// ------------------- sum ------------------
START_TEST(sum_diag_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    diag_init(&A, (i - 0.123) * 7.9 / (i + 11.88));
    diag_init(&B, (i - 0.123) * 7.9 / (i + 11.88));
    diag_init(&res_matr, 2 * (i - 0.123) * 7.9 / (i + 11.88));
    res = s21_sum_matrix(&A, &B, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(sum_index_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    index_init(&A);
    zero_init(&B);
    index_init(&res_matr);
    res = s21_sum_matrix(&A, &B, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(sum_random_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  // double a, b, r;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        // a = rand() / RAND_MAX * (5.123456 - 1.987654) + 1.987654;
        // b = rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
        // r  = a + b;
        A.matrix[j][k] =
            pow(-1, rand()) * rand() / RAND_MAX * (5.123456 - 1.987654) +
            1.987654;
        B.matrix[j][k] =
            pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) +
            3.654321;
        res_matr.matrix[j][k] = A.matrix[j][k] + B.matrix[j][k];
      }

    res = s21_sum_matrix(&A, &B, &C);
    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(sum_nm_err) {
  matrix_t A, B, C;
  int res;
  s21_create_matrix(4, 3, &A);
  s21_create_matrix(4, 1, &B);
  zero_init(&A);
  zero_init(&B);
  res = s21_sum_matrix(&A, &B, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

// ------------------- sub ------------------
START_TEST(sub_diag_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    diag_init(&A, (i - 0.123) * 7.9 / (i + 11.88));
    diag_init(&B, (i - 0.123) * 7.9 / (i + 11.88));
    zero_init(&res_matr);
    res = s21_sub_matrix(&A, &B, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(sub_index_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    index_init(&A);
    zero_init(&B);
    index_init(&res_matr);
    res = s21_sub_matrix(&A, &B, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(sub_random_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        A.matrix[j][k] =
            pow(-1, rand()) * rand() / RAND_MAX * (5.123456 - 1.987654) +
            1.987654;
        B.matrix[j][k] =
            pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) +
            3.654321;
        res_matr.matrix[j][k] = A.matrix[j][k] - B.matrix[j][k];
      }

    res = s21_sub_matrix(&A, &B, &C);
    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(sub_nm_err) {
  matrix_t A, B, C;
  int res;
  s21_create_matrix(1, 3, &A);
  s21_create_matrix(1, 5, &B);
  zero_init(&A);
  zero_init(&B);
  res = s21_sub_matrix(&A, &B, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

// ------------------- mult_num ------------------
START_TEST(mult_num_diag_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    diag_init(&A, (i - 0.123) * 7.9 / (i + 11.88));
    a = pow(-1, rand()) * rand() / RAND_MAX * (9.99 - 1.11) + 1.11;
    diag_init(&res_matr, a * (i - 0.123) * 7.9 / (i + 11.88));
    res = s21_mult_number(&A, a, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(mult_num_index_ok) {
  matrix_t A, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    index_init(&A);
    index_init(&res_matr);
    res = s21_mult_number(&A, 1, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(mult_num_random_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        A.matrix[j][k] =
            pow(-1, rand()) * rand() / RAND_MAX * (5.123 - 1.987) + 1.987;
        // a = pow(-1, rand()) * rand() / RAND_MAX * (9.99 - 1.11) + 1.11;
        // a = 1.11;
        a = rand() / RAND_MAX * (9.99 - 1.11) + 1.11;
        res_matr.matrix[j][k] = a * A.matrix[j][k];
      }

    res = s21_mult_number(&A, 1.11, &C);
    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

// ------------------- mult_matr ------------------
START_TEST(mult_matr_diag_ok) {
  matrix_t A, B, C, res_matr;
  double a, b;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    a = rand() / RAND_MAX * (5.123456 - 1.987654) + 1.987654;
    b = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    diag_init(&B, b);
    diag_init(&res_matr, a * b);
    res = s21_mult_matrix(&A, &B, &C);

    ck_assert_int_eq(res, code_OK);
    for (int n = 0; n < i; n++)
      for (int m = 0; m < i; m++) {
        ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(mult_matr_index_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    s21_create_matrix(i, i, &res_matr);
    index_init(&A);
    diag_init(&B, 1.0);
    index_init(&res_matr);
    res = s21_mult_matrix(&A, &B, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(mult_matr_simple_ok) {
  matrix_t A, B, C, res_matr;
  int res;
  s21_create_matrix(3, 2, &A);
  A.matrix[0][0] = 1;  // {1, 2}
  A.matrix[0][1] = 2;  // {2, 3}
  A.matrix[1][0] = 2;  // {3, 4}
  A.matrix[1][1] = 3;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 4;
  s21_create_matrix(2, 1, &B);
  B.matrix[0][0] = 4;  // {4}
  B.matrix[1][0] = 5;  // {5}
  s21_create_matrix(3, 1, &res_matr);
  res_matr.matrix[0][0] = 14;  // {14}
  res_matr.matrix[1][0] = 23;  // {23}
  res_matr.matrix[2][0] = 32;  // {32}
  res = s21_mult_matrix(&A, &B, &C);
  ck_assert_int_eq(res, code_OK);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < 1; m++) {
      ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(mult_matr_nm_err) {
  matrix_t A, B, C;
  int res;
  s21_create_matrix(1, 3, &A);
  s21_create_matrix(1, 5, &B);
  zero_init(&A);
  zero_init(&B);
  res = s21_mult_matrix(&A, &B, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

// ------------------- trans ------------------
START_TEST(trans_diag_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    diag_init(&res_matr, a);
    res = s21_transpose(&A, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(trans_up_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    up_init(&A, a);
    low_init(&res_matr, a);
    res = s21_transpose(&A, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(trans_low_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    low_init(&A, a);
    up_init(&res_matr, a);
    res = s21_transpose(&A, &C);

    ck_assert_int_eq(res, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(trans_random_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 1; i <= 7; i++)
    for (int j = 1; j <= 8; j++) {
      s21_create_matrix(i, j, &A);
      s21_create_matrix(j, i, &res_matr);
      for (int k = 0; k < i; k++)
        for (int l = 0; l < j; l++) {
          a = pow(-1, rand()) * rand() / RAND_MAX * (5.123 - 1.987) + 1.987;
          A.matrix[k][l] = a;
          res_matr.matrix[l][k] = a;
        }

      res = s21_transpose(&A, &C);
      ck_assert_int_eq(res, code_OK);
      for (int k = 0; k < j; k++)
        for (int l = 0; l < i; l++) {
          ck_assert_double_eq_tol(C.matrix[k][l], res_matr.matrix[k][l], eps);
        }
      s21_remove_matrix(&A);
      s21_remove_matrix(&C);
      s21_remove_matrix(&res_matr);
    }
}
END_TEST

// ------------------- minr ------------------
START_TEST(minr_diag_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  for (int i = 2; i < 5; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    diag_init(&res_matr, pow(a, i - 1));
    res = s21_calc_complements(&A, &C);

    ck_assert_int_eq(res, code_OK);
    for (int n = 0; n < i; n++)
      for (int m = 0; m < i; m++) {
        ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(minr_simple_ok_1) {
  matrix_t A, C, res_matr;
  int res;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;  // {1, 2, 3}
  A.matrix[0][1] = 2;  // {0, 4, 2}
  A.matrix[0][2] = 3;  // {5, 2, 1}
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  s21_create_matrix(3, 3, &res_matr);
  res_matr.matrix[0][0] = 0;    // {0, 10, -20}
  res_matr.matrix[0][1] = 10;   // {4, -14, 8}
  res_matr.matrix[0][2] = -20;  // {-8, -2, 4}
  res_matr.matrix[1][0] = 4;
  res_matr.matrix[1][1] = -14;
  res_matr.matrix[1][2] = 8;
  res_matr.matrix[2][0] = -8;
  res_matr.matrix[2][1] = -2;
  res_matr.matrix[2][2] = 4;
  res = s21_calc_complements(&A, &C);
  ck_assert_int_eq(res, code_OK);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++) {
      ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(minr_simple_ok_2) {
  matrix_t A, C, res_matr;
  int res;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 0.1;  // {0.1, 2.3, 1.2}
  A.matrix[0][1] = 2.3;  // {3.4, 4.5, 5.6}
  A.matrix[0][2] = 1.2;  // {7.8, 6.7, 8.9}
  A.matrix[1][0] = 3.4;
  A.matrix[1][1] = 4.5;
  A.matrix[1][2] = 5.6;
  A.matrix[2][0] = 7.8;
  A.matrix[2][1] = 6.7;
  A.matrix[2][2] = 8.9;
  s21_create_matrix(3, 3, &res_matr);
  res_matr.matrix[0][0] = 2.53;    // {2.53,  13.42, -12.32}
  res_matr.matrix[0][1] = 13.42;   // {-12.43, -8.47,  17.27}
  res_matr.matrix[0][2] = -12.32;  // {7.48,  3.52,  -7.37}
  res_matr.matrix[1][0] = -12.43;
  res_matr.matrix[1][1] = -8.47;
  res_matr.matrix[1][2] = 17.27;
  res_matr.matrix[2][0] = 7.48;
  res_matr.matrix[2][1] = 3.52;
  res_matr.matrix[2][2] = -7.37;
  res = s21_calc_complements(&A, &C);
  ck_assert_int_eq(res, code_OK);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++) {
      ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(minr_random_ok) {
  matrix_t A, C, res_matr;
  double a;
  int res;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &res_matr);
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      a = rand() / RAND_MAX * (5.123 - 1.987) + 1.987;
      A.matrix[k][l] = a;
      res_matr.matrix[1 - l][1 - k] = pow(-1, k + l) * a;
    }

  res = s21_calc_complements(&A, &C);
  ck_assert_int_eq(res, code_OK);
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      ck_assert_double_eq_tol(C.matrix[k][l], res_matr.matrix[k][l], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(minr_nm_err) {
  matrix_t A, C;
  int res;
  s21_create_matrix(2, 3, &A);
  zero_init(&A);
  res = s21_calc_complements(&A, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(minr_nm_1_err) {
  matrix_t A, C;
  int res;
  s21_create_matrix(1, 1, &A);
  zero_init(&A);
  res = s21_calc_complements(&A, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
}
END_TEST

// ------------------- det ------------------
START_TEST(det_diag_ok) {
  matrix_t A;
  double a, res;
  int error;
  for (int i = 1; i < 5; i++) {
    s21_create_matrix(i, i, &A);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    error = s21_determinant(&A, &res);

    ck_assert_int_eq(error, code_OK);
    ck_assert_double_eq_tol(res, pow(a, i), eps);
    s21_remove_matrix(&A);
  }
}
END_TEST

START_TEST(det_up_ok) {
  matrix_t A;
  double a, res;
  int error;
  for (int i = 1; i < 5; i++) {
    s21_create_matrix(i, i, &A);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    error = s21_determinant(&A, &res);

    ck_assert_int_eq(error, code_OK);
    ck_assert_double_eq_tol(res, pow(a, i), eps);
    s21_remove_matrix(&A);
  }
}
END_TEST

START_TEST(det_low_ok) {
  matrix_t A;
  double a, res;
  int error;
  for (int i = 1; i < 5; i++) {
    s21_create_matrix(i, i, &A);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    error = s21_determinant(&A, &res);

    ck_assert_int_eq(error, code_OK);
    ck_assert_double_eq_tol(res, pow(a, i), eps);
    s21_remove_matrix(&A);
  }
}
END_TEST

START_TEST(det_simple_ok_1) {
  matrix_t A;
  double res = -40.0, det;
  int error;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1.0;  // {1, 2, 3}
  A.matrix[0][1] = 2.0;  // {0, 4, 2}
  A.matrix[0][2] = 3.0;  // {5, 2, 1}
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 4.0;
  A.matrix[1][2] = 2.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = 2.0;
  A.matrix[2][2] = 1.0;
  error = s21_determinant(&A, &det);
  ck_assert_int_eq(error, code_OK);
  ck_assert_double_eq_tol(res, det, eps);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(det_simple_ok_2) {
  matrix_t A;
  double res = 16.335, det;
  int error;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 0.1;  // {0.1, 2.3, 1.2}
  A.matrix[0][1] = 2.3;  // {3.4, 4.5, 5.6}
  A.matrix[0][2] = 1.2;  // {7.8, 6.7, 8.9}
  A.matrix[1][0] = 3.4;
  A.matrix[1][1] = 4.5;
  A.matrix[1][2] = 5.6;
  A.matrix[2][0] = 7.8;
  A.matrix[2][1] = 6.7;
  A.matrix[2][2] = 8.9;
  error = s21_determinant(&A, &det);
  ck_assert_int_eq(error, code_OK);
  ck_assert_double_eq_tol(res, det, eps);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(det_simple_ok_3) {
  matrix_t A;
  double res = 0.0, x = 0.0, det;
  int error;
  int n = 3;
  s21_create_matrix(n, n, &A);
  for (int i = 0; i < n; i++)      // {1, 2, 3}
    for (int j = 0; j < n; j++) {  // {4, 5, 6}
      A.matrix[i][j] = ++x;        // {7, 8, 9}
    }
  error = s21_determinant(&A, &det);
  ck_assert_int_eq(error, code_OK);
  ck_assert_double_eq_tol(res, det, eps);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(det_nm_err) {
  matrix_t A;
  double det;
  int error;
  s21_create_matrix(2, 3, &A);
  zero_init(&A);
  error = s21_determinant(&A, &det);
  ck_assert_int_eq(error, code_ERR_CULC);
  s21_remove_matrix(&A);
}
END_TEST

// ------------------- inv ------------------
START_TEST(inv_diag_ok) {
  matrix_t A, C, res_matr;
  double a;
  int error;
  for (int i = 1; i < 5; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &res_matr);
    a = pow(-1, rand()) * rand() / RAND_MAX * (7.987654 - 3.654321) + 3.654321;
    diag_init(&A, a);
    diag_init(&res_matr, 1 / a);
    error = s21_inverse_matrix(&A, &C);

    ck_assert_int_eq(error, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(C.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&C);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(inv_simple_ok_1) {
  matrix_t A, C, res_matr;
  int res;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;  // {1, 2, 3}
  A.matrix[0][1] = 2;  // {0, 4, 2}
  A.matrix[0][2] = 3;  // {5, 2, 1}
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  s21_create_matrix(3, 3, &res_matr);
  res_matr.matrix[0][0] = 0;     // {  0,   -0.1,  0.2}
  res_matr.matrix[0][1] = -0.1;  // {-0.25, 0.35,	0.05}
  res_matr.matrix[0][2] = 0.2;   // { 0.5,  -0.2, -0.1}
  res_matr.matrix[1][0] = -0.25;
  res_matr.matrix[1][1] = 0.35;
  res_matr.matrix[1][2] = 0.05;
  res_matr.matrix[2][0] = 0.5;
  res_matr.matrix[2][1] = -0.2;
  res_matr.matrix[2][2] = -0.1;
  res = s21_inverse_matrix(&A, &C);
  ck_assert_int_eq(res, code_OK);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++) {
      ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(inv_simple_ok_2) {
  matrix_t A, C, res_matr;
  int res;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 2;  // {2,  5,  7}
  A.matrix[0][1] = 5;  // {6,  3,  4}
  A.matrix[0][2] = 7;  // {5, -2, -3}
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  s21_create_matrix(3, 3, &res_matr);
  res_matr.matrix[0][0] = 1;   // {1,   -1,   1}
  res_matr.matrix[0][1] = -1;  // {-38, 41, -34}
  res_matr.matrix[0][2] = 1;   // {27, -29,  24}
  res_matr.matrix[1][0] = -38;
  res_matr.matrix[1][1] = 41;
  res_matr.matrix[1][2] = -34;
  res_matr.matrix[2][0] = 27;
  res_matr.matrix[2][1] = -29;
  res_matr.matrix[2][2] = 24;
  res = s21_inverse_matrix(&A, &C);
  ck_assert_int_eq(res, code_OK);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++) {
      ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(inv_big_ok_1) {
  matrix_t A, C, res_matr;
  int res;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 2;  // {2,  5,  7}
  A.matrix[0][1] = 5;  // {6,  3,  4}
  A.matrix[0][2] = 7;  // {5, -2, -3}
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  s21_create_matrix(3, 3, &res_matr);
  res_matr.matrix[0][0] = 1;   // {1,   -1,   1}
  res_matr.matrix[0][1] = -1;  // {-38, 41, -34}
  res_matr.matrix[0][2] = 1;   // {27, -29,  24}
  res_matr.matrix[1][0] = -38;
  res_matr.matrix[1][1] = 41;
  res_matr.matrix[1][2] = -34;
  res_matr.matrix[2][0] = 27;
  res_matr.matrix[2][1] = -29;
  res_matr.matrix[2][2] = 24;
  res = s21_inverse_matrix(&A, &C);
  ck_assert_int_eq(res, code_OK);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++) {
      ck_assert_double_eq_tol(C.matrix[n][m], res_matr.matrix[n][m], eps);
    }
  s21_remove_matrix(&A);
  s21_remove_matrix(&C);
  s21_remove_matrix(&res_matr);
}
END_TEST

START_TEST(inv_nm_err) {
  matrix_t A, C;
  int res;
  s21_create_matrix(2, 3, &A);
  zero_init(&A);
  res = s21_inverse_matrix(&A, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(inv_det0_err) {
  matrix_t A, C;
  int res;
  s21_create_matrix(3, 3, &A);
  zero_init(&A);
  res = s21_inverse_matrix(&A, &C);
  ck_assert_int_eq(res, code_ERR_CULC);
  s21_remove_matrix(&A);
}
END_TEST

// ------------------- others ------------------
START_TEST(sum_sub) {
  matrix_t A, B, sum, res_matr;
  int res1, res2;
  for (int i = 1; i < 10; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &B);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        A.matrix[j][k] =
            pow(-1, j + k) * rand() / RAND_MAX * (5.123456 - 1.987654) +
            1.987654;
        B.matrix[j][k] =
            pow(-1, j + k) * rand() / RAND_MAX * (7.987654 - 3.654321) +
            3.654321;
      }

    res1 = s21_sum_matrix(&A, &B, &sum);
    res2 = s21_sub_matrix(&sum, &B, &res_matr);
    ck_assert_int_eq(res1, code_OK);
    ck_assert_int_eq(res2, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(A.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&sum);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

START_TEST(inv_mult) {
  matrix_t A, E, inv_matr, res_matr;
  int res1, res2;
  for (int i = 1; i < 5; i++) {
    s21_create_matrix(i, i, &A);
    s21_create_matrix(i, i, &E);
    diag_init(&E, 1.0);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        A.matrix[j][k] =
            pow(-1, j + k) * rand() / RAND_MAX * (5.123456 - 1.987654) +
            1.987654;
      }
    res1 = s21_inverse_matrix(&A, &inv_matr);
    res2 = s21_mult_matrix(&inv_matr, &A, &res_matr);
    ck_assert_int_eq(res1, code_OK);
    ck_assert_int_eq(res2, code_OK);
    for (int j = 0; j < i; j++)
      for (int k = 0; k < i; k++) {
        ck_assert_double_eq_tol(E.matrix[j][k], res_matr.matrix[j][k], eps);
      }
    s21_remove_matrix(&A);
    s21_remove_matrix(&E);
    s21_remove_matrix(&inv_matr);
    s21_remove_matrix(&res_matr);
  }
}
END_TEST

Suite *s21_matr(void) {
  Suite *s;
  TCase *tc_core;
  s = suite_create("s21_matrix");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, create_ok);
  tcase_add_test(tc_core, create_nm_err);

  tcase_add_test(tc_core, eq_diag_success);
  tcase_add_test(tc_core, eq_index_success);
  tcase_add_test(tc_core, eq_nm_failure);
  tcase_add_test(tc_core, eq_elem_failure);

  tcase_add_test(tc_core, sum_diag_ok);
  tcase_add_test(tc_core, sum_index_ok);
  tcase_add_test(tc_core, sum_random_ok);
  tcase_add_test(tc_core, sum_nm_err);

  tcase_add_test(tc_core, sub_diag_ok);
  tcase_add_test(tc_core, sub_index_ok);
  tcase_add_test(tc_core, sub_random_ok);
  tcase_add_test(tc_core, sub_nm_err);

  tcase_add_test(tc_core, mult_num_diag_ok);
  tcase_add_test(tc_core, mult_num_index_ok);
  tcase_add_test(tc_core, mult_num_random_ok);

  tcase_add_test(tc_core, mult_matr_diag_ok);
  tcase_add_test(tc_core, mult_matr_index_ok);
  tcase_add_test(tc_core, mult_matr_simple_ok);
  tcase_add_test(tc_core, mult_matr_nm_err);

  tcase_add_test(tc_core, trans_diag_ok);
  tcase_add_test(tc_core, trans_up_ok);
  tcase_add_test(tc_core, trans_low_ok);
  tcase_add_test(tc_core, trans_random_ok);

  tcase_add_test(tc_core, minr_diag_ok);
  tcase_add_test(tc_core, minr_simple_ok_1);
  tcase_add_test(tc_core, minr_simple_ok_2);
  tcase_add_test(tc_core, minr_random_ok);
  tcase_add_test(tc_core, minr_nm_err);
  tcase_add_test(tc_core, minr_nm_1_err);

  tcase_add_test(tc_core, det_diag_ok);
  tcase_add_test(tc_core, det_up_ok);
  tcase_add_test(tc_core, det_low_ok);
  tcase_add_test(tc_core, det_simple_ok_1);
  tcase_add_test(tc_core, det_simple_ok_2);
  tcase_add_test(tc_core, det_simple_ok_3);
  tcase_add_test(tc_core, det_nm_err);

  tcase_add_test(tc_core, inv_diag_ok);
  tcase_add_test(tc_core, inv_simple_ok_1);
  tcase_add_test(tc_core, inv_simple_ok_2);
  tcase_add_test(tc_core, inv_big_ok_1);
  tcase_add_test(tc_core, inv_nm_err);
  tcase_add_test(tc_core, inv_det0_err);

  tcase_add_test(tc_core, sum_sub);
  tcase_add_test(tc_core, inv_mult);

  suite_add_tcase(s, tc_core);
  return s;
}

int main(void) {
  int no_failed;
  Suite *s;
  SRunner *sr;

  s = s21_matr();
  sr = srunner_create(s);

  // srunner_set_fork_status(sr, CK_NOFORK);

  srunner_run_all(sr, CK_NORMAL);
  no_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
