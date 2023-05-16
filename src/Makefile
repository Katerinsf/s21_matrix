FLAGS = -Wall -Werror -Wextra -std=c11
FUNC_C = helper.c s21_create.c s21_remove.c s21_eq.c s21_sum.c s21_sub.c \
		 s21_mult_num.c s21_mult_matr.c s21_transpose.c s21_calc.c s21_det.c s21_inv.c
FUNC_O = helper.o s21_create.o s21_remove.o s21_eq.o s21_sum.o s21_sub.o \
		 s21_mult_num.o s21_mult_matr.o s21_transpose.o s21_calc.o s21_det.o s21_inv.o
VAR = $(shell pkg-config --cflags --libs check)

all: s21_matrix.a

rebuild: clean all

clean:
	rm -rf *.o out *.info *.gcda *.gcno *.gcov *.gch *.out *.a *.txt test

s21_matrix.o:
	gcc $(FLAGS) -c $(FUNC_C)

s21_matrix.a: s21_matrix.o
	ar rc s21_matrix.a $(FUNC_O)
	ranlib s21_matrix.a 

test: clean s21_matrix.a
	gcc $(FLAGS) -fprofile-arcs -ftest-coverage -o test $(FUNC_C) s21_test.c -lcheck -lm -lpthread
	./test


gcov_report: test
	gcov s21_test.c
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory out

style:
	cp ../materials/linters/.clang-format ./
	clang-format -n *.c

valgrind:
	valgrind --leak-check=full -s --track-origins=yes --log-file=valgrind-out.txt ./test
# valgrind --leak-check=full -s --track-origins=yes ./test

# main: main.c
# gcc $(FLAGS) main.c $(FUNC_C) 
# # gcc $(FLAGS) main.c $$FUNC.c -lm 
# ./a.out

# verter:
# cd ../materials/build
# sh run.sh

#valgrind: s21_matrix.a
# $(CC) s21_test.c s21_matrix.a -lcheck -o test_leaks
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test_leaks
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./test



