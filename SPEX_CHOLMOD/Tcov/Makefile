#-------------------------------------------------------------------------------
# SPEX_CHOLMOD/Tcov/Makefile: compile and run SPEX_CHOLMOD test coverage
#-------------------------------------------------------------------------------

# SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
# Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
# SPEX_CHOLMOD/License for the license.

#-------------------------------------------------------------------------------

default: run

include ../../SuiteSparse_config/SuiteSparse_config.mk

# Linux test coverage
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter -I../Include -I../Source -I../Demo \
        -I../../SuiteSparse_config \
        -DSPEX_GMP_LIST_INIT=2 -DSPEX_CHOLMOD_TCOV -DSPEX_DEBUG

LDLIBS += -lm -lgmp -lmpfr -lsuitesparseconfig

# run all statement coverage tests, and then check for 100% coverage
run: runtests
	./covall

all: test 

CS = tcov_malloc_test.o \
    spex_backward_sub.o             spex_get_scattered_v.o \
    SPEX_calloc.o                   SPEX_gmp.o \
    spex_cppu.o \
    SPEX_create_default_options.o   SPEX_initialize.o \
    spex_create_mpq_array.o         SPEX_initialize_expert.o \
    spex_create_mpz_array.o         spex_insert_new_entry.o \
    spex_cumsum.o                   spex_ipge.o \
    spex_delete_mpq_array.o         SPEX_LUU.o \
    spex_delete_mpz_array.o         SPEX_malloc.o \
    SPEX_matrix_alloc.o \
    spex_dppu1.o                    SPEX_matrix_free.o \
    spex_dppu2.o                    SPEX_realloc.o \
    spex_finalize_and_insert_vk.o   SPEX_solve.o \
    SPEX_finalize.o                 spex_triangular_solve.o \
    spex_find_next_nz.o             SPEX_vector_alloc.o \
    spex_forward_sub.o              SPEX_vector_free.o \
    SPEX_free.o                     SPEX_vector_realloc.o \
    spex_get_nnz_pattern.o          spex_verify.o \
#spex_delete_scattered_vector.o  spex_create_scattered_vector.o  

$(CS): ../Include/SPEX_CHOLMOD.h ../Source/spex_internal.h tcov_malloc_test.h 

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_%.c:
	- ln -s ../Source/$@

spex_%.c:
	- ln -s ../Source/$@

#SLIPLU.c:
#	- ln -s ../Demo/$@

#demos.c:
#	- ln -s ../Demo/$@

test: $(CS) tcov_test.c
	$(CC) $(LDFLAGS) tcov_test.c $(CFLAGS) -o tcov_test $(CS) $(LDLIBS)

#SLIPLU: $(CS) SLIPLU.c demos.c
#	$(CC) $(LDFLAGS) SLIPLU.c demos.c $(CFLAGS) -o SLIPLU $(CS) $(LDLIBS)

# run all statement coverage tests
runtests: all
	- ./tcov_test

# To run with valgrind:
V = valgrind #--leak-check=full

# run all statement coverage tests but with valgrind
vtests: all
	- $(V) ./tcov_test

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) tcov_test *.out *.a cov.sort out
	- $(RM) -r SPEX_*.c spex_*.c *.dSYM

purge: distclean

distclean: clean

