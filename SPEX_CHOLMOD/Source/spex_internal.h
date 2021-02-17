//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_internal: include file for internal use in SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_CHOLMOD.h instead.

#ifndef SPEX_CHOLMOD_INTERNAL_H
#define SPEX_CHOLMOD_INTERNAL_H

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-value"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------C Libraries------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Standard C libraries
#include <setjmp.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
// debugging
//------------------------------------------------------------------------------

#ifdef SPEX_DEBUG

#ifdef MATLAB_MEX_FILE

#define ASSERT(x)                                                             \
{                                                                             \
    if (!(x))                                                                 \
    {                                                                         \
        mexErrMsgTxt ("assertion failed: %s line %d\n", __FILE__, __LINE__) ; \
    }                                                                         \
}

#else

#include <assert.h>
#define ASSERT(x) assert (x)

#endif

#else

// debugging disabled
#define ASSERT(x)

#endif

//------------------------------------------------------------------------------
//-------------------------Common Macros----------------------------------------
//------------------------------------------------------------------------------

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

#define SPEX_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SPEX_FLIP(i) (-(i)-2)
#define SPEX_UNFLIP(i) (((i) < 0) ? SPEX_FLIP(i) : (i))
#define SPEX_MARKED(Ap,j) (Ap [j] < 0)
#define SPEX_MARK(Ap,j) { Ap [j] = SPEX_FLIP (Ap [j]) ; }

// SPEX_CHECK(method) is a macro that calls a SPEX_CHOLMOD method and checks the
// status; if a failure occurs, it frees all allocated workspace and returns
// the error status to the caller.  To use SPEX_CHECK, the #include'ing file
// must declare a SPEX_info info, and must define SPEX_FREE_ALL as a macro that
// frees all workspace if an error occurs. The method can be a scalar as well,
// so that SPEX_CHECK(info) works.

// the default is to free nothing
#ifndef SPEX_FREE_ALL
#define SPEX_FREE_ALL
#endif

#define SPEX_CHECK(method)      \
{                               \
    info = (method) ;           \
    if (info != SPEX_OK)        \
    {                           \
        SPEX_FREE_ALL ;         \
        return (info) ;         \
    }                           \
}

#include "SPEX_CHOLMOD.h"
#ifdef SPEX_CHOLMOD_TCOV
    /* include this header to make macro GOTCHA available */
    #include "../Tcov/tcov_malloc_test.h"
#endif
#ifndef GOTCHA
    #define GOTCHA
#endif

//------------------------------------------------------------------------------
// printing control
//------------------------------------------------------------------------------

// SPEX_CHOLMOD uses SuiteSparse_config.printf_func instead of a mere call to printf
// (the default function is printf, or mexPrintf when in MATLAB).  If this
// function pointer is NULL, no printing is done.

#define SPEX_PRINTF(...)                                    \
{                                                           \
    if (SuiteSparse_config.printf_func != NULL)             \
    {                                                       \
        SuiteSparse_config.printf_func (__VA_ARGS__) ;      \
    }                                                       \
}

#define SPEX_PR1(...) { if (pr >= 1) SPEX_PRINTF (__VA_ARGS__) }
#define SPEX_PR2(...) { if (pr >= 2) SPEX_PRINTF (__VA_ARGS__) }
#define SPEX_PR3(...) { if (pr >= 3) SPEX_PRINTF (__VA_ARGS__) }

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------functions for GMP wrapper----------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// uncomment this to print memory debugging info
// #define SPEX_GMP_MEMORY_DEBUG

#ifdef SPEX_GMP_MEMORY_DEBUG
void spex_gmp_dump ( void ) ;
#endif

extern int64_t spex_gmp_ntrials ;

#ifndef SPEX_GMP_LIST_INIT
// A size of 32 ensures that the list never needs to be increased in size.
// The test coverage suite in SPEX_CHOLMOD/Tcov reduces this initial size to
// exercise the code, in SPEX_CHOLMOD/Tcov/Makefile.
#define SPEX_GMP_LIST_INIT 32
#endif


bool spex_gmp_init (void) ;

void spex_gmp_finalize (void) ;

void *spex_gmp_allocate (size_t size) ;

void spex_gmp_free (void *p, size_t size) ;

void *spex_gmp_reallocate (void *p_old, size_t old_size, size_t new_size );

void spex_gmp_failure (int status) ;


// Defines printing to be done
#define SPEX_DEFAULT_PRINT_LEVEL 0

// MPFR precision used (quad is default)
#define SPEX_DEFAULT_PRECISION 128

//------------------------------------------------------------------------------
// Type of MPFR rounding used.
//------------------------------------------------------------------------------

// The MPFR library utilizes an internal rounding scheme. The options are
//  MPFR_RNDN: round to nearest (roundTiesToEven in IEEE 754-2008),
//  MPFR_RNDZ: round toward zero (roundTowardZero in IEEE 754-2008),
//  MPFR_RNDU: round toward plus infinity (roundTowardPositive in
//             IEEE 754-2008),
//  MPFR_RNDD: round toward minus infinity (roundTowardNegative in
//             IEEE 754-2008),
//  MPFR_RNDA: round away from zero.
//  MPFR_RNDF: faithful rounding. This is not stable
//
// SPEX_CHOLMOD utilizes MPFR_RNDN by default.

#define SPEX_DEFAULT_MPFR_ROUND MPFR_RNDN

//------------------------------------------------------------------------------
// Macros to utilize the default if option is NULL
//------------------------------------------------------------------------------

#define SPEX_OPTION(option,parameter,default_value) \
    ((option == NULL) ? (default_value) : (option->parameter))

#define SPEX_OPTION_PREC(option) \
    SPEX_OPTION (option, prec, SPEX_DEFAULT_PRECISION)

#define SPEX_OPTION_PRINT_LEVEL(option) \
    SPEX_OPTION (option, print_level, SPEX_DEFAULT_PRINT_LEVEL)

#define SPEX_OPTION_ROUND(option) \
    SPEX_OPTION (option, round, SPEX_DEFAULT_MPFR_ROUND)

//------------------------------------------------------------------------------
// Field access macros for MPZ/MPQ/MPFR struct
//------------------------------------------------------------------------------

// (similar definition in gmp-impl.h and mpfr-impl.h)

#define SPEX_MPZ_SIZ(x)   ((x)->_mp_size)
#define SPEX_MPZ_PTR(x)   ((x)->_mp_d)
#define SPEX_MPZ_ALLOC(x) ((x)->_mp_alloc)
#define SPEX_MPQ_NUM(x)   mpq_numref(x)
#define SPEX_MPQ_DEN(x)   mpq_denref(x)
#define SPEX_MPFR_MANT(x) ((x)->_mpfr_d)
#define SPEX_MPFR_EXP(x)  ((x)->_mpfr_exp)
#define SPEX_MPFR_PREC(x) ((x)->_mpfr_prec)
#define SPEX_MPFR_SIGN(x) ((x)->_mpfr_sign)

/*re-define but same result: */
#define SPEX_MPFR_REAL_PTR(x) (&((x)->_mpfr_d[-1]))

/* Invalid exponent value (to track bugs...) */
#define SPEX_MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Macros to set the pointer in mpz_t/mpq_t/mpfr_t variable to NULL. It is best
 * practice to call these macros immediately after mpz_t/mpq_t/mpfr_t variable
 * is declared, and before the mp*_init function is called. It would help to
 * prevent error when SPEX_MP*_CLEAR is called before the variable is
 * successfully initialized.
 */

#define SPEX_MPZ_SET_NULL(x)                \
    SPEX_MPZ_PTR(x) = NULL;                 \
    SPEX_MPZ_SIZ(x) = 0;                    \
    SPEX_MPZ_ALLOC(x) = 0;

#define SPEX_MPQ_SET_NULL(x)                     \
    SPEX_MPZ_PTR(SPEX_MPQ_NUM(x)) = NULL;        \
    SPEX_MPZ_SIZ(SPEX_MPQ_NUM(x)) = 0;           \
    SPEX_MPZ_ALLOC(SPEX_MPQ_NUM(x)) = 0;         \
    SPEX_MPZ_PTR(SPEX_MPQ_DEN(x)) = NULL;        \
    SPEX_MPZ_SIZ(SPEX_MPQ_DEN(x)) = 0;           \
    SPEX_MPZ_ALLOC(SPEX_MPQ_DEN(x)) = 0;

#define SPEX_MPFR_SET_NULL(x)               \
    SPEX_MPFR_MANT(x) = NULL;               \
    SPEX_MPFR_PREC(x) = 0;                  \
    SPEX_MPFR_SIGN(x) = 1;                  \
    SPEX_MPFR_EXP(x) = SPEX_MPFR_EXP_INVALID;

/* GMP does not give a mechanism to tell a user when an mpz, mpq, or mpfr
 * item has been cleared; thus, if mp*_clear is called on an object that
 * has already been cleared, gmp will crash. It is also not possible to
 * set a mp*_t = NULL. Thus, this mechanism modifies the internal GMP
 * size of entries to avoid crashing in the case that a mp*_t is cleared
 * multiple times.
 */

#define SPEX_MPZ_CLEAR(x)                        \
{                                                \
    if ((x) != NULL && SPEX_MPZ_PTR(x) != NULL)  \
    {                                            \
        mpz_clear(x);                            \
        SPEX_MPZ_SET_NULL(x);                    \
    }                                            \
}

#define SPEX_MPQ_CLEAR(x)                   \
{                                           \
    SPEX_MPZ_CLEAR(SPEX_MPQ_NUM(x));        \
    SPEX_MPZ_CLEAR(SPEX_MPQ_DEN(x));        \
}

#define SPEX_MPFR_CLEAR(x)                        \
{                                                 \
    if ((x) != NULL && SPEX_MPFR_MANT(x) != NULL) \
    {                                             \
        mpfr_clear(x);                            \
        SPEX_MPFR_SET_NULL(x);                    \
    }                                             \
}


// ============================================================================
//                           Internal Functions
// ============================================================================

// check if SPEX_initialize* has been called
bool spex_initialized ( void ) ;        // true if called, false if not
void spex_set_initialized (bool s) ;    // set global initialzed flag to s

// typecast a double value to int64, accounting for Infs and Nans
static inline int64_t spex_cast_double_to_int64 (double x)
{
    if (isnan (x))
    {
        return (0) ;
    }
    else if (x > INT64_MAX)
    {
        return (INT64_MAX) ;
    }
    else if (x < INT64_MIN)
    {
        return (INT64_MIN) ;
    }
    else
    {
        return ((int64_t) (x)) ;
    }
}

/*
SPEX_info spex_cast_array
(
    void *Y,                // output array, of size n
    SPEX_type ytype,        // type of Y
    void *X,                // input array, of size n
    SPEX_type xtype,        // type of X
    int64_t n,              // size of Y and X
    mpq_t y_scale,          // scale factor applied if y is mpz_t
    mpq_t x_scale,          // scale factor applied if x is mpz_t
    const SPEX_options *option
) ;

SPEX_info spex_cast_matrix
(
    SPEX_matrix **Y_handle,     // nz-by-1 dense matrix to create
    SPEX_type Y_type,           // type of Y
    SPEX_matrix *A,             // matrix with nz entries
    const SPEX_options *option
) ;

// (void *) pointer to the values of A.  A must be non-NULL with a valid type
#define SPEX_X(A)                                                           \
    ((A->type == SPEX_MPZ  ) ? (void *) A->x.mpz   :                        \
    ((A->type == SPEX_MPQ  ) ? (void *) A->x.mpq   :                        \
    ((A->type == SPEX_MPFR ) ? (void *) A->x.mpfr  :                        \
    ((A->type == SPEX_INT64) ? (void *) A->x.int64 : (void *) A->x.fp64))))


// return an error if A->kind (csc, triplet, dense) is wrong
#define SPEX_REQUIRE_KIND(A,required_kind) \
    if (A == NULL || A->kind != required_kind) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_KIND(A,required_kind) \
    ASSERT (A != NULL && A->kind == required_kind)

// return an error if A->type (mpz, mpq, mpfr, int64, or double) is wrong
#define SPEX_REQUIRE_TYPE(A,required_type) \
    if (A == NULL || A->type != required_type) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_TYPE(A,required_type) \
    ASSERT (A != NULL && A->type == required_type)

// return an error if A->kind or A->type is wrong
#define SPEX_REQUIRE(A,required_kind,required_type)     \
    SPEX_REQUIRE_KIND (A,required_kind) ;               \
    SPEX_REQUIRE_TYPE (A,required_type) ;

#define ASSERT_MATRIX(A,required_kind,required_type)    \
    ASSERT_KIND (A,required_kind) ;                     \
    ASSERT_TYPE (A,required_type) ;
*/

mpz_t* spex_create_mpz_array
(
    int64_t n            // size of the array
);

void spex_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int64_t n       // Size of x
);

mpq_t* spex_create_mpq_array
(
    int64_t n            // size of the array
);

void spex_delete_mpq_array
(
    mpq_t **x,      // mpq array to be deleted
    int64_t n       // Size of x
);

SPEX_info spex_cumsum
(
    int64_t *p,          // vector to store the sum of c
    int64_t *c,          // vector which is summed
    int64_t n            // size of c
);

/*
typedef struct
{
    int64_t n;    // number of entries
    int64_t nz;   // number of nonzeros
    int64_t *i;   // array of size nzmax that contains the column/row indices
                  // of nnz
    //bool i_shallow;// if true, i is shallow, and original i should be updated
                  // when this struct is free'd
    mpz_t *x;     // array of size n that contains the values of each entry
} spex_scattered_vector;

spex_scattered_vector* spex_create_scattered_vector
(
    const int64_t n             // number of entries in v
);

void spex_delete_scattered_vector
(
    spex_scattered_vector **sv  // scattered vector to be deleted
);*/
#define spex_scattered_vector SPEX_vector
#define spex_scattered_vector_alloc SPEX_vector_alloc
#define spex_scattered_vector_free SPEX_vector_free

SPEX_info spex_get_scattered_v
(
    spex_scattered_vector **sv_handle,// output vector in scattered form
    SPEX_vector *v,              // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const bool keep_v            // indicate if the mpz values should be kept
);

SPEX_info spex_cppu
(
    SPEX_matrix *L,  // matrix L
    SPEX_matrix *U,  // matrix U
    mpq_t *S,        // array of size 3*n that stores pending scales
    mpz_t *d,        // array of size n that stores the unscaled pivot
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    const mpq_t vk_scale,// scale factor for newly inserted column vk, which
                     // should be in col k of L in the last iteration when used.
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    const int64_t *P,// row permutation
    const int64_t *P_inv,// inverse of row permutation
    int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                       // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
);

SPEX_info spex_dppu1
(
    SPEX_matrix *L,  // matrix L
    SPEX_matrix *U,  // matrix U
    mpq_t *S,        // array of size 3*n that stores pending scales
    mpz_t *d,        // array of size n that stores the unscaled pivot
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    const mpq_t vk_scale,// scale factor for newly inserted column vk, which
                     // should be in col k of L in the last iteration when used.
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    int64_t *Ldiag,  // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                       // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
);


SPEX_info spex_dppu2
(
    SPEX_matrix *L,  // matrix L
    SPEX_matrix *U,  // matrix U
    mpq_t *S,        // array of size 3*n that stores pending scales
    mpz_t *d,        // array of size n that stores the unscaled pivot
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    const mpq_t vk_scale,// scale factor for newly inserted column vk, which
                     // should be in col k of L in the last iteration when used.
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    int64_t *Ldiag,  // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                       // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
);

SPEX_info spex_finalize_and_insert_vk
(
    spex_scattered_vector *vk_dense, //scattered version of the solution for
                      // LDx=v using the first k-1 columns of L
    int64_t *h,       // history vector for vk_dense
    SPEX_matrix *U,   // matrix U
    SPEX_matrix *L,   // matrix L
    mpq_t *S,         // array of size 3*n that stores pending scales
    mpz_t *d,         // array of unscaled pivots
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const mpz_t *sd,  // array of scaled pivots
    const int64_t *Q, // the column permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t k,  // the column index in L that vk_dense will be inserted
    const int64_t diag,// the index of entry in vk_dense that will be diagonal
    const mpq_t one   // a const mpq number, just to avoid constantly alloc
);

SPEX_info spex_find_next_nz
(
    int64_t *next,                  // the col/row index of next nnz
    spex_scattered_vector *Ak_dense,// the scattered vector
    int64_t *perm_inv,              // inverse of permutation
    int64_t k
);

SPEX_info spex_get_nnz_pattern    // find the nnz pattern of L and U
(
    // OUTPUT:
    int64_t **Ldiag,              // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    int64_t **Lr_offdiag,         // Lr_offdiag[k] gives the column index of the
                                  // last off-diagonal nnz in k-th row of L.
                                  // -1 if no off diagonal entry
    int64_t **Uci,                // the row index for col-wise nnz pattern of U
    int64_t **Ucp,                // col pointers for col-wise pattern of U
    int64_t **Ucx,                // find the value of k-th entry as
                                  // U->v[Uci[k]]->x[Ucx[k]]
    // INPUT:
    const SPEX_matrix *L,         // the target matrix L
    const SPEX_matrix *U,         // the target matrix U
    const int64_t *P,             // row permutation
    const SPEX_options *option     // command option
);

SPEX_info spex_insert_new_entry
(
    mpz_t vi,          // the entry to be inserted as i-th entry of v1
    SPEX_vector *v1,   // the vector that would add new entry
    mpq_t S1,          // pending scale for v1
    const SPEX_vector *v2,// the other vector that is in same frame as v1
    mpq_t S2,          // pending scale for v2
    mpq_t S3,          // pending scale for frame that holds v1 and v2
    mpz_t d,           // the unscale pivot in frame of v1 and v2
    const int64_t i,   // the index of vi when inserted to v1
    const int64_t v2_diag,// the pointer to the diagonal entry in v2
    const mpq_t one    // a constant mpq number, just to avoid constantly alloc
);

SPEX_info spex_ipge // perform IPGE on x based on v
(
    spex_scattered_vector *sv_x,// array of size n for x in the scattered form.
                    // x could be dense by setting sv_x->i = NULL.
    mpq_t x_scale,  // pending scale for x
    int64_t *h,     // history vector for x, x[i] was last updated in the
                    // SPEX_FLIP(h[i])-th iteration
    int64_t *prev,  // prev is the index of the found previous entry of the last
                    // one (i.e., 2nd last entry) in v(perm). update if !prev
    const SPEX_vector *v,// the vector that contains the j-th pivot used to
                    // compute x in the j-th IPGE iteration
    const int64_t *perm, // permutation
    const int64_t *perm_inv, // inverse of permutation
    const mpz_t *sd,// array of scaled pivots
    const mpq_t v_scale1, // the first pending scale for v
    const mpq_t v_scale2, // the second pending scale for v
    const mpz_t prev_unscaled_diag,// the (j-1)-th unscaled diagonal entry
    const int64_t diag_j,// x[diag_j] is the entry in v with index perm[j]
    const int64_t j // column index of v in L
);

SPEX_info spex_triangular_solve // perform REF triangular solve for LDx=v
(
    spex_scattered_vector *sv_x,// the scattered version of solution for LDx=v,
                        // using the first k-1 columns of L
    mpq_t x_scale,      // pending scale for x
    int64_t *h,         // history vector for x
    int64_t *last_update,// the number of finished IPGE iterations, which is
                        // also the number of columns in L used last time
    int64_t *i_2ndlast, // i_2ndlast is the index of the found last nnz entry
                        // of x[P] less than n, this could be NULL if not needed
    const int64_t k,    // compute x up to k-th IPGE iteration, that is, using
                        // the first k-1 columns of L
    const SPEX_matrix *L,// matrix L
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const mpq_t *S,     // the pending scale factor matrix
    const mpz_t *sd,    // array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
);

SPEX_info spex_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    mpq_t x_scale,      // pending scale for x, initially 1
    int64_t *h,         // history vector for x
    const SPEX_matrix *L,// matrix L
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const mpq_t *S,     // the pending scale factor matrix
    const mpz_t *sd,    // array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
);

SPEX_info spex_backward_sub  // performs sparse REF backward substitution
(
    SPEX_vector *x,         // right hand side vector
    mpq_t x_scale,          // pending scale for x, applying this x_scale to
                            // x will result in a non-integer value
    const SPEX_matrix *U,   // input upper triangular matrix
    const mpq_t *S,         // the pending scale factor matrix
    const mpz_t *sd,        // array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv    // inverse of column permutation
);

SPEX_info spex_verify
(   
    bool *correct,         // indicate if the verification is passed
    const SPEX_matrix *L,  // lower triangular matrix
    const SPEX_matrix *U,  // upper triangular matrix
    const SPEX_matrix *A,  // Input matrix
    int64_t *h,            // history vector
    const mpz_t *sd,       // array of scaled pivots
    const mpq_t *S,        // the pending scale factor matrix
    const int64_t *P,      // row permutation
    const int64_t *P_inv,  // inverse of row permutation
    const int64_t *Q,      // column permutation
    const int64_t *Q_inv,  // inverse of column permutation
    const int64_t *Ldiag,  // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const SPEX_options *option// command options
);

#endif

